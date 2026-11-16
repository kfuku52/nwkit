"""Optional external MCMCTree reference, without an R dependency.

PAML's calibration intervals are soft priors. This adapter deliberately keeps
that distinction from the native engine's hard chronology constraints.
"""

import hashlib
import os
import shutil
import signal
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib.resources import files
from pathlib import Path
from threading import Event

import numpy as np
import pandas as pd

from nwkit.clade_index import CladeIndex
from nwkit.radte_model import DatingFit
from nwkit.radte_sequence import AA_CODES, DNA_CODES, read_alignment_sequences
from nwkit.time_tree import (
    mcmctree_posterior_from_table,
    paml_node_mapping,
    read_mcmctree_table,
)

PAML_MODELS = {"jc69": 0, "hky": 4, "lg": 2}
DNA_LETTERS = str.maketrans("", "", "".join(DNA_CODES))
AA_LETTERS = str.maketrans("", "", "".join(AA_CODES))
PAML_DNA_TRANSLATION = str.maketrans("UX", "TN")
PAML_AA_TRANSLATION = str.maketrans("BJZ", "XXX")


def paml_tree(c):
    aliases = {str(n.name): f"g{i}" for i, n in enumerate(c.gene.leaves(), start=1)}
    by_id = c.events.set_index("gene_clade_id").to_dict("index")
    index = CladeIndex(c.gene)
    members: dict[str, list] = {}
    for node in c.gene.traverse(strategy="preorder"):
        rec = by_id[index.clade_id_for_node(node)]
        if rec["event_type"] == "speciation":
            members.setdefault(rec["shared_age_id"], []).append(node)
    mirrors = {
        key: i
        for i, (key, nodes) in enumerate(
            ((key, nodes) for key, nodes in members.items() if len(nodes) > 1), start=1
        )
    }
    annotations = {}
    for node in c.gene.traverse(strategy="preorder"):
        if node.is_leaf:
            continue
        rec = by_id[index.clade_id_for_node(node)]
        group = rec["shared_age_id"]
        pieces = []
        if group in mirrors:
            pieces.append(f"#{mirrors[group]}")
        selected = node is c.gene or rec["event_type"] == "speciation"
        driver = group not in mirrors or node is members[group][0]
        if selected and driver:
            lower, upper = rec["age_min"] / c.scale, rec["age_max"] / c.scale
            if lower == upper:
                epsilon = max(1e-8, abs(lower) * 1e-6)
                lower, upper = max(0, lower - epsilon), upper + epsilon
            if node is c.gene:
                pieces.append(f">{lower:.15g}<{upper:.15g}")
            else:
                pieces.append(f"B{{{lower:.15g},{upper:.15g}}}")
        if pieces:
            annotation = " ".join(pieces)
            if node is not c.gene and selected and driver:
                annotation = "[" + annotation + "]"
            annotations[node] = annotation
    tokens, stack = [], [c.gene]
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            tokens.append(item)
        elif item.is_leaf:
            tokens.append(aliases[str(item.name)])
        else:
            stack.extend(
                [
                    annotations.get(item, ""),
                    ")",
                    item.children[1],
                    ",",
                    item.children[0],
                    "(",
                ]
            )
    return f"{len(aliases)} 1\n" + "".join(tokens) + ";\n", aliases, len(mirrors)


def paml_alignment(sequences, aliases, alphabet):
    protein = alphabet == "protein"
    letters = AA_LETTERS if protein else DNA_LETTERS
    # PAML uses N for unknown DNA bases and X for unknown protein residues.
    translation = PAML_AA_TRANSLATION if protein else PAML_DNA_TRANSLATION
    length = len(next(iter(sequences.values())))
    rows = [f"{len(sequences)} {length}\n"]
    for name, seq in sequences.items():
        invalid = seq.translate(letters)
        if invalid:
            raise ValueError(
                f"Unsupported {alphabet} alignment character: {invalid[0]}"
            )
        rows.append(f"{aliases[name]}  {seq.translate(translation)}\n")
    return "".join(rows)


def paml_sequence_settings(args, sequences):
    model = args.substitution_model
    if model is None:
        model = (
            "lg"
            if any(seq.translate(DNA_LETTERS) for seq in sequences.values())
            else "hky"
        )
    if model not in PAML_MODELS:
        raise ValueError("MCMCTree supports JC69 or HKY for DNA and LG for protein.")
    if model == "lg":
        if args.mcmctree_likelihood == "exact":
            raise ValueError(
                "MCMCTree requires --mcmctree-likelihood approximate for protein alignments."
            )
        return model, "approximate"
    return model, args.mcmctree_likelihood or "exact"


def paml_control(args, mirrors, seed, *, model, usedata):
    values = dict(
        seed=seed,
        seqfile="input.phy",
        treefile="input.tre",
        outfile="out.txt",
        mcmcfile="mcmc.txt",
        ndata=1,
        seqtype=2 if model == "lg" else 0,
        usedata="2 in.BV" if usedata == 2 else usedata,
        clock=args.mcmctree_clock or 2,
        model=PAML_MODELS[model],
        alpha=0,
        ncatG=4,
        cleandata=0,
        BDparas="1 1 0.1 M",
        rgene_gamma=args.mcmctree_rate_prior or "2 20 1",
        sigma2_gamma=args.mcmctree_variance_prior or "1 10 1",
        kappa_gamma="6 2",
        alpha_gamma="1 1",
        burnin=args.mcmctree_burnin if args.mcmctree_burnin is not None else 2000,
        sampfreq=args.mcmctree_sampfreq or 10,
        nsample=args.mcmctree_samples or 20000,
        duplication=int(mirrors > 0),
        print=1,
    )
    if model == "lg" and usedata == 3:
        values["aaRatefile"] = "lg.dat"
    return "".join(f"{key} = {value}\n" for key, value in values.items())


def chain_diagnostics(chains):
    """Split R-hat and conservative per-chain initial-positive-sequence ESS.

    These are basic trace diagnostics, not a guarantee of convergence. They
    are reported on ages, rate parameters, and log likelihood when available.
    """
    count = min(len(chain) for chain in chains)
    if len(chains) < 2 or count < 20:
        return dict(
            status="insufficient-chains-or-samples", max_split_rhat=None, min_ess=None
        )
    columns = [col for col in chains[0] if col != "Gen"]
    values = np.array([chain[columns].iloc[:count].to_numpy(float) for chain in chains])
    half = count // 2
    split = np.concatenate([values[:, :half], values[:, -half:]], axis=0)
    within = np.var(split, axis=1, ddof=1).mean(axis=0)
    between = half * np.var(split.mean(axis=1), axis=0, ddof=1)
    variance = (half - 1) / half * within + between / half
    valid = within > 0
    rhat = np.sqrt(
        np.divide(variance, within, out=np.full_like(within, np.inf), where=valid)
    )
    per_chain_ess = []
    for chain in values:
        centered = chain - chain.mean(axis=0)
        transform = np.fft.rfft(centered, n=2 * count, axis=0)
        acov = np.fft.irfft(transform * transform.conj(), n=2 * count, axis=0)[:count]
        acov /= np.arange(count, 0, -1)[:, None]
        autocorrelation = np.divide(
            acov, acov[:1], out=np.zeros_like(acov), where=acov[:1] > 0
        )
        pair = (
            autocorrelation[1 : 2 * ((count - 1) // 2) + 1]
            .reshape(-1, 2, len(columns))
            .sum(axis=1)
        )
        active = np.cumprod(pair > 0, axis=0)
        tau = 1 + 2 * np.sum(
            np.minimum.accumulate(np.maximum(pair, 0), axis=0) * active, axis=0
        )
        per_chain_ess.append(count / np.maximum(1, tau))
    ess = np.sum(per_chain_ess, axis=0)
    maximum = float(np.max(rhat))
    return dict(
        status="passed-basic-diagnostics"
        if maximum < 1.05 and float(ess.min()) >= 200
        else "convergence-not-established",
        max_split_rhat=maximum if np.isfinite(maximum) else None,
        min_ess=float(ess.min()),
        columns=columns,
    )


def _wait_for_paml(process, timeout, cancelled):
    if cancelled is None:
        return process.wait(timeout=timeout)
    deadline = time.monotonic() + timeout
    while not cancelled.is_set():
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise subprocess.TimeoutExpired(process.args, timeout)
        try:
            return process.wait(timeout=min(remaining, 0.25))
        except subprocess.TimeoutExpired:
            continue
    raise ValueError("MCMCTree chain cancelled.")


def execute_paml(executable, work, args, *, environment=None, cancelled=None):
    with (
        open(work / "stdout.txt", "w", encoding="utf-8") as stdout,
        open(work / "stderr.txt", "w", encoding="utf-8") as stderr,
        subprocess.Popen(
            [executable, "mcmctree.ctl"],
            cwd=work,
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            env=environment,
            start_new_session=os.name == "posix",
        ) as process,
    ):
        try:
            status = _wait_for_paml(process, args.mcmctree_timeout or 1800, cancelled)
            if status:
                raise ValueError(f"MCMCTree {work.name} exited with status {status}.")
        except BaseException:
            # During approximation MCMCTree launches CODEML or BASEML itself.
            if os.name == "posix":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                process.kill()
            process.wait()
            raise


def prepare_paml_approximation(
    directory, executable, args, tree, alignment, mirrors, *, model
):
    environment = os.environ.copy()
    environment["PATH"] = (
        str(Path(executable).parent) + os.pathsep + environment.get("PATH", "")
    )
    engine = "codeml" if model == "lg" else "baseml"
    if shutil.which(engine, path=environment["PATH"]) is None:
        raise ValueError(
            f"Approximate MCMCTree requires {engine.upper()} beside MCMCTree or on PATH."
        )
    work = directory / "approximation"
    work.mkdir()
    (work / "input.tre").write_text(tree, encoding="utf-8")
    (work / "input.phy").write_text(alignment, encoding="utf-8")
    if model == "lg":
        lg = files("nwkit").joinpath("data_model/lg.txt").read_text(encoding="utf-8")
        matrix = "\n".join(line for line in lg.splitlines() if not line.startswith("#"))
        (work / "lg.dat").write_text(matrix + "\n", encoding="utf-8")
    control = paml_control(args, mirrors, max(1, args.seed), model=model, usedata=3)
    (work / "mcmctree.ctl").write_text(control, encoding="utf-8")
    execute_paml(executable, work, args, environment=environment)
    summary = work / "out.BV"
    if not summary.is_file() or summary.stat().st_size == 0:
        raise ValueError("MCMCTree did not produce a nonempty out.BV approximation.")
    return summary


def _read_chain(work, args, gene):
    table = read_mcmctree_table(work / "mcmc.txt")
    expected = args.mcmctree_samples or 20000
    thin = args.mcmctree_sampfreq or 10
    required = ["Gen", "mu", "lnL"]
    if (args.mcmctree_clock or 2) > 1:
        required.append("sigma2")
    if any(column not in table for column in required):
        raise ValueError(
            "MCMCTree trace is missing required generation or model columns."
        )
    # Recent PAML versions emit a generation-1 row before the thinned series.
    # With sampling frequency one that row is already part of the series.
    initial_record = int(len(table) == expected + 1 and table.Gen.iloc[0] == 1)
    retained = table.iloc[initial_record:]
    if len(retained) != expected or not np.array_equal(
        retained.Gen, np.arange(1, expected + 1) * thin
    ):
        raise ValueError("MCMCTree did not produce the requested sample generations.")
    values = retained[required].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (retained.mu < 0).any():
        raise ValueError("MCMCTree trace requires finite values and nonnegative rates.")
    if "sigma2" in retained and (retained.sigma2 < 0).any():
        raise ValueError("MCMCTree trace contains negative rate variances.")
    return retained, mcmctree_posterior_from_table(table, gene, burnin=initial_record)


def _run_chains(executable, works, args):
    workers = min(len(works), args.mcmctree_workers or 2, os.cpu_count() or 1)
    if workers == 1:
        for work in works:
            execute_paml(executable, work, args)
        return
    with ThreadPoolExecutor(max_workers=workers) as pool:
        cancelled = Event()
        pending = [
            pool.submit(execute_paml, executable, work, args, cancelled=cancelled)
            for work in works
        ]
        try:
            for finished in as_completed(pending):
                finished.result()
        finally:
            cancelled.set()
            for future in pending:
                future.cancel()


def _sample_ages(c, posterior, group_by_node, node_ids, samples):
    seen = set()
    for node in c.gene.leaves():
        samples[:, group_by_node[node]] = 0
    for j, node_id in enumerate(posterior.node_ids):
        group = group_by_node[node_ids[node_id]]
        values = posterior.values[:, j]
        if group in seen and not np.allclose(
            samples[:, group], values, rtol=1e-7, atol=1e-8
        ):
            raise ValueError("MCMCTree mirror ages differ within posterior samples.")
        samples[:, group] = values
        seen.add(group)


def run_mcmctree(c, args):
    binary = args.mcmctree_bin or "mcmctree"
    executable = shutil.which(binary)
    if executable is None:
        raise ValueError(f"MCMCTree executable not found: {binary}")
    executable = os.path.join(
        os.path.realpath(os.path.dirname(executable)), os.path.basename(executable)
    )
    tree, aliases, mirrors = paml_tree(c)
    sequences = read_alignment_sequences(args.alignment, aliases)
    model, likelihood = paml_sequence_settings(args, sequences)
    alignment = paml_alignment(
        sequences, aliases, "protein" if model == "lg" else "dna"
    )
    del sequences
    parent = Path(args.out_prefix).absolute().parent
    parent.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix=".nwkit-mcmctree-", dir=parent))
    chain_tables = []
    node_ids = paml_node_mapping(c.gene)
    group_by_node = dict(zip(c.nodes, c.group_by_node, strict=True))
    count = args.mcmctree_chains or 2
    expected = args.mcmctree_samples or 20000
    try:
        summary = (
            prepare_paml_approximation(
                directory, executable, args, tree, alignment, mirrors, model=model
            )
            if likelihood == "approximate"
            else None
        )
        works = []
        for chain in range(count):
            work = directory / f"chain-{chain + 1}"
            work.mkdir()
            (work / "input.tre").write_text(tree, encoding="utf-8")
            (work / "input.phy").write_text(alignment, encoding="utf-8")
            (work / "mcmctree.ctl").write_text(
                paml_control(
                    args,
                    mirrors,
                    max(1, args.seed) + chain,
                    model=model,
                    usedata=2 if summary is not None else 1,
                ),
                encoding="utf-8",
            )
            if summary is not None:
                shutil.copyfile(summary, work / "in.BV")
            works.append(work)
        _run_chains(executable, works, args)
        samples = np.full((count * expected, len(c.groups)), np.nan)
        for chain, work in enumerate(works):
            table, posterior = _read_chain(work, args, c.gene)
            _sample_ages(
                c,
                posterior,
                group_by_node,
                node_ids,
                samples[chain * expected : (chain + 1) * expected],
            )
            chain_tables.append(table)
    except Exception as exc:
        raise ValueError(
            f"MCMCTree reference failed; inspect {directory}: {exc}"
        ) from exc
    sampled = np.isfinite(samples).all(axis=0)
    ages, lo, hi = (
        c.initial.copy(),
        np.full(len(c.groups), np.nan),
        np.full(len(c.groups), np.nan),
    )
    ages[sampled] = samples[:, sampled].mean(axis=0)
    tail = (1 - args.interval_level) / 2
    lo[sampled], hi[sampled] = np.quantile(
        samples[:, sampled], [tail, 1 - tail], axis=0
    )
    if np.any(c.durations(ages) <= 0):
        raise ValueError("MCMCTree posterior mean has nonpositive gene branches.")
    diagnostics = chain_diagnostics(chain_tables)
    trace = pd.concat(
        [table.assign(chain=i + 1) for i, table in enumerate(chain_tables)],
        ignore_index=True,
    )
    mean_rate = trace.mu.mean() / c.scale
    if not np.isfinite(mean_rate) or mean_rate <= 0:
        raise ValueError(
            f"MCMCTree reference failed; inspect {directory}: invalid mean rate."
        )
    fit = DatingFit(
        ages,
        np.full(len(c.edges), np.nan),
        float(np.log(mean_rate)),
        float(np.sqrt(trace.sigma2.mean())) if "sigma2" in trace else 0.0,
        -float(trace.lnL.mean()),
        np.array([]),
        diagnostics["status"] == "passed-basic-diagnostics",
        [],
        [diagnostics["status"], "PAML-soft-priors-not-native-hard-bounds"],
        lo,
        hi,
        "mcmc-equal-tail",
        samples,
    )
    metadata = dict(
        work_directory=str(directory),
        executable=executable,
        executable_sha256=hashlib.sha256(Path(executable).read_bytes()).hexdigest(),
        likelihood=likelihood,
        clock=args.mcmctree_clock or 2,
        model=model,
        time_scale=c.scale,
        mirror_groups=mirrors,
        diagnostics=diagnostics,
        leaf_aliases=aliases,
    )
    return fit, trace, metadata
