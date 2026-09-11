"""Optional external MCMCTree reference, without an R dependency.

PAML's calibration intervals are soft priors. This adapter deliberately keeps
that distinction from the native engine's hard chronology constraints.
"""

import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit.clade_index import CladeIndex
from nwkit.fasta import parse_fasta
from nwkit.radte_model import DatingFit
from nwkit.radte_sequence import read_alignment
from nwkit.time_tree import paml_node_mapping, read_mcmctree_posterior


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


def paml_alignment(path, aliases):
    read_alignment(path, list(aliases), "dna")
    with open(path, encoding="utf-8") as handle:
        records = parse_fasta(handle)
    sequences = {
        r.name: "".join("".join(r.raw.splitlines()[1:]).split())
        .upper()
        .replace("U", "T")
        for r in records
    }
    length = len(next(iter(sequences.values())))
    return f"{len(sequences)} {length}\n" + "".join(
        f"{aliases[name]}  {seq}\n" for name, seq in sequences.items()
    )


def paml_control(args, mirrors, seed):
    model = args.substitution_model or "hky"
    if model not in {"jc69", "hky"}:
        raise ValueError("The MCMCTree reference currently supports DNA JC69 or HKY.")
    values = dict(
        seed=seed,
        seqfile="input.phy",
        treefile="input.tre",
        outfile="out.txt",
        mcmcfile="mcmc.txt",
        ndata=1,
        seqtype=0,
        usedata="2 in.BV" if args.mcmctree_likelihood == "approximate" else 1,
        clock=args.mcmctree_clock or 2,
        model=0 if model == "jc69" else 4,
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


def execute_paml(executable, work, args, *, environment=None):
    with (
        open(work / "stdout.txt", "w", encoding="utf-8") as stdout,
        open(work / "stderr.txt", "w", encoding="utf-8") as stderr,
    ):
        completed = subprocess.run(
            [executable, "mcmctree.ctl"],
            cwd=work,
            stdout=stdout,
            stderr=stderr,
            timeout=args.mcmctree_timeout or 1800,
            check=False,
            env=environment,
        )
    if completed.returncode != 0:
        raise ValueError(f"MCMCTree exited with status {completed.returncode}.")


def prepare_paml_approximation(directory, executable, args, tree, alignment, mirrors):
    if args.mcmctree_likelihood != "approximate":
        return None
    environment = os.environ.copy()
    environment["PATH"] = (
        str(Path(executable).parent) + os.pathsep + environment.get("PATH", "")
    )
    if shutil.which("baseml", path=environment["PATH"]) is None:
        raise ValueError(
            "Approximate MCMCTree requires BASEML beside MCMCTree or on PATH."
        )
    work = directory / "approximation"
    work.mkdir()
    (work / "input.tre").write_text(tree, encoding="utf-8")
    (work / "input.phy").write_text(alignment, encoding="utf-8")
    control = paml_control(args, mirrors, max(1, args.seed)).replace(
        "usedata = 2 in.BV", "usedata = 3"
    )
    (work / "mcmctree.ctl").write_text(control, encoding="utf-8")
    execute_paml(executable, work, args, environment=environment)
    summary = work / "out.BV"
    if not summary.is_file() or summary.stat().st_size == 0:
        raise ValueError("MCMCTree did not produce a nonempty out.BV approximation.")
    return summary


def run_mcmctree(c, args):
    binary = args.mcmctree_bin or "mcmctree"
    executable = shutil.which(binary)
    if executable is None:
        raise ValueError(f"MCMCTree executable not found: {binary}")
    executable = os.path.join(
        os.path.realpath(os.path.dirname(executable)), os.path.basename(executable)
    )
    parent = Path(args.out_prefix).absolute().parent
    parent.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix=".nwkit-mcmctree-", dir=parent))
    tree, aliases, mirrors = paml_tree(c)
    alignment = paml_alignment(args.alignment, aliases)
    chain_tables, age_samples = [], []
    node_ids = paml_node_mapping(c.gene)
    group_by_node = dict(zip(c.nodes, c.group_by_node, strict=True))
    count = args.mcmctree_chains or 2
    try:
        summary = prepare_paml_approximation(
            directory, executable, args, tree, alignment, mirrors
        )
        for chain in range(count):
            work = directory / f"chain-{chain + 1}"
            work.mkdir()
            (work / "input.tre").write_text(tree, encoding="utf-8")
            (work / "input.phy").write_text(alignment, encoding="utf-8")
            (work / "mcmctree.ctl").write_text(
                paml_control(args, mirrors, max(1, args.seed) + chain), encoding="utf-8"
            )
            if summary is not None:
                shutil.copyfile(summary, work / "in.BV")
            execute_paml(executable, work, args)
            table = pd.read_csv(work / "mcmc.txt", sep=r"\s+")
            expected = args.mcmctree_samples or 20000
            thin = args.mcmctree_sampfreq or 10
            initial_record = (
                len(table) == expected + 1
                and "Gen" in table
                and table.Gen.iloc[0] == 1
                and np.array_equal(
                    table.Gen.iloc[1:], np.arange(1, expected + 1) * thin
                )
            )
            # PAML 4.10.10 emits an extra initial state (Gen=1) before its
            # regularly thinned series. Keep the regular series for ESS/CI.
            if initial_record:
                table = table.iloc[1:].reset_index(drop=True)
            if len(table) != expected:
                raise ValueError(
                    "MCMCTree did not produce the requested number of samples."
                )
            posterior = read_mcmctree_posterior(
                str(work / "mcmc.txt"), c.gene, burnin=int(initial_record)
            )
            chain_tables.append(table)
            samples = np.full((len(table), len(c.groups)), np.nan)
            for node in c.gene.leaves():
                samples[:, group_by_node[node]] = 0
            for j, node_id in enumerate(posterior.node_ids):
                group = group_by_node[node_ids[node_id]]
                values = posterior.values[:, j]
                if np.isfinite(samples[:, group]).all() and not np.allclose(
                    samples[:, group], values, rtol=1e-7, atol=1e-8
                ):
                    raise ValueError(
                        "MCMCTree mirror ages differ within posterior samples."
                    )
                samples[:, group] = values
            age_samples.append(samples)
    except Exception as exc:
        raise ValueError(
            f"MCMCTree reference failed; inspect {directory}: {exc}"
        ) from exc
    samples = np.concatenate(age_samples)
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
    fit = DatingFit(
        ages,
        np.full(len(c.edges), np.nan),
        float(np.log(trace.mu.mean() / c.scale)),
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
        likelihood=args.mcmctree_likelihood or "exact",
        clock=args.mcmctree_clock or 2,
        model=args.substitution_model or "hky",
        time_scale=c.scale,
        mirror_groups=mirrors,
        diagnostics=diagnostics,
        leaf_aliases=aliases,
    )
    return fit, trace, metadata
