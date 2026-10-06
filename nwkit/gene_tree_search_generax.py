"""Independent per-topology GeneRax EVAL fits on one unchanged alignment."""

import gzip
import hashlib
import math
import os
import shlex
import signal
import subprocess
from dataclasses import dataclass
from decimal import Decimal
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

from nwkit.gene_tree_search_model import topology_key
from nwkit.reconcile import _validate_rooted_binary_tree
from nwkit.util import copy_tree_iteratively, read_tree, write_tree

# RAxML-NG inline nucleotide matrices; an amino-acid N is asparagine.
DNA_MODELS = frozenset(
    "JC K80 F81 HKY TN93EF TN93 K81 K81UF TPM2 TPM2UF TPM3 TPM3UF "
    "TIM1 TIM1UF TIM2 TIM2UF TIM3 TIM3UF TVMEF TVM SYM GTR".split()
)


def tree_text(tree, *, declaration=True):
    handle = StringIO()
    tree = copy_tree_iteratively(tree)
    for node in tree.traverse():
        if not node.is_leaf:
            node.name = None
    write_tree(tree, SimpleNamespace(outfile=handle), 1, quiet=True, props=[])
    text = handle.getvalue()
    if not declaration and text.startswith("[&R]"):
        text = text[4:]
    return text


def read_alignment(path, expected_names, *, subst_model=None):
    opener = gzip.open if str(path).endswith(".gz") else open
    records: dict[str, str] = {}
    current = None
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:]
                if (
                    not current
                    or current in records
                    or any(c.isspace() for c in current)
                ):
                    raise ValueError(
                        "Alignment FASTA requires unique, nonempty exact tip IDs without descriptions."
                    )
                records[current] = ""
            elif current is None:
                raise ValueError("Alignment must be FASTA.")
            else:
                if any(c.isspace() for c in line):
                    raise ValueError(
                        "Whitespace within alignment sequences is unsupported."
                    )
                records[current] += line.upper()
    if set(records) != set(expected_names):
        raise ValueError(
            f"Alignment/tree tip mismatch: missing={sorted(set(expected_names) - records.keys())}, extra={sorted(records.keys() - set(expected_names))}"
        )
    lengths = {len(s) for s in records.values()}
    if len(lengths) != 1 or 0 in lengths:
        raise ValueError("Alignment sequences must have one common nonzero length.")
    matrix = (subst_model or "").split("+", 1)[0].split("{", 1)[0].upper()
    missing = {"-", "?", "X"}
    if matrix in DNA_MODELS or matrix.startswith("DNA"):
        missing.add("N")
    if any(not set(s) - missing for s in records.values()):
        raise ValueError("Alignment contains a sequence with no observed characters.")
    return records


@dataclass(frozen=True)
class Evaluation:
    sequence_log_likelihood: float
    reconciliation_log_likelihood: float
    optimized_tree: str
    rounding_bound: float = 0.0

    @property
    def joint(self):
        return self.sequence_log_likelihood + self.reconciliation_log_likelihood


def _parse_result(prefix, candidate, *, root_policy):
    directory = prefix / "results" / candidate.id
    stats = directory / "stats.txt"
    numbers = stats.read_text().splitlines()[0].split()
    if len(numbers) != 2:
        raise ValueError(f"Unexpected GeneRax stats format: {stats}")
    seq, rec = map(float, numbers)
    if not all(math.isfinite(v) for v in (seq, rec)):
        raise ValueError(f"Nonfinite GeneRax likelihood: {stats}")
    tree = read_tree(directory / "geneTree.newick", 1, True, quiet=True)
    _validate_rooted_binary_tree(tree, "GeneRax fitted tree")
    if any(n.dist is not None and n.dist < 0 for n in tree.traverse()):
        raise ValueError(f"Negative GeneRax fitted branch length for {candidate.id}.")
    if topology_key(tree, rooted=root_policy == "keep") != topology_key(
        candidate.tree, rooted=root_policy == "keep"
    ):
        raise ValueError(f"GeneRax EVAL changed topology/tips for {candidate.id}.")
    # GeneRax stats may be rounded. The displayed decimal place gives a
    # conservative bound even when trailing zeroes have been omitted.
    bounds = []
    for value in numbers:
        exponent = Decimal(value).as_tuple().exponent
        if not isinstance(exponent, int):
            raise ValueError(f"Nonfinite GeneRax likelihood: {stats}")
        bounds.append(10.0**exponent / 2)
    rounding_bound = sum(bounds, 0.0)
    return Evaluation(seq, rec, tree_text(tree), rounding_bound)


def _run_generax(argv, *, cwd, log, timeout):
    # MPI launchers can leave workers alive if only the launcher is killed.
    # A private POSIX group also stops workers if the launcher fails first.
    with subprocess.Popen(
        argv,
        cwd=cwd,
        stdin=subprocess.DEVNULL,
        stdout=log,
        stderr=subprocess.STDOUT,
        start_new_session=os.name == "posix",
    ) as process:
        try:
            returncode = process.wait(timeout=timeout)
            if returncode:
                raise subprocess.CalledProcessError(returncode, argv)
        except BaseException:
            if os.name == "posix":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                process.kill()
            process.wait()
            raise


def _evaluate_round(
    candidates,
    workdir,
    aligned,
    species,
    mapfile,
    *,
    command,
    subst_model,
    rec_model,
    root_policy,
    seed,
    timeout,
    starting_text,
):
    workdir.mkdir()
    families = ["[FAMILIES]"]
    for candidate in candidates:
        tree_path = workdir / f"{candidate.id}.nwk"
        tree_path.write_text(starting_text[candidate.id].removeprefix("[&R]"))
        families.extend(
            (
                f"- {candidate.id}",
                f"starting_gene_tree = {tree_path.name}",
                f"alignment = {os.path.relpath(aligned, workdir)}",
                f"mapping = {os.path.relpath(mapfile, workdir)}",
                f"subst_model = {subst_model}",
            )
        )
    family_path = workdir / "families.txt"
    family_path.write_text("\n".join(families) + "\n")
    prefix = workdir / "eval"
    argv = [
        *shlex.split(command),
        "--families",
        family_path.name,
        "--species-tree",
        os.path.relpath(species, workdir),
        "--prefix",
        prefix.name,
        "--strategy",
        "EVAL",
        "--rec-model",
        rec_model,
        "--per-family-rates",
        "--rec-weight",
        "1.0",
        "--seed",
        str(seed),
    ]
    if root_policy == "keep":
        argv.append("--enforce-gene-tree-root")
    with (workdir / "generax.log").open("w") as log:
        _run_generax(
            argv,
            cwd=workdir,
            log=log,
            timeout=timeout,
        )
    results = {
        c.id: _parse_result(prefix, c, root_policy=root_policy) for c in candidates
    }
    return results, {
        "command": argv,
        "cwd": str(workdir),
        "log": str(workdir / "generax.log"),
        "likelihoods": {
            name: {
                "sequence": fit.sequence_log_likelihood,
                "reconciliation": fit.reconciliation_log_likelihood,
                "rounding_bound": fit.rounding_bound,
            }
            for name, fit in results.items()
        },
    }


def evaluate_candidates(
    candidates,
    species_tree,
    mapping,
    alignment,
    workdir,
    *,
    command="generax",
    subst_model,
    rec_model="UndatedDL",
    root_policy="optimize",
    seed=12345,
    timeout=3600,
    rounds=2,
):
    workdir = Path(workdir).resolve()
    # Never reuse a directory containing old scores or overwrite an input there.
    names = tuple(sorted(mapping))
    records = read_alignment(alignment, names, subst_model=subst_model)
    if len(names) < 3:
        raise ValueError("GeneRax evaluation requires at least three gene tips.")
    if not subst_model or any(c.isspace() for c in subst_model):
        raise ValueError(
            "--subst-model must be an explicit GeneRax model without whitespace."
        )
    matrix = subst_model.split("+", 1)[0].split("{", 1)[0]
    if os.path.isfile(subst_model) or "/" in matrix or "\\" in matrix:
        raise ValueError(
            "--subst-model requires an inline model, not a model-file path."
        )
    if any(
        any(c.isspace() or c in "(),:;[]" for c in label)
        for label in (*mapping, *species_tree.leaf_names())
    ):
        raise ValueError(
            "GeneRax labels cannot contain whitespace, commas, colons, semicolons, parentheses or brackets."
        )
    if not shlex.split(command) or rounds < 1:
        raise ValueError(
            "GeneRax command must not be empty and EVAL rounds must be positive."
        )
    if "#" in subst_model:
        raise ValueError("GeneRax family-file models cannot contain '#'.")
    workdir.mkdir(parents=True, exist_ok=False)
    aligned = workdir / "alignment.fasta"
    aligned.write_text("".join(f">{n}\n{records[n]}\n" for n in names))
    species = workdir / "species.nwk"
    species.write_text(tree_text(species_tree, declaration=False))
    mapfile = workdir / "mapping.tsv"
    mapfile.write_text("".join(f"{n}\t{mapping[n]}\n" for n in names))
    starting_text = {c.id: tree_text(c.tree, declaration=False) for c in candidates}
    history = []
    best_evaluations: dict[str, Evaluation] = {}
    best_rounds = {}
    for number in range(1, rounds + 1):
        evaluations, record = _evaluate_round(
            candidates,
            workdir / f"round-{number:03d}",
            aligned,
            species,
            mapfile,
            command=command,
            subst_model=subst_model,
            rec_model=rec_model,
            root_policy=root_policy,
            seed=seed,
            timeout=timeout,
            starting_text=starting_text,
        )
        history.append(record)
        for name, fit in evaluations.items():
            previous = best_evaluations.get(name)
            if previous is None or fit.joint > previous.joint:
                best_evaluations[name] = fit
                best_rounds[name] = number
        starting_text = {
            name: fit.optimized_tree for name, fit in best_evaluations.items()
        }
    metadata = {
        "rounds": history,
        "best_round_per_topology": best_rounds,
        "alignment_sha256": hashlib.sha256(aligned.read_bytes()).hexdigest(),
        "alignment_tips": len(records),
        "alignment_sites": len(next(iter(records.values()))),
        "substitution_model": subst_model,
        "reconciliation_model": rec_model,
        "root_policy": root_policy,
        "per_family_rates": True,
        "rec_weight": 1.0,
    }
    return best_evaluations, metadata
