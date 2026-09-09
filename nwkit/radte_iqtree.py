"""IQ-TREE likelihood/derivative adapter; reconciliation and dating stay in NWKIT.

Uses a persistent IQ-TREE worker by default; IQ2MC exports supply the initial
unclocked prefit and the explicit subprocess mode. MCMCTree is never launched.
Both protocols verify lengths and map bipartitions to the reconciled rooted tree.
"""

import gzip
import re
import shutil
import subprocess
import tempfile
from collections import OrderedDict
from pathlib import Path

import numpy as np
from ete4 import Tree

from nwkit.fasta import parse_fasta
from nwkit.radte_sequence import read_alignment

BASES = {
    "JC": "dna",
    "HKY": "dna",
    "GTR": "dna",
    "F81": "dna",
    "POISSON": "protein",
    "LG": "protein",
    "WAG": "protein",
    "JTT": "protein",
    "GY": "codon",
    "MG": "codon",
    "ECMK07": "codon",
    "ECMREST": "codon",
}
ALIASES = {"jc69": "JC", "gy94": "GY", "lg-f": "LG+F", "ecmrest": "ECMrest"}
TOKEN = re.compile(r"([A-Za-z0-9]+)(?:\{([0-9.eE,+-]+)\})?")


def model_tokens(model):
    tokens = TOKEN.findall(model)
    reconstructed = "+".join(
        name + ("{" + values + "}" if values else "") for name, values in tokens
    )
    if (
        reconstructed.upper() != model.upper()
        or not tokens
        or tokens[0][0].upper() not in BASES
    ):
        raise ValueError(
            "Unsupported IQ-TREE dating model; use a supported reversible single model with +F/+I/+Gk/+Rk."
        )
    base, parameters = tokens[0]
    if parameters:
        allowed_counts = {"GY": {1, 2}, "MG": {1}, "HKY": {1}, "GTR": {6}}
        if len(parameters.split(",")) not in allowed_counts.get(base.upper(), set()):
            raise ValueError("Unsupported IQ-TREE base-model parameter list.")
    seen = set()
    for name, parameters in tokens[1:]:
        name = name.upper()
        kind = "F" if name in {"F", "FQ", "F1X4", "F3X4"} else name[:1]
        if (
            kind not in {"F", "I", "G", "R"}
            or kind in seen
            or (kind in {"G", "R"} and not re.fullmatch(r"[GR][2-9][0-9]*", name))
            or (kind == "I" and name != "I")
            or (kind == "F" and name not in {"F", "FQ", "F1X4", "F3X4"})
        ):
            raise ValueError("Unsupported or repeated IQ-TREE model modifier: " + name)
        if kind == "F" and parameters:
            raise ValueError(
                "Explicit frequency vectors are not supported by this IQ-TREE adapter."
            )
        seen.add(kind)
    if {"G", "R"} <= seen:
        raise ValueError("Choose Gamma or FreeRate, not both.")
    return tokens


def requested_model(args):
    if getattr(args, "iqtree_model", None):
        controls = (
            "kappa",
            "omega",
            "codon_frequencies",
            "gtr_exchangeabilities",
            "gamma_shape",
            "gamma_categories",
        )
        if any(getattr(args, key, None) is not None for key in controls):
            raise ValueError(
                "--iqtree-model specifies the complete model; omit separate frequency/rate/model-parameter controls."
            )
        model_tokens(args.iqtree_model)
        return args.iqtree_model
    from nwkit.radte_sequence_fit import default_sequence_model

    base = args.substitution_model or default_sequence_model(args.alignment)
    model = ALIASES.get(base, base.upper())
    if base == "gy94":
        if (args.omega is None) != (args.kappa is None):
            raise ValueError(
                "IQ-TREE GY94 requires both --omega and --kappa fixed, or neither."
            )
        if args.omega is not None:
            model += f"{{{args.omega},{args.kappa}}}"
    elif args.omega is not None:
        raise ValueError("--omega requires GY94.")
    if args.kappa is not None and base != "gy94":
        if base != "hky":
            raise ValueError("--kappa applies only to HKY or GY94.")
        model += f"{{{args.kappa}}}"
    if args.gtr_exchangeabilities is not None:
        if base != "gtr":
            raise ValueError("--gtr-exchangeabilities requires GTR.")
        model += "{" + args.gtr_exchangeabilities + "}"
    codon = base in {"gy94", "ecmk07", "ecmrest"}
    if codon:
        freq = args.codon_frequencies or ("f3x4" if base == "gy94" else "model")
        if freq != "model":
            model += "+" + freq.upper()
    elif args.codon_frequencies is not None or args.genetic_code is not None:
        raise ValueError("Codon controls require a codon model.")
    categories = args.gamma_categories if args.gamma_categories is not None else 4
    if categories < 1:
        raise ValueError("Gamma categories must be positive.")
    if categories > 1:
        model += f"+G{categories}"
        if args.gamma_shape is not None:
            model += f"{{{args.gamma_shape}}}"
    model_tokens(model)
    return model


def checkpoint(path):
    """Read only the scalar/array fields needed from IQ-TREE's checkpoint."""
    result = {}
    section = ""
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(" ") and ": " in line:
                key, value = line.lstrip().rstrip("\n").split(": ", 1)
                result[section + "." + key] = value
            elif line.rstrip().endswith(":"):
                section = line.strip()[:-1]
    return result


def freeze_model(model, data):
    tokens = model_tokens(model)
    base, values = tokens[0]
    upper = base.upper()
    if upper == "GY":
        # Even a user-fixed omega alone leaves kappa free in IQ-TREE.
        values = data["ModelCodon.omega"] + "," + data["ModelCodon.kappa"]
    if not values:
        if upper in {"GY", "MG"}:
            values = data["ModelCodon.omega"]
            if upper == "GY":
                values += "," + data["ModelCodon.kappa"]
        elif upper in {"GTR", "HKY"}:
            rates = [float(x) for x in data["ModelDNA.rates"].split(",")]
            values = (
                ",".join(format(x, ".12g") for x in rates)
                if upper == "GTR"
                else format(rates[1] / rates[0], ".12g")
            )

    def rate_value(field):
        matches = [
            value
            for key, value in data.items()
            if key.startswith("Rate") and key.endswith("." + field)
        ]
        if len(matches) != 1:
            raise ValueError("Unsupported IQ-TREE rate checkpoint: " + field)
        return matches[0]

    frozen = base + ("{" + values + "}" if values else "")
    for name, values in tokens[1:]:
        if not values and name.upper().startswith("G"):
            values = rate_value("gamma_shape")
        elif not values and name.upper().startswith("R"):
            props = rate_value("prop").split(",")
            rates = rate_value("rates").split(",")
            values = ",".join(
                x.strip() for pair in zip(props, rates, strict=True) for x in pair
            )
        elif not values and name.upper() == "I":
            values = rate_value("p_invar")
        frozen += "+" + name + ("{" + values + "}" if values else "")
    return frozen


def split_key(node, names: frozenset[str]):
    side = frozenset(str(n.name) for n in node.leaves())
    other: frozenset[str] = names - side
    first, second = tuple(sorted(side)), tuple(sorted(other))
    return min((first, second), key=lambda value: (len(value), value))


def read_export(prefix, *, single_tip_root=False):
    export = Path(str(prefix) + ".mcmctree.hessian")
    if not export.is_file():
        raise ValueError(
            "IQ-TREE did not produce its required IQ2MC derivative export."
        )
    lines = [line.strip() for line in export.read_text().splitlines() if line.strip()]
    if len(lines) < 6:
        raise ValueError("Incomplete IQ-TREE derivative export.")
    tree = Tree(lines[1], parser=1)
    nodes = [node for node in tree.traverse("preorder") if node is not tree]
    if single_tip_root:
        # IQ2MC puts the first root-tip edge last in its derivative vectors,
        # while leaving that tip first in the accompanying Newick.
        nodes = nodes[1:] + nodes[:1]
    count = 2 * int(lines[0]) - 3
    lengths = np.fromstring(lines[2], sep=" ")
    gradient = np.fromstring(lines[3], sep=" ")
    hessian = np.array([np.fromstring(line, sep=" ") for line in lines[5:]])
    if (
        len(nodes) != count
        or lengths.shape != (count,)
        or gradient.shape != (count,)
        or hessian.shape != (count, count)
        or lines[4] != "Hessian"
    ):
        raise ValueError("Unsupported IQ-TREE Hessian export layout.")
    if not np.allclose(lengths, [n.dist for n in nodes], rtol=1e-5, atol=1e-9):
        raise ValueError("IQ-TREE Hessian branch order does not match its tree.")
    for label, array in (
        ("lengths", lengths),
        ("gradient", gradient),
        ("Hessian", hessian),
    ):
        if not np.isfinite(array).all():
            raise ValueError(
                f"Nonfinite IQ-TREE {label} export at lengths {lengths.tolist()}; tree {lines[1]}."
            )
    data = checkpoint(str(prefix) + ".ckp.gz")
    score = float(data["CandidateSet.0"].split()[0])
    return tree, nodes, lengths, -score, -gradient, -hessian, data


class IQTreeLikelihood:
    def __init__(
        self,
        chronology,
        alignment,
        model,
        *,
        executable="iqtree",
        threads=1,
        seed=1,
        genetic_code=1,
        mode="persistent",
    ):
        if genetic_code != 1:
            raise ValueError("IQ-TREE dating currently supports genetic code 1 only.")
        if threads < 1:
            raise ValueError("IQ-TREE threads must be positive.")
        self.executable = shutil.which(executable)
        if not self.executable:
            raise ValueError("IQ-TREE executable not found: " + executable)
        if mode not in {"persistent", "subprocess"}:
            raise ValueError("Unknown IQ-TREE mode: " + mode)
        self.mode = mode
        self.session = None
        self.chronology = chronology
        self.edges = chronology.edges
        self.edge_id = {node: i for i, node in enumerate(self.edges)}
        self.model = model
        self.alphabet = BASES[model_tokens(model)[0][0].upper()]
        self.names = [str(n.name) for n in chronology.gene.leaves()]
        if len(self.names) < 3:
            raise ValueError(
                "IQ-TREE sequence dating requires at least three gene sequences."
            )
        read_alignment(alignment, self.names, self.alphabet)
        with open(alignment, encoding="utf-8") as handle:
            records = parse_fasta(handle)
        self.sequences = {
            r.name: "".join("".join(r.raw.splitlines()[1:]).split()).upper()
            for r in records
        }
        self.aliases = {name: f"T{i}" for i, name in enumerate(self.names)}
        self.all_aliases = frozenset(self.aliases.values())
        self.threads, self.seed = threads, seed
        self._temporary = tempfile.TemporaryDirectory(prefix="nwkit-iqtree-")
        self.directory = Path(self._temporary.name)
        self.alignment = self.directory / "alignment.fa"
        self.alignment.write_text(
            "".join(
                f">{self.aliases[name]}\n{self.sequences[name]}\n"
                for name in self.names
            )
        )
        self.cache = OrderedDict()
        self.evaluations = 0
        self.version = self._command(["--version"]).splitlines()[0]
        self.initial_lengths = np.array([node.dist for node in self.edges])
        fitted = self._evaluate(self.initial_lengths, model, fixed=False)
        try:
            self.frozen_model = freeze_model(model, fitted[-1])
        except KeyError as exc:
            raise ValueError(
                "IQ-TREE checkpoint lacks required model parameters."
            ) from exc
        self.initial_lengths = fitted[0]
        try:
            if self.mode == "persistent":
                self._start_session()
                checked = self._session_evaluate(self.initial_lengths)
            else:
                checked = self._evaluate(
                    self.initial_lengths, self.frozen_model, fixed=True
                )
        except BaseException:
            self.close()
            raise
        if abs(checked[1] - fitted[1]) > 1e-4:
            self.close()
            raise ValueError(
                "IQ-TREE frozen-model likelihood does not reproduce the fitted model."
            )
        self.prefit_nll = checked[1]
        self.fitted_checkpoint = fitted[-1]

    def _start_session(self):
        from nwkit.radte_iqtree_session import IQTreeSession

        treepath = self.directory / "session.nwk"
        treepath.write_text(self._newick(self.initial_lengths))
        self.session = IQTreeSession(
            [
                self.executable,
                "--likelihood-session",
                "-s",
                str(self.alignment),
                "--seqtype",
                {"dna": "DNA", "protein": "AA", "codon": "CODON"}[self.alphabet],
                "-m",
                self.frozen_model,
                "-te",
                str(treepath),
                "-blfix",
                "--prefix",
                str(self.directory / "session"),
                "-T",
                str(self.threads),
                "--seed",
                str(self.seed),
                "-keep-ident",
                "--safe",
                "-prec",
                "50",
                "-blmin",
                "1e-12",
            ]
        )

        def canonical(side):
            other = self.all_aliases - side
            return min(
                (tuple(sorted(side)), tuple(sorted(other))), key=lambda x: (len(x), x)
            )

        keys = [canonical(side) for side in self.session.splits]
        expected = [
            canonical(frozenset(self.aliases[str(n.name)] for n in e.leaves()))
            for e in self.edges
        ]
        if (
            len(keys) != len(set(keys))
            or set(keys) != set(expected)
            or any(not side <= self.all_aliases for side in self.session.splits)
        ):
            raise ValueError(
                "IQ-TREE session changed the fixed topology or sequence set."
            )
        self.session_mapping = np.array([keys.index(key) for key in expected])
        combined = np.bincount(self.session_mapping, weights=self.initial_lengths)
        if not np.allclose(combined, self.session.lengths, rtol=1e-7, atol=1e-12):
            raise ValueError("IQ-TREE session changed fixed initial branch lengths.")

    def _session_evaluate(self, lengths, *, second_derivatives=True):
        assert self.session is not None
        combined = np.bincount(self.session_mapping, weights=lengths)
        nll, gradient, diagonal = self.session.evaluate(
            combined, second_derivatives=second_derivatives
        )
        self.evaluations += 1
        # Only diagonal curvature is exported; full observed curvature remains
        # finite differences of exact scores in build_quadratic().
        return lengths.copy(), nll, gradient, diagonal, self.session_mapping, {}

    def close(self):
        if self.session is not None:
            self.session.close()

    def _command(self, args):
        result = subprocess.run(
            [self.executable, *args],
            text=True,
            capture_output=True,
            timeout=600,
            check=False,
        )
        if result.returncode:
            raise ValueError(
                "IQ-TREE failed: " + (result.stdout + result.stderr)[-6000:]
            )
        return result.stdout

    def _newick(self, lengths):
        parts: dict[object, str] = {}
        for node in self.chronology.gene.traverse("postorder"):
            value = (
                self.aliases[str(node.name)]
                if node.is_leaf
                else "(" + ",".join(parts[c] for c in node.children) + ")"
            )
            if node in self.edge_id:
                value += ":" + format(lengths[self.edge_id[node]], ".17g")
            parts[node] = value
        return parts[self.chronology.gene] + ";"

    def _evaluate(self, lengths, model, *, fixed):
        with tempfile.TemporaryDirectory(dir=self.directory) as scratch:
            prefix = Path(scratch) / "run"
            treepath = Path(scratch) / "tree.nwk"
            treepath.write_text(self._newick(lengths))
            command = [
                "-s",
                str(self.alignment),
                "--seqtype",
                {"dna": "DNA", "protein": "AA", "codon": "CODON"}[self.alphabet],
                "-m",
                model,
                "-te",
                str(treepath),
                "--dating",
                "mcmctree",
                "--prefix",
                str(prefix),
                "-T",
                str(self.threads),
                "--seed",
                str(self.seed),
                "-keep-ident",
                "--safe",
                # IQ-TREE serializes internal trees with fixed decimal places.
                # Preserve the tiny positive lengths explored by the clock fit.
                "-prec",
                "50",
                "-blmin",
                "1e-12",
            ]
            if fixed:
                command += ["-blfix"]
            self._command(command)
            tree, nodes, exported, nll, gradient, hessian, data = read_export(
                prefix, single_tip_root=self.chronology.gene.children[0].is_leaf
            )
        keys = [split_key(node, self.all_aliases) for node in nodes]
        source = Tree(self._newick(lengths), parser=1)
        source_keys = [
            split_key(node, self.all_aliases)
            for node in source.traverse("preorder")
            if node is not source
        ]
        # Chronology.edges need not be preorder.
        source_by_clade = {
            frozenset(n.name for n in node.leaves()): split_key(node, self.all_aliases)
            for node in source.traverse()
            if node is not source
        }
        mapping = np.array(
            [
                keys.index(
                    source_by_clade[
                        frozenset(self.aliases[str(n.name)] for n in edge.leaves())
                    ]
                )
                for edge in self.edges
            ]
        )
        if (
            len(keys) != len(set(keys))
            or set(source_keys) != set(keys)
            or set(n.name for n in tree.leaves()) != self.all_aliases
        ):
            raise ValueError("IQ-TREE changed the sequence set or fixed topology.")
        combined = np.bincount(mapping, weights=lengths)
        # Newick has more precision except when tiny lengths round to zero.
        exported = np.array(
            [
                node.dist if node.dist > 1e-8 else length
                for node, length in zip(nodes, exported, strict=True)
            ]
        )
        if fixed and not np.allclose(combined, exported, rtol=1e-7, atol=1e-9):
            raise ValueError("IQ-TREE changed fixed branch lengths.")
        fitted_lengths = exported[mapping] * lengths / combined[mapping]
        self.evaluations += 1
        return fitted_lengths, nll, gradient, hessian, mapping, data

    def evaluate(self, lengths, *, second_derivatives=True):
        lengths = np.asarray(lengths, dtype=float)
        if not np.isfinite(lengths).all() or np.any(lengths <= 0):
            raise ValueError(
                "Sequence likelihood requires finite positive branch lengths."
            )
        key = lengths.tobytes()
        if key not in self.cache or (second_derivatives and self.cache[key][3] is None):
            self.cache[key] = (
                self._session_evaluate(lengths, second_derivatives=second_derivatives)
                if self.mode == "persistent"
                else self._evaluate(lengths, self.frozen_model, fixed=True)
            )
            if len(self.cache) > 32:
                self.cache.popitem(last=False)
        return self.cache[key]

    def value_gradient(self, lengths):
        _, nll, gradient, _, mapping, _ = self.evaluate(
            lengths, second_derivatives=False
        )
        return nll, gradient[mapping]

    def bootstrap(self, rng):
        width = 3 if self.alphabet == "codon" else 1
        count = len(self.sequences[self.names[0]]) // width
        columns = rng.integers(0, count, count)
        path = self.directory / "bootstrap.fa"
        path.write_text(
            "".join(
                f">{name}\n"
                + "".join(
                    self.sequences[name][i * width : (i + 1) * width] for i in columns
                )
                + "\n"
                for name in self.names
            )
        )
        return IQTreeLikelihood(
            self.chronology,
            path,
            self.model,
            executable=self.executable,
            threads=self.threads,
            seed=self.seed,
            mode=self.mode,
        )


def prepare_iqtree(c, args):
    from nwkit.radte import _hash_file

    model = requested_model(args)
    exact = IQTreeLikelihood(
        c,
        args.alignment,
        model,
        executable=getattr(args, "iqtree_executable", None) or "iqtree",
        threads=args.iqtree_threads
        if getattr(args, "iqtree_threads", None) is not None
        else 1,
        seed=args.seed,
        mode=getattr(args, "iqtree_mode", None) or "persistent",
        genetic_code=args.genetic_code if args.genetic_code is not None else 1,
    )
    return exact, dict(
        engine="iqtree",
        model=model,
        frozen_model=exact.frozen_model,
        iqtree_version=exact.version,
        iqtree_mode=exact.mode,
        iqtree_session_protocol=1 if exact.mode == "persistent" else None,
        iqtree_binary_sha256=_hash_file(exact.executable),
        derivative_method="iqtree-score; finite-difference-log-length-hessian",
        alignment_sites=len(exact.sequences[exact.names[0]])
        // (3 if exact.alphabet == "codon" else 1),
        genetic_code=args.genetic_code if args.genetic_code is not None else 1,
        branch_length_unit="iqtree_substitutions_per_site",
        prefit=dict(
            status="estimated-unclocked-conditional-model", objective=exact.prefit_nll
        ),
        fitted_parameters={
            key: value
            for key, value in exact.fitted_checkpoint.items()
            if key.startswith(("Model", "Rate"))
        },
    )
