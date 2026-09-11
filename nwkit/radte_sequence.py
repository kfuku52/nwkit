"""Fixed-topology reversible sequence likelihood and analytic branch gradients.

The rooted pruning representation automatically treats the two root edges as
one identifiable unrooted length. No independent root-split observation enters
the sequence likelihood. Site patterns and equal-probability gamma categories
are integrated exactly (conditional on the chosen substitution parameters).
"""

from dataclasses import dataclass
from functools import cached_property
from importlib.resources import files
from typing import Any

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import gamma

from nwkit.fasta import parse_fasta
from nwkit.radte_codon import CODON_MODELS, CODONS, codon_matrix, encode_codons

DNA_STATES = "ACGT"
DNA_CODES = dict(zip(DNA_STATES, DNA_STATES, strict=True))
DNA_CODES.update(
    R="AG",
    Y="CT",
    S="CG",
    W="AT",
    K="GT",
    M="AC",
    B="CGT",
    D="AGT",
    H="ACT",
    V="ACG",
    N="ACGT",
    X="ACGT",
    **{"?": "ACGT", "-": "ACGT", "U": "T"},
)
AA_STATES = "ARNDCQEGHILKMFPSTWYV"


def read_alignment(path, names, alphabet="dna"):
    with open(path, encoding="utf-8") as handle:
        records = parse_fasta(handle)
    seqs = {}
    for record in records:
        if record.name in seqs:
            raise ValueError(f"Duplicate alignment sequence: {record.name}")
        seqs[record.name] = (
            "".join(record.raw.splitlines()[1:]).replace(" ", "").upper()
        )
    if set(seqs) != set(names):
        raise ValueError("Alignment names must exactly match the gene-tree tips.")
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1 or not next(iter(lengths)):
        raise ValueError("Alignment must have a positive, equal sequence length.")
    if alphabet == "codon":
        return encode_codons(seqs, names, DNA_CODES)
    states = DNA_STATES if alphabet == "dna" else AA_STATES
    codes = (
        DNA_CODES
        if alphabet == "dna"
        else {
            **dict(zip(states, states, strict=True)),
            "B": "DN",
            "Z": "EQ",
            "J": "IL",
            "X": states,
            "?": states,
            "-": states,
        }
    )
    matrix = np.empty((len(names), next(iter(lengths))), dtype=np.uint32)
    state_id = {s: i for i, s in enumerate(states)}
    for i, name in enumerate(names):
        for j, code in enumerate(seqs[name]):
            if code not in codes:
                raise ValueError(f"Unsupported {alphabet} alignment character: {code}")
            matrix[i, j] = sum(1 << state_id[s] for s in codes[code])
    return matrix, states


def gamma_rates(shape, categories):
    if categories < 1 or not np.isfinite(shape) or shape <= 0:
        raise ValueError("Gamma shape and category count must be positive.")
    if categories == 1:
        return np.ones(1)
    cutpoints = gamma.ppf(
        np.arange(categories + 1) / categories, a=shape, scale=1 / shape
    )
    moments = gamma.cdf(cutpoints, a=shape + 1, scale=1 / shape)
    rates = categories * np.diff(moments)
    if np.any(rates <= 0):
        raise ValueError("Gamma categories underflowed; increase gamma shape.")
    return rates / rates.mean()


def substitution_matrix(
    model, patterns, weights, states, kappa=2.0, exchangeabilities=None
):
    count = len(states)
    frequencies = np.full(count, 0.5)
    for i in range(count):
        frequencies[i] += np.sum(weights * np.sum(patterns == (1 << i), axis=0))
    frequencies /= frequencies.sum()
    if model in {"jc69", "poisson"}:
        frequencies[:] = 1 / count
    exchange = np.ones((count, count))
    if model == "hky":
        if not np.isfinite(kappa) or kappa <= 0:
            raise ValueError("HKY kappa must be positive.")
        for a, b in [(0, 2), (1, 3)]:
            exchange[a, b] = exchange[b, a] = kappa
    if model == "gtr":
        exchangeabilities = (
            np.ones(6) if exchangeabilities is None else np.asarray(exchangeabilities)
        )
        if (
            exchangeabilities.shape != (6,)
            or not np.isfinite(exchangeabilities).all()
            or np.any(exchangeabilities <= 0)
        ):
            raise ValueError(
                "GTR requires six positive exchangeabilities in AC, AG, AT, CG, CT, GT order."
            )
        exchange[np.triu_indices(4, 1)] = exchangeabilities
        exchange[np.tril_indices(4, -1)] = exchange.T[np.tril_indices(4, -1)]
    if model in {"lg", "lg-f"}:
        text = files("nwkit").joinpath("data_model/lg.txt").read_text(encoding="utf-8")
        values = np.fromstring(
            " ".join(line for line in text.splitlines() if not line.startswith("#")),
            sep=" ",
        )
        if len(values) != 210:
            raise ValueError("Invalid packaged LG matrix.")
        exchange[np.tril_indices(20, -1)] = values[:190]
        exchange[np.triu_indices(20, 1)] = exchange.T[np.triu_indices(20, 1)]
        if model == "lg":
            frequencies = values[190:] / values[190:].sum()
    matrix = exchange * frequencies[None, :]
    np.fill_diagonal(matrix, 0)
    np.fill_diagonal(matrix, -matrix.sum(axis=1))
    matrix /= -(frequencies @ np.diag(matrix))
    return matrix, frequencies


class SequenceLikelihood:
    def __init__(
        self,
        chronology,
        alignment,
        *,
        model="hky",
        kappa=2.0,
        gamma_shape=1.0,
        gamma_categories=4,
        matrix=None,
        exchangeabilities=None,
        omega=0.5,
        codon_frequencies=None,
        genetic_code=1,
    ):
        self.chronology = chronology
        self.fit_settings: dict[str, bool] = {}
        self.model = model
        if (
            model
            not in {"jc69", "hky", "gtr", "f81", "poisson", "lg", "lg-f"} | CODON_MODELS
        ):
            raise ValueError("Unsupported sequence substitution model.")
        self.omega = omega
        self.genetic_code = genetic_code
        self.codon_frequencies = codon_frequencies or (
            "f3x4" if model == "gy94" else "model"
        )
        if model in CODON_MODELS and genetic_code != 1:
            raise ValueError("Codon models currently require standard genetic code 1.")
        self.names = [str(n.name) for n in chronology.gene.leaves()]
        alphabet = "protein" if model in {"poisson", "lg", "lg-f"} else "dna"
        if model in CODON_MODELS:
            alphabet = "codon"
        if matrix is None:
            matrix, self.states = read_alignment(alignment, self.names, alphabet)
        else:
            self.states = (
                CODONS
                if alphabet == "codon"
                else (AA_STATES if alphabet == "protein" else DNA_STATES)
            )
        self.raw_matrix = matrix
        self.patterns, self.counts = np.unique(matrix, axis=1, return_counts=True)
        self.rates = gamma_rates(gamma_shape, gamma_categories)
        self.kappa, self.gamma_shape, self.gamma_categories = (
            kappa,
            gamma_shape,
            gamma_categories,
        )
        self.exchangeabilities = (
            np.ones(6)
            if exchangeabilities is None
            else np.asarray(exchangeabilities, dtype=float).copy()
        )
        self.update_parameters()
        self.edges = chronology.edges
        self.edge_id = {n: i for i, n in enumerate(self.edges)}
        self.tip_data = {}
        for i, node in enumerate(chronology.gene.leaves()):
            self.tip_data[node] = (
                (
                    self.patterns[i, :, None]
                    >> np.arange(len(self.states), dtype=self.patterns.dtype)
                )
                & 1
            ).astype(float)
        self.postorder = list(chronology.gene.traverse(strategy="postorder"))
        self.preorder = list(chronology.gene.traverse())

    def update_parameters(self):
        self.rates = gamma_rates(self.gamma_shape, self.gamma_categories)
        if self.model in CODON_MODELS:
            self.q, self.pi = codon_matrix(
                self.model,
                self.patterns,
                self.counts,
                self.kappa,
                self.omega,
                self.codon_frequencies,
            )
        else:
            self.q, self.pi = substitution_matrix(
                self.model,
                self.patterns,
                self.counts,
                self.states,
                self.kappa,
                self.exchangeabilities,
            )
        rootpi = np.sqrt(self.pi)
        symmetric = rootpi[:, None] * self.q / rootpi[None, :]
        self.eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
        self.eigenvalues = np.minimum(self.eigenvalues, 0.0)
        # Every irreducible rate generator has an exact stationary eigenvalue.
        # A tiny negative eigensolver residual otherwise destroys probability
        # mass on long branches, including valid optimizer trial parameters.
        self.eigenvalues[-1] = 0.0
        self.left = eigenvectors / rootpi[:, None]
        self.right = eigenvectors.T * rootpi[None, :]

    def transition(self, length, rate):
        powers = np.exp(self.eigenvalues * length * rate)
        # expm1 preserves off-diagonal probabilities on very short branches.
        p = (
            np.eye(len(self.pi))
            + (self.left * np.expm1(self.eigenvalues * length * rate)) @ self.right
        )
        dp = (self.left * (powers * self.eigenvalues * rate)) @ self.right
        return np.maximum(p, 0), dp

    def _category(self, lengths, rate):
        partials, messages, transition, derivatives, scales = {}, {}, {}, {}, {}
        for node in self.postorder:
            if node.is_leaf:
                partials[node] = self.tip_data[node]
                scales[node] = np.zeros(len(self.counts))
            else:
                value = np.ones_like(next(iter(self.tip_data.values())))
                scale = np.zeros(len(self.counts))
                for child in node.children:
                    p, dp = self.transition(lengths[self.edge_id[child]], rate)
                    transition[child], derivatives[child] = p, dp
                    messages[child] = partials[child] @ p.T
                    value *= messages[child]
                    scale += scales[child]
                norm = value.max(axis=1)
                if np.any(norm <= 0):
                    raise ValueError(
                        "Sequence likelihood underflow: incompatible zero-length branches."
                    )
                partials[node] = value / norm[:, None]
                scales[node] = scale + np.log(norm)
        root = self.chronology.gene
        log_likelihood = np.log(partials[root] @ self.pi) + scales[root]
        outside = {root: np.broadcast_to(self.pi, partials[root].shape)}
        gradient = np.zeros((len(self.edges), len(self.counts)))
        for node in self.preorder:
            if node.is_leaf:
                continue
            first, second = node.children
            for child, sibling in [(first, second), (second, first)]:
                context = outside[node] * messages[sibling]
                denominator = np.sum(context * messages[child], axis=1)
                numerator = np.sum(
                    context * (partials[child] @ derivatives[child].T), axis=1
                )
                gradient[self.edge_id[child]] = numerator / denominator
                value = context @ transition[child]
                outside[child] = value / value.max(axis=1)[:, None]
        return log_likelihood, gradient

    def value_gradient(self, lengths):
        if not np.isfinite(lengths).all() or np.any(lengths <= 0):
            raise ValueError(
                "Sequence likelihood requires finite positive branch lengths."
            )
        category_values, gradients = zip(
            *(self._category(lengths, r) for r in self.rates), strict=True
        )
        category_values = np.array(category_values)
        total = logsumexp(category_values, axis=0)
        posterior = np.exp(category_values - total)
        gradient = np.einsum(
            "kp,kep,p->e", posterior, np.asarray(gradients), self.counts
        )
        return -float((total - np.log(len(self.rates))) @ self.counts), -gradient

    def bootstrap(self, rng):
        columns = rng.integers(0, self.raw_matrix.shape[1], self.raw_matrix.shape[1])
        result = SequenceLikelihood(
            self.chronology,
            None,
            model=self.model,
            omega=self.omega,
            codon_frequencies=self.codon_frequencies,
            genetic_code=self.genetic_code,
            kappa=self.kappa,
            gamma_shape=self.gamma_shape,
            gamma_categories=self.gamma_categories,
            matrix=self.raw_matrix[:, columns],
            exchangeabilities=self.exchangeabilities,
        )
        result.fit_settings = getattr(self, "fit_settings", {})
        return result


@dataclass
class QuadraticLikelihood:
    """Taylor likelihood in identifiable log unrooted lengths.

    The root pair is added before taking logs. This is essential: a rooted
    branch Hessian has an unidentifiable root-split direction.
    """

    center: np.ndarray
    gradient: np.ndarray
    hessian: np.ndarray
    nll: float
    mapping: np.ndarray
    exact: Any = None

    def value_gradient(self, lengths):
        combined = np.bincount(
            self.mapping, weights=lengths, minlength=len(self.center)
        )
        difference = np.log(combined) - self.center
        gradient = self.gradient + self.hessian @ difference
        value = (
            self.nll
            + self.gradient @ difference
            + 0.5 * difference @ self.hessian @ difference
        )
        return float(value), (gradient / combined)[self.mapping]

    @cached_property
    def covariance(self):
        return np.linalg.inv(self.hessian)

    def check(self, lengths, tolerance=0.1):
        approx_value, approx_gradient = self.value_gradient(lengths)
        exact_value, exact_gradient = self.exact.value_gradient(lengths)
        error = abs(approx_value - exact_value)
        score = np.bincount(
            self.mapping, weights=(approx_gradient - exact_gradient) * lengths
        )
        # Score error in one-standard-error likelihood coordinates. Raw score
        # differences depend on alignment length and do not measure local fit.
        gradient_error = float(np.sqrt(max(0, score @ self.covariance @ score)))
        return (
            error <= tolerance and gradient_error <= tolerance,
            float(error),
            gradient_error,
        )


def build_quadratic(likelihood, initial_lengths, maxiter=1000):
    count = len(initial_lengths)
    if count > 600:
        raise ValueError(
            "Quadratic likelihood exceeds the 600-edge dense Hessian limit."
        )
    root_children = likelihood.chronology.gene.children
    root_ids = [likelihood.edge_id[n] for n in root_children]
    mapping = []
    next_id = 1
    for i in range(count):
        if i in root_ids:
            mapping.append(0)
        else:
            mapping.append(next_id)
            next_id += 1
    mapping = np.array(mapping)
    combined = np.bincount(mapping, weights=initial_lengths)
    fractions = np.ones(count)
    fractions[root_ids] = initial_lengths[root_ids] / combined[0]

    def objective(x):
        lengths = np.exp(x[mapping]) * fractions
        nll, gradient = likelihood.value_gradient(lengths)
        return nll, np.bincount(mapping, weights=gradient * lengths)

    result = minimize(
        objective,
        np.log(combined),
        jac=True,
        method="L-BFGS-B",
        bounds=[(-18.0, 4.0)] * len(combined),
        options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-5},
    )
    if not result.success or np.max(np.abs(result.jac)) > 0.01:
        raise ValueError(
            "Unclocked branch likelihood did not converge to an interior optimum."
        )
    if np.any(result.x < -17.9) or np.any(result.x > 3.9):
        raise ValueError(
            "Branch likelihood lies on a numerical bound; use exact likelihood."
        )
    hessian = np.empty((len(combined), len(combined)))
    for j in range(len(combined)):
        step = np.zeros(len(combined))
        step[j] = 1e-4
        hessian[:, j] = (
            objective(result.x + step)[1] - objective(result.x - step)[1]
        ) / 2e-4
    hessian = (hessian + hessian.T) / 2
    eigen = np.linalg.eigvalsh(hessian)
    if eigen[0] <= max(1e-8, eigen[-1] * 1e-10):
        raise ValueError(
            "Branch likelihood Hessian is singular/nonpositive; use exact likelihood."
        )
    return QuadraticLikelihood(
        result.x, result.jac, hessian, float(result.fun), mapping, likelihood
    )
