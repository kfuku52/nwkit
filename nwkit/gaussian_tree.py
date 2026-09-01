"""Shared linear-Gaussian processes on rooted trees.

The process layer deliberately separates biological model semantics from the
consumer.  ASR can retain a flat or stationary root and condition every node,
whereas regression can use a fixed root, profile fixed effects, and request a
normalized tip covariance.  Both consumers nevertheless use the same branch
transitions and innovation variances.
"""

import math
from dataclasses import dataclass, replace
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import sparse

from nwkit.clade_index import LcaIndex
from nwkit.sparse_laplace import SparseCovarianceModel


@dataclass(frozen=True, slots=True)
class GaussianTransition:
    """One branch transition ``child = slope * parent + intercept + error``."""

    slope: float
    intercept: float
    variance: float

    def __post_init__(self) -> None:
        values = (float(self.slope), float(self.intercept), float(self.variance))
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Gaussian branch transitions must be finite.")
        if values[2] < 0.0:
            raise ValueError("Gaussian branch variances must be non-negative.")
        object.__setattr__(self, "slope", values[0])
        object.__setattr__(self, "intercept", values[1])
        object.__setattr__(self, "variance", values[2])

    def scaled_variance(self, scale: float) -> "GaussianTransition":
        """Return the same affine transition with its innovation variance scaled."""

        return replace(self, variance=self.variance * scale)


@dataclass(frozen=True, slots=True)
class GaussianRootPrior:
    """Explicit root treatment for a tree process.

    ``fixed`` has zero variance, ``flat`` is improper, and ``gaussian`` or
    ``stationary`` have a proper finite variance.  ``stationary`` is retained
    as a distinct label so model metadata cannot silently lose that semantic.
    """

    mode: str
    mean: float = 0.0
    variance: float | None = 0.0

    def __post_init__(self) -> None:
        if self.mode not in {"fixed", "flat", "gaussian", "stationary"}:
            raise ValueError("Unsupported Gaussian root-prior mode: {}.".format(self.mode))
        mean = float(self.mean)
        if not math.isfinite(mean):
            raise ValueError("Gaussian root-prior means must be finite.")
        object.__setattr__(self, "mean", mean)
        if self.mode == "flat":
            if self.variance is not None:
                raise ValueError("A flat Gaussian root prior must not have a variance.")
            return
        if self.variance is None:
            raise ValueError("A proper Gaussian root prior requires a variance.")
        variance = float(self.variance)
        if not math.isfinite(variance) or variance < 0.0:
            raise ValueError("Gaussian root-prior variances must be non-negative and finite.")
        if self.mode == "fixed" and variance != 0.0:
            raise ValueError("A fixed Gaussian root prior must have zero variance.")
        if self.mode in {"gaussian", "stationary"} and variance <= 0.0:
            raise ValueError("A proper Gaussian root prior must have positive variance.")
        object.__setattr__(self, "variance", variance)

    @property
    def is_proper(self) -> bool:
        return self.mode != "flat"

    def scaled_variance(self, scale: float) -> "GaussianRootPrior":
        if self.mode in {"flat", "fixed"}:
            return self
        assert self.variance is not None
        return replace(self, variance=self.variance * scale)


@dataclass(frozen=True)
class GaussianTreeProcess:
    """A rooted scalar linear-Gaussian process in original tree-node space."""

    tree: Any
    transitions: Mapping[Any, GaussianTransition]
    root: GaussianRootPrior
    model: str
    parameter: float | None = None

    def __post_init__(self) -> None:
        expected = {node for node in self.tree.traverse() if not node.is_root}
        transitions = dict(self.transitions)
        missing = expected - set(transitions)
        extra = set(transitions) - expected
        if missing or extra:
            raise ValueError("Gaussian tree transitions must cover every non-root node exactly.")
        object.__setattr__(self, "transitions", transitions)

    def scaled_variance(self, scale: float) -> "GaussianTreeProcess":
        """Scale every stochastic variance while preserving affine means."""

        scale = float(scale)
        if not math.isfinite(scale) or scale < 0.0:
            raise ValueError("Gaussian process variance scales must be non-negative and finite.")
        if scale == 0.0 and self.root.mode in {"gaussian", "stationary"}:
            raise ValueError("A proper Gaussian root prior cannot be scaled to zero variance.")
        return replace(
            self,
            transitions={
                node: transition.scaled_variance(scale)
                for node, transition in self.transitions.items()
            },
            root=self.root.scaled_variance(scale),
        )

    def marginal_moments(self) -> tuple[dict[Any, float], dict[Any, float]]:
        """Return unconditional means and variances for a proper-root process."""

        if not self.root.is_proper:
            raise ValueError("A flat-root process has no finite unconditional moments.")
        assert self.root.variance is not None
        means = {self.tree: self.root.mean}
        variances = {self.tree: self.root.variance}
        for node in self.tree.traverse(strategy="preorder"):
            if node.is_root:
                continue
            transition = self.transitions[node]
            means[node] = transition.slope * means[node.up] + transition.intercept
            variances[node] = (
                transition.slope * transition.slope * variances[node.up]
                + transition.variance
            )
            if not math.isfinite(means[node]) or not math.isfinite(variances[node]):
                raise ValueError("Gaussian tree moments exceed floating-point range.")
        return means, variances

    def covariance(self, nodes: Sequence[Any]) -> np.ndarray:
        """Return a covariance matrix for arbitrary nodes under a proper root."""

        nodes = tuple(nodes)
        _, variances = self.marginal_moments()
        known = set(variances)
        if any(node not in known for node in nodes):
            raise ValueError("Gaussian covariance requested a node outside the process tree.")
        ancestor_loadings: dict[Any, dict[Any, float]] = {}
        for node in nodes:
            current = node
            loading = 1.0
            by_ancestor = {current: loading}
            while not current.is_root:
                loading *= self.transitions[current].slope
                current = current.up
                by_ancestor[current] = loading
            ancestor_loadings[node] = by_ancestor
        lca = LcaIndex(self.tree)
        covariance = np.empty((len(nodes), len(nodes)), dtype=float)
        for first, first_node in enumerate(nodes):
            for second in range(first, len(nodes)):
                second_node = nodes[second]
                ancestor = lca.common_ancestor(first_node, second_node)
                value = (
                    ancestor_loadings[first_node][ancestor]
                    * ancestor_loadings[second_node][ancestor]
                    * variances[ancestor]
                )
                covariance[first, second] = value
                covariance[second, first] = value
        if not np.isfinite(covariance).all():
            raise ValueError("Gaussian tree covariance exceeds floating-point range.")
        return covariance

    def ordered_leaves(self, leaf_names: Sequence[str]) -> tuple[Any, ...]:
        names = tuple(str(name) for name in leaf_names)
        if len(set(names)) != len(names):
            raise ValueError("Gaussian covariance requested duplicated tree-tip names.")
        leaf_by_name = {str(leaf.name): leaf for leaf in self.tree.leaves()}
        missing = sorted(set(names) - set(leaf_by_name))
        if missing:
            raise ValueError(
                "Gaussian covariance requested absent tree tips: {}.".format(
                    ", ".join(missing)
                )
            )
        return tuple(leaf_by_name[name] for name in names)

    def tip_covariance(self, leaf_names: Sequence[str]) -> np.ndarray:
        """Return an unnormalized covariance in the requested tip order."""

        return self.covariance(self.ordered_leaves(leaf_names))

    def sparse_tip_model(
        self,
        leaf_names: Sequence[str],
        *,
        normalize: bool = True,
    ) -> SparseCovarianceModel:
        """Return a sparse latent representation of the requested tip covariance."""

        if not self.root.is_proper:
            raise ValueError("A flat-root process has no proper sparse tip covariance.")
        leaves = self.ordered_leaves(leaf_names)
        required: set[Any] = set()
        for leaf in leaves:
            node = leaf
            while node is not None and node not in required:
                required.add(node)
                node = node.up

        assert self.root.variance is not None
        state_by_node: dict[Any, int] = {self.tree: -1}
        state_scale_by_node: dict[Any, float] = {self.tree: 0.0}
        marginal_variance = {self.tree: self.root.variance}
        state_parent: list[int] = []
        state_transition: list[float] = []
        state_innovation: list[float] = []
        if self.root.variance > 0.0:
            state_by_node[self.tree] = 0
            state_scale_by_node[self.tree] = 1.0
            state_parent.append(-1)
            state_transition.append(0.0)
            state_innovation.append(self.root.variance)

        for node in self.tree.traverse(strategy="preorder"):
            if node.is_root or node not in required:
                continue
            transition = self.transitions[node]
            parent_state = state_by_node[node.up]
            parent_scale = state_scale_by_node[node.up]
            marginal_variance[node] = (
                transition.slope * transition.slope * marginal_variance[node.up]
                + transition.variance
            )
            if transition.variance == 0.0:
                state_by_node[node] = parent_state
                state_scale_by_node[node] = transition.slope * parent_scale
                continue
            state_by_node[node] = len(state_parent)
            state_scale_by_node[node] = 1.0
            state_parent.append(parent_state)
            state_transition.append(
                transition.slope * parent_scale if parent_state >= 0 else 0.0
            )
            state_innovation.append(transition.variance)

        tip_variances = np.asarray([marginal_variance[leaf] for leaf in leaves])
        normalization = float(np.mean(tip_variances)) if normalize else 1.0
        if not math.isfinite(normalization) or normalization <= 0.0:
            raise ValueError("Sparse Gaussian tree covariance has zero tip variance.")
        normalized_innovation = np.asarray(state_innovation, dtype=float) / normalization
        if not len(normalized_innovation) or np.any(normalized_innovation <= 0.0):
            raise ValueError("Sparse Gaussian innovations must be positive.")

        rows: list[int] = []
        columns: list[int] = []
        values: list[float] = []

        def add(row: int, column: int, value: float) -> None:
            rows.append(row)
            columns.append(column)
            values.append(value)

        for child, (parent, coefficient, innovation) in enumerate(
            zip(
                state_parent,
                state_transition,
                normalized_innovation,
                strict=True,
            )
        ):
            inverse = 1.0 / innovation
            add(child, child, inverse)
            if parent >= 0:
                add(parent, parent, coefficient * coefficient * inverse)
                add(child, parent, -coefficient * inverse)
                add(parent, child, -coefficient * inverse)
        n_states = len(state_parent)
        precision = sparse.coo_matrix(
            (values, (rows, columns)), shape=(n_states, n_states)
        ).tocsc()
        tip_states = np.asarray([state_by_node[leaf] for leaf in leaves], dtype=int)
        if np.any(tip_states < 0):
            raise ValueError("Every requested tip must have positive evolutionary variance.")
        tip_scales = np.asarray([state_scale_by_node[leaf] for leaf in leaves], dtype=float)
        tip_loading = sparse.csr_matrix(
            (
                tip_scales,
                (np.arange(len(leaves), dtype=int), tip_states),
            ),
            shape=(len(leaves), n_states),
        )
        return SparseCovarianceModel(
            precision=precision,
            tip_loading=tip_loading,
            logdet_covariance=float(np.sum(np.log(normalized_innovation))),
            sampling_parent=np.asarray(state_parent, dtype=int),
            sampling_transition=np.asarray(state_transition, dtype=float),
            sampling_variance=normalized_innovation,
            covariance_scale=normalization,
        )


def brownian_transition(variance: float) -> GaussianTransition:
    """Return a unit-slope Brownian branch transition."""

    return GaussianTransition(1.0, 0.0, variance)


def ou_transition(
    length: float,
    alpha: float,
    variance_rate: float = 1.0,
    optimum: float = 0.0,
) -> tuple[GaussianTransition, float]:
    """Return an OU branch transition and its stationary root variance."""

    alpha = float(alpha)
    variance_rate = float(variance_rate)
    if not all(math.isfinite(value) for value in (alpha, variance_rate)):
        raise ValueError("OU transition parameters must be finite.")
    if alpha <= 0.0 or variance_rate < 0.0:
        raise ValueError("OU requires non-negative variance and positive alpha.")
    root_variance = variance_rate / (2.0 * alpha)
    if not math.isfinite(root_variance):
        raise ValueError("The stationary OU variance exceeds floating-point range.")
    return (
        ou_transition_from_root_variance(length, alpha, root_variance, optimum),
        root_variance,
    )


def ou_transition_from_root_variance(
    length: float,
    alpha: float,
    root_variance: float,
    optimum: float = 0.0,
) -> GaussianTransition:
    """Return an OU transition parameterized by stationary state variance."""

    length = float(length)
    alpha = float(alpha)
    root_variance = float(root_variance)
    optimum = float(optimum)
    if not all(math.isfinite(value) for value in (length, alpha, root_variance, optimum)):
        raise ValueError("OU transition parameters must be finite.")
    if length < 0.0 or alpha <= 0.0 or root_variance < 0.0:
        raise ValueError(
            "OU requires non-negative lengths and root variance, and positive alpha."
        )
    exponent = alpha * length
    if not math.isfinite(exponent):
        raise ValueError("OU alpha multiplied by a branch length exceeds floating-point range.")
    decay = math.exp(-exponent)
    innovation = root_variance * (-math.expm1(-2.0 * exponent))
    if root_variance > 0.0 and length > 0.0 and innovation <= 0.0:
        raise ValueError("A positive OU branch variance underflowed.")
    intercept = (-math.expm1(-exponent)) * optimum
    return GaussianTransition(decay, intercept, innovation)


def exponential_rate_edge_variance(start: float, length: float, rate: float) -> float:
    """Integrate ``exp(rate * time)`` over one branch without cancellation."""

    start = float(start)
    length = float(length)
    rate = float(rate)
    if not all(math.isfinite(value) for value in (start, length, rate)):
        raise ValueError("Rate-change branch parameters must be finite.")
    if start < 0.0 or length < 0.0:
        raise ValueError("Rate-change models require non-negative times.")
    if length == 0.0 or rate == 0.0:
        return length
    start_exponent = rate * start
    edge_exponent = rate * length
    if not math.isfinite(start_exponent) or not math.isfinite(edge_exponent):
        raise ValueError("Rate change multiplied by tree depth exceeds floating-point range.")
    try:
        variance = math.exp(start_exponent) * math.expm1(edge_exponent) / rate
    except OverflowError as exc:
        raise ValueError("A rate-change branch variance exceeds floating-point range.") from exc
    if not math.isfinite(variance) or variance <= 0.0:
        raise ValueError("A positive rate-change branch variance is not finite.")
    return variance
