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
            raise ValueError(
                "Unsupported Gaussian root-prior mode: {}.".format(self.mode)
            )
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
            raise ValueError(
                "Gaussian root-prior variances must be non-negative and finite."
            )
        if self.mode == "fixed" and variance != 0.0:
            raise ValueError("A fixed Gaussian root prior must have zero variance.")
        if self.mode in {"gaussian", "stationary"} and variance <= 0.0:
            raise ValueError(
                "A proper Gaussian root prior must have positive variance."
            )
        object.__setattr__(self, "variance", variance)

    @property
    def is_proper(self) -> bool:
        return self.mode != "flat"

    def scaled_variance(self, scale: float) -> "GaussianRootPrior":
        if self.mode in {"flat", "fixed"}:
            return self
        assert self.variance is not None
        return replace(self, variance=self.variance * scale)


@dataclass
class _SparseTreeState:
    state_by_node: dict[Any, int]
    state_scale_by_node: dict[Any, float]
    marginal_variance: dict[Any, float]
    parents: list[int]
    coefficients: list[float]
    innovations: list[float]


def _required_ancestry(tree, leaves, all_leaves_requested):
    if all_leaves_requested:
        return None
    required = set()
    for leaf in leaves:
        node = leaf
        while node is not None and node not in required:
            required.add(node)
            node = node.up
    return required


def _initial_sparse_state(process, required):
    assert process.root.variance is not None
    state = _SparseTreeState(
        state_by_node={process.tree: -1},
        state_scale_by_node={process.tree: 0.0},
        marginal_variance={process.tree: process.root.variance},
        parents=[],
        coefficients=[],
        innovations=[],
    )
    if process.root.variance > 0.0:
        state.state_by_node[process.tree] = 0
        state.state_scale_by_node[process.tree] = 1.0
        state.parents.append(-1)
        state.coefficients.append(0.0)
        state.innovations.append(process.root.variance)
    for node in process.tree.traverse(strategy="preorder"):
        if node.is_root or (required is not None and node not in required):
            continue
        transition = process.transitions[node]
        parent_state = state.state_by_node[node.up]
        parent_scale = state.state_scale_by_node[node.up]
        state.marginal_variance[node] = (
            transition.slope * transition.slope * state.marginal_variance[node.up]
            + transition.variance
        )
        if transition.variance == 0.0:
            state.state_by_node[node] = parent_state
            state.state_scale_by_node[node] = transition.slope * parent_scale
            continue
        state.state_by_node[node] = len(state.parents)
        state.state_scale_by_node[node] = 1.0
        state.parents.append(parent_state)
        state.coefficients.append(
            transition.slope * parent_scale if parent_state >= 0 else 0.0
        )
        state.innovations.append(transition.variance)
    return state


def _tip_variance_normalization(state, leaves, normalize):
    tip_variances = np.asarray(
        [state.marginal_variance[leaf] for leaf in leaves], dtype=float
    )
    if (
        not len(tip_variances)
        or np.any(~np.isfinite(tip_variances))
        or np.any(tip_variances < 0.0)
    ):
        raise ValueError("Sparse Gaussian tree covariance has invalid tip variance.")
    maximum = float(np.max(tip_variances))
    if maximum <= 0.0:
        raise ValueError("Sparse Gaussian tree covariance has zero tip variance.")
    normalization = (
        maximum * float(np.mean(tip_variances / maximum)) if normalize else 1.0
    )
    if not math.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("Sparse Gaussian tree covariance has zero tip variance.")
    return normalization


def _build_sparse_precision(parents, coefficients, innovations):
    parent_values = np.asarray(parents, dtype=int)
    coefficient_values = np.asarray(coefficients, dtype=float)
    innovation_values = np.asarray(innovations, dtype=float)
    size = len(parent_values)
    if (
        len(coefficient_values) != size
        or len(innovation_values) != size
        or np.any(~np.isfinite(coefficient_values))
        or np.any(~np.isfinite(innovation_values))
        or np.any(innovation_values < 1.0 / np.finfo(float).max)
    ):
        return None
    inverse = 1.0 / innovation_values
    children = np.arange(size, dtype=int)
    has_parent = parent_values >= 0
    child_with_parent = children[has_parent]
    parents_with_child = parent_values[has_parent]
    coefficients_with_parent = coefficient_values[has_parent]
    inverse_with_parent = inverse[has_parent]
    parent_diagonal = (
        coefficients_with_parent * coefficients_with_parent * inverse_with_parent
    )
    cross = -coefficients_with_parent * inverse_with_parent
    if any(
        np.any(~np.isfinite(values)) for values in (inverse, parent_diagonal, cross)
    ):
        return None
    rows = np.concatenate(
        (children, parents_with_child, child_with_parent, parents_with_child)
    )
    columns = np.concatenate(
        (children, parents_with_child, parents_with_child, child_with_parent)
    )
    values = np.concatenate((inverse, parent_diagonal, cross, cross))
    result = sparse.coo_matrix((values, (rows, columns)), shape=(size, size)).tocsc()
    return result if np.isfinite(result.data).all() else None


def _normal_sparse_precision(state, normalization):
    with np.errstate(over="ignore", under="ignore", divide="ignore"):
        innovations = np.asarray(state.innovations, dtype=float) / normalization
    representable = (
        len(innovations) > 0
        and np.all(np.isfinite(innovations))
        and np.all(innovations > 0.0)
        and np.all(np.isfinite(state.coefficients))
        and np.all(np.isfinite(list(state.state_scale_by_node.values())))
    )
    precision = (
        _build_sparse_precision(state.parents, state.coefficients, innovations)
        if representable
        else None
    )
    return precision, innovations


def _normalized_standard_deviation(variance, normalization):
    log_sd = 0.5 * (math.log(float(variance)) - math.log(normalization))
    try:
        value = math.exp(log_sd)
    except OverflowError as exc:
        raise ValueError(
            "Sparse Gaussian covariance exceeds floating-point dynamic range."
        ) from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(
            "Sparse Gaussian covariance exceeds floating-point dynamic range."
        )
    return value


def _rescaled_sparse_state(process, required, normalization):
    assert process.root.variance is not None
    state = _SparseTreeState(
        state_by_node={process.tree: -1},
        state_scale_by_node={process.tree: 0.0},
        marginal_variance={},
        parents=[],
        coefficients=[],
        innovations=[],
    )
    if process.root.variance > 0.0:
        state.state_by_node[process.tree] = 0
        state.state_scale_by_node[process.tree] = _normalized_standard_deviation(
            process.root.variance, normalization
        )
        state.parents.append(-1)
        state.coefficients.append(0.0)
        state.innovations.append(1.0)
    for node in process.tree.traverse(strategy="preorder"):
        if node.is_root or (required is not None and node not in required):
            continue
        transition = process.transitions[node]
        parent_state = state.state_by_node[node.up]
        parent_scale = state.state_scale_by_node[node.up]
        if transition.variance == 0.0:
            state.state_by_node[node] = parent_state
            state.state_scale_by_node[node] = transition.slope * parent_scale
            continue
        child_scale = _normalized_standard_deviation(transition.variance, normalization)
        coefficient = (
            transition.slope * parent_scale / child_scale if parent_state >= 0 else 0.0
        )
        state.state_by_node[node] = len(state.parents)
        state.state_scale_by_node[node] = child_scale
        state.parents.append(parent_state)
        state.coefficients.append(coefficient)
        state.innovations.append(1.0)
    return state


def _sparse_covariance_model(state, leaves, precision, innovations, normalization):
    tip_states = np.asarray([state.state_by_node[leaf] for leaf in leaves], dtype=int)
    if np.any(tip_states < 0):
        raise ValueError(
            "Every requested tip must have positive evolutionary variance."
        )
    tip_scales = np.asarray(
        [state.state_scale_by_node[leaf] for leaf in leaves], dtype=float
    )
    if np.any(~np.isfinite(tip_scales)):
        raise ValueError(
            "Sparse Gaussian covariance exceeds floating-point dynamic range."
        )
    tip_loading = sparse.csr_matrix(
        (tip_scales, (np.arange(len(leaves), dtype=int), tip_states)),
        shape=(len(leaves), len(state.parents)),
    )
    logdet_covariance = float(np.sum(np.log(innovations)))
    if not math.isfinite(logdet_covariance):
        raise ValueError(
            "Sparse Gaussian covariance log-determinant exceeds floating-point range."
        )
    return SparseCovarianceModel(
        precision=precision,
        tip_loading=tip_loading,
        logdet_covariance=logdet_covariance,
        sampling_parent=np.asarray(state.parents, dtype=int),
        sampling_transition=np.asarray(state.coefficients, dtype=float),
        sampling_variance=innovations,
        covariance_scale=normalization,
    )


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
            raise ValueError(
                "Gaussian tree transitions must cover every non-root node exactly."
            )
        object.__setattr__(self, "transitions", transitions)

    def scaled_variance(self, scale: float) -> "GaussianTreeProcess":
        """Scale every stochastic variance while preserving affine means."""

        scale = float(scale)
        if not math.isfinite(scale) or scale < 0.0:
            raise ValueError(
                "Gaussian process variance scales must be non-negative and finite."
            )
        if scale == 0.0 and self.root.mode in {"gaussian", "stationary"}:
            raise ValueError(
                "A proper Gaussian root prior cannot be scaled to zero variance."
            )
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
        if not self.root.is_proper:
            raise ValueError("A flat-root process has no finite unconditional moments.")
        lca = LcaIndex(self.tree)
        if any(node not in lca.index_by_node for node in nodes):
            raise ValueError(
                "Gaussian covariance requested a node outside the process tree."
            )
        requested = [lca.index_by_node[node] for node in nodes]
        assert self.root.variance is not None
        variances = np.empty(len(lca.nodes), dtype=float)
        cumulative_loading = np.ones(len(lca.nodes), dtype=float)
        path_log_abs = np.zeros(len(lca.nodes), dtype=float)
        path_zero_count = np.zeros(len(lca.nodes), dtype=np.int64)
        path_negative_parity = np.zeros(len(lca.nodes), dtype=bool)
        variances[0] = self.root.variance
        unit_loadings = True
        direct_loadings = True
        for index, node in enumerate(lca.nodes[1:], start=1):
            parent = lca.ancestors[0][index]
            transition = self.transitions[node]
            slope = transition.slope
            variances[index] = slope * slope * variances[parent] + transition.variance
            cumulative_loading[index] = cumulative_loading[parent] * slope
            unit_loadings = unit_loadings and slope == 1.0
            direct_loadings = direct_loadings and (
                cumulative_loading[index] != 0.0
                and math.isfinite(float(cumulative_loading[index]))
            )
            path_log_abs[index] = path_log_abs[parent]
            path_zero_count[index] = path_zero_count[parent]
            path_negative_parity[index] = path_negative_parity[parent]
            if slope == 0.0:
                path_zero_count[index] += 1
            else:
                path_log_abs[index] += math.log(abs(slope))
                path_negative_parity[index] ^= slope < 0.0
        if not np.isfinite(variances).all():
            raise ValueError("Gaussian tree moments exceed floating-point range.")

        def stable_covariance(first, second, ancestor):
            ancestor_variance = float(variances[ancestor])
            if ancestor_variance == 0.0:
                return 0.0
            if (
                path_zero_count[first] > path_zero_count[ancestor]
                or path_zero_count[second] > path_zero_count[ancestor]
            ):
                return 0.0
            log_value = (
                math.log(ancestor_variance)
                + path_log_abs[first]
                + path_log_abs[second]
                - 2.0 * path_log_abs[ancestor]
            )
            try:
                value = math.exp(log_value)
            except OverflowError:
                return math.inf
            negative = path_negative_parity[first] ^ path_negative_parity[second]
            return -value if negative else value

        covariance = np.empty((len(nodes), len(nodes)), dtype=float)
        common_ancestor = lca.common_ancestor_indices
        if unit_loadings:
            for first, first_index in enumerate(requested):
                covariance[first, first] = variances[first_index]
                for second in range(first + 1, len(nodes)):
                    ancestor = common_ancestor(first_index, requested[second])
                    value = variances[ancestor]
                    covariance[first, second] = value
                    covariance[second, first] = value
        elif direct_loadings:
            with np.errstate(over="ignore", under="ignore", divide="ignore"):
                ancestor_factor = variances / (cumulative_loading * cumulative_loading)
            direct_loadings = bool(np.isfinite(ancestor_factor).all())
            if direct_loadings:
                for first, first_index in enumerate(requested):
                    covariance[first, first] = variances[first_index]
                    first_loading = cumulative_loading[first_index]
                    for second in range(first + 1, len(nodes)):
                        second_index = requested[second]
                        ancestor = common_ancestor(first_index, second_index)
                        value = (
                            first_loading
                            * cumulative_loading[second_index]
                            * ancestor_factor[ancestor]
                        )
                        if not math.isfinite(float(value)) or (
                            value == 0.0 and variances[ancestor] > 0.0
                        ):
                            value = stable_covariance(
                                first_index, second_index, ancestor
                            )
                        covariance[first, second] = value
                        covariance[second, first] = value
        if not unit_loadings and not direct_loadings:
            for first, first_index in enumerate(requested):
                covariance[first, first] = variances[first_index]
                for second in range(first + 1, len(nodes)):
                    second_index = requested[second]
                    ancestor = common_ancestor(first_index, second_index)
                    value = stable_covariance(first_index, second_index, ancestor)
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
        leaves = tuple(leaf_by_name[name] for name in names)
        required = _required_ancestry(
            self.tree, leaves, len(leaves) == len(leaf_by_name)
        )
        state = _initial_sparse_state(self, required)
        normalization = _tip_variance_normalization(state, leaves, normalize)
        precision, normalized_innovation = _normal_sparse_precision(
            state, normalization
        )
        if precision is None:
            # A single global normalization can make a positive innovation
            # underflow to zero (and its precision overflow).  Re-express each
            # stochastic node in units of its own normalized standard deviation;
            # all innovation variances then equal one while tip loadings preserve
            # the exact normalized covariance.
            state = _rescaled_sparse_state(self, required, normalization)
            normalized_innovation = np.asarray(state.innovations, dtype=float)
            precision = _build_sparse_precision(
                state.parents, state.coefficients, normalized_innovation
            )
            if precision is None or not np.all(
                np.isfinite(list(state.state_scale_by_node.values()))
            ):
                raise ValueError(
                    "Sparse Gaussian covariance exceeds floating-point dynamic range."
                )
        return _sparse_covariance_model(
            state, leaves, precision, normalized_innovation, normalization
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
    if not all(
        math.isfinite(value) for value in (length, alpha, root_variance, optimum)
    ):
        raise ValueError("OU transition parameters must be finite.")
    if length < 0.0 or alpha <= 0.0 or root_variance < 0.0:
        raise ValueError(
            "OU requires non-negative lengths and root variance, and positive alpha."
        )
    exponent = alpha * length
    if not math.isfinite(exponent):
        raise ValueError(
            "OU alpha multiplied by a branch length exceeds floating-point range."
        )
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
        raise ValueError(
            "Rate change multiplied by tree depth exceeds floating-point range."
        )
    try:
        variance = math.exp(start_exponent) * math.expm1(edge_exponent) / rate
    except OverflowError as exc:
        raise ValueError(
            "A rate-change branch variance exceeds floating-point range."
        ) from exc
    if not math.isfinite(variance) or variance <= 0.0:
        raise ValueError("A positive rate-change branch variance is not finite.")
    return variance
