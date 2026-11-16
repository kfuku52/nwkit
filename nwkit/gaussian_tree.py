"""Shared linear-Gaussian processes on rooted trees.

The process layer deliberately separates biological model semantics from the
consumer.  ASR can retain a flat or stationary root and condition every node,
whereas regression can use a fixed root, profile fixed effects, and request a
normalized tip covariance.  Both consumers nevertheless use the same branch
transitions and innovation variances.
"""

import math
from dataclasses import dataclass, replace
from fractions import Fraction
from sys import float_info
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import sparse

from nwkit.clade_index import LcaIndex
from nwkit.sparse_laplace import SparseCovarianceModel

_MIN_NORMAL = float(np.finfo(float).tiny)


def _scaled_gaussian_variance(slope, variance):
    squared = slope * slope
    if _MIN_NORMAL <= squared < math.inf:
        return squared * variance
    return slope * (slope * variance)


def _recover_affine_mean(slope, mean, intercept, result):
    # A finite opposing intercept can recover an overflowing product.
    if not math.isfinite(mean) or math.copysign(1.0, result) == math.copysign(
        1.0, intercept
    ):
        return result
    try:
        return float(Fraction(slope) * Fraction(mean) + Fraction(intercept))
    except OverflowError:
        return result


_COVARIANCE_WORKSPACE_BYTES = 4 * 1024**2


def _propagated_variance(slope, variance, innovation):
    squared = slope * slope
    if squared < float_info.min or not math.isfinite(squared):
        # Applying the slope to the variance first keeps the intermediate
        # between the starting and ending magnitudes.
        if variance == 0.0:
            return innovation
        return slope * (slope * variance) + innovation
    return squared * variance + innovation


def _fill_covariance_blocks(
    covariance, lca, requested, variances, loading, factors, stable
):
    """Every requested pair is written at its LCA, with bounded product tiles."""
    local: list[list[int]] = [[] for _ in lca.nodes]
    for output, index in enumerate(requested):
        local[index].append(output)
    descendants: dict[int, np.ndarray] = {}
    output_loading = loading[np.asarray(requested, dtype=int)]
    for ancestor in range(len(lca.nodes) - 1, -1, -1):
        own = np.asarray(local[ancestor], dtype=int)
        groups = [own] if len(own) else []
        for child in lca.nodes[ancestor].children:
            group = descendants.pop(lca.index_by_node[child])
            if len(group):
                groups.append(group)
        merged = np.empty(sum(map(len, groups)), dtype=int)
        if len(own):
            own_value = variances[ancestor]
            if factors is not None:
                # Repeated requests for one node are off-diagonal pairs in
                # the original API: retain its loading-product arithmetic.
                with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                    own_value = (
                        loading[ancestor] * loading[ancestor] * factors[ancestor]
                    )
                if not np.isfinite(own_value) or (
                    own_value == 0 and variances[ancestor] > 0
                ):
                    own_value = stable(ancestor, ancestor, ancestor)
            covariance[own[:, None], own] = own_value
            covariance[own, own] = variances[ancestor]
        # Unit-loading blocks start in a zero-filled matrix. Still assign
        # negative zero explicitly so the pairwise and block paths agree.
        skip_cross_blocks = (
            factors is None
            and variances[ancestor] == 0.0
            and not np.signbit(variances[ancestor])
        )
        start = 0
        for group in groups:
            if start and not skip_cross_blocks:
                # Tile both axes, allowing for the product and finite/zero masks.
                max_entries = max(1, _COVARIANCE_WORKSPACE_BYTES // 24)
                columns = min(len(group), max_entries)
                rows = max(1, min(start, max_entries // columns))
                for second in range(0, len(group), columns):
                    right = group[second : second + columns]
                    for first in range(0, start, rows):
                        left = merged[first : min(start, first + rows)]
                        if factors is None:
                            covariance[left[:, None], right] = variances[ancestor]
                            covariance[right[:, None], left] = variances[ancestor]
                            continue
                        with np.errstate(
                            over="ignore", under="ignore", invalid="ignore"
                        ):
                            block = (
                                output_loading[left, None]
                                * output_loading[right][None, :]
                            )
                            block *= factors[ancestor]
                        bad = ~np.isfinite(block) | (
                            (block == 0) & (variances[ancestor] > 0)
                        )
                        for row, column in zip(*np.nonzero(bad), strict=True):
                            one, two = sorted((int(left[row]), int(right[column])))
                            block[row, column] = stable(
                                requested[one], requested[two], ancestor
                            )
                        covariance[left[:, None], right] = block
                        covariance[right[:, None], left] = block.T
                        del block, bad
            merged[start : start + len(group)] = group
            start += len(group)
        descendants[ancestor] = merged


@dataclass(frozen=True, slots=True)
class GaussianTransition:
    """One branch transition ``child = slope * parent + intercept + error``."""

    slope: float
    intercept: float
    variance: float

    def __post_init__(self) -> None:
        slope = float(self.slope)
        intercept = float(self.intercept)
        variance = float(self.variance)
        if not (
            math.isfinite(slope)
            and math.isfinite(intercept)
            and math.isfinite(variance)
        ):
            raise ValueError("Gaussian branch transitions must be finite.")
        if variance < 0.0:
            raise ValueError("Gaussian branch variances must be non-negative.")
        object.__setattr__(self, "slope", slope)
        object.__setattr__(self, "intercept", intercept)
        object.__setattr__(self, "variance", variance)

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
        state.marginal_variance[node] = _propagated_variance(
            transition.slope, state.marginal_variance[node.up], transition.variance
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
        or np.any(innovation_values <= 0.0)
    ):
        return None
    with np.errstate(over="ignore", under="ignore"):
        inverse = 1.0 / innovation_values
    if not np.isfinite(inverse).all():
        return None
    children = np.arange(size, dtype=int)
    has_parent = parent_values >= 0
    child_with_parent = children[has_parent]
    parents_with_child = parent_values[has_parent]
    coefficients_with_parent = coefficient_values[has_parent]
    inverse_with_parent = inverse[has_parent]
    parent_diagonal, cross = _sparse_precision_branch_terms(
        coefficients_with_parent, inverse_with_parent
    )
    if not np.isfinite(parent_diagonal).all() or not np.isfinite(cross).all():
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


def _sparse_precision_branch_terms(coefficients, inverse):
    with np.errstate(over="ignore", under="ignore"):
        cross = -coefficients * inverse
        squared = coefficients * coefficients
        parent_diagonal = squared * inverse
        extreme = (squared < _MIN_NORMAL) | ~np.isfinite(squared)
        parent_diagonal[extreme] = -coefficients[extreme] * cross[extreme]
    return parent_diagonal, cross


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


def _fill_direct_covariance(
    covariance, requested, lca, ancestor_factors, loadings, variances
):
    ancestor_at = lca.indexed_ancestors()
    count = len(requested)
    batch_size = max(1, 16_384 // count)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        for start in range(0, count, batch_size):
            stop = min(count, start + batch_size)
            first = requested[start:stop, None]
            ancestors = ancestor_at(first, requested[None, :])
            values = ancestor_factors[ancestors]
            if loadings is not None:
                values = loadings[first] * loadings[requested[None, :]] * values
                same_node = first == requested[None, :]
                values = np.where(same_node, variances[ancestors], values)
                if not np.isfinite(values).all():
                    return False
                if np.any(
                    (values == 0.0) & (ancestor_factors[ancestors] > 0.0) & ~same_node
                ):
                    return False
            covariance[start:stop, :] = values
    return True


def _finite_loading_factors(squares, factors, variances):
    # Squared loadings and their reciprocally scaled variances are only
    # intermediates. Subnormal intermediates can corrupt normal covariances.
    return bool(
        np.isfinite(factors).all()
        and np.isfinite(squares).all()
        and np.all(squares >= _MIN_NORMAL)
        and np.all((factors >= _MIN_NORMAL) | (variances == 0.0))
    )


def _gaussian_covariance_paths(lca, transitions, root_variance):
    variances = [root_variance]
    loadings = [1.0]
    unit_loadings = True
    direct_loadings = True
    for index, node in enumerate(lca.nodes[1:], start=1):
        parent = lca.ancestors[0][index]
        transition = transitions[node]
        slope = transition.slope
        variance = _scaled_gaussian_variance(slope, variances[parent])
        variances.append(variance + transition.variance)
        loading = loadings[parent] * slope
        loadings.append(loading)
        unit_loadings = unit_loadings and slope == 1.0
        direct_loadings = direct_loadings and loading != 0.0 and math.isfinite(loading)
    values = np.asarray(variances)
    if not np.isfinite(values).all():
        raise ValueError("Gaussian tree moments exceed floating-point range.")
    return values, np.asarray(loadings), unit_loadings, direct_loadings


def _log_covariance_paths(lca, transitions):
    log_abs = np.zeros(len(lca.nodes), dtype=float)
    zero_count = np.zeros(len(lca.nodes), dtype=np.int64)
    negative_parity = np.zeros(len(lca.nodes), dtype=bool)
    for index, node in enumerate(lca.nodes[1:], start=1):
        parent = lca.ancestors[0][index]
        slope = transitions[node].slope
        log_abs[index] = log_abs[parent]
        zero_count[index] = zero_count[parent]
        negative_parity[index] = negative_parity[parent]
        if slope == 0.0:
            zero_count[index] += 1
        else:
            log_abs[index] += math.log(abs(slope))
            negative_parity[index] ^= slope < 0.0
    return log_abs, zero_count, negative_parity


def _fill_affine_covariance(covariance, requested, lca, factors, loadings, variances):
    if len(requested) >= 8:
        return _fill_direct_covariance(
            covariance, np.asarray(requested), lca, factors, loadings, variances
        )
    for row, first in enumerate(requested):
        covariance[row, row] = variances[first]
        for column in range(row + 1, len(requested)):
            second = requested[column]
            ancestor = lca.common_ancestor_indices(first, second)
            factor = float(factors[ancestor])
            value = float(loadings[first]) * float(loadings[second]) * factor
            if not math.isfinite(value) or (value == 0.0 and factor > 0.0):
                return False
            covariance[row, column] = covariance[column, row] = value
    return True


def _fill_log_covariance(
    covariance, requested, lca, variances, log_loadings, zeros, negative
):
    ancestor_at = lca.indexed_ancestors()
    count = len(requested)
    batch_size = max(1, 16_384 // count)
    second = requested[None, :]
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        log_variances = np.log(variances)
        for start in range(0, count, batch_size):
            stop = min(count, start + batch_size)
            first = requested[start:stop, None]
            ancestors = ancestor_at(first, second)
            log_values = (
                log_variances[ancestors]
                + log_loadings[first]
                + log_loadings[second]
                - 2.0 * log_loadings[ancestors]
            )
            values = np.exp(log_values)
            values[negative[first] ^ negative[second]] *= -1.0
            disconnected = (zeros[first] > zeros[ancestors]) | (
                zeros[second] > zeros[ancestors]
            )
            values[disconnected | (variances[ancestors] == 0.0)] = 0.0
            covariance[start:stop] = np.where(
                first == second, variances[ancestors], values
            )
    # A symmetric pair must use one summation order for its two entries.
    for row in range(count):
        covariance[row + 1 :, row] = covariance[row, row + 1 :]


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
        if transitions.keys() != expected:
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
            parent_mean = means[node.up]
            mean = transition.slope * parent_mean + transition.intercept
            if not math.isfinite(mean):
                mean = _recover_affine_mean(
                    transition.slope, parent_mean, transition.intercept, mean
                )
            means[node] = mean
            variances[node] = (
                _scaled_gaussian_variance(transition.slope, variances[node.up])
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
        variances, cumulative_loading, unit_loadings, direct_loadings = (
            _gaussian_covariance_paths(lca, self.transitions, self.root.variance)
        )

        covariance = np.empty((len(nodes), len(nodes)), dtype=float)
        common_ancestor = lca.common_ancestor_indices
        if len(requested) >= 32 and unit_loadings:
            _fill_direct_covariance(
                covariance, np.asarray(requested), lca, variances, None, variances
            )
            return covariance
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
                loading_squares = cumulative_loading * cumulative_loading
                ancestor_factor = variances / loading_squares
            direct_loadings = _finite_loading_factors(
                loading_squares, ancestor_factor, variances
            )
            if direct_loadings and requested:
                if _fill_affine_covariance(
                    covariance,
                    requested,
                    lca,
                    ancestor_factor,
                    cumulative_loading,
                    variances,
                ):
                    return covariance
                direct_loadings = False
        if requested and not unit_loadings and not direct_loadings:
            _fill_log_covariance(
                covariance,
                np.asarray(requested, dtype=int),
                lca,
                variances,
                *_log_covariance_paths(lca, self.transitions),
            )
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


def _product_ratio(*factors: float, denominator: float = 1.0) -> float:
    """Multiply finite nonnegative factors without intermediate range loss."""
    if any(value == 0.0 for value in factors):
        return 0.0
    mantissa, exponent = 1.0, 0
    for value in factors:
        part, power = math.frexp(value)
        mantissa *= part
        exponent += power
    divisor, power = math.frexp(denominator)
    try:
        return math.ldexp(mantissa / divisor, exponent - power)
    except OverflowError as exc:
        raise ValueError("Gaussian transition overflows floating-point range.") from exc


def _mean_at_tiny_variance_ratio(mean1, mean2, variance1, variance2):
    """Retain a mean contribution even when its weight is subnormal."""
    # The precise weight rounds to one throughout this range.
    if variance2 < variance1:
        mean1, mean2, variance1, variance2 = mean2, mean1, variance2, variance1
    contribution = _product_ratio(abs(mean2), variance1, denominator=variance2)
    return mean1 + math.copysign(contribution, mean2)


def _stationary_ou_variance(alpha, variance_rate):
    divisor = 2.0 * alpha
    if math.isfinite(divisor):
        return variance_rate / divisor
    return (variance_rate / alpha) / 2


def ou_stationary_variance(alpha: float, variance_rate: float) -> float:
    """Convert OU diffusion to stationary variance for positive alpha."""

    denominator = 2.0 * alpha
    if math.isfinite(denominator):
        return variance_rate / denominator
    return 0.5 * (variance_rate / alpha)


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
    root_variance = ou_stationary_variance(alpha, variance_rate)
    if not math.isfinite(root_variance):
        raise ValueError("The stationary OU variance exceeds floating-point range.")
    return (
        _ou_transition_from_root_variance(
            length, alpha, root_variance, optimum, variance_rate
        ),
        root_variance,
    )


def ou_branch_transition(
    length: float, alpha: float, variance_rate: float = 1.0, optimum: float = 0.0
) -> GaussianTransition:
    """Return the OU transition without requiring a representable stationary root."""

    if not all(
        math.isfinite(value) for value in (length, alpha, variance_rate, optimum)
    ):
        raise ValueError("OU transition parameters must be finite.")
    if length < 0.0 or alpha <= 0.0 or variance_rate < 0.0:
        raise ValueError(
            "OU requires non-negative lengths and diffusion, and positive alpha."
        )
    root_variance = ou_stationary_variance(alpha, variance_rate)
    if 0.0 < root_variance < math.inf or variance_rate == 0.0:
        return ou_transition_from_root_variance(length, alpha, root_variance, optimum)
    if length == 0.0:
        return GaussianTransition(1.0, 0.0, 0.0)
    exponent = alpha * length
    doubled = 2.0 * exponent
    # expm1(-x)/(-x) tends to one as x approaches zero. Compute the finite
    # interval innovation directly before considering a stationary limit.
    if math.isfinite(doubled):
        correction = -math.expm1(-doubled) / doubled if doubled else 1.0
        innovation = variance_rate * (length * correction)
    else:
        innovation = root_variance
    if not math.isfinite(innovation):
        raise ValueError("The OU branch variance exceeds floating-point range.")
    if innovation == 0.0:
        raise ValueError("A positive OU branch variance underflowed.")
    intercept = (
        (-math.expm1(-exponent)) * optimum
        if exponent >= 2.0**-1022
        else _small_ou_increment(optimum, alpha, length, 0)
    )
    return GaussianTransition(math.exp(-exponent), intercept, innovation)


def ou_transition_from_root_variance(
    length: float,
    alpha: float,
    root_variance: float,
    optimum: float = 0.0,
) -> GaussianTransition:
    """Return an OU transition parameterized by stationary state variance."""
    return _ou_transition_from_root_variance(length, alpha, root_variance, optimum)


def _ou_transition_from_root_variance(
    length, alpha, root_variance, optimum, variance_rate=None
):
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
    # An overflowing positive exponent is already past the stationary limit
    # in floating point. exp and expm1 give the limiting transition at infinity.
    decay = math.exp(-exponent)
    fraction = -math.expm1(-2.0 * exponent)
    attenuation = -math.expm1(-exponent)
    innovation = root_variance * fraction
    intercept = attenuation * optimum
    if 0.0 < length and exponent < float_info.min:
        innovation = (
            _product_ratio(2.0, alpha, length, root_variance)
            if variance_rate is None
            else _product_ratio(length, variance_rate)
        )
        intercept = math.copysign(_product_ratio(alpha, length, abs(optimum)), optimum)
    elif root_variance < float_info.min and variance_rate is not None:
        innovation = _product_ratio(variance_rate, 0.5, fraction, denominator=alpha)
    has_diffusion = (
        root_variance > 0.0 if variance_rate is None else variance_rate > 0.0
    )
    if has_diffusion and length > 0.0 and innovation <= 0.0:
        raise ValueError("A positive OU branch variance underflowed.")
    return GaussianTransition(decay, intercept, innovation)


def _small_ou_increment(value, alpha, length, extra_exponent):
    mantissa_value, exponent_value = math.frexp(value)
    mantissa_alpha, exponent_alpha = math.frexp(alpha)
    mantissa_length, exponent_length = math.frexp(length)
    return math.ldexp(
        mantissa_value * mantissa_alpha * mantissa_length,
        exponent_value + exponent_alpha + exponent_length + extra_exponent,
    )


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
    # Use the direct integral while its numerator is comfortably representable.
    # On extreme scales the numerator can round to zero or infinity even when
    # dividing by the rate would give a normal branch variance.
    if abs(start_exponent) <= 700.0 and 1e-8 <= abs(edge_exponent) <= 700.0:
        numerator = math.exp(start_exponent) * math.expm1(edge_exponent)
        if math.isfinite(numerator) and abs(numerator) >= 2.0**-1022:
            variance = numerator / rate
            if math.isfinite(variance) and variance > 0.0:
                return variance
    if abs(edge_exponent) < 1e-8:
        factor = 1.0 + 0.5 * edge_exponent
        if abs(start_exponent) <= 700.0:
            variance = (length * factor) * math.exp(start_exponent)
            if math.isfinite(variance) and variance > 0.0:
                return variance
        log_variance = start_exponent + math.log(length) + math.log(factor)
    elif edge_exponent > 1.0:
        log_variance = (
            start_exponent
            + edge_exponent
            + math.log1p(-math.exp(-edge_exponent))
            - math.log(rate)
        )
    else:
        log_variance = (
            start_exponent
            + math.log(abs(math.expm1(edge_exponent)))
            - math.log(abs(rate))
        )
    try:
        variance = math.exp(log_variance)
    except OverflowError as exc:
        raise ValueError(
            "A rate-change branch variance exceeds floating-point range."
        ) from exc
    if not math.isfinite(variance) or variance <= 0.0:
        raise ValueError("A positive rate-change branch variance is not finite.")
    return variance


def _rate_edge_variance_from_logs(start_exponent, edge_exponent, rate):
    if edge_exponent > 1:
        log_variance = (start_exponent + edge_exponent) + math.log1p(
            -math.exp(-edge_exponent)
        )
    else:
        log_variance = start_exponent + math.log(abs(math.expm1(edge_exponent)))
    return _checked_rate_exp(log_variance - math.log(abs(rate)))


def _checked_rate_exp(log_variance):
    try:
        return math.exp(log_variance)
    except OverflowError as exc:
        raise ValueError(
            "A rate-change branch variance exceeds floating-point range."
        ) from exc
