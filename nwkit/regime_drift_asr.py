"""Regime-specific Brownian diffusion and directional drift ASR."""

import math
from dataclasses import dataclass
from typing import Mapping

import numpy as np

from nwkit.compiled_tree import CompiledTree
from nwkit.continuous_asr import _ldexp
from nwkit.gaussian_inference import (
    condition_gaussian_tree,
    gaussian_tree_likelihood,
)
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
)
from nwkit.optimization import MultistartResult, deterministic_multistart
from nwkit.rooting_state import require_rooted
from nwkit.util import validate_unique_named_leaves


@dataclass(frozen=True, slots=True)
class RegimeDriftFit:
    sigma2_by_regime: Mapping[str, float]
    sigma2_estimated: bool
    drift_by_regime: Mapping[str, float]
    drift_estimated: bool
    restricted_log_likelihood: float | None
    num_observed: int
    num_effective_observations: int
    num_observed_positions: int
    residual_df: int
    regimes: tuple[str, ...]
    root_regime: str
    regime_map_source: str
    regime_parameters_source: str | None
    fit_status: str
    optimizer_success: bool
    optimizer_message: str
    optimizer_starts: int
    optimizer_converged_starts: int
    optimizer_failed_starts: int
    sigma2_bounds_by_regime: Mapping[str, tuple[float, float]] | None


def _finite(value, label, *, positive=False):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be numeric, not boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be finite.") from exc
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "positive and finite" if positive else "finite"
        raise ValueError(f"{label} must be {qualifier}.")
    return result


def build_regime_drift_process(
    tree,
    regime_assignment,
    *,
    sigma2_by_regime,
    drift_by_regime,
):
    """Build ``child = parent + drift[r] t + N(0, sigma2[r] t)``."""

    regimes = regime_assignment.regimes
    if set(sigma2_by_regime) != set(regimes) or set(drift_by_regime) != set(regimes):
        raise ValueError("BMS-DRIFT parameters must exactly match the regime map.")
    sigma = {
        regime: _finite(
            sigma2_by_regime[regime], f"sigma2 for '{regime}'", positive=True
        )
        for regime in regimes
    }
    drift = {
        regime: _finite(drift_by_regime[regime], f"drift for '{regime}'")
        for regime in regimes
    }
    transitions = {}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            continue
        length = _finite(node.dist, "BMS-DRIFT branch lengths")
        if length < 0.0:
            raise ValueError("BMS-DRIFT branch lengths must be non-negative.")
        regime = regime_assignment.by_node[node]
        transitions[node] = GaussianTransition(
            1.0,
            drift[regime] * length,
            sigma[regime] * length,
        )
    return GaussianTreeProcess(
        tree,
        transitions,
        GaussianRootPrior("flat", 0.0, None),
        "bms-drift",
    )


def _path_exposures(tree, assignment):
    regimes = assignment.regimes
    index = {regime: position for position, regime in enumerate(regimes)}
    exposures = {tree: np.zeros(len(regimes), dtype=float)}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            continue
        values = exposures[node.up].copy()
        values[index[assignment.by_node[node]]] += float(node.dist)
        exposures[node] = values
    return exposures


def _drift_design_rank(observed_leaves, exposures, varying):
    if not len(varying):
        return 1
    matrix = np.column_stack(
        (
            np.ones(len(observed_leaves), dtype=float),
            np.vstack([exposures[leaf] for leaf in observed_leaves])[:, varying],
        )
    )
    return _regime_design_rank(matrix)


def _regime_design_rank(matrix):
    scales = np.max(np.abs(matrix), axis=0, initial=0.0)
    nonzero = scales > 0.0
    if not np.any(nonzero):
        return 0
    normalized = matrix[:, nonzero] / scales[nonzero]
    return int(np.linalg.matrix_rank(normalized / np.linalg.norm(normalized, axis=0)))


def _observed_split_regime_lengths(tree, observed_leaves, assignment):
    """Sum branch lengths by the bipartition they induce on observed tips."""
    regime_index = {regime: index for index, regime in enumerate(assignment.regimes)}
    observed = set(observed_leaves)
    num_observed = len(observed)
    # Observed tips in depth-first order form a contiguous interval per clade.
    descendant_counts = {}
    rows = {}
    seen = 0
    for node in tree.traverse(strategy="postorder"):
        local = int(node in observed)
        seen += local
        count = local + sum(descendant_counts.pop(child) for child in node.children)
        descendant_counts[node] = count
        if node.is_root or count in {0, num_observed}:
            continue
        # Orient each split away from the first observed tip. The complement
        # of any interval that contains the first tip is also contiguous.
        start = seen - count
        if start == 0:
            split = (seen, num_observed)
        else:
            split = (start, seen)
        row = rows.get(split)
        if row is None:
            row = np.zeros(len(regime_index), dtype=float)
            rows[split] = row
        row[regime_index[assignment.by_node[node]]] += float(node.dist)
    return np.asarray(list(rows.values())).reshape(-1, len(regime_index))


def _sigma_design_rank(tree, observed_leaves, assignment):
    # The pairwise distance vector determines the lengths of all distinct
    # splits of the observed subtree, so both designs have the same nullspace.
    return _regime_design_rank(
        _observed_split_regime_lengths(tree, observed_leaves, assignment)
    )


def _drift_fit_scales(tree, observed_values):
    positive_times = [
        float(node.dist)
        for node in tree.traverse()
        if not node.is_root and float(node.dist) > 0.0
    ]
    _, time_exponent = math.frexp(max(positive_times, default=1.0))
    _, trait_exponent = math.frexp(float(np.max(np.abs(observed_values))))
    normalized = np.ldexp(observed_values, -trait_exponent)
    variance = max(float(np.var(normalized)), np.finfo(float).eps)
    mean_time = (
        float(np.mean(np.ldexp(positive_times, -time_exponent)))
        if positive_times
        else 1.0
    )
    return math.log(variance / mean_time), trait_exponent, time_exponent


def compute_regime_drift_marginals(
    tree,
    values_by_leaf,
    regime_assignment,
    *,
    sigma2=None,
    drift=None,
    regime_parameters=None,
    regime_parameters_source=None,
    standard_errors=None,
    _tree_validated=False,
):
    """Fit and condition regime-specific Brownian rates and drifts."""

    if not _tree_validated:
        validate_unique_named_leaves(tree, option_name="--infile", context=" for 'asr'")
        require_rooted(tree, "ASR requires a rooted tree.")
    regimes = regime_assignment.regimes
    fixed_sigma: dict[str, float] | None
    fixed_drift: dict[str, float] | None
    if regime_parameters is not None:
        if sigma2 is not None or drift is not None:
            raise ValueError(
                "--regime-parameters cannot be combined with --sigma2 or --drift."
            )
        if set(regime_parameters) != set(regimes):
            raise ValueError(
                "BMS-DRIFT fixed parameters must exactly match the regime map."
            )
        fixed_sigma = {
            regime: _finite(
                regime_parameters[regime]["sigma2"],
                f"sigma2 for '{regime}'",
                positive=True,
            )
            for regime in regimes
        }
        fixed_drift = {
            regime: _finite(regime_parameters[regime]["drift"], f"drift for '{regime}'")
            for regime in regimes
        }
    else:
        fixed_sigma = (
            None
            if sigma2 is None
            else dict.fromkeys(regimes, _finite(sigma2, "--sigma2", positive=True))
        )
        fixed_drift = (
            None if drift is None else dict.fromkeys(regimes, _finite(drift, "--drift"))
        )
    observed_leaves = [
        leaf for leaf in tree.leaves() if values_by_leaf.get(str(leaf.name)) is not None
    ]
    observed_values = np.asarray(
        [float(values_by_leaf[str(leaf.name)]) for leaf in observed_leaves], dtype=float
    )
    if len(observed_values) < 2 or np.any(~np.isfinite(observed_values)):
        raise ValueError("BMS-DRIFT requires at least two finite observed tips.")
    varying = np.arange(len(regimes), dtype=int)
    if fixed_drift is None:
        exposures = _path_exposures(tree, regime_assignment)
        if _drift_design_rank(observed_leaves, exposures, varying) < len(regimes) + 1:
            raise ValueError(
                "BMS-DRIFT regime drifts are not identifiable from the observed "
                "root-to-tip regime-time design; fix drifts or revise the regime map."
            )
    if fixed_sigma is None and _sigma_design_rank(
        tree, observed_leaves, regime_assignment
    ) < len(regimes):
        raise ValueError(
            "BMS-DRIFT regime rates are not identifiable from observed pairwise "
            "path-regime lengths; fix rates or revise the regime map."
        )
    free_count = (0 if fixed_sigma is not None else len(regimes)) + (
        0 if fixed_drift is not None else len(regimes)
    )
    if len(observed_values) - 1 <= free_count:
        raise ValueError(
            "BMS-DRIFT has too many free parameters for the observed tips; "
            "fix rates/drifts or supply more data."
        )
    initial_log_rate, trait_exponent, time_exponent = _drift_fit_scales(
        tree, observed_values
    )
    drift_exponent = trait_exponent - time_exponent
    names: list[tuple[str, str]] = []
    initial: list[float] = []
    bounds: list[tuple[float | None, float | None]] = []
    if fixed_sigma is None:
        initial_log_rate += (2 * trait_exponent - time_exponent) * math.log(2.0)
        min_log_rate = math.log(math.ulp(0.0))
        max_log_rate = math.log(np.finfo(float).max)
        initial_log_rate = min(max_log_rate, max(min_log_rate, initial_log_rate))
        log_sigma_bounds = (
            max(initial_log_rate - math.log(1e8), min_log_rate),
            min(initial_log_rate + math.log(1e8), max_log_rate),
        )
        sigma_bounds = tuple(math.exp(bound) for bound in log_sigma_bounds)
        for regime in regimes:
            names.append(("sigma2", regime))
            initial.append(initial_log_rate)
            bounds.append(log_sigma_bounds)
    if fixed_drift is None:
        for regime in regimes:
            names.append(("drift", regime))
            initial.append(0.0)
            bounds.append((None, None))
    compiled = CompiledTree.from_tree(tree)

    def unpack(parameters):
        supplied = {
            name: float(value) for name, value in zip(names, parameters, strict=True)
        }
        sigma_map = (
            dict(fixed_sigma)
            if fixed_sigma is not None
            else {regime: math.exp(supplied[("sigma2", regime)]) for regime in regimes}
        )
        drift_map = (
            dict(fixed_drift)
            if fixed_drift is not None
            else {
                regime: _ldexp(
                    supplied[("drift", regime)], drift_exponent, "BMS-DRIFT trend"
                )
                for regime in regimes
            }
        )
        return sigma_map, drift_map

    def process_for(parameters):
        sigma_map, drift_map = unpack(parameters)
        return build_regime_drift_process(
            tree,
            regime_assignment,
            sigma2_by_regime=sigma_map,
            drift_by_regime=drift_map,
        )

    def objective(parameters):
        try:
            return -gaussian_tree_likelihood(
                process_for(parameters),
                values_by_leaf,
                standard_errors=standard_errors,
                compiled_tree=compiled,
            ).log_likelihood
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100

    if names:
        optimized = deterministic_multistart(objective, initial, bounds, maxiter=1200)
    else:
        fixed_parameters = np.asarray(initial, dtype=float)
        optimized = MultistartResult(
            fixed_parameters, 0.0, True, "all parameters fixed", 0, 0, 0
        )
    conditioned = condition_gaussian_tree(
        process_for(optimized.x),
        values_by_leaf,
        standard_errors=standard_errors,
        compiled_tree=compiled,
    )
    fitted_sigma, fitted_drift = unpack(optimized.x)
    residual_df = conditioned.likelihood_rank
    likelihood = (
        conditioned.log_likelihood if residual_df == len(observed_values) - 1 else None
    )
    statuses = []
    if likelihood is None:
        statuses.append("singular_support")
    if fixed_sigma is None:
        tolerance = 1e-5
        if any(
            value <= sigma_bounds[0] * (1.0 + tolerance)
            for value in fitted_sigma.values()
        ):
            statuses.append("sigma2_lower_boundary")
        if any(
            value >= sigma_bounds[1] * (1.0 - tolerance)
            for value in fitted_sigma.values()
        ):
            statuses.append("sigma2_upper_boundary")
    fit = RegimeDriftFit(
        sigma2_by_regime=fitted_sigma,
        sigma2_estimated=fixed_sigma is None,
        drift_by_regime=fitted_drift,
        drift_estimated=fixed_drift is None,
        restricted_log_likelihood=likelihood,
        num_observed=conditioned.num_observed,
        num_effective_observations=conditioned.num_observed,
        num_observed_positions=conditioned.num_observed_positions,
        residual_df=residual_df,
        regimes=regimes,
        root_regime=regime_assignment.root_regime,
        regime_map_source=regime_assignment.source,
        regime_parameters_source=regime_parameters_source,
        fit_status="+".join(dict.fromkeys(statuses)) if statuses else "ok",
        optimizer_success=optimized.success,
        optimizer_message=optimized.message,
        optimizer_starts=optimized.starts,
        optimizer_converged_starts=optimized.converged_starts,
        optimizer_failed_starts=optimized.failed_starts,
        sigma2_bounds_by_regime=(
            None if fixed_sigma is not None else dict.fromkeys(regimes, sigma_bounds)
        ),
    )
    return conditioned.marginals, fit
