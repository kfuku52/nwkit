"""Linear-Gaussian ancestral reconstruction for branch-regime BM and OU."""

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import minimize

from nwkit.continuous_asr import (
    GaussianMarginal,
    _finite_number,
    _ldexp,
    compute_bm_marginals,
)
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
    brownian_transition,
    ou_transition_from_root_variance,
)
from nwkit.ou_asr import (
    _combine,
    _Factor,
    _log_normal_density,
    _prepare_data,
)

_LOG_2 = math.log(2.0)


_Transition = GaussianTransition


@dataclass(frozen=True)
class BrownianRegimeFit:
    sigma2_by_regime: dict[str, float]
    sigma2_estimated: bool
    sigma2_bounds_by_regime: dict[str, tuple[float, float]] | None
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


@dataclass(frozen=True)
class OrnsteinUhlenbeckRegimeFit:
    alpha: float
    alpha_estimated: bool
    theta_by_regime: dict[str, float]
    theta_estimated: bool
    sigma2: float
    sigma2_estimated: bool
    root_variance: float
    log_likelihood: float
    num_observed: int
    num_effective_observations: int
    num_observed_positions: int
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
    alpha_bounds: tuple[float, float]
    root_variance_bounds: tuple[float, float] | None


def _regime_by_group(tree, data, assignment):
    regimes = [assignment.root_regime] * len(data.parents)
    seen = set()
    for node in tree.traverse(strategy="preorder"):
        group = data.node_groups[node]
        if node.is_root or group in seen:
            continue
        seen.add(group)
        regimes[group] = assignment.by_node[node]
    return tuple(regimes)


def _propagate_up(factor, transition):
    if factor.marginal is None:
        return factor
    marginal = factor.marginal
    variance = marginal.variance + transition.variance
    if not math.isfinite(variance) or variance < 0.0:
        raise ValueError("A regime-model pruning variance is not finite.")
    slope = transition.slope
    if slope == 0.0:
        if variance == 0.0:
            if marginal.mean != transition.intercept:
                raise ValueError(
                    "Conflicting exact observations have zero likelihood under the regime model."
                )
            return _Factor(None, factor.log_weight)
        return _Factor(
            None,
            factor.log_weight
            + _log_normal_density(marginal.mean, transition.intercept, variance),
        )
    parent_mean = (marginal.mean - transition.intercept) / slope
    parent_variance = variance / slope / slope
    if not math.isfinite(parent_mean) or not math.isfinite(parent_variance):
        raise ValueError("A regime-model pruning message exceeds floating-point range.")
    return _Factor(
        GaussianMarginal(parent_mean, parent_variance),
        factor.log_weight - math.log(abs(slope)),
    )


def _predict_down(factor, transition):
    if factor.marginal is None:
        return factor
    marginal = factor.marginal
    mean = transition.slope * marginal.mean + transition.intercept
    variance = (
        transition.slope * transition.slope * marginal.variance + transition.variance
    )
    if not math.isfinite(mean) or not math.isfinite(variance) or variance < 0.0:
        raise ValueError("A regime-model smoothing message is not finite.")
    return _Factor(GaussianMarginal(mean, variance), factor.log_weight)


def _prune(data, transitions, root_prior):
    inside = list(data.local)
    upward = [_Factor(None) for _ in data.parents]
    for group in range(len(data.parents) - 1, -1, -1):
        factor = inside[group]
        for child in data.children[group]:
            factor = _combine(factor, upward[child])
        inside[group] = factor
        if data.parents[group] >= 0:
            upward[group] = _propagate_up(factor, transitions[group])
    if root_prior is None:
        if inside[0].marginal is None:
            raise ValueError("A flat root prior requires informative observations.")
        log_likelihood = inside[0].log_weight
    else:
        posterior_root = _combine(_Factor(root_prior), inside[0])
        log_likelihood = posterior_root.log_weight
    if not math.isfinite(log_likelihood):
        raise ValueError(
            "The observations have zero likelihood under the regime model."
        )
    return tuple(inside), tuple(upward), log_likelihood


def _smooth(data, transitions, root_prior, upward):
    outside = [_Factor(None) for _ in data.parents]
    outside[0] = _Factor(root_prior)
    posterior = [None for _ in data.parents]
    for group in range(len(data.parents)):
        base = _combine(outside[group], data.local[group])
        children = data.children[group]
        prefix = [base]
        for child in children:
            prefix.append(_combine(prefix[-1], upward[child]))
        posterior[group] = prefix[-1].marginal
        suffix = _Factor(None)
        for child_index in range(len(children) - 1, -1, -1):
            child = children[child_index]
            cavity = _combine(prefix[child_index], suffix)
            outside[child] = _predict_down(cavity, transitions[child])
            suffix = _combine(upward[child], suffix)
    if any(item is None for item in posterior):
        raise ValueError("Failed to calculate regime-model ancestral marginals.")
    return tuple(posterior)


def _restore_posterior(data, smoothed):
    restored = tuple(
        GaussianMarginal(
            data.trait_center
            + _ldexp(item.mean, data.trait_exponent, "Ancestral mean"),
            _ldexp(item.variance, 2 * data.trait_exponent, "Ancestral variance"),
        )
        for item in smoothed
    )
    posterior = {
        node: GaussianMarginal(data.exact_values[group], 0.0)
        if group in data.exact_values
        else restored[group]
        for node, group in data.node_groups.items()
    }
    if any(
        not math.isfinite(item.mean)
        or not math.isfinite(item.variance)
        or item.variance < 0.0
        for item in posterior.values()
    ):
        raise ValueError(
            "A regime-model ancestral marginal is outside floating-point range."
        )
    return posterior


def _rate_scale(data):
    values = np.asarray(data.observed_values, dtype=float)
    errors = np.asarray(data.observed_errors, dtype=float)
    positive_lengths = [value for value in data.lengths if value > 0.0]
    mean_length = float(np.mean(positive_lengths)) if positive_lengths else 1.0
    variance = float(np.var(values)) + float(np.mean(errors * errors))
    floor = max(1.0, float(np.max(np.abs(values)))) ** 2 * np.finfo(float).eps
    return max(variance / mean_length, floor / mean_length)


def _bms_regime_design_rank(tree, data, regime_assignment):
    observed_groups = {
        group for group, local in enumerate(data.local) if local.marginal is not None
    }
    representative: dict[int, Any] = {}
    for leaf in tree.leaves():
        group = data.node_groups[leaf]
        if group in observed_groups:
            representative.setdefault(group, leaf)
    ordered_groups = sorted(observed_groups)
    regimes = regime_assignment.regimes
    regime_index = {regime: index for index, regime in enumerate(regimes)}
    positive_lengths = [
        float(node.dist)
        for node in tree.traverse()
        if not node.is_root and float(node.dist) > 0.0
    ]
    length_scale = max(positive_lengths, default=1.0)
    rng = np.random.default_rng(1729)
    signatures = []
    for _ in range(max(8, 4 * len(regimes))):
        weights = rng.standard_normal(len(ordered_groups))
        weights -= np.mean(weights)
        norm = float(np.linalg.norm(weights))
        if norm == 0.0:
            continue
        weights /= norm
        local = {
            representative[group]: float(weight)
            for group, weight in zip(ordered_groups, weights, strict=True)
        }
        descendant: dict[Any, float] = {}
        for node in tree.traverse(strategy="postorder"):
            descendant[node] = local.get(node, 0.0) + math.fsum(
                descendant[child] for child in node.children
            )
        row = np.zeros(len(regimes), dtype=float)
        for node in tree.traverse():
            if node.is_root or float(node.dist) == 0.0:
                continue
            index = regime_index[regime_assignment.by_node[node]]
            row[index] += float(node.dist) / length_scale * descendant[node] ** 2
        signatures.append(row)
    matrix = np.vstack(signatures)
    norms = np.linalg.norm(matrix, axis=0)
    threshold = 1024.0 * np.finfo(float).eps * max(1.0, float(np.max(norms)))
    if np.any(norms <= threshold):
        return int(np.sum(norms > threshold))
    return int(np.linalg.matrix_rank(matrix / norms))


def _unique_starts(starts):
    unique = []
    seen = set()
    for start in starts:
        key = tuple(float(value) for value in start)
        if key not in seen:
            seen.add(key)
            unique.append(np.asarray(start, dtype=float))
    return unique


def _positive_multistart(objective, initial, lower, upper):
    initial_log = np.log(initial)
    lower_log = np.log(lower)
    upper_log = np.log(upper)
    starts = [initial_log]
    starts.extend(
        lower_log + fraction * (upper_log - lower_log)
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0)
    )
    if len(initial) > 1:
        starts.extend(
            [
                np.where(np.arange(len(initial)) % 2, lower_log, upper_log),
                np.where(np.arange(len(initial)) % 2, upper_log, lower_log),
            ]
        )
    unique = _unique_starts(starts)
    bounds = list(zip(lower_log, upper_log, strict=True))
    candidates = []
    failures = []
    for start in unique:
        result = minimize(objective, start, method="L-BFGS-B", bounds=bounds)
        messages = [str(result.message)]
        success = (
            bool(result.success)
            and math.isfinite(float(result.fun))
            and float(result.fun) < 1e99
            and np.all(np.isfinite(result.x))
        )
        if not success:
            fallback = minimize(
                objective,
                start,
                method="Powell",
                bounds=bounds,
                options={"xtol": 1e-8, "ftol": 1e-10, "maxiter": 800},
            )
            messages.append(f"Powell fallback: {fallback.message}")
            if (
                bool(fallback.success)
                and math.isfinite(float(fallback.fun))
                and float(fallback.fun) < 1e99
                and np.all(np.isfinite(fallback.x))
            ):
                result, success = fallback, True
        if success:
            candidates.append((float(result.fun), result, "; ".join(messages)))
        else:
            failures.append("; ".join(messages))
    if not candidates:
        raise ValueError(
            "Regime-model parameter optimization failed: "
            + ("; ".join(dict.fromkeys(failures)) or "no finite fit")
        )
    selected = min(candidates, key=lambda item: item[0])
    return selected, len(unique), len(candidates), len(failures)


def _fixed_bms_marginals(
    tree, values_by_leaf, regime_assignment, sigma2_by_regime, standard_errors
):
    transitions = {}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            continue
        transformed_length = sigma2_by_regime[regime_assignment.by_node[node]] * float(
            node.dist
        )
        if not math.isfinite(transformed_length) or transformed_length < 0.0:
            raise ValueError(
                "A BMS rate × branch length is outside floating-point range."
            )
        transitions[node] = brownian_transition(transformed_length)
    process = GaussianTreeProcess(
        tree=tree,
        transitions=transitions,
        root=GaussianRootPrior("flat", 0.0, None),
        model="bms",
    )
    posterior, bm_fit = compute_bm_marginals(
        tree,
        values_by_leaf,
        sigma2=1.0,
        standard_errors=standard_errors,
        _tree_validated=True,
        _process=process,
    )
    return posterior, bm_fit


def compute_bms_marginals(
    tree,
    values_by_leaf,
    regime_assignment,
    *,
    sigma2=None,
    sigma2_by_regime=None,
    regime_parameters_source=None,
    standard_errors=None,
    _tree_validated=False,
):
    """Return flat-root marginals under regime-specific Brownian rates."""

    if sigma2 is not None and sigma2_by_regime is not None:
        raise ValueError("--sigma2 cannot be combined with --regime-parameters.")
    fixed = sigma2 is not None or sigma2_by_regime is not None
    fixed_scale = 0.0
    if sigma2 is not None:
        sigma2 = _finite_number(sigma2, "--sigma2", nonnegative=True)
        sigma2_by_regime = dict.fromkeys(regime_assignment.regimes, sigma2)
    elif sigma2_by_regime is not None:
        if set(sigma2_by_regime) != set(regime_assignment.regimes):
            raise ValueError(
                "BMS fixed-rate regimes must exactly match the regime assignment."
            )
        sigma2_by_regime = {
            regime: _finite_number(
                sigma2_by_regime[regime],
                f"sigma2 for regime '{regime}'",
                nonnegative=True,
            )
            for regime in regime_assignment.regimes
        }
    if sigma2_by_regime:
        fixed_scale = math.sqrt(max(sigma2_by_regime.values(), default=0.0))
    data = _prepare_data(
        tree,
        values_by_leaf,
        standard_errors,
        process_sd=fixed_scale,
        tree_validated=_tree_validated,
    )
    regimes_by_group = _regime_by_group(tree, data, regime_assignment)

    if fixed:
        assert sigma2_by_regime is not None
        posterior, bm_fit = _fixed_bms_marginals(
            tree,
            values_by_leaf,
            regime_assignment,
            sigma2_by_regime,
            standard_errors,
        )
        singular_support = (
            bm_fit.num_effective_observations < data.num_effective_observations
        )
        fit = BrownianRegimeFit(
            sigma2_by_regime=dict(sigma2_by_regime),
            sigma2_estimated=False,
            sigma2_bounds_by_regime=None,
            restricted_log_likelihood=(
                None if singular_support else bm_fit.restricted_log_likelihood
            ),
            num_observed=bm_fit.num_observed,
            num_effective_observations=data.num_effective_observations,
            num_observed_positions=data.num_observed_positions,
            residual_df=data.num_effective_observations - 1,
            regimes=regime_assignment.regimes,
            root_regime=regime_assignment.root_regime,
            regime_map_source=regime_assignment.source,
            regime_parameters_source=regime_parameters_source,
            fit_status=(
                "singular_zero_boundary" if singular_support else bm_fit.fit_status
            ),
            optimizer_success=True,
            optimizer_message="all regime rates fixed",
            optimizer_starts=0,
            optimizer_converged_starts=0,
            optimizer_failed_starts=0,
        )
        return posterior, fit

    def transitions_for(scaled_rates):
        return tuple(
            _Transition(1.0, 0.0, scaled_rates[regimes_by_group[group]] * length)
            for group, length in enumerate(data.lengths)
        )

    if data.num_observed_positions - 1 < len(regime_assignment.regimes):
        raise ValueError(
            "BMS regime rates are not separately identifiable: the number of "
            "regimes must not exceed distinct observed positions minus one."
        )
    if _bms_regime_design_rank(tree, data, regime_assignment) < len(
        regime_assignment.regimes
    ):
        raise ValueError(
            "BMS regime rates are not separately identifiable from the "
            "observed branch-regime covariance design; fix rates or revise the regime map."
        )
    scale = _rate_scale(data)
    lower = np.full(len(regime_assignment.regimes), scale * 1e-8)
    upper = np.full(len(regime_assignment.regimes), scale * 1e8)
    initial = np.full(len(regime_assignment.regimes), scale)

    def objective(log_rates):
        try:
            values = np.exp(log_rates)
            rates = dict(zip(regime_assignment.regimes, values, strict=True))
            return -_prune(data, transitions_for(rates), None)[2]
        except (ValueError, ArithmeticError):
            return 1e100

    (objective_value, selected, message), starts, converged, failed = (
        _positive_multistart(objective, initial, lower, upper)
    )
    del objective_value
    values = np.exp(selected.x)
    scaled_rates = dict(zip(regime_assignment.regimes, values, strict=True))
    _, upward, log_likelihood = _prune(data, transitions_for(scaled_rates), None)
    fitted = {
        regime: _ldexp(value, 2 * data.trait_exponent, "Fitted BMS sigma2")
        for regime, value in scaled_rates.items()
    }
    bounds_by_regime = {
        regime: (
            _ldexp(lower[index], 2 * data.trait_exponent, "BMS sigma2 bound"),
            _ldexp(upper[index], 2 * data.trait_exponent, "BMS sigma2 bound"),
        )
        for index, regime in enumerate(regime_assignment.regimes)
    }
    tolerance = 1e-5
    lower_hit = any(values <= lower * (1.0 + tolerance))
    upper_hit = any(values >= upper * (1.0 - tolerance))
    statuses = []
    if lower_hit:
        statuses.append("sigma2_lower_boundary")
    if upper_hit:
        statuses.append("sigma2_upper_boundary")
    fit_status = "+".join(statuses) if statuses else "ok"
    optimizer_success = bool(selected.success)
    optimizer_message = (
        f"deterministic multistart: {converged}/{starts} starts converged; {message}"
    )
    transitions = transitions_for(scaled_rates)
    smoothed = _smooth(data, transitions, None, upward)
    posterior = _restore_posterior(data, smoothed)
    residual_df = data.num_effective_observations - 1
    restored_log_likelihood = (
        log_likelihood - residual_df * data.trait_exponent * _LOG_2
    )
    fit = BrownianRegimeFit(
        sigma2_by_regime=fitted,
        sigma2_estimated=True,
        sigma2_bounds_by_regime=bounds_by_regime,
        restricted_log_likelihood=restored_log_likelihood,
        num_observed=data.num_observed,
        num_effective_observations=data.num_effective_observations,
        num_observed_positions=data.num_observed_positions,
        residual_df=residual_df,
        regimes=regime_assignment.regimes,
        root_regime=regime_assignment.root_regime,
        regime_map_source=regime_assignment.source,
        regime_parameters_source=regime_parameters_source,
        fit_status=fit_status,
        optimizer_success=optimizer_success,
        optimizer_message=optimizer_message,
        optimizer_starts=starts,
        optimizer_converged_starts=converged,
        optimizer_failed_starts=failed,
    )
    return posterior, fit


def _ou_transitions(data, regimes_by_group, alpha, root_variance, theta_by_regime):
    transitions = []
    for group, length in enumerate(data.lengths):
        theta = theta_by_regime[regimes_by_group[group]]
        transitions.append(
            ou_transition_from_root_variance(length, alpha, root_variance, theta)
        )
    return tuple(transitions)


def _oum_theta_design_rank(tree, data, regime_assignment, alpha):
    regime_index = {
        regime: index for index, regime in enumerate(regime_assignment.regimes)
    }
    coefficients: dict[Any, np.ndarray] = {}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            value = np.zeros(len(regime_index), dtype=float)
            value[regime_index[regime_assignment.root_regime]] = 1.0
        else:
            decay = math.exp(-alpha * float(node.dist))
            value = decay * coefficients[node.up]
            value = value.copy()
            value[regime_index[regime_assignment.by_node[node]]] += 1.0 - decay
        coefficients[node] = value
    by_group: dict[int, np.ndarray] = {}
    for node, group in data.node_groups.items():
        by_group.setdefault(group, coefficients[node])
    matrix = np.vstack(
        [
            by_group[group]
            for group, local in enumerate(data.local)
            if local.marginal is not None
        ]
    )
    norms = np.linalg.norm(matrix, axis=0)
    threshold = 1024.0 * np.finfo(float).eps * max(1.0, float(np.max(norms)))
    if np.any(norms <= threshold):
        return int(np.sum(norms > threshold))
    normalized = matrix / norms
    return int(np.linalg.matrix_rank(normalized))


def compute_oum_marginals(
    tree,
    values_by_leaf,
    regime_assignment,
    *,
    alpha=None,
    sigma2=None,
    theta=None,
    theta_by_regime=None,
    regime_parameters_source=None,
    alpha_bounds=None,
    standard_errors=None,
    _tree_validated=False,
):
    """Return stationary-root OUM marginals via the shared Gaussian engine."""

    if theta is not None and theta_by_regime is not None:
        raise ValueError("--theta cannot be combined with --regime-parameters.")
    fixed_parameters = None
    if theta_by_regime is not None:
        if set(theta_by_regime) != set(regime_assignment.regimes):
            raise ValueError(
                "OUM fixed-optimum regimes must exactly match the regime assignment."
            )
        fixed_parameters = {
            regime: {"theta": value} for regime, value in theta_by_regime.items()
        }
    from nwkit.regime_gaussian_asr import compute_regime_ou_marginals

    posterior, shared_fit = compute_regime_ou_marginals(
        tree,
        values_by_leaf,
        regime_assignment,
        model="OUM",
        alpha=alpha,
        sigma2=sigma2,
        theta=theta,
        regime_parameters=fixed_parameters,
        regime_parameters_source=regime_parameters_source,
        alpha_bounds=alpha_bounds,
        standard_errors=standard_errors,
        _tree_validated=_tree_validated,
    )
    posterior = {
        node: GaussianMarginal(marginal.mean, marginal.variance)
        for node, marginal in posterior.items()
    }
    assert shared_fit.alpha is not None
    assert shared_fit.sigma2 is not None
    fit = OrnsteinUhlenbeckRegimeFit(
        alpha=shared_fit.alpha,
        alpha_estimated=shared_fit.alpha_estimated,
        theta_by_regime=dict(shared_fit.theta_by_regime),
        theta_estimated=shared_fit.theta_estimated,
        sigma2=shared_fit.sigma2,
        sigma2_estimated=shared_fit.sigma2_estimated,
        root_variance=shared_fit.root_variance,
        log_likelihood=shared_fit.log_likelihood,
        num_observed=shared_fit.num_observed,
        num_effective_observations=shared_fit.num_effective_observations,
        num_observed_positions=shared_fit.num_observed_positions,
        regimes=shared_fit.regimes,
        root_regime=shared_fit.root_regime,
        regime_map_source=shared_fit.regime_map_source,
        regime_parameters_source=shared_fit.regime_parameters_source,
        fit_status=shared_fit.fit_status,
        optimizer_success=shared_fit.optimizer_success,
        optimizer_message=shared_fit.optimizer_message,
        optimizer_starts=shared_fit.optimizer_starts,
        optimizer_converged_starts=shared_fit.optimizer_converged_starts,
        optimizer_failed_starts=shared_fit.optimizer_failed_starts,
        alpha_bounds=shared_fit.alpha_bounds,
        root_variance_bounds=shared_fit.root_variance_bounds,
    )
    return posterior, fit
