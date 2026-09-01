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
    default_alpha_bounds,
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
        transition.slope * transition.slope * marginal.variance
        + transition.variance
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
        raise ValueError("The observations have zero likelihood under the regime model.")
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
        raise ValueError("A regime-model ancestral marginal is outside floating-point range.")
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
        group
        for group, local in enumerate(data.local)
        if local.marginal is not None
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
            row[index] += (
                float(node.dist) / length_scale * descendant[node] ** 2
            )
        signatures.append(row)
    matrix = np.vstack(signatures)
    norms = np.linalg.norm(matrix, axis=0)
    threshold = 1024.0 * np.finfo(float).eps * max(1.0, float(np.max(norms)))
    if np.any(norms <= threshold):
        return int(np.sum(norms > threshold))
    return int(np.linalg.matrix_rank(matrix / norms))


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
    unique = []
    seen = set()
    for start in starts:
        key = tuple(float(value) for value in start)
        if key not in seen:
            seen.add(key)
            unique.append(np.asarray(start, dtype=float))
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
        transformed_length = (
            sigma2_by_regime[regime_assignment.by_node[node]] * float(node.dist)
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
    if (
        _bms_regime_design_rank(tree, data, regime_assignment)
        < len(regime_assignment.regimes)
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
    _, upward, log_likelihood = _prune(
        data, transitions_for(scaled_rates), None
    )
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
    """Return stationary-root OU marginals with regime-specific optima."""

    if theta is not None and theta_by_regime is not None:
        raise ValueError("--theta cannot be combined with --regime-parameters.")
    bounds = default_alpha_bounds(tree) if alpha_bounds is None else alpha_bounds
    bounds = float(bounds[0]), float(bounds[1])
    if (
        not all(math.isfinite(value) and value > 0.0 for value in bounds)
        or bounds[0] >= bounds[1]
    ):
        raise ValueError("OUM alpha bounds must be two increasing positive finite values.")
    fixed_alpha = None if alpha is None else _finite_number(alpha, "--alpha")
    fixed_sigma2 = None if sigma2 is None else _finite_number(sigma2, "--sigma2")
    if fixed_alpha is not None and fixed_alpha <= 0.0:
        raise ValueError("--alpha must be strictly positive for stationary OUM.")
    if fixed_sigma2 is not None and fixed_sigma2 <= 0.0:
        raise ValueError("--sigma2 must be strictly positive for stationary OUM.")
    if theta is not None:
        theta = _finite_number(theta, "--theta")
        theta_by_regime = dict.fromkeys(regime_assignment.regimes, theta)
    elif theta_by_regime is not None:
        if set(theta_by_regime) != set(regime_assignment.regimes):
            raise ValueError(
                "OUM fixed-optimum regimes must exactly match the regime assignment."
            )
        theta_by_regime = {
            regime: _finite_number(theta_by_regime[regime], f"theta for '{regime}'")
            for regime in regime_assignment.regimes
        }
    theta_estimated = theta_by_regime is None
    process_sd = 0.0
    if fixed_alpha and fixed_sigma2:
        process_sd = math.sqrt(fixed_sigma2 / (2.0 * fixed_alpha))
    data = _prepare_data(
        tree,
        values_by_leaf,
        standard_errors,
        process_sd=process_sd,
        tree_validated=_tree_validated,
    )
    regimes_by_group = _regime_by_group(tree, data, regime_assignment)
    if theta_estimated:
        design_alpha = (
            fixed_alpha
            if fixed_alpha is not None
            else math.exp(
                0.5 * (math.log(bounds[0]) + math.log(bounds[1]))
            )
        )
        if (
            _oum_theta_design_rank(tree, data, regime_assignment, design_alpha)
            < len(regime_assignment.regimes)
        ):
            raise ValueError(
                "OUM regime optima are not separately identifiable from the "
                "observed branch-regime design; fix theta values or revise the regime map."
            )
    scaled_fixed_sigma2 = (
        None
        if fixed_sigma2 is None
        else _ldexp(fixed_sigma2, -2 * data.trait_exponent, "--sigma2")
    )
    scaled_fixed_theta = None
    if theta_by_regime is not None:
        scaled_fixed_theta = {}
        for regime, value in theta_by_regime.items():
            difference = value - data.trait_center
            if not math.isfinite(difference):
                raise ValueError("OUM theta minus the trait center is not finite.")
            scaled_fixed_theta[regime] = _ldexp(
                difference, -data.trait_exponent, "OUM theta"
            )
    free_names = []
    initial = []
    parameter_bounds: list[tuple[float | None, float | None]] = []
    if fixed_alpha is None:
        free_names.append("log_alpha")
        initial.append(0.5 * (math.log(bounds[0]) + math.log(bounds[1])))
        parameter_bounds.append((math.log(bounds[0]), math.log(bounds[1])))
    values = np.asarray(data.observed_values, dtype=float)
    errors = np.asarray(data.observed_errors, dtype=float)
    variance_scale = max(
        float(np.var(values)) + float(np.mean(errors * errors)),
        float(np.ptp(values)) ** 2 if len(values) > 1 else 0.0,
        np.finfo(float).eps,
    )
    root_variance_bounds = (variance_scale / 1e10, variance_scale * 1e10)
    if fixed_sigma2 is None:
        free_names.append("log_root_variance")
        initial.append(math.log(max(float(np.var(values)), variance_scale * 1e-3)))
        parameter_bounds.append(
            (math.log(root_variance_bounds[0]), math.log(root_variance_bounds[1]))
        )
    if theta_estimated:
        global_mean = float(np.mean(values))
        for regime in regime_assignment.regimes:
            free_names.append(f"theta:{regime}")
            initial.append(global_mean)
            parameter_bounds.append((None, None))
    if data.num_observed_positions < len(free_names):
        raise ValueError(
            "OUM free parameters are not separately identifiable: fix alpha, "
            "sigma2, or regime theta values, or supply more observed positions."
        )

    def unpack(parameters):
        supplied = dict(zip(free_names, parameters, strict=True))
        current_alpha = (
            math.exp(supplied["log_alpha"])
            if fixed_alpha is None
            else fixed_alpha
        )
        if fixed_sigma2 is None:
            root_variance = math.exp(supplied["log_root_variance"])
        else:
            assert scaled_fixed_sigma2 is not None
            root_variance = scaled_fixed_sigma2 / (2.0 * current_alpha)
        current_sigma2 = 2.0 * current_alpha * root_variance
        current_theta = (
            {regime: supplied[f"theta:{regime}"] for regime in regime_assignment.regimes}
            if theta_estimated
            else scaled_fixed_theta
        )
        return current_alpha, current_sigma2, root_variance, current_theta

    def objective(parameters):
        try:
            current_alpha, _, root_variance, current_theta = unpack(parameters)
            transitions = _ou_transitions(
                data,
                regimes_by_group,
                current_alpha,
                root_variance,
                current_theta,
            )
            root_prior = GaussianMarginal(
                current_theta[regime_assignment.root_regime], root_variance
            )
            return -_prune(data, transitions, root_prior)[2]
        except (ValueError, ArithmeticError):
            return 1e100

    starts = [np.asarray(initial, dtype=float)]
    candidates = []
    failures = []
    if free_names:
        covariance_indices = [
            index
            for index, name in enumerate(free_names)
            if name in {"log_alpha", "log_root_variance"}
        ]
        for fraction in (0.1, 0.5, 0.9):
            start = np.asarray(initial, dtype=float)
            for index in covariance_indices:
                lower, upper = parameter_bounds[index]
                assert lower is not None and upper is not None
                start[index] = lower + fraction * (upper - lower)
            starts.append(start)
        for start in starts:
            result = minimize(
                objective,
                start,
                method="L-BFGS-B",
                bounds=parameter_bounds,
                options={"maxiter": 1000, "ftol": 1e-12},
            )
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
                    bounds=parameter_bounds,
                    options={"maxiter": 1500, "xtol": 1e-8, "ftol": 1e-10},
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
                "Failed to estimate OUM parameters: "
                + ("; ".join(dict.fromkeys(failures)) or "no finite fit")
            )
        _, selected, selected_message = min(candidates, key=lambda item: item[0])
        selected_parameters = selected.x
        selected_success = bool(selected.success)
    else:
        selected_parameters = np.asarray([], dtype=float)
        selected_message = "all parameters fixed"
        selected_success = True
        candidates.append((objective(selected_parameters), None, selected_message))
    fitted_alpha, fitted_sigma2, root_variance, fitted_theta = unpack(
        selected_parameters
    )
    if theta_estimated and (
        _oum_theta_design_rank(tree, data, regime_assignment, fitted_alpha)
        < len(regime_assignment.regimes)
    ):
        raise ValueError(
            "OUM regime optima are not separately identifiable at the fitted alpha."
        )
    transitions = _ou_transitions(
        data, regimes_by_group, fitted_alpha, root_variance, fitted_theta
    )
    root_prior = GaussianMarginal(
        fitted_theta[regime_assignment.root_regime], root_variance
    )
    _, upward, log_likelihood = _prune(data, transitions, root_prior)
    smoothed = _smooth(data, transitions, root_prior, upward)
    posterior = _restore_posterior(data, smoothed)
    restored_theta = {
        regime: data.trait_center
        + _ldexp(value, data.trait_exponent, "Fitted OUM theta")
        for regime, value in fitted_theta.items()
    }
    restored_sigma2 = _ldexp(
        fitted_sigma2, 2 * data.trait_exponent, "Fitted OUM sigma2"
    )
    restored_root_variance = _ldexp(
        root_variance, 2 * data.trait_exponent, "Fitted OUM root variance"
    )
    restored_bounds = (
        None
        if fixed_sigma2 is not None
        else tuple(
            _ldexp(value, 2 * data.trait_exponent, "OUM root variance bound")
            for value in root_variance_bounds
        )
    )
    statuses = []
    tolerance = 1e-5
    if fixed_alpha is None:
        if fitted_alpha <= bounds[0] * (1.0 + tolerance):
            statuses.append("alpha_lower_boundary")
        elif fitted_alpha >= bounds[1] * (1.0 - tolerance):
            statuses.append("alpha_upper_boundary")
    if fixed_sigma2 is None:
        if root_variance <= root_variance_bounds[0] * (1.0 + tolerance):
            statuses.append("root_variance_lower_boundary")
        elif root_variance >= root_variance_bounds[1] * (1.0 - tolerance):
            statuses.append("root_variance_upper_boundary")
    fit = OrnsteinUhlenbeckRegimeFit(
        alpha=fitted_alpha,
        alpha_estimated=fixed_alpha is None,
        theta_by_regime=restored_theta,
        theta_estimated=theta_estimated,
        sigma2=restored_sigma2,
        sigma2_estimated=fixed_sigma2 is None,
        root_variance=restored_root_variance,
        log_likelihood=log_likelihood
        - data.likelihood_dimensions * data.trait_exponent * _LOG_2,
        num_observed=data.num_observed,
        num_effective_observations=data.num_effective_observations,
        num_observed_positions=data.num_observed_positions,
        regimes=regime_assignment.regimes,
        root_regime=regime_assignment.root_regime,
        regime_map_source=regime_assignment.source,
        regime_parameters_source=regime_parameters_source,
        fit_status="+".join(statuses) if statuses else "ok",
        optimizer_success=selected_success,
        optimizer_message=(
            f"deterministic multistart: {len(candidates)}/{len(starts)} starts "
            f"converged; {selected_message}"
        ),
        optimizer_starts=len(starts),
        optimizer_converged_starts=len(candidates),
        optimizer_failed_starts=len(failures),
        alpha_bounds=bounds,
        root_variance_bounds=restored_bounds,
    )
    return posterior, fit
