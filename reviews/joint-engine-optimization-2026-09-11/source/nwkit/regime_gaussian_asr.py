"""Shared affine-Gaussian ASR for OU regime models.

The four Hansen/OUwie-style model names differ only in which branch-regime
parameters vary: every model has regime-specific optima, OUMA additionally
varies attraction, OUMV varies diffusion, and OUMVA varies both.  All of them
are conditioned by :mod:`nwkit.gaussian_inference`.
"""

import math
from dataclasses import dataclass
from typing import Mapping

import numpy as np

from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_inference import (
    condition_gaussian_tree,
    gaussian_process_parameter_rank,
    gaussian_tree_likelihood,
)
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTreeProcess,
    ou_transition,
)
from nwkit.optimization import MultistartResult, deterministic_multistart
from nwkit.ou_asr import compute_ou_marginals, default_alpha_bounds
from nwkit.rooting_state import require_rooted
from nwkit.util import validate_unique_named_leaves

REGIME_OU_MODELS = ("OUM", "OUMA", "OUMV", "OUMVA")


@dataclass(frozen=True, slots=True)
class RegimeOUFit:
    model: str
    alpha: float | None
    alpha_by_regime: Mapping[str, float]
    alpha_estimated: bool
    theta_by_regime: Mapping[str, float]
    theta_estimated: bool
    sigma2: float | None
    sigma2_by_regime: Mapping[str, float]
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
    root_variance_bounds: tuple[float, float] | None = None


@dataclass(frozen=True, slots=True)
class _RegimeOUSetup:
    model: str
    regimes: tuple[str, ...]
    alpha_bounds: tuple[float, float]
    sigma_bounds: tuple[float, float]
    root_variance_bounds: tuple[float, float]
    initial_alpha: float
    initial_sigma: float
    initial_root_variance: float
    observed_mean: float
    num_observed: int
    varying_alpha: bool
    varying_sigma: bool
    fixed_alpha: Mapping[str, float] | None
    fixed_sigma: Mapping[str, float] | None
    fixed_theta: Mapping[str, float] | None


def regime_parameter_columns(model):
    """Return columns required by a complete fixed regime-parameter table."""

    model = str(model).upper()
    if model not in REGIME_OU_MODELS:
        raise ValueError(f"Unsupported OU regime model: {model}.")
    columns = ["theta"]
    if model in {"OUMA", "OUMVA"}:
        columns.append("alpha")
    if model in {"OUMV", "OUMVA"}:
        columns.append("sigma2")
    return tuple(columns)


def _positive(value, label):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be numeric, not boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be positive and finite.") from exc
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{label} must be positive and finite.")
    return result


def _finite(value, label):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be numeric, not boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be finite.") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite.")
    return result


def build_regime_ou_process(
    tree,
    regime_assignment,
    *,
    alpha_by_regime,
    sigma2_by_regime,
    theta_by_regime,
):
    """Build one stationary-root OU process with edge-specific regimes."""

    regimes = regime_assignment.regimes
    expected = set(regimes)
    supplied = {
        "alpha": set(alpha_by_regime),
        "sigma2": set(sigma2_by_regime),
        "theta": set(theta_by_regime),
    }
    mismatch = [name for name, values in supplied.items() if values != expected]
    if mismatch:
        raise ValueError(
            "Regime process parameters must exactly match the regime map: "
            + ", ".join(mismatch)
        )
    alpha_values = {
        regime: _positive(alpha_by_regime[regime], f"alpha for '{regime}'")
        for regime in regimes
    }
    sigma_values = {
        regime: _positive(sigma2_by_regime[regime], f"sigma2 for '{regime}'")
        for regime in regimes
    }
    theta_values = {
        regime: _finite(theta_by_regime[regime], f"theta for '{regime}'")
        for regime in regimes
    }
    transitions = {}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            continue
        regime = regime_assignment.by_node[node]
        try:
            length = float(node.dist)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("OU regime branch lengths must be numeric.") from exc
        if not math.isfinite(length) or length < 0.0:
            raise ValueError(
                "OU regime branch lengths must be non-negative and finite."
            )
        transition, _ = ou_transition(
            length,
            alpha_values[regime],
            sigma_values[regime],
            theta_values[regime],
        )
        transitions[node] = transition
    root_regime = regime_assignment.root_regime
    root_variance = sigma_values[root_regime] / (2.0 * alpha_values[root_regime])
    if not math.isfinite(root_variance) or root_variance <= 0.0:
        raise ValueError("The OU regime root variance is outside floating-point range.")
    return GaussianTreeProcess(
        tree=tree,
        transitions=transitions,
        root=GaussianRootPrior("stationary", theta_values[root_regime], root_variance),
        model="regime-ou",
    )


def _validated_fixed(parameters, regimes, columns):
    if parameters is None:
        return None
    if set(parameters) != set(regimes):
        raise ValueError(
            "Fixed OU regime parameters must exactly match the regime assignment."
        )
    validated = {}
    for regime in regimes:
        values = parameters[regime]
        if set(values) != set(columns):
            raise ValueError(
                "Fixed OU regime parameter columns must be: " + ", ".join(columns)
            )
        validated[regime] = {
            name: (
                _positive(values[name], f"{name} for '{regime}'")
                if name in {"alpha", "sigma2"}
                else _finite(values[name], f"{name} for '{regime}'")
            )
            for name in columns
        }
    return validated


def _validate_regime_parameter_conflicts(fixed_table, columns, alpha, sigma2, theta):
    if fixed_table is None:
        return
    conflicts = []
    if theta is not None:
        conflicts.append("--theta")
    if "alpha" in columns and alpha is not None:
        conflicts.append("--alpha")
    if "sigma2" in columns and sigma2 is not None:
        conflicts.append("--sigma2")
    if conflicts:
        raise ValueError(
            "--regime-parameters cannot be combined with " + ", ".join(conflicts) + "."
        )


def _regime_fixed_maps(
    regimes, fixed_table, varying_alpha, varying_sigma, alpha, sigma2, theta
):
    alpha_map = None
    sigma_map = None
    theta_map = None
    if fixed_table is not None:
        theta_map = {regime: fixed_table[regime]["theta"] for regime in regimes}
        if varying_alpha:
            alpha_map = {regime: fixed_table[regime]["alpha"] for regime in regimes}
        if varying_sigma:
            sigma_map = {regime: fixed_table[regime]["sigma2"] for regime in regimes}
    if alpha is not None:
        alpha_map = dict.fromkeys(regimes, _positive(alpha, "--alpha"))
    if sigma2 is not None:
        sigma_map = dict.fromkeys(regimes, _positive(sigma2, "--sigma2"))
    if theta is not None:
        theta_map = dict.fromkeys(regimes, _finite(theta, "--theta"))
    return alpha_map, sigma_map, theta_map


def _informative_edge_regimes(tree, values_by_leaf, regime_assignment):
    informative = {}
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            informative[node] = values_by_leaf.get(str(node.name)) is not None
        else:
            informative[node] = any(informative[child] for child in node.children)
    return {
        regime_assignment.by_node[node]
        for node in tree.traverse()
        if not node.is_root and float(node.dist) > 0.0 and informative[node]
    }


def _validate_regime_ou_representation(
    tree,
    values_by_leaf,
    regime_assignment,
    *,
    varying_alpha,
    varying_sigma,
    fixed_alpha,
    fixed_sigma,
    fixed_theta,
):
    represented = _informative_edge_regimes(tree, values_by_leaf, regime_assignment)
    root_regime = regime_assignment.root_regime
    required = set()
    if fixed_theta is None:
        required.update(regime_assignment.regimes)
    if varying_alpha and fixed_alpha is None:
        required.update(regime_assignment.regimes)
    if varying_sigma and fixed_sigma is None:
        required.update(regime_assignment.regimes)
    missing = sorted(required - represented - {root_regime})
    if missing:
        raise ValueError(
            "Free OU parameters require each non-root regime to occur on a "
            "positive-length branch ancestral to an observation; unrepresented "
            "regime(s): " + ", ".join(missing)
        )
    if (
        varying_alpha
        and varying_sigma
        and fixed_alpha is None
        and fixed_sigma is None
        and root_regime not in represented
    ):
        raise ValueError(
            "Root-regime alpha and sigma2 are confounded when that regime occurs "
            "only at the root; fix one component or map the root regime to an "
            "informative positive-length branch."
        )


def _regime_ou_scales(tree, observed, alpha_bounds):
    observed_mean = float(np.mean(observed))
    trait_scale = max(
        math.sqrt(float(np.var(observed))),
        float(np.ptp(observed)) / 4.0,
        abs(observed_mean) * math.sqrt(np.finfo(float).eps),
        math.sqrt(np.finfo(float).eps),
    )
    positive_lengths = [
        float(node.dist)
        for node in tree.traverse()
        if not node.is_root and float(node.dist) > 0.0
    ]
    mean_length = float(np.mean(positive_lengths)) if positive_lengths else 1.0
    initial_alpha = math.sqrt(alpha_bounds[0] * alpha_bounds[1])
    initial_sigma = max(
        2.0 * initial_alpha * trait_scale * trait_scale, np.finfo(float).tiny
    )
    sigma_bounds = (
        max(initial_sigma * 1e-8, np.finfo(float).tiny),
        min(initial_sigma * 1e8, np.finfo(float).max / max(1.0, mean_length)),
    )
    variance_scale = max(
        float(np.var(observed)),
        float(np.ptp(observed)) ** 2 if len(observed) > 1 else 0.0,
        np.finfo(float).eps,
    )
    initial_root_variance = max(float(np.var(observed)), variance_scale * 1e-3)
    root_variance_bounds = (variance_scale / 1e10, variance_scale * 1e10)
    return (
        observed_mean,
        initial_alpha,
        initial_sigma,
        sigma_bounds,
        initial_root_variance,
        root_variance_bounds,
    )


def _prepare_regime_ou_setup(
    tree,
    values_by_leaf,
    regime_assignment,
    model,
    alpha,
    sigma2,
    theta,
    regime_parameters,
    alpha_bounds,
    tree_validated,
):
    model = str(model).upper()
    columns = regime_parameter_columns(model)
    regimes = tuple(regime_assignment.regimes)
    if not tree_validated:
        validate_unique_named_leaves(tree, option_name="--infile", context=" for 'asr'")
        require_rooted(tree, "ASR requires a rooted tree.")
    if not regimes:
        raise ValueError("An OU regime model requires at least one regime.")
    fixed_table = _validated_fixed(regime_parameters, regimes, columns)
    _validate_regime_parameter_conflicts(fixed_table, columns, alpha, sigma2, theta)
    bounds = default_alpha_bounds(tree) if alpha_bounds is None else alpha_bounds
    bounds = float(bounds[0]), float(bounds[1])
    if (
        not all(math.isfinite(value) and value > 0.0 for value in bounds)
        or bounds[0] >= bounds[1]
    ):
        raise ValueError("OU regime alpha bounds must be increasing and positive.")
    observed = np.asarray(
        [float(value) for value in values_by_leaf.values() if value is not None],
        dtype=float,
    )
    if not len(observed) or not np.all(np.isfinite(observed)):
        raise ValueError("OU regime ASR requires finite observed tip values.")
    (
        observed_mean,
        initial_alpha,
        initial_sigma,
        sigma_bounds,
        initial_root_variance,
        root_variance_bounds,
    ) = _regime_ou_scales(tree, observed, bounds)
    varying_alpha = model in {"OUMA", "OUMVA"}
    varying_sigma = model in {"OUMV", "OUMVA"}
    fixed_alpha, fixed_sigma, fixed_theta = _regime_fixed_maps(
        regimes,
        fixed_table,
        varying_alpha,
        varying_sigma,
        alpha,
        sigma2,
        theta,
    )
    _validate_regime_ou_representation(
        tree,
        values_by_leaf,
        regime_assignment,
        varying_alpha=varying_alpha,
        varying_sigma=varying_sigma,
        fixed_alpha=fixed_alpha,
        fixed_sigma=fixed_sigma,
        fixed_theta=fixed_theta,
    )
    return _RegimeOUSetup(
        model=model,
        regimes=regimes,
        alpha_bounds=bounds,
        sigma_bounds=sigma_bounds,
        root_variance_bounds=root_variance_bounds,
        initial_alpha=initial_alpha,
        initial_sigma=initial_sigma,
        initial_root_variance=initial_root_variance,
        observed_mean=observed_mean,
        num_observed=len(observed),
        varying_alpha=varying_alpha,
        varying_sigma=varying_sigma,
        fixed_alpha=fixed_alpha,
        fixed_sigma=fixed_sigma,
        fixed_theta=fixed_theta,
    )


def _regime_ou_parameter_vectors(tree, values_by_leaf, regime_assignment, setup):
    names: list[tuple[str, str]] = []
    initial: list[float] = []
    parameter_bounds: list[tuple[float | None, float | None]] = []
    if setup.fixed_alpha is None:
        keys = setup.regimes if setup.varying_alpha else ("__shared__",)
        for regime in keys:
            names.append(("alpha", regime))
            initial.append(math.log(setup.initial_alpha))
            parameter_bounds.append(
                (
                    math.log(setup.alpha_bounds[0]),
                    math.log(setup.alpha_bounds[1]),
                )
            )
    if setup.fixed_sigma is None:
        if setup.model == "OUM":
            names.append(("root_variance", "__shared__"))
            initial.append(math.log(setup.initial_root_variance))
            parameter_bounds.append(
                (
                    math.log(setup.root_variance_bounds[0]),
                    math.log(setup.root_variance_bounds[1]),
                )
            )
        else:
            keys = setup.regimes if setup.varying_sigma else ("__shared__",)
            for regime in keys:
                names.append(("sigma2", regime))
                initial.append(math.log(setup.initial_sigma))
                parameter_bounds.append(
                    (
                        math.log(setup.sigma_bounds[0]),
                        math.log(setup.sigma_bounds[1]),
                    )
                )
    if setup.fixed_theta is None:
        _add_regime_theta_parameters(
            tree,
            values_by_leaf,
            regime_assignment,
            setup,
            names,
            initial,
            parameter_bounds,
        )
    return names, initial, parameter_bounds


def _add_regime_theta_parameters(
    tree,
    values_by_leaf,
    regime_assignment,
    setup,
    names,
    initial,
    parameter_bounds,
):
    for regime in setup.regimes:
        names.append(("theta", regime))
        regime_tips = [
            float(values_by_leaf[str(leaf.name)])
            for leaf in tree.leaves()
            if values_by_leaf.get(str(leaf.name)) is not None
            and regime_assignment.by_node[leaf] == regime
        ]
        initial.append(
            float(np.mean(regime_tips)) if regime_tips else setup.observed_mean
        )
        parameter_bounds.append((None, None))


def _unpack_regime_ou_parameters(parameters, names, setup):
    supplied = {
        name: float(value) for name, value in zip(names, parameters, strict=True)
    }
    if setup.fixed_alpha is not None:
        alpha_map = dict(setup.fixed_alpha)
    elif setup.varying_alpha:
        alpha_map = {
            regime: math.exp(supplied[("alpha", regime)]) for regime in setup.regimes
        }
    else:
        alpha_map = dict.fromkeys(
            setup.regimes, math.exp(supplied[("alpha", "__shared__")])
        )
    if setup.fixed_sigma is not None:
        sigma_map = dict(setup.fixed_sigma)
    elif setup.model == "OUM":
        root_variance = math.exp(supplied[("root_variance", "__shared__")])
        shared_alpha = next(iter(alpha_map.values()))
        sigma_map = dict.fromkeys(setup.regimes, 2.0 * shared_alpha * root_variance)
    elif setup.varying_sigma:
        sigma_map = {
            regime: math.exp(supplied[("sigma2", regime)]) for regime in setup.regimes
        }
    else:
        sigma_map = dict.fromkeys(
            setup.regimes, math.exp(supplied[("sigma2", "__shared__")])
        )
    theta_map = (
        dict(setup.fixed_theta)
        if setup.fixed_theta is not None
        else {regime: supplied[("theta", regime)] for regime in setup.regimes}
    )
    return alpha_map, sigma_map, theta_map


def _regime_ou_boundary_status(setup, fitted_alpha, fitted_sigma):
    tolerance = 1e-5
    statuses = []
    checks = [(setup.fixed_alpha, fitted_alpha, setup.alpha_bounds, "alpha")]
    if setup.model != "OUM":
        checks.append((setup.fixed_sigma, fitted_sigma, setup.sigma_bounds, "sigma2"))
    for fixed, values, bounds, label in checks:
        if fixed is not None:
            continue
        if any(value <= bounds[0] * (1.0 + tolerance) for value in values.values()):
            statuses.append(f"{label}_lower_boundary")
        if any(value >= bounds[1] * (1.0 - tolerance) for value in values.values()):
            statuses.append(f"{label}_upper_boundary")
    if setup.model == "OUM" and setup.fixed_sigma is None:
        root_variance = next(iter(fitted_sigma.values())) / (
            2.0 * next(iter(fitted_alpha.values()))
        )
        if root_variance <= setup.root_variance_bounds[0] * (1.0 + tolerance):
            statuses.append("root_variance_lower_boundary")
        elif root_variance >= setup.root_variance_bounds[1] * (1.0 - tolerance):
            statuses.append("root_variance_upper_boundary")
    return "+".join(dict.fromkeys(statuses)) if statuses else "ok"


def _regime_theta_design_rank(
    tree, values_by_leaf, regime_assignment, alpha_map, sigma_map
):
    observed_leaves = [
        leaf for leaf in tree.leaves() if values_by_leaf.get(str(leaf.name)) is not None
    ]
    columns = []
    for focal_regime in regime_assignment.regimes:
        theta_map = {
            regime: float(regime == focal_regime)
            for regime in regime_assignment.regimes
        }
        process = build_regime_ou_process(
            tree,
            regime_assignment,
            alpha_by_regime=alpha_map,
            sigma2_by_regime=sigma_map,
            theta_by_regime=theta_map,
        )
        means, _ = process.marginal_moments()
        columns.append([means[leaf] for leaf in observed_leaves])
    matrix = np.asarray(columns, dtype=float).T
    norms = np.linalg.norm(matrix, axis=0)
    threshold = (
        1024.0 * np.finfo(float).eps * max(1.0, float(np.max(norms, initial=0.0)))
    )
    nonzero = norms > threshold
    if not np.any(nonzero):
        return 0
    return int(np.linalg.matrix_rank(matrix[:, nonzero] / norms[nonzero]))


def _validate_regime_theta_design(
    tree, values_by_leaf, regime_assignment, setup, alpha_map, sigma_map
):
    if setup.fixed_theta is not None:
        return
    if _regime_theta_design_rank(
        tree, values_by_leaf, regime_assignment, alpha_map, sigma_map
    ) < len(setup.regimes):
        raise ValueError(
            "OU regime optima are not separately identifiable from the observed "
            "branch-regime design; fix theta values or revise the regime map."
        )


def _single_regime_ou_reduction(
    tree,
    values_by_leaf,
    regime_assignment,
    setup,
    *,
    regime_parameters_source,
    standard_errors,
):
    """Use the canonical OU1 fitter when every regime parameter is shared."""

    regime = setup.regimes[0]

    def fixed_value(values):
        return None if values is None else values[regime]

    marginals, ordinary = compute_ou_marginals(
        tree,
        values_by_leaf,
        alpha=fixed_value(setup.fixed_alpha),
        sigma2=fixed_value(setup.fixed_sigma),
        theta=fixed_value(setup.fixed_theta),
        alpha_bounds=setup.alpha_bounds,
        standard_errors=standard_errors,
        _tree_validated=True,
    )
    alpha_by_regime = {regime: ordinary.alpha}
    sigma_by_regime = {regime: ordinary.sigma2}
    theta_by_regime = {regime: ordinary.theta}
    return marginals, RegimeOUFit(
        model=setup.model,
        alpha=None if setup.varying_alpha else ordinary.alpha,
        alpha_by_regime=alpha_by_regime,
        alpha_estimated=ordinary.alpha_estimated,
        theta_by_regime=theta_by_regime,
        theta_estimated=ordinary.theta_estimated,
        sigma2=None if setup.varying_sigma else ordinary.sigma2,
        sigma2_by_regime=sigma_by_regime,
        sigma2_estimated=ordinary.sigma2_estimated,
        root_variance=ordinary.root_variance,
        log_likelihood=ordinary.log_likelihood,
        num_observed=ordinary.num_observed,
        num_effective_observations=ordinary.num_effective_observations,
        num_observed_positions=ordinary.num_observed_positions,
        regimes=setup.regimes,
        root_regime=regime_assignment.root_regime,
        regime_map_source=regime_assignment.source,
        regime_parameters_source=regime_parameters_source,
        fit_status=ordinary.fit_status,
        optimizer_success=ordinary.optimizer_success,
        optimizer_message=(
            "single-regime reduction to stationary OU; " + ordinary.optimizer_message
        ),
        optimizer_starts=ordinary.optimizer_starts,
        optimizer_converged_starts=ordinary.optimizer_converged_starts,
        optimizer_failed_starts=ordinary.optimizer_failed_starts,
        alpha_bounds=ordinary.alpha_bounds,
        root_variance_bounds=ordinary.root_variance_bounds,
    )


def compute_regime_ou_marginals(
    tree,
    values_by_leaf,
    regime_assignment,
    *,
    model="OUM",
    alpha=None,
    sigma2=None,
    theta=None,
    regime_parameters=None,
    regime_parameters_source=None,
    alpha_bounds=None,
    standard_errors=None,
    _tree_validated=False,
):
    """Fit and condition OUM, OUMA, OUMV, or OUMVA."""

    setup = _prepare_regime_ou_setup(
        tree,
        values_by_leaf,
        regime_assignment,
        model,
        alpha,
        sigma2,
        theta,
        regime_parameters,
        alpha_bounds,
        _tree_validated,
    )
    if len(setup.regimes) == 1:
        return _single_regime_ou_reduction(
            tree,
            values_by_leaf,
            regime_assignment,
            setup,
            regime_parameters_source=regime_parameters_source,
            standard_errors=standard_errors,
        )
    names, initial, parameter_bounds = _regime_ou_parameter_vectors(
        tree, values_by_leaf, regime_assignment, setup
    )
    if setup.num_observed <= len(names):
        raise ValueError(
            f"{setup.model} has {len(names)} free parameters but only "
            f"{setup.num_observed} "
            "observed tips; fix parameters or supply more data."
        )
    compiled = CompiledTree.from_tree(tree)

    def unpack(parameters):
        return _unpack_regime_ou_parameters(parameters, names, setup)

    initial_alpha, initial_sigma, _ = unpack(initial)
    _validate_regime_theta_design(
        tree,
        values_by_leaf,
        regime_assignment,
        setup,
        initial_alpha,
        initial_sigma,
    )

    def process_for(parameters):
        alpha_map, sigma_map, theta_map = unpack(parameters)
        return build_regime_ou_process(
            tree,
            regime_assignment,
            alpha_by_regime=alpha_map,
            sigma2_by_regime=sigma_map,
            theta_by_regime=theta_map,
        )

    observed_nodes = [
        leaf for leaf in tree.leaves() if values_by_leaf.get(str(leaf.name)) is not None
    ]
    if gaussian_process_parameter_rank(
        process_for, initial, observed_nodes, parameter_bounds
    ) < len(names):
        raise ValueError(
            f"{setup.model} free parameters are not jointly identifiable from "
            "the observed mean/covariance design; fix parameters or revise the "
            "regime map."
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
        optimized = deterministic_multistart(
            objective,
            initial,
            parameter_bounds,
            maxiter=1200,
        )
        process = process_for(optimized.x)
        conditioned = condition_gaussian_tree(
            process,
            values_by_leaf,
            standard_errors=standard_errors,
            compiled_tree=compiled,
        )
    else:
        fixed_parameters = np.asarray(initial, dtype=float)
        process = process_for(fixed_parameters)
        conditioned = condition_gaussian_tree(
            process,
            values_by_leaf,
            standard_errors=standard_errors,
            compiled_tree=compiled,
        )
        optimized = MultistartResult(
            fixed_parameters,
            -conditioned.log_likelihood,
            True,
            "all parameters fixed",
            0,
            0,
            0,
        )
    fitted_alpha, fitted_sigma, fitted_theta = unpack(optimized.x)
    _validate_regime_theta_design(
        tree,
        values_by_leaf,
        regime_assignment,
        setup,
        fitted_alpha,
        fitted_sigma,
    )
    shared_alpha = None if setup.varying_alpha else next(iter(fitted_alpha.values()))
    shared_sigma = None if setup.varying_sigma else next(iter(fitted_sigma.values()))
    fit = RegimeOUFit(
        model=setup.model,
        alpha=shared_alpha,
        alpha_by_regime=fitted_alpha,
        alpha_estimated=setup.fixed_alpha is None,
        theta_by_regime=fitted_theta,
        theta_estimated=setup.fixed_theta is None,
        sigma2=shared_sigma,
        sigma2_by_regime=fitted_sigma,
        sigma2_estimated=setup.fixed_sigma is None,
        root_variance=float(process.root.variance),
        log_likelihood=conditioned.log_likelihood,
        num_observed=conditioned.num_observed,
        num_effective_observations=conditioned.likelihood_rank,
        num_observed_positions=conditioned.num_observed_positions,
        regimes=setup.regimes,
        root_regime=regime_assignment.root_regime,
        regime_map_source=regime_assignment.source,
        regime_parameters_source=regime_parameters_source,
        fit_status=_regime_ou_boundary_status(setup, fitted_alpha, fitted_sigma),
        optimizer_success=optimized.success,
        optimizer_message=optimized.message,
        optimizer_starts=optimized.starts,
        optimizer_converged_starts=optimized.converged_starts,
        optimizer_failed_starts=optimized.failed_starts,
        alpha_bounds=setup.alpha_bounds,
        root_variance_bounds=(
            setup.root_variance_bounds
            if setup.model == "OUM" and setup.fixed_sigma is None
            else None
        ),
    )
    return conditioned.marginals, fit
