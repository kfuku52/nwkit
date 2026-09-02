"""Flat-root ASR for Brownian processes with transformed evolutionary time."""

import math
from dataclasses import dataclass

import numpy as np

from nwkit.continuous_asr import compute_bm_marginals
from nwkit.evolution import (
    build_evolutionary_process,
    encoded_evolution_parameter,
    optimization_parameterization,
    parameter_near_boundary,
    validate_evolution_parameter,
)
from nwkit.model_specs import evolution_model_spec
from nwkit.optimization import ScalarFitCache, global_bounded_scalar_minimize

SUPPORTED_TRANSFORMED_ASR_MODELS = ("lambda", "kappa", "delta", "eb", "acdc")


@dataclass(frozen=True)
class TransformedBrownianFit:
    evolution_model: str
    evolution_parameter_name: str
    evolution_parameter: float
    evolution_parameter_estimated: bool
    evolution_parameter_bounds: tuple[float, float]
    sigma2: float
    sigma2_estimated: bool
    restricted_log_likelihood: float | None
    num_observed: int
    num_effective_observations: int
    residual_df: int
    fit_status: str
    optimizer_success: bool
    optimizer_message: str
    optimizer_grid_evaluations: int
    profile_ci_level: float | None = None
    evolution_parameter_ci_lower: float | None = None
    evolution_parameter_ci_upper: float | None = None
    profile_ci_lower_boundary_limited: bool = False
    profile_ci_upper_boundary_limited: bool = False

    @property
    def eb_rate(self):
        return self.evolution_parameter

    @property
    def eb_rate_estimated(self):
        return self.evolution_parameter_estimated

    @property
    def eb_rate_bounds(self):
        return self.evolution_parameter_bounds


def default_parameter_bounds(tree, model):
    model = str(model).lower()
    if model not in SUPPORTED_TRANSFORMED_ASR_MODELS:
        raise ValueError(f"Unsupported transformed ASR model: {model}.")
    encoded_bounds, decode = optimization_parameterization(tree, model, allow_zero=True)
    values = float(decode(encoded_bounds[0])), float(decode(encoded_bounds[1]))
    return min(values), max(values)


def parse_parameter_bounds(value, tree, model):
    if value in (None, ""):
        return default_parameter_bounds(tree, model)
    items = [item.strip() for item in str(value).split(",")]
    if len(items) != 2:
        raise ValueError("'--evolution-parameter-bounds' must contain two values.")
    try:
        lower, upper = (float(item) for item in items)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            "'--evolution-parameter-bounds' values must be numeric."
        ) from exc
    if not all(math.isfinite(item) for item in (lower, upper)) or lower >= upper:
        raise ValueError(
            "'--evolution-parameter-bounds' must contain increasing finite values."
        )
    validate_evolution_parameter(model, lower)
    validate_evolution_parameter(model, upper)
    return lower, upper


def _process(tree, model, parameter):
    return build_evolutionary_process(
        tree,
        model=model,
        parameter=parameter,
        root_mode="flat",
        allow_zero=True,
    )


def compute_transformed_bm_marginals(
    tree,
    values_by_leaf,
    *,
    model,
    sigma2=None,
    evolution_parameter=None,
    evolution_parameter_bounds=None,
    profile_ci_level=None,
    standard_errors=None,
    _tree_validated=False,
):
    """Fit and condition a transformed-time Brownian process."""

    model = str(model).lower()
    if model not in SUPPORTED_TRANSFORMED_ASR_MODELS:
        raise ValueError(f"Unsupported transformed ASR model: {model}.")
    bounds = (
        default_parameter_bounds(tree, model)
        if evolution_parameter_bounds is None
        else (
            float(evolution_parameter_bounds[0]),
            float(evolution_parameter_bounds[1]),
        )
    )
    if not all(math.isfinite(value) for value in bounds) or bounds[0] >= bounds[1]:
        raise ValueError("Evolution-parameter bounds must be increasing and finite.")
    validate_evolution_parameter(model, bounds[0])
    validate_evolution_parameter(model, bounds[1])

    def fit_at(parameter):
        process = _process(tree, model, parameter)
        return compute_bm_marginals(
            tree,
            values_by_leaf,
            sigma2=sigma2,
            standard_errors=standard_errors,
            _tree_validated=_tree_validated,
            _process=process,
        )

    if evolution_parameter is not None:
        fitted_parameter = validate_evolution_parameter(model, evolution_parameter)
        assert fitted_parameter is not None
        posterior, brownian_fit = fit_at(fitted_parameter)
        parameter_estimated = False
        optimizer_success = True
        optimizer_message = "evolution parameter fixed"
        evaluations = 1
    else:
        if model == "delta":
            # Ultrametricity is a structural DELTA requirement, not an invalid
            # optimizer candidate.  Probe once so the useful model-level error
            # is not converted into an all-nonfinite search failure.
            _process(tree, model, 1.0)
        encoded_bounds = (
            encoded_evolution_parameter(model, bounds[0]),
            encoded_evolution_parameter(model, bounds[1]),
        )
        if encoded_bounds[0] > encoded_bounds[1]:
            encoded_bounds = encoded_bounds[::-1]
        _, decode = optimization_parameterization(tree, model, allow_zero=True)
        evaluated = {}

        def candidate(parameter):
            result = fit_at(parameter)
            likelihood = result[1].restricted_log_likelihood
            if likelihood is None or not math.isfinite(likelihood):
                raise ValueError("The transformed process has no ordinary likelihood.")
            return result

        cache = ScalarFitCache(
            candidate,
            lambda result: -float(result[1].restricted_log_likelihood),
            invalid_errors=(ValueError, ArithmeticError),
        )

        def objective(encoded):
            parameter = float(decode(float(encoded)))
            value = cache.objective(parameter)
            evaluated[parameter] = value
            return value

        optimized = global_bounded_scalar_minimize(
            objective, encoded_bounds, grid_size=25
        )
        if not math.isfinite(float(optimized.fun)):
            raise ValueError(f"Failed to find a finite {model.upper()} likelihood.")
        fitted_parameter = float(decode(float(optimized.x)))
        posterior, brownian_fit = fit_at(fitted_parameter)
        finite_objectives = np.asarray(
            [value for value in evaluated.values() if math.isfinite(value)], dtype=float
        )
        if len(finite_objectives) > 1 and np.ptp(finite_objectives) <= 1e-10 * max(
            1.0, abs(float(np.min(finite_objectives)))
        ):
            raise ValueError(
                f"{model.upper()} parameter is not identifiable on this tree and "
                "observation pattern; fix --evolution-parameter."
            )
        parameter_estimated = True
        optimizer_success = bool(optimized.success)
        optimizer_message = str(optimized.message)
        evaluations = len(evaluated)

    statuses = []
    if parameter_estimated and parameter_near_boundary(
        tree,
        model,
        fitted_parameter,
        allow_zero=True,
        parameter_bounds=bounds,
    ):
        statuses.append("evolution_parameter_boundary")
    if brownian_fit.fit_status != "ok":
        statuses.append(brownian_fit.fit_status)
    spec = evolution_model_spec(model)
    profile_interval = None
    if profile_ci_level is not None:
        if not parameter_estimated:
            raise ValueError(
                "Profile-likelihood intervals require an estimated evolution parameter."
            )
        level = float(profile_ci_level)
        if not 0.0 < level < 1.0:
            raise ValueError("Profile-likelihood level must be between zero and one.")
        from nwkit.asr_diagnostics import profile_likelihood_interval

        def profile(parameter):
            result = fit_at(float(parameter))[1].restricted_log_likelihood
            return -math.inf if result is None else float(result)

        profile_interval = profile_likelihood_interval(
            profile,
            fitted_parameter,
            bounds,
            level=level,
        )
    fit = TransformedBrownianFit(
        evolution_model=model,
        evolution_parameter_name=spec.parameter_name or "parameter",
        evolution_parameter=fitted_parameter,
        evolution_parameter_estimated=parameter_estimated,
        evolution_parameter_bounds=bounds,
        sigma2=brownian_fit.sigma2,
        sigma2_estimated=brownian_fit.sigma2_estimated,
        restricted_log_likelihood=brownian_fit.restricted_log_likelihood,
        num_observed=brownian_fit.num_observed,
        num_effective_observations=brownian_fit.num_effective_observations,
        residual_df=brownian_fit.residual_df,
        fit_status="+".join(statuses) if statuses else "ok",
        optimizer_success=optimizer_success,
        optimizer_message=optimizer_message,
        optimizer_grid_evaluations=evaluations,
        profile_ci_level=None if profile_ci_level is None else float(profile_ci_level),
        evolution_parameter_ci_lower=(
            None if profile_interval is None else profile_interval.lower
        ),
        evolution_parameter_ci_upper=(
            None if profile_interval is None else profile_interval.upper
        ),
        profile_ci_lower_boundary_limited=(
            False
            if profile_interval is None
            else profile_interval.lower_boundary_limited
        ),
        profile_ci_upper_boundary_limited=(
            False
            if profile_interval is None
            else profile_interval.upper_boundary_limited
        ),
    )
    return posterior, fit
