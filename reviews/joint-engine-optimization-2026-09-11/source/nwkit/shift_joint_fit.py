"""Joint ML covariance/mean fitting for native OU shift layouts."""

import math
from dataclasses import replace
from functools import cached_property, lru_cache

import numpy as np
from scipy.optimize import minimize

from nwkit.gaussian_whitening import TreeWhitening
from nwkit.optimization import global_bounded_scalar_minimize
from nwkit.shift_joint_dense import DenseJointContext, dense_eligible
from nwkit.shift_joint_model import evaluate_joint, evaluate_separable
from nwkit.shift_native_identifiability import covariance_identifiability
from nwkit.shift_native_model import (
    NativeTraitFit,
    covariance_geometry,
    restore_trait_fit,
)
from nwkit.vector_whitening import VectorObservationPlan


class JointFitContext:
    """Bounded reuse of shared-alpha tree factors within one search/data object."""

    def __init__(self, data, root_model):
        self.data, self.root_model = data, root_model
        self.factor = lru_cache(maxsize=32)(self._factor)

    @cached_property
    def plan(self):
        return VectorObservationPlan.build(
            self.data.tree.compiled, np.isfinite(self.data.values)
        )

    @cached_property
    def dense(self):
        return (
            DenseJointContext(self.data, self.root_model)
            if dense_eligible(self.data)
            else None
        )

    def _factor(self, alpha):
        slopes, innovations, root = covariance_geometry(
            self.data.tree, alpha, 1.0, self.root_model
        )
        return TreeWhitening.build(
            self.data.tree.compiled,
            self.data.tree.compiled.leaf_indices,
            slopes,
            innovations,
            np.zeros(len(self.data.values)),
            root_variance=root,
        )


def joint_residual_rank_excluded(data, layout, arguments):
    """A separable submodel with singular residual covariance has unbounded ML.

    Free trait-specific alpha contains the shared-alpha submodel. Fixed unequal
    alpha and fixed positive observation errors are not covered by this test.
    """
    options = arguments.get("options")
    if options is None or options.trait_covariance != "full":
        return False
    alpha = arguments.get("alpha_height")
    alpha = None if alpha is None else np.asarray(alpha)
    noise = arguments.get("measurement_variance")
    return bool(
        arguments.get("process_variance") is None
        and (noise is None or not np.any(noise))
        and not np.any(data.variances)
        and np.isfinite(data.values).all()
        and (alpha is None or np.all(alpha == alpha.flat[0]))
        and len(data.values) - len(layout.groups) < len(data.trait_names)
    )


def _correlation(parameters, p, diagonal):
    if diagonal:
        return np.eye(p)
    matrix = np.eye(p)
    matrix[np.tril_indices(p, -1)] = parameters
    matrix /= np.linalg.norm(matrix, axis=1)[:, None]
    return matrix @ matrix.T


def _correlation_parameters(covariance):
    scale = np.sqrt(np.diag(covariance))
    correlation = covariance / scale[:, None] / scale[None]
    factor = np.linalg.cholesky(correlation)
    factor /= np.diag(factor)[:, None]
    return factor[np.tril_indices(len(scale), -1)]


def _initial_covariance(data, layout, alpha):
    residual = np.zeros_like(data.values)
    for j in range(len(data.trait_names)):
        mask = np.isfinite(data.values[:, j])
        design = layout.design(data.tree, float(alpha[j]))[mask]
        beta = np.linalg.lstsq(design, data.values[mask, j], rcond=None)[0]
        residual[mask, j] = data.values[mask, j] - design @ beta
    covariance = residual.T @ residual / len(residual)
    # Regularization is used only to provide an interior starting point, never
    # in the objective or reported estimate.
    diagonal = np.maximum(np.diag(covariance), 0.01)
    return 0.8 * covariance + 0.2 * np.diag(diagonal)


def _richardson_gradient(function, x, bounds):
    """Independent finite differences with fourth-order interior accuracy.

    Variances can be small enough that a fixed central step's O(h^2) bias
    exceeds the stationarity tolerance. Keep the tolerance; cancel that bias.
    One-sided three-point stencils handle constrained boundary coordinates.
    """
    center = function(x)
    gradient = np.empty(len(x))
    for j, (lower, upper) in enumerate(bounds):
        h = 1e-5 * max(1.0, abs(x[j]))
        left_room = math.inf if lower is None else x[j] - lower
        right_room = math.inf if upper is None else upper - x[j]
        direction = 0
        if min(left_room, right_room) < h:
            direction = 1 if right_room >= left_room else -1
            h = min(h, max(left_room, right_room) / 2)

        def estimate(step, j=j, direction=direction):
            left, right = x.copy(), x.copy()
            if direction:
                left[j] += direction * step
                right[j] += direction * 2 * step
                return (
                    direction
                    * (-3 * center + 4 * function(left) - function(right))
                    / (2 * step)
                )
            left[j] -= step
            right[j] += step
            return (function(right) - function(left)) / (2 * step)

        gradient[j] = (4 * estimate(h / 2) - estimate(h)) / 3
    return gradient


def _projected_gradient(function, x, bounds, *, accurate=False):
    gradient = (
        _richardson_gradient(function, x, bounds) if accurate else np.empty(len(x))
    )
    for j, (lower, upper) in enumerate(bounds):
        if not accurate:
            h = 1e-5 * max(1.0, abs(x[j]))
            left, right = x.copy(), x.copy()
            left[j] = max(lower, x[j] - h) if lower is not None else x[j] - h
            right[j] = min(upper, x[j] + h) if upper is not None else x[j] + h
            gradient[j] = (function(right) - function(left)) / (right[j] - left[j])
        if lower is not None and x[j] <= lower + 1e-7 and gradient[j] > 0:
            gradient[j] = 0
        if upper is not None and x[j] >= upper - 1e-7 and gradient[j] < 0:
            gradient[j] = 0
    return float(np.max(np.abs(gradient))) if len(gradient) else 0.0


def _parameter_score(
    x,
    scores,
    p,
    alpha_count,
    process_count,
    correlation_count,
    noise_count,
    fixed_process,
):
    """Chain physical covariance/alpha scores through the optimizer coordinates."""
    alpha_score, covariance_score, noise_score = scores
    result = np.zeros_like(x)
    if alpha_count:
        result[:alpha_count] = (
            np.sum(alpha_score) if alpha_count == 1 else alpha_score
        ) * np.exp(x[:alpha_count])
    offset = alpha_count
    std = x[offset : offset + p] if process_count else np.sqrt(fixed_process)
    corr_offset = offset + process_count
    raw = np.eye(p)
    indices = np.tril_indices(p, -1)
    if correlation_count:
        raw[indices] = x[corr_offset : corr_offset + correlation_count]
    norms = np.linalg.norm(raw, axis=1)
    factor = raw / norms[:, None]
    corr = factor @ factor.T
    if process_count:
        result[offset : offset + p] = 2 * (covariance_score * corr) @ std
    if correlation_count:
        dfactor = 2 * (covariance_score * std[:, None] * std[None]) @ factor
        draw = (dfactor - np.sum(dfactor * factor, axis=1)[:, None] * factor) / norms[
            :, None
        ]
        result[corr_offset : corr_offset + correlation_count] = draw[indices]
    if noise_count:
        result[-noise_count:] = noise_score
    return result


def _general_fit(
    data, layout, options, fixed_alpha, fixed_process, fixed_noise, context
):
    p = len(data.trait_names)
    diagonal = options.trait_covariance == "diagonal"
    alpha_count = (
        0 if fixed_alpha is not None else (1 if options.alpha_model == "shared" else p)
    )
    correlation_count = 0 if diagonal else p * (p - 1) // 2
    process_count = p if fixed_process is None else 0
    noise_count = p if options.estimate_measurement_error else 0
    dense = context.dense if options.joint_engine == "auto" else None
    alpha_bounds = tuple(np.log(options.alpha_height_bounds))
    std_bound = math.sqrt(options.variance_bounds[1])
    bounds = (
        [alpha_bounds] * alpha_count
        + [(0.0, std_bound)] * process_count
        + [(None, None)] * correlation_count
        + [(0.0, options.variance_bounds[1] if dense is not None else std_bound)]
        * noise_count
    )
    evaluations, invalid = 0, 0

    def decode(x):
        offset = alpha_count
        alpha = (
            np.broadcast_to(np.exp(x[:offset]), (p,)).copy()
            if alpha_count
            else fixed_alpha
        )
        std = x[offset : offset + p] if process_count else np.sqrt(fixed_process)
        offset += process_count
        corr = _correlation(x[offset : offset + correlation_count], p, diagonal)
        offset += correlation_count
        noise = fixed_noise
        if noise_count:
            noise = x[offset : offset + p]
            if dense is None:
                noise = noise**2
        return alpha, std[:, None] * corr * std[None], noise

    def evaluate(x):
        nonlocal evaluations
        evaluations += 1
        alpha, covariance, noise = decode(x)
        if dense is not None:
            return dense.evaluate(layout, alpha, covariance, noise)
        return evaluate_joint(
            data,
            layout,
            alpha,
            covariance,
            noise,
            root_model=options.root_model,
            _plan=context.plan,
        )

    def objective(x):
        nonlocal invalid
        try:
            return -evaluate(x).log_likelihood
        except (ValueError, ArithmeticError, np.linalg.LinAlgError):
            invalid += 1
            return math.inf

    def objective_score(x):
        nonlocal evaluations, invalid
        assert dense is not None
        evaluations += 1
        try:
            alpha, covariance, noise = decode(x)
            fit, scores = dense.evaluate(
                layout, alpha, covariance, noise, gradient=True
            )
            gradient = _parameter_score(
                x,
                scores,
                p,
                alpha_count,
                process_count,
                correlation_count,
                noise_count,
                fixed_process,
            )
            if not np.isfinite(gradient).all():
                raise ValueError("Nonfinite analytic joint likelihood score.")
            return -fit.log_likelihood, gradient
        except (ValueError, ArithmeticError, np.linalg.LinAlgError):
            invalid += 1
            return math.inf, np.zeros_like(x)

    function = objective_score if dense is not None else objective
    starts, accepted = [], []
    for index in range(options.optimizer_starts):
        seed_alpha = (
            np.full(p, (0.3, 1.0, 3.0)[index % 3])
            if fixed_alpha is None
            else fixed_alpha
        )
        initial = _initial_covariance(data, layout, seed_alpha)
        noise_fraction = (0.0, 0.25, 0.65)[index % 3] if noise_count else 0.0
        x = []
        if alpha_count:
            x.extend(np.log(seed_alpha[:alpha_count]))
        if process_count:
            x.extend(np.sqrt(np.diag(initial) * (1 - noise_fraction)))
        if correlation_count:
            x.extend(
                _correlation_parameters(initial)
                if index % 2 == 0
                else np.zeros(correlation_count)
            )
        if noise_count:
            initial_noise = np.diag(initial) * noise_fraction
            # Squared standard deviations have a zero derivative at zero even
            # when increasing variance improves the objective. Direct variance
            # coordinates preserve the one-sided score for analytic gradients.
            x.extend(initial_noise if dense is not None else np.sqrt(initial_noise))
        x = np.asarray(x)
        if len(x):
            result = minimize(
                function,
                x,
                jac=True if dense is not None else None,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": options.maxiter, "ftol": 1e-12, "gtol": 1e-6},
            )
            score = (
                _projected_gradient(
                    objective, result.x, bounds, accurate=dense is not None
                )
                if np.isfinite(result.fun)
                else math.inf
            )
            if result.success and score > 1e-3:
                result = minimize(
                    function,
                    result.x,
                    method="L-BFGS-B",
                    jac=True if dense is not None else "3-point",
                    bounds=bounds,
                    options={"maxiter": options.maxiter, "ftol": 1e-14, "gtol": 1e-7},
                )
                score = (
                    _projected_gradient(
                        objective, result.x, bounds, accurate=dense is not None
                    )
                    if np.isfinite(result.fun)
                    else math.inf
                )
            success = bool(result.success and np.isfinite(result.fun) and score <= 1e-3)
            message = str(result.message)
            optimum = result.x
        else:
            score, success, message, optimum = (
                0.0,
                math.isfinite(objective(x)),
                "Fixed covariance",
                x,
            )
        starts.append(
            {
                "success": success,
                "message": message,
                "projected_score": score if math.isfinite(score) else None,
            }
        )
        if success:
            accepted.append(evaluate(optimum))
    if not accepted:
        raise ValueError(
            "Joint covariance optimizer failed checked convergence: " + str(starts)
        )
    best = max(accepted, key=lambda fit: fit.log_likelihood)
    return best, {
        "evaluations": evaluations,
        "invalid_evaluations": invalid,
        "starts": starts,
        "success": True,
        "message": "Checked multistart joint ML",
    }


def _profile_fit(data, layout, options, fixed_alpha, context):
    scores = {}
    best = None
    invalid = 0

    def objective(logalpha):
        nonlocal best, invalid
        alpha = math.exp(logalpha)
        if alpha not in scores:
            try:
                fit = evaluate_separable(
                    data,
                    layout,
                    alpha,
                    root_model=options.root_model,
                    diagonal=options.trait_covariance == "diagonal",
                    factor=context.factor(alpha),
                )
                scores[alpha] = -fit.log_likelihood
                if best is None or fit.log_likelihood > best.log_likelihood:
                    best = fit
            except (ValueError, ArithmeticError, np.linalg.LinAlgError):
                invalid += 1
                scores[alpha] = math.inf
        return scores[alpha]

    if fixed_alpha is not None:
        alpha = float(fixed_alpha[0])
        fit = evaluate_separable(
            data,
            layout,
            alpha,
            root_model=options.root_model,
            diagonal=options.trait_covariance == "diagonal",
            factor=context.factor(alpha),
        )
        return fit, {
            "evaluations": 1,
            "invalid_evaluations": 0,
            "success": True,
            "message": "Exact covariance profile at fixed alpha",
        }
    result = global_bounded_scalar_minimize(
        objective, np.log(options.alpha_height_bounds)
    )
    if not result.success or best is None:
        raise ValueError("Joint shared-alpha profile failed: " + str(result.message))
    return best, {
        "evaluations": len(scores),
        "invalid_evaluations": invalid,
        "success": True,
        "message": result.message,
    }


def fit_joint_layout(
    data,
    layout,
    *,
    options,
    alpha_height=None,
    process_variance=None,
    measurement_variance=None,
    _joint_context=None,
):
    """Fit full covariance, or a shared-alpha diagonal covariance model.

    Supplied variances use original units. All finite alpha coordinates and the
    common Brownian/independent limits are checked. Mixed finite/infinite alpha
    limits are not part of the joint model. Mixed finite/zero boundary faces in
    estimated trait-specific alpha are not exhaustively enumerated.
    """
    options.validate()
    if joint_residual_rank_excluded(
        data,
        layout,
        {
            "options": options,
            "alpha_height": alpha_height,
            "process_variance": process_variance,
            "measurement_variance": measurement_variance,
        },
    ):
        raise ValueError(
            "Insufficient residual rank for a finite full covariance ML estimate."
        )
    p = len(data.trait_names)
    if options.estimate_measurement_error and measurement_variance is not None:
        raise ValueError(
            "Additional measurement variance cannot be fixed and estimated together."
        )

    def vector(value, label, variance=False):
        if value is None:
            return None
        array = np.broadcast_to(np.asarray(value, float), (p,)).copy()
        if (
            np.isnan(array).any()
            or np.any(array < 0)
            or (variance and not np.isfinite(array).all())
        ):
            raise ValueError(f"{label} requires nonnegative values.")
        if variance:
            original = array.copy()
            array = array / data.scales / data.scales
            if not np.isfinite(array).all() or np.any((original > 0) & (array == 0)):
                raise ValueError(
                    "Fixed variance is not representable in normalized units."
                )
        return array

    alpha = vector(alpha_height, "Alpha")
    process = vector(process_variance, "Process variance", True)
    noise = vector(measurement_variance, "Measurement variance", True)
    noise = np.zeros(p) if noise is None else noise
    if (
        alpha is not None
        and options.alpha_model == "shared"
        and np.any(alpha != alpha[0])
    ):
        raise ValueError("Shared alpha requires the same fixed value for all traits.")
    if alpha is not None and np.isinf(alpha).any() and not np.isinf(alpha).all():
        raise ValueError(
            "Joint covariance does not define mixed finite/infinite alpha limits."
        )
    context = _joint_context or JointFitContext(data, options.root_model)
    if context.data is not data or context.root_model != options.root_model:
        raise ValueError(
            "Joint fit cache belongs to a different dataset or root model."
        )
    joint, records = _fit_joint_modes(
        data, layout, options, alpha, process, noise, context
    )
    return _joint_result(data, layout, options, alpha, process, joint, records)


def _fit_joint_modes(data, layout, options, alpha, process, noise, context):
    p = len(data.trait_names)
    alphas = [alpha]
    if alpha is None:
        if options.root_model == "OUfixedRoot":
            alphas.append(np.zeros(p))
        alphas.append(np.full(p, math.inf))
    candidates, records = [], []
    for candidate_alpha in alphas:
        mode = "finite" if candidate_alpha is None else str(float(candidate_alpha[0]))
        local_options = options
        # At the common independent limit a free process covariance and extra
        # diagonal noise cannot be separated. Use the canonical all-process
        # allocation and report the ambiguity, rather than optimize a ridge.
        if (
            candidate_alpha is not None
            and np.isinf(candidate_alpha).all()
            and process is None
        ):
            local_options = replace(options, estimate_measurement_error=False)
        shared = (
            options.alpha_model == "shared"
            if candidate_alpha is None
            else np.all(candidate_alpha == candidate_alpha[0])
        )
        separable = (
            options.joint_engine == "auto"
            and shared
            and process is None
            and not local_options.estimate_measurement_error
            and not np.any(noise)
            and not np.any(data.variances)
            and np.isfinite(data.values).all()
        )
        try:
            if candidate_alpha is not None:
                for j in range(p):
                    matrix = layout.design(data.tree, float(candidate_alpha[j]))[
                        np.isfinite(data.values[:, j])
                    ]
                    norms = np.linalg.norm(matrix, axis=0)
                    if (
                        np.any(norms == 0)
                        or np.linalg.matrix_rank(matrix / norms) < matrix.shape[1]
                    ):
                        raise ValueError("The observed mean design is rank deficient.")
            if separable:
                fit, diagnostic = _profile_fit(
                    data, layout, local_options, candidate_alpha, context
                )
            else:
                fit, diagnostic = _general_fit(
                    data,
                    layout,
                    local_options,
                    candidate_alpha,
                    process,
                    noise,
                    context,
                )
            records.append(
                {
                    "alpha_mode": mode,
                    "success": True,
                    "log_likelihood": fit.log_likelihood,
                    "optimizer": diagnostic,
                }
            )
            candidates.append(fit)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            records.append({"alpha_mode": mode, "success": False, "message": str(exc)})
    if not candidates:
        raise ValueError("Joint native fixed-layout fitting failed: " + str(records))
    maximum = max(fit.log_likelihood for fit in candidates)
    tied = [fit for fit in candidates if fit.log_likelihood >= maximum - 1e-9]
    joint = min(
        tied,
        key=lambda fit: (
            np.all(np.isfinite(fit.alpha_height)) and np.any(fit.alpha_height > 0),
            float(np.sum(fit.alpha_height)),
        ),
    )
    return joint, records


def _joint_result(data, layout, options, alpha, process, joint, records):
    p = len(data.trait_names)
    evaluations = sum(r.get("optimizer", {}).get("evaluations", 0) for r in records)
    variance_boundary = bool(
        (
            process is None
            and np.any(
                np.diag(joint.covariance_coordinate)
                >= options.variance_bounds[1] * (1 - 1e-5)
            )
        )
        or (
            options.estimate_measurement_error
            and np.any(
                joint.measurement_variance >= options.variance_bounds[1] * (1 - 1e-5)
            )
        )
    )
    ambiguity = bool(
        options.estimate_measurement_error
        and process is None
        and np.isinf(joint.alpha_height).all()
    )
    diagnostics = {
        "evaluations": evaluations,
        "invalid_evaluations": sum(
            r.get("optimizer", {}).get("invalid_evaluations", 0) for r in records
        ),
        "alpha_candidates": records,
        "continuous_globally_optimal": False,
        "complete_alpha_modes": all(r["success"] for r in records)
        and not (alpha is None and options.alpha_model == "trait-specific" and p > 1),
        "evaluated_alpha_modes_succeeded": all(r["success"] for r in records),
        "alpha_boundary_scope": "shared_limits_only; mixed_trait_specific_boundary_faces_not_enumerated",
        "alpha_height_bounds": list(options.alpha_height_bounds),
        "normalized_variance_bounds": [0.0, options.variance_bounds[1]],
        "alpha_estimated": alpha is None,
        "process_variance_estimated": process is None,
        "measurement_variance_estimated": options.estimate_measurement_error,
        "nuisance_variance_at_numerical_bound": variance_boundary,
        "equivalent_modes_disagree_on_variance_decomposition": ambiguity,
        "likelihood_scope": "joint_across_traits",
    }
    fits, traits = [], []
    q = len(layout.groups)
    for j in range(p):
        fit = NativeTraitFit(
            float(joint.alpha_height[j]),
            float(joint.process_tip_covariance[j, j]),
            float(joint.measurement_variance[j]),
            joint.coefficients[:, j].copy(),
            joint.coefficient_covariance[
                j * q : (j + 1) * q, j * q : (j + 1) * q
            ].copy(),
            joint.predicted[:, j].copy(),
            joint.log_likelihood,
            joint.quadratic,
            int(np.isfinite(data.values[:, j]).sum()),
            q,
            options.root_model,
        )
        record = restore_trait_fit(data, j, fit)
        record["log_likelihood"] = None
        record["likelihood_scope"] = "joint_across_traits; no additive_trait_likelihood"
        record["optimizer"] = dict(
            diagnostics, evaluations=evaluations if j == 0 else 0
        )
        identification = covariance_identifiability(
            data,
            j,
            fit,
            alpha_free=alpha is None,
            process_free=process is None,
            noise_free=options.estimate_measurement_error,
        )
        identification["scope"] = (
            "conservative_marginal_diagnostic; joint_covariance_identifiability_not_certified"
        )
        identification["between_mode_variance_ambiguity"] = ambiguity
        if ambiguity:
            identification["variance_decomposition_supported"] = False
        record["identifiability_diagnostic"] = identification
        record["covariance_components_identifiable"] = identification[
            "variance_decomposition_supported"
        ]
        fits.append(fit)
        traits.append(record)
    covariance_parameters = (
        p * (p + 1) // 2 if options.trait_covariance == "full" else p
    )
    if process is not None:
        covariance_parameters -= p
    alpha_parameters = (
        0 if alpha is not None else (1 if options.alpha_model == "shared" else p)
    )
    scale = data.scales[:, None] * data.scales[None]
    joint_metadata = {
        "alpha_model": options.alpha_model,
        "process_tip_covariance": (joint.process_tip_covariance * scale).tolist(),
        "diffusion_covariance": None
        if joint.diffusion_covariance is None
        else (joint.diffusion_covariance * scale / data.tree.height).tolist(),
        "measurement_covariance": np.diag(
            joint.measurement_variance * data.scales**2
        ).tolist(),
        "coefficient_order": "trait_major_then_scaled_regime_coefficients",
        "coefficient_covariance": (
            joint.coefficient_covariance
            * np.repeat(data.scales, q)[:, None]
            * np.repeat(data.scales, q)[None]
        ).tolist(),
        "covariance_parameter_count": covariance_parameters,
        "alpha_parameter_count": alpha_parameters,
        "measurement_parameter_count": p if options.estimate_measurement_error else 0,
        "engine": joint.engine,
        "optimizer": diagnostics,
        "identifiability": "marginal_diagnostics_only; full_joint_rank_not_certified",
        "independent_limit_variance_allocation": "all_process_nonidentifiable"
        if ambiguity
        else None,
    }
    adjustment = sum(
        int(np.isfinite(data.values[:, j]).sum()) * math.log(data.scales[j])
        for j in range(p)
    )
    return {
        "layout": layout,
        "fits": fits,
        "traits": traits,
        "joint_fit": joint,
        "joint_covariance": joint_metadata,
        "log_likelihood": joint.log_likelihood - adjustment,
        "root_model": options.root_model,
        "trait_covariance": options.trait_covariance,
    }
