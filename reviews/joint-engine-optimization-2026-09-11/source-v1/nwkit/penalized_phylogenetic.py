"""Exploratory elastic-net regression with fixed phylogenetic covariance shape.

Gaussian loss is GLS squared error / (2*n). Non-Gaussian loss is the
Laplace negative log likelihood / n. The latter reuses the GLMM mode solvers.
This module deliberately performs no coefficient hypothesis tests.
"""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import minimize
from scipy.special import expit

from nwkit.phylogenetic_glmm import _random_mode, _scalar_random_mode

FAMILIES = ("gaussian", "binomial", "poisson", "negative-binomial")


@dataclass
class SelectionFit:
    coefficients: np.ndarray  # original predictor units, intercept first
    standardized_coefficients: np.ndarray
    center: np.ndarray
    scale: np.ndarray
    active: np.ndarray
    family: str
    strength: float
    l1_ratio: float
    covariance: np.ndarray
    random_mode: np.ndarray
    random_variance: float
    dispersion: float | None
    objective: float
    projected_gradient: float
    boundary_warning: str
    iterations: int

    def predict(self, x, cross_covariance=None):
        """Plug-in mean; optionally condition the latent effect on training tips.

        Neither form integrates latent uncertainty. The caller records which
        estimand it evaluates; held-out responses are never used here.
        """
        x = np.asarray(x, dtype=float)
        if x.ndim != 2 or x.shape[1] != len(self.center) or not np.isfinite(x).all():
            raise ValueError(
                "Prediction design must be finite and match fitted columns."
            )
        # Use the fitted coordinates to avoid cancellation from a large origin.
        linear = (
            self.standardized_coefficients[0]
            + ((x[:, self.active] - self.center[self.active]) / self.scale[self.active])
            @ self.standardized_coefficients[1:][self.active]
        )
        if cross_covariance is not None:
            linear = linear + np.asarray(cross_covariance) @ np.linalg.solve(
                self.covariance, self.random_mode
            )
        if self.family == "gaussian":
            return linear
        if self.family == "binomial":
            return expit(linear)
        with np.errstate(over="raise"):
            return np.exp(linear)


def validate_inputs(x, y, covariance, family):
    x, y, covariance = map(lambda a: np.asarray(a, dtype=float), (x, y, covariance))
    if family not in FAMILIES:
        raise ValueError(f"Unsupported selection family: {family}.")
    if x.ndim != 2 or y.ndim != 1 or len(y) != len(x) or len(y) < 4:
        raise ValueError(
            "Selection needs a 2D design, 1D response, and at least four training tips."
        )
    if not all(np.isfinite(a).all() for a in (x, y, covariance)):
        raise ValueError(
            "Selection inputs must be finite; missing values require explicit upstream handling."
        )
    if covariance.shape != (len(y), len(y)) or not np.allclose(
        covariance, covariance.T, rtol=1e-10, atol=0
    ):
        raise ValueError("Selection covariance must be symmetric and tip-matched.")
    np.linalg.cholesky(covariance)
    if np.unique(y).size < 2:
        raise ValueError("Every training fold must contain a variable response.")
    if family == "binomial" and not np.isin(y, [0, 1]).all():
        raise ValueError("Binomial selection requires numeric 0/1 responses.")
    if family in ("poisson", "negative-binomial") and (
        (y < 0).any() or (y != np.floor(y)).any()
    ):
        raise ValueError("Count selection requires non-negative integer responses.")
    return x, y, covariance


class LaplaceLoss:
    """Marginal loss and analytic fixed-linear gradient, using GLMM solvers."""

    def __init__(self, y, covariance, family):
        self.y = y
        self.family = family
        self.inverse = cho_solve(cho_factor(covariance), np.eye(len(y)))
        self.logdet = np.linalg.slogdet(covariance)[1]

    def state(self, linear, nuisance):
        variance = np.exp(nuisance[0])
        dispersion = np.exp(nuisance[1]) if self.family == "negative-binomial" else None
        precision = self.inverse / variance
        n = len(self.y)
        if self.family == "binomial":
            counts = np.column_stack([self.y, 1 - self.y])
            mode, joint, weights = _random_mode(
                counts,
                linear[:, None],
                precision,
                family="binomial",
                thresholds=np.empty(0),
            )
            weights = np.asarray(weights).reshape(n)
            mu = expit(linear + mode)
            gradient = mu - self.y
            third = weights * (1 - 2 * mu)
        else:
            mode, joint, weights = _scalar_random_mode(
                self.y,
                linear,
                precision,
                np.arange(n),
                family=self.family,
                dispersion=dispersion,
                zero_probability=None,
                offset=np.zeros(n),
                trials=None,
                censor_lower=None,
                censor_upper=None,
            )
            if self.family == "poisson":
                mu = np.exp(linear + mode)
                gradient = mu - self.y
                third = weights
            else:
                assert dispersion is not None
                size = 1 / dispersion
                fraction = expit(linear + mode - np.log(size))
                gradient = (size + self.y) * fraction - self.y
                third = weights * (1 - 2 * fraction)
        hessian = precision + np.diag(weights)
        factor = cho_factor(hessian)
        inverse = cho_solve(factor, np.eye(n))
        correction = np.diag(inverse) * third
        # Implicit derivative of the random mode: du/deta = -H^-1 W.
        eta_gradient = gradient + 0.5 * (correction - weights * (inverse @ correction))
        logdet_hessian = 2 * np.log(np.diag(factor[0])).sum()
        value = (joint + 0.5 * (self.logdet + n * nuisance[0] + logdet_hessian)) / n
        return float(value), eta_gradient / n, mode

    def value_gradient(self, linear, nuisance):
        value, gradient, _mode = self.state(linear, nuisance)
        nuisance_gradient = np.empty(len(nuisance))
        # Only one/two variance parameters use finite differences, never p columns.
        for j in range(len(nuisance)):
            step = np.zeros(len(nuisance))
            step[j] = 1e-4
            nuisance_gradient[j] = (
                self.state(linear, nuisance + step)[0]
                - self.state(linear, nuisance - step)[0]
            ) / 2e-4
        return value, gradient, nuisance_gradient


def _projected_gradient(parameters, gradient, bounds):
    out = gradient.copy()
    for i, (lo, hi) in enumerate(bounds):
        if lo is not None and parameters[i] <= lo + 1e-9:
            out[i] = min(0, out[i])
        if hi is not None and parameters[i] >= hi - 1e-9:
            out[i] = max(0, out[i])
    return float(np.max(np.abs(out)))


def _minimize_checked(trial_objective, objective, parameters, bounds, maxiter):
    result = minimize(
        trial_objective,
        parameters,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-6, "maxls": 40},
    )
    value, gradient = objective(result.x)
    kkt = _projected_gradient(result.x, gradient, bounds)
    iterations = int(result.nit)
    # L-BFGS can stop on negligible objective change while its stored curvature
    # still gives a poor step (especially near a variance boundary). Restart
    # from that state with fresh curvature; keep the same objective, bounds,
    # total iteration budget, and independent stationarity requirement.
    for _ in range(2):
        if (
            not result.success
            or not np.isfinite(kkt)
            or kkt <= 2e-4
            or iterations >= maxiter
        ):
            break
        retry = minimize(
            trial_objective,
            result.x,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={
                "maxiter": maxiter - iterations,
                "ftol": 1e-12,
                "gtol": 1e-6,
                "maxls": 40,
            },
        )
        retry_value, retry_gradient = objective(retry.x)
        iterations += int(retry.nit)
        if not np.isfinite(retry_value) or retry_value > value:
            break
        result, value, gradient = retry, retry_value, retry_gradient
        kkt = _projected_gradient(result.x, gradient, bounds)
    if (
        not result.success
        or not np.isfinite(value)
        or not np.isfinite(kkt)
        or kkt > 2e-4
    ):
        raise RuntimeError(
            f"Elastic-net optimizer did not converge (projected gradient={kkt:g}): {result.message}"
        )
    return result, value, kkt, iterations


def _gaussian_coordinates(
    design, y, covariance, free, penalized, strength, ratio, beta, maxiter
):
    """Solve the convex GLS elastic net without line-search stopping criteria.

    Project out the unpenalized space after whitening, then minimize each
    penalized coordinate exactly. Storage is O(n*p), including for p > n.
    """
    lower = np.linalg.cholesky(covariance)
    whitened = solve_triangular(lower, design, lower=True) / np.sqrt(len(y))
    response = solve_triangular(lower, y, lower=True) / np.sqrt(len(y))
    q, r = np.linalg.qr(whitened[:, free], mode="reduced")
    columns = whitened[:, penalized]
    projected = columns - q @ (q.T @ columns)
    target = response - q @ (q.T @ response)
    diagonal = np.sum(projected * projected, axis=0)
    lasso, ridge = strength * ratio, strength * (1 - ratio)
    values = beta[penalized].copy()
    residual = target - projected @ values
    kkt = np.inf
    for iteration in range(1, maxiter + 1):  # noqa: B007 - returned as iteration count
        for j in range(len(values)):
            score = projected[:, j] @ residual + diagonal[j] * values[j]
            numerator = np.sign(score) * max(abs(score) - lasso, 0)
            updated = numerator / (diagonal[j] + ridge) if numerator else 0.0
            residual -= projected[:, j] * (updated - values[j])
            values[j] = updated
        # Recompute residuals to avoid accumulated coordinate-update roundoff.
        residual = target - projected @ values
        gradient = -projected.T @ residual + ridge * values
        violation = np.where(
            values != 0,
            np.abs(gradient + lasso * np.sign(values)),
            np.maximum(np.abs(gradient) - lasso, 0),
        )
        kkt = float(np.max(violation, initial=0))
        if kkt <= 5e-7:
            break
    beta[penalized] = values
    beta[free] = solve_triangular(r, q.T @ (response - columns @ values))
    residual = whitened @ beta - response
    gradient = whitened.T @ residual
    gradient[penalized] += ridge * values
    violation = np.where(
        values != 0,
        np.abs(gradient[penalized] + lasso * np.sign(values)),
        np.maximum(np.abs(gradient[penalized]) - lasso, 0),
    )
    kkt = max(
        float(np.max(violation, initial=0)), float(np.max(np.abs(gradient[free])))
    )
    value = float(
        0.5 * residual @ residual
        + lasso * np.abs(values).sum()
        + 0.5 * ridge * (values @ values)
    )
    if not np.isfinite(value) or not np.isfinite(kkt) or kkt > 1e-6:
        raise RuntimeError(
            f"Gaussian elastic-net coordinates did not converge (projected gradient={kkt:g})."
        )
    return beta, value, kkt, iteration


def _validate_penalty_options(
    maxiter, family, strength, l1_ratio, unpenalized, columns
):
    if maxiter is None:
        maxiter = 20000 if family == "gaussian" else 2000
    if not isinstance(maxiter, (int, np.integer)) or maxiter < 1:
        raise ValueError("maxiter must be a positive integer.")
    if (
        not np.isfinite(strength)
        or strength <= 0
        or not np.isfinite(l1_ratio)
        or not 0 < l1_ratio <= 1
    ):
        raise ValueError("strength must be positive; l1_ratio must be in (0, 1].")
    unpenalized = tuple(unpenalized)
    if len(set(unpenalized)) != len(unpenalized) or any(
        i < 0 or i >= columns for i in unpenalized
    ):
        raise ValueError("Invalid or duplicate unpenalized column index.")
    return maxiter, unpenalized


def fit_elastic_net(
    x,
    y,
    covariance,
    *,
    family="gaussian",
    strength=0.1,
    l1_ratio=0.5,
    unpenalized=(),
    initial=None,
    maxiter=None,
):
    """Fit one penalty point; all preprocessing is learned from these rows only."""
    x, y, covariance = validate_inputs(x, y, covariance, family)
    maxiter, unpenalized = _validate_penalty_options(
        maxiter, family, strength, l1_ratio, unpenalized, x.shape[1]
    )
    # Scale deviations before squaring: std(x) under/overflows for valid units.
    # Exact constants alone are excluded; an absolute cutoff changes selection
    # when the same predictor is expressed in different units.
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        center = (x / len(x)).sum(axis=0)
        deviations = x - center
        magnitude = np.max(np.abs(deviations), axis=0)
        active = np.any(x != x[0], axis=0)
        safe_magnitude = np.where(active, magnitude, 1)
        scale = safe_magnitude * np.sqrt(
            np.mean((deviations / safe_magnitude) ** 2, axis=0)
        )
    if any(not active[i] for i in unpenalized):
        raise ValueError("An unpenalized predictor is constant in a training fold.")
    scale = np.where(active, scale, 1)
    standardized = deviations / scale
    standardized[:, ~active] = 0
    design = np.column_stack([np.ones(len(y)), standardized])
    free = np.array([0] + [i + 1 for i in unpenalized])
    penalized = np.array(
        [i + 1 for i in range(x.shape[1]) if active[i] and i not in unpenalized],
        dtype=int,
    )
    if len(free) >= len(y) or np.linalg.matrix_rank(design[:, free]) != len(free):
        raise ValueError(
            "Unpenalized design must have full rank and fewer columns than training tips."
        )
    f, p = len(free), len(penalized)
    nuisance_size = 0 if family == "gaussian" else 1 + (family == "negative-binomial")
    parameters = np.zeros(f + 2 * p + nuisance_size)
    parameters[0] = (
        np.mean(y)
        if family == "gaussian"
        else (
            np.log(np.mean(y) / (1 - np.mean(y)))
            if family == "binomial"
            else np.log(np.mean(y))
        )
    )
    if nuisance_size:
        parameters[-nuisance_size:] = -1
    if initial is not None:
        beta = initial.standardized_coefficients
        parameters[:f] = beta[free]
        parameters[f : f + p] = np.maximum(beta[penalized], 0)
        parameters[f + p : f + 2 * p] = np.maximum(-beta[penalized], 0)
        if nuisance_size:
            parameters[-nuisance_size] = np.log(initial.random_variance)
        if family == "negative-binomial":
            parameters[-1] = np.log(initial.dispersion)
    bounds: list[tuple[float | None, float | None]] = [(None, None)] * f + [
        (0, None)
    ] * (2 * p)
    bounds += [(-12.0, 6.0)] * nuisance_size
    gaussian_factor = cho_factor(covariance) if family == "gaussian" else None
    laplace = None if family == "gaussian" else LaplaceLoss(y, covariance, family)

    def decode(theta):
        beta = np.zeros(design.shape[1])
        beta[free] = theta[:f]
        beta[penalized] = theta[f : f + p] - theta[f + p : f + 2 * p]
        return beta

    def objective(theta):
        beta = decode(theta)
        linear = design @ beta
        if family == "gaussian":
            residual = linear - y
            eta_gradient = cho_solve(gaussian_factor, residual) / len(y)
            value = 0.5 * residual @ eta_gradient
            nuisance_gradient = np.empty(0)
        else:
            assert laplace is not None
            value, eta_gradient, nuisance_gradient = laplace.value_gradient(
                linear, theta[-nuisance_size:]
            )
        beta_gradient = design.T @ eta_gradient
        ridge = strength * (1 - l1_ratio)
        lasso = strength * l1_ratio
        value += lasso * np.sum(theta[f : f + 2 * p]) + 0.5 * ridge * np.sum(
            beta[penalized] ** 2
        )
        penalized_gradient = beta_gradient[penalized] + ridge * beta[penalized]
        gradient = np.concatenate(
            [
                beta_gradient[free],
                penalized_gradient + lasso,
                -penalized_gradient + lasso,
                nuisance_gradient,
            ]
        )
        return float(value), gradient

    def trial_objective(theta):
        # The line search may probe a state whose latent mode cannot be solved.
        # Reject that trial, then require a valid final state and KKT residual.
        try:
            return objective(theta)
        except (RuntimeError, FloatingPointError, np.linalg.LinAlgError):
            return 1e100, np.zeros_like(theta)

    if family == "gaussian":
        # Quasi-Newton steps locate the solution quickly; exact coordinate
        # minimization supplies the independent stationarity refinement rather
        # than trusting a machine-dependent relative-objective stopping flag.
        seed_budget = maxiter // 2
        seed_iterations = 0
        start_beta = decode(parameters)
        if seed_budget:
            seed = minimize(
                objective,
                parameters,
                method="L-BFGS-B",
                jac=True,
                bounds=bounds,
                options={
                    "maxiter": seed_budget,
                    "ftol": 1e-12,
                    "gtol": 1e-6,
                    "maxls": 40,
                },
            )
            seed_iterations = int(seed.nit)
            if np.isfinite(seed.fun) and np.isfinite(seed.x).all():
                start_beta = decode(seed.x)
        beta, value, kkt, iterations = _gaussian_coordinates(
            design,
            y,
            covariance,
            free,
            penalized,
            strength,
            l1_ratio,
            start_beta,
            maxiter - seed_iterations,
        )
        iterations += seed_iterations
    else:
        result, value, kkt, iterations = _minimize_checked(
            trial_objective, objective, parameters, bounds, maxiter
        )
        beta = decode(result.x)
    coefficients = beta.copy()
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        coefficients[1:] /= scale
        coefficients[0] -= center @ coefficients[1:]
    warning = ""
    if family == "gaussian":
        mode = y - design @ beta
        variance = float(mode @ cho_solve(gaussian_factor, mode) / len(y))
        dispersion = None
        if l1_ratio == 1 and len(penalized):
            score = design[:, penalized].T @ cho_solve(gaussian_factor, -mode) / len(y)
            equicorrelated = (
                np.abs(np.abs(score) - strength) <= 1e-5 * max(1, strength)
            ) | (beta[penalized] != 0)
            indices = np.r_[free, penalized[equicorrelated]]
            if np.linalg.matrix_rank(design[:, indices]) < len(indices):
                warning = "lasso_equicorrelation_design_rank_deficient"
    else:
        nuisance = result.x[-nuisance_size:]
        assert laplace is not None
        _, _, mode = laplace.state(design @ beta, nuisance)
        variance = float(np.exp(nuisance[0]))
        dispersion = (
            float(np.exp(nuisance[1])) if family == "negative-binomial" else None
        )
        if np.any(nuisance < -11.99) or np.any(nuisance > 5.99):
            warning = "variance_or_dispersion_at_optimization_boundary"
    return SelectionFit(
        coefficients,
        beta,
        center,
        scale,
        active,
        family,
        strength,
        l1_ratio,
        covariance,
        mode,
        variance,
        dispersion,
        value,
        kkt,
        warning,
        iterations,
    )


def prediction_loss(y, prediction, family):
    """Per-tip squared error (Gaussian/count) or log loss (binary)."""
    y, prediction = np.asarray(y), np.asarray(prediction)
    if not np.isfinite(prediction).all():
        raise ValueError("Non-finite held-out prediction.")
    if family == "binomial":
        prob = np.clip(prediction, np.finfo(float).eps, 1 - np.finfo(float).eps)
        return -(y * np.log(prob) + (1 - y) * np.log1p(-prob))
    with np.errstate(over="raise", invalid="raise"):
        return (y - prediction) ** 2
