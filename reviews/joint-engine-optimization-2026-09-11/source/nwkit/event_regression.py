"""Event-average estimating equations, separate from covariance likelihoods."""

import numpy as np

from nwkit.gaussian import center_response, project_covariance
from nwkit.measurement_error import (
    _gmrf_uncertainty_update,
    _predictor_error_covariance,
    _structured_uncertainty_loading,
)
from nwkit.sparse_laplace import (
    ContinuousPredictorUncertainty,
    GmrfPredictorUncertainty,
    GroupedPredictorUncertainty,
    JointPredictorUncertainty,
    SparseLatentModel,
    prepare_sparse_latent_sampler,
)


def resolve_event_weighting(event_weighting=None, regression_estimand=None):
    """Resolve the explicit estimand and the legacy event-weighting alias."""
    aliases = {"event-average": "event", "common": "contrast"}
    if regression_estimand is not None and regression_estimand not in aliases:
        raise ValueError("regression_estimand must be event-average or common.")
    if event_weighting is not None and event_weighting not in {"event", "contrast"}:
        raise ValueError("event_weighting must be event or contrast.")
    if regression_estimand is not None:
        resolved = aliases[regression_estimand]
        if event_weighting is not None and event_weighting != resolved:
            raise ValueError("--regression-estimand conflicts with --event-weighting.")
        return resolved
    return event_weighting or "event"


def event_average_operator(design, groups):
    """Return L for beta_E = L y, with total loss weight one per event.

    The design is in raw contrast units. Evolutionary variances do not enter
    the weights: inverse-variance weighting would change the estimand under
    heterogeneous effects. The operator depends only on the fixed design and
    event membership, not fitted nuisance parameters or the response.
    """
    design = np.asarray(design, dtype=float)
    groups = np.asarray(groups)
    if design.ndim != 2 or groups.shape != (len(design),):
        raise ValueError("Event membership must align with design rows.")
    if not np.isfinite(design).all():
        raise ValueError("Event-average design must be finite.")
    _, inverse, counts = np.unique(groups, return_inverse=True, return_counts=True)
    weights = 1.0 / counts[inverse]
    weighted = design * np.sqrt(weights)[:, None]
    if np.linalg.matrix_rank(weighted) != design.shape[1]:
        raise ValueError("Event-average design is rank deficient.")
    # Use a QR solve rather than normal equations to avoid squaring condition
    # numbers when the predictor units differ substantially.
    q, r = np.linalg.qr(weighted, mode="reduced")
    return np.linalg.solve(r, q.T) * np.sqrt(weights)[None, :]


def event_average_estimate(response, design, groups):
    """Compute the event-average coefficient and its linear operator."""
    response, offset = center_response(response, design)
    operator = event_average_operator(design, groups)
    return operator @ response + offset, operator


def event_average_fit(response, design, groups, nuisance_fit):
    """Attach event-average estimates to an independently fitted Gaussian C.

    L C L' is the sampling covariance of this estimator, not inverse GLS
    information. With fitted C it is a model-based plug-in approximation.
    Because L has no nuisance parameters, its first-order derivative with
    respect to those parameters is zero. The likelihood and random-effect
    state remain those of the auxiliary common-coefficient fit.
    """
    beta, operator = event_average_estimate(response, design, groups)
    covariance = nuisance_fit["covariance"]
    factor = nuisance_fit["cholesky"]
    if "covariance_for_beta" in nuisance_fit:
        covariance, factor = nuisance_fit["covariance_for_beta"](beta)
    return {
        **nuisance_fit,
        "beta": beta,
        "beta_covariance": project_covariance(covariance, operator),
        "covariance": covariance,
        "cholesky": factor,
        "event_operator": operator,
    }


def prepare_shared_response_sampler(
    beta, component_variances, component_factors, uncertainties, columns, n_observations
):
    """Sample only event/lineage/predictor effects, without tip residual noise.

    Tip-space shape bootstrap already draws evolutionary and sampling noise.
    Adding only these shared terms preserves the original sampling covariance
    between selected contrasts, other contrasts, and the ancestral mean.
    """
    loadings = [
        np.sqrt(component_variances[name]) * factor
        for name, factor in component_factors.items()
    ]
    samplers = []
    for uncertainty, selected in zip(uncertainties, columns, strict=True):
        if isinstance(uncertainty, GmrfPredictorUncertainty):
            loading, precision, precision_factor = _gmrf_uncertainty_update(
                beta, selected, uncertainty
            )
            model = SparseLatentModel(precision, loading, precision_factor, 0.0, {}, {})
            samplers.append(prepare_sparse_latent_sampler(model))
        elif isinstance(
            uncertainty,
            (
                ContinuousPredictorUncertainty,
                GroupedPredictorUncertainty,
                JointPredictorUncertainty,
            ),
        ):
            loadings.append(
                _structured_uncertainty_loading(beta, selected, uncertainty)
            )
        else:
            covariance = _predictor_error_covariance(
                beta, uncertainty, selected, n_observations
            )
            eigenvalues, vectors = np.linalg.eigh(covariance)
            tolerance = (
                np.finfo(float).eps
                * max(1.0, float(np.max(np.abs(eigenvalues))))
                * n_observations
                * 100
            )
            if np.min(eigenvalues) < -tolerance:
                raise ValueError(
                    "Predictor sampling covariance is not positive semidefinite."
                )
            loadings.append(vectors * np.sqrt(np.maximum(eigenvalues, 0.0)))

    def draw(rng):
        result = np.zeros(n_observations)
        for loading in loadings:
            result += np.asarray(
                loading @ rng.standard_normal(loading.shape[1])
            ).reshape(-1)
        for sampler in samplers:
            result += sampler.sample(rng)
        return result

    return draw
