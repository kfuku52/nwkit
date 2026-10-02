"""Candidate scanning and plug-in, search-wide parametric bootstrap for counts."""

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from scipy.optimize import minimize

from nwkit.wgd_count_model import CountLikelihood, MultiplicationEvent


@dataclass(frozen=True)
class CountFit:
    rates: np.ndarray
    root_mean: float
    log_likelihood: float
    num_parameters: int
    max_count: int
    state_error: float
    converged: bool
    boundary: bool
    message: str
    event: MultiplicationEvent | None = None

    @property
    def aic(self):
        return 2 * self.num_parameters - 2 * self.log_likelihood


@dataclass(frozen=True)
class ScanCandidate:
    node: int
    event_fit: CountFit
    burst_fit: CountFit
    improvement: float
    burst_aic_difference: float
    p_value: float | None = None
    p_value_mc_se: float | None = None


@dataclass(frozen=True)
class CountScan:
    background: CountFit
    candidates: tuple[ScanCandidate, ...]
    bootstrap_statistics: tuple[float, ...] = ()
    calibration: str = "not-run"
    fractions: tuple[float, ...] = (0.25, 0.5, 0.75)
    multiplicity: int = 2


def _objective_gradient(objective, parameters, bounds, *, step_factor=1.0):
    """Bound-aware Richardson derivatives of scalar or per-family objectives.

    Subtract per-family values before summing to avoid cancellation. Cancelling
    the three-point stencil's leading error is important near root mean one,
    where a fixed log-parameter step can otherwise have substantial h^2 bias.
    """
    gradient = np.empty_like(parameters)
    value = np.asarray(objective(parameters))
    for column, (low, high) in enumerate(bounds):
        step = (
            step_factor
            * np.cbrt(np.finfo(float).eps)
            * max(1.0, abs(parameters[column]))
        )
        direction = 1 if parameters[column] + 2 * step <= high else -1
        central = low <= parameters[column] - step and parameters[column] + step <= high
        estimates = []
        # Use the same stencil direction at both scales, including near bounds.
        for increment in (step, step / 2):
            first, second = parameters.copy(), parameters.copy()
            if central:
                first[column] += increment
                second[column] -= increment
                difference = np.asarray(objective(first)) - objective(second)
                estimate = np.sum(difference) / (2 * increment)
            else:
                first[column] += direction * increment
                second[column] += direction * 2 * increment
                difference = 4 * (np.asarray(objective(first)) - value) - (
                    np.asarray(objective(second)) - value
                )
                estimate = direction * np.sum(difference) / (2 * increment)
            estimates.append(estimate)
        gradient[column] = (4 * estimates[1] - estimates[0]) / 3
    return gradient


def _projected_gradient_norm(parameters, gradient, bounds):
    lower, upper = np.asarray(bounds).T
    return float(
        np.max(np.abs(parameters - np.clip(parameters - gradient, lower, upper)))
    )


def _fit_at_bound(
    model,
    max_count,
    *,
    node=None,
    fraction=0.5,
    multiplicity=2,
    initial=None,
    max_iterations=200,
):
    positive_lengths = [t for t in model.tree.lengths[1:] if t > 0]
    if not positive_lengths:
        raise ValueError("Count fitting needs at least one positive branch length.")
    time_unit = float(np.median(positive_lengths))
    groups = model.num_groups
    observed_max = float(np.nanmax(model.counts))
    root_upper = max(10.0, observed_max * 2)
    bounds = [(-16.0, np.log(20.0))] * (2 * groups) + [(0.0, np.log(root_upper))]
    if node is not None:
        bounds.append((0.0, 1.0))

    def unpack(parameters):
        rates = np.exp(parameters[: 2 * groups]).reshape(groups, 2) / time_unit
        root_mean = float(np.exp(parameters[2 * groups]))
        event = (
            None
            if node is None
            else MultiplicationEvent(
                node, float(parameters[-1]), fraction, multiplicity
            )
        )
        return rates, root_mean, event

    def objective(parameters):
        rates, root_mean, event = unpack(parameters)
        value = model.log_likelihood(rates, root_mean, max_count, event)
        if not np.isfinite(value):
            raise ValueError("Count optimizer encountered a nonfinite likelihood.")
        return -value

    def family_objective(parameters):
        rates, root_mean, event = unpack(parameters)
        values = model.family_log_likelihoods(rates, root_mean, max_count, event)
        if not np.all(np.isfinite(values)):
            raise ValueError("Count optimizer encountered a nonfinite likelihood.")
        return -values

    starts = []
    for rate in (0.05, 0.5):
        parameters = np.full(2 * groups + 1, np.log(rate))
        parameters[-1] = np.log(1.1)
        if node is not None:
            parameters = np.append(parameters, 0.5)
        starts.append(parameters)
    if initial is not None:
        rates = initial.rates
        if len(rates) != groups:
            raise ValueError("Initial fit rate groups do not match the count model.")
        parameters = np.append(
            np.log(rates.ravel() * time_unit), np.log(initial.root_mean)
        )
        if node is not None:
            parameters = np.append(
                parameters, initial.event.retention if initial.event else 0.0
            )
        starts.insert(0, parameters)
    results = [
        minimize(
            objective,
            start,
            method="L-BFGS-B",
            jac="3-point",
            bounds=bounds,
            options={"maxiter": max_iterations, "ftol": 0.0, "gtol": 1e-5},
        )
        for start in starts
    ]
    valid, messages = [], []

    def verify(result, label=""):
        # Validation uses twice the polishing step, an independently evaluated
        # derivative estimate rather than the optimizer's reported Jacobian.
        gradient = _objective_gradient(
            family_objective, result.x, bounds, step_factor=2.0
        )
        projected = _projected_gradient_norm(result.x, gradient, bounds)
        messages.append(f"{label}{result.message}; projected gradient={projected:.6g}")
        if result.success and np.isfinite(result.fun) and projected <= 1e-5:
            valid.append(result)

    for result in results:
        verify(result)
    if not valid:
        # Function stagnation can precede true stationarity when the initial
        # three-point Jacobian is biased. Retain that failed attempt and use
        # the remaining per-start iteration budget with higher-order gradients.
        for result in results:
            remaining = max_iterations - int(result.nit)
            if not result.success or remaining <= 0:
                continue
            polished = minimize(
                objective,
                result.x,
                method="L-BFGS-B",
                jac=lambda parameters: _objective_gradient(
                    family_objective, parameters, bounds
                ),
                bounds=bounds,
                options={"maxiter": remaining, "ftol": 0.0, "gtol": 1e-5},
            )
            verify(polished, "Higher-order polish: ")
    if not valid:
        raise ValueError(f"Count optimizer did not converge: {'; '.join(messages)}")
    best = min(valid, key=lambda result: result.fun)
    rates, root_mean, event = unpack(best.x)
    boundary = bool(
        any(
            min(abs(value - low), abs(value - high)) < 1e-4
            for value, (low, high) in zip(
                best.x[: 2 * groups], bounds[: 2 * groups], strict=True
            )
        )
        or abs(best.x[2 * groups] - bounds[2 * groups][1]) < 1e-4
    )
    return CountFit(
        rates,
        root_mean,
        -float(best.fun),
        len(best.x),
        max_count,
        float("inf"),
        False,
        boundary,
        str(best.message),
        event,
    )


def fit_counts(
    model: CountLikelihood,
    *,
    node: int | None = None,
    fraction: float = 0.5,
    multiplicity: int = 2,
    initial: CountFit | None = None,
    max_states: int = 256,
    state_tolerance: float = 1e-7,
    max_iterations: int = 200,
) -> CountFit:
    """Fit and independently check per-family truncation error at twice the bound."""
    if not np.isfinite(state_tolerance) or state_tolerance <= 0:
        raise ValueError("State tolerance must be finite and positive.")
    bound = max(8, int(np.nanmax(model.counts)) + 2)
    if initial is not None:
        bound = max(bound, initial.max_count)
    if bound * 2 > max_states:
        raise ValueError(
            "--max-states must allow at least twice the initial count bound."
        )
    while 2 * bound <= max_states:
        fit = _fit_at_bound(
            model,
            bound,
            node=node,
            fraction=fraction,
            multiplicity=multiplicity,
            initial=initial,
            max_iterations=max_iterations,
        )
        error = model.convergence(fit.rates, fit.root_mean, bound, fit.event)
        if error <= state_tolerance:
            return replace(fit, state_error=error, converged=True)
        initial = fit
        bound *= 2
    raise ValueError(
        f"State truncation did not converge within --max-states={max_states}; "
        "increase the bound or inspect high-copy families and rate boundaries."
    )


def _burst_model(model, node):
    groups = list(model.branch_groups)
    groups[node] = model.num_groups
    # If the original group had one branch, compact the remaining group IDs.
    remap = {value: index for index, value in enumerate(sorted(set(groups[1:])))}
    assignments = tuple(remap.get(value, 0) for value in groups)
    return CountLikelihood(
        model.tree,
        model.counts,
        detection=model.detection,
        rate_scales=model.rate_scales,
        branch_groups=assignments,
        ascertainment=model.ascertainment,
    )


def scan_counts(
    model: CountLikelihood,
    *,
    nodes: tuple[int, ...] | None = None,
    fractions: tuple[float, ...] = (0.25, 0.5, 0.75),
    multiplicity: int = 2,
    max_states: int = 256,
    state_tolerance: float = 1e-7,
    max_iterations: int = 200,
    background: CountFit | None = None,
    compare_bursts: bool = True,
) -> CountScan:
    """Scan one event at a time; multiple simultaneous events are not inferred."""
    if not fractions or any(not np.isfinite(f) or not 0 <= f <= 1 for f in fractions):
        raise ValueError("Event fractions need a nonempty finite grid in [0, 1].")
    if nodes is None:
        nodes = tuple(
            node for node, time in enumerate(model.tree.lengths[1:], 1) if time > 0
        )
    if (
        not nodes
        or len(set(nodes)) != len(nodes)
        or any(
            isinstance(node, bool)
            or not isinstance(node, int)
            or not 1 <= node < len(model.tree.parents)
            or model.tree.lengths[node] <= 0
            for node in nodes
        )
    ):
        raise ValueError(
            "Candidate nodes must uniquely identify positive-length non-root branches."
        )
    options: dict[str, Any] = dict(
        max_states=max_states,
        state_tolerance=state_tolerance,
        max_iterations=max_iterations,
    )
    background = background or fit_counts(model, **options)
    candidates = []
    for node in nodes:
        fits = [
            fit_counts(
                model,
                node=node,
                fraction=fraction,
                multiplicity=multiplicity,
                initial=background,
                **options,
            )
            for fraction in fractions
        ]
        fit = max(fits, key=lambda result: result.log_likelihood)
        # q=0 is exactly the fitted null even if a local optimizer misses it.
        if fit.log_likelihood < background.log_likelihood:
            fit = replace(
                background,
                event=MultiplicationEvent(node, 0, fractions[0], multiplicity),
                num_parameters=background.num_parameters + 1,
            )
        if compare_bursts:
            burst_model = _burst_model(model, node)
            burst_rates = np.zeros((burst_model.num_groups, 2))
            for branch in range(1, len(model.tree.parents)):
                burst_rates[burst_model.branch_groups[branch]] = background.rates[
                    model.branch_groups[branch]
                ]
            burst_initial = replace(background, rates=burst_rates)
            burst = fit_counts(burst_model, initial=burst_initial, **options)
        else:
            burst = background
        candidates.append(
            ScanCandidate(
                node,
                fit,
                burst,
                max(0.0, 2 * (fit.log_likelihood - background.log_likelihood)),
                burst.aic - fit.aic,
            )
        )
    candidates.sort(
        key=lambda candidate: (
            -candidate.improvement,
            model.tree.clade_ids[candidate.node],
        )
    )
    return CountScan(
        background,
        tuple(candidates),
        fractions=tuple(fractions),
        multiplicity=multiplicity,
    )


def calibrate_scan(
    model: CountLikelihood,
    observed: CountScan,
    draws: int,
    seed: int,
    *,
    fractions: tuple[float, ...] | None = None,
    multiplicity: int | None = None,
    max_states: int = 256,
    state_tolerance: float = 1e-7,
    max_iterations: int = 200,
    progress=None,
) -> CountScan:
    """Plug-in null bootstrap of the maximum statistic across searched branches.

    This handles zero-retention and finite candidate-grid selection. It is not
    a supremum-nuisance test and does not guarantee control under misspecified
    rates, ascertainment, annotation errors, or simultaneous genome events.
    """
    if isinstance(draws, bool) or not isinstance(draws, int) or draws < 1:
        raise ValueError("Bootstrap draws must be a positive integer.")
    if not isinstance(seed, int) or seed < 0:
        raise ValueError("Seed must be a nonnegative integer.")
    fractions = observed.fractions if fractions is None else tuple(fractions)
    multiplicity = observed.multiplicity if multiplicity is None else multiplicity
    if fractions != observed.fractions or multiplicity != observed.multiplicity:
        raise ValueError(
            "Bootstrap fractions and multiplicity must match the observed search."
        )
    if not observed.candidates:
        raise ValueError("Bootstrap needs a nonempty observed candidate search.")
    rng = np.random.default_rng(seed)
    statistics = []
    for draw in range(draws):
        counts = model.simulate(
            observed.background.rates, observed.background.root_mean, rng
        )
        simulated = CountLikelihood(
            model.tree,
            counts,
            detection=model.detection,
            rate_scales=model.rate_scales,
            branch_groups=model.branch_groups,
            ascertainment=model.ascertainment,
        )
        result = scan_counts(
            simulated,
            nodes=tuple(candidate.node for candidate in observed.candidates),
            fractions=fractions,
            multiplicity=multiplicity,
            max_states=max_states,
            state_tolerance=state_tolerance,
            max_iterations=max_iterations,
            compare_bursts=False,
        )
        statistics.append(max(candidate.improvement for candidate in result.candidates))
        if progress is not None:
            progress(draw + 1, draws)
    candidates = []
    for candidate in observed.candidates:
        exceedances = sum(value >= candidate.improvement - 1e-9 for value in statistics)
        p_value = (1 + exceedances) / (1 + draws)
        mc_se = float(np.sqrt(p_value * (1 - p_value) / (draws + 1)))
        candidates.append(replace(candidate, p_value=p_value, p_value_mc_se=mc_se))
    return replace(
        observed,
        candidates=tuple(candidates),
        bootstrap_statistics=tuple(statistics),
        calibration="plugin-parametric-bootstrap-search-maximum",
    )
