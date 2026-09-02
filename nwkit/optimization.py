"""Shared, deterministic optimization helpers for phylogenetic models."""

import math
import warnings
from dataclasses import dataclass
from functools import lru_cache
from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize, minimize_scalar


class FitResourceError(ValueError):
    """A requested fit exceeds a resource limit, not an invalid model candidate."""


class ScalarFitCache:
    """Retain scalar objectives for every trial and at most two complete fits.

    Optimizers revisit trial points, but retaining every covariance and its
    factorization makes storage grow with the search length. Evicted fits are
    rebuilt only when their full results are requested (usually at selection).
    Resource failures must escape the search instead of biasing its optimum.
    """

    def __init__(self, fit, score, *, invalid_errors=(ValueError,)):
        self._fit = fit
        self._score = score
        self._invalid_errors = invalid_errors
        self._objectives = {}
        self.candidate = lru_cache(maxsize=2)(self._candidate)

    def _candidate(self, value):
        try:
            candidate = self._fit(value)
        except FitResourceError:
            raise
        except self._invalid_errors:
            candidate = None
        self._objectives[value] = (
            float("inf") if candidate is None else float(self._score(candidate))
        )
        return candidate

    def objective(self, value):
        if value not in self._objectives:
            self.candidate(value)
        return self._objectives[value]

    def best(self, values):
        """Select using cheap scalar scores, then retain only the winning fit."""
        try:
            finite = [value for value in values if math.isfinite(self.objective(value))]
            if not finite:
                return None
            return self.candidate(min(finite, key=self.objective))
        finally:
            self.candidate.cache_clear()


def _bounded_scalar_minimize(function, bounds):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in scalar subtract",
            category=RuntimeWarning,
            module=r"scipy\.optimize\._optimize",
        )
        return minimize_scalar(
            function,
            bounds=bounds,
            method="bounded",
            options={"xatol": 1e-7, "maxiter": 500},
        )


def global_bounded_scalar_minimize(function, bounds, *, grid_size=17):
    """Search a bounded one-dimensional objective without assuming unimodality."""

    lower, upper = (float(bounds[0]), float(bounds[1]))
    if not all(math.isfinite(value) for value in (lower, upper)) or lower >= upper:
        raise ValueError("Scalar optimizer bounds must be increasing and finite.")
    if grid_size < 3:
        raise ValueError("A global scalar search requires at least three grid points.")
    grid = np.linspace(lower, upper, grid_size)
    grid_values = np.asarray([float(function(value)) for value in grid], dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(grid_values))
    if not len(finite_indices):
        return SimpleNamespace(
            x=(lower + upper) / 2.0,
            fun=float("inf"),
            success=False,
            message="global grid search found no finite objective",
            grid_evaluations=grid_size,
            converged_refinements=0,
        )

    local_indices = []
    for index in finite_indices:
        if index in {0, grid_size - 1}:
            continue
        if (
            grid_values[index] <= grid_values[index - 1]
            and grid_values[index] <= grid_values[index + 1]
        ):
            local_indices.append(int(index))
    ranked_indices = sorted(
        (int(index) for index in finite_indices if 0 < index < grid_size - 1),
        key=lambda index: grid_values[index],
    )
    refinement_indices = []
    for index in local_indices + ranked_indices:
        if index not in refinement_indices:
            refinement_indices.append(index)
        if len(refinement_indices) == 4:
            break

    candidates = [
        SimpleNamespace(
            x=float(grid[index]),
            fun=float(grid_values[index]),
            success=index in {0, grid_size - 1},
            message="global grid point",
        )
        for index in finite_indices
    ]
    successful_refinements = 0
    for index in refinement_indices:
        refined = _bounded_scalar_minimize(
            function,
            (float(grid[index - 1]), float(grid[index + 1])),
        )
        if math.isfinite(float(refined.fun)):
            candidates.append(refined)
            successful_refinements += int(bool(refined.success))
    best = min(candidates, key=lambda candidate: float(candidate.fun))
    best_is_boundary = math.isclose(float(best.x), lower) or math.isclose(
        float(best.x), upper
    )
    return SimpleNamespace(
        x=float(best.x),
        fun=float(best.fun),
        success=bool(best_is_boundary or successful_refinements),
        message=(
            "global grid search ({} points; {} successful local refinement(s))"
        ).format(grid_size, successful_refinements),
        grid_evaluations=grid_size,
        converged_refinements=successful_refinements,
    )


@dataclass(frozen=True, slots=True)
class MultistartResult:
    x: np.ndarray
    fun: float
    success: bool
    message: str
    starts: int
    converged_starts: int
    failed_starts: int


def _unique_vectors(vectors):
    unique = []
    seen = set()
    for vector in vectors:
        array = np.asarray(vector, dtype=float)
        key = tuple(float(value) for value in array)
        if key not in seen:
            unique.append(array)
            seen.add(key)
    return unique


def deterministic_multistart(
    objective,
    initial,
    bounds,
    *,
    maxiter=1000,
    fallback=True,
    additional_starts=(),
    minimizer=minimize,
    ftol=1e-12,
):
    """Run reproducible L-BFGS-B starts with an optional Powell fallback."""

    initial = np.asarray(initial, dtype=float)
    bounds = tuple((lower, upper) for lower, upper in bounds)
    if initial.ndim != 1 or len(initial) != len(bounds):
        raise ValueError("Optimizer initial values and bounds have different sizes.")
    if not len(initial):
        value = float(objective(initial))
        if not math.isfinite(value):
            raise ValueError("The fixed-parameter likelihood is not finite.")
        return MultistartResult(initial, value, True, "all parameters fixed", 0, 0, 0)

    starts = [initial]
    finite_indices = [
        index
        for index, (lower, upper) in enumerate(bounds)
        if lower is not None
        and upper is not None
        and math.isfinite(float(lower))
        and math.isfinite(float(upper))
    ]
    for fraction in (0.1, 0.25, 0.5, 0.75, 0.9):
        start = initial.copy()
        for index in finite_indices:
            lower, upper = bounds[index]
            assert lower is not None and upper is not None
            start[index] = float(lower) + fraction * (float(upper) - float(lower))
        starts.append(start)
    if len(finite_indices) > 1:
        for reverse in (False, True):
            start = initial.copy()
            fractions = np.linspace(0.2, 0.8, len(finite_indices))
            if reverse:
                fractions = fractions[::-1]
            for fraction, index in zip(fractions, finite_indices, strict=True):
                lower, upper = bounds[index]
                assert lower is not None and upper is not None
                start[index] = float(lower) + float(fraction) * (
                    float(upper) - float(lower)
                )
            starts.append(start)
    starts.extend(additional_starts)
    starts = _unique_vectors(starts)

    candidates = []
    failures = []
    for start in starts:
        options = {"maxiter": maxiter}
        if ftol is not None:
            options["ftol"] = ftol
        result = minimizer(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options=options,
        )
        messages = [str(result.message)]
        success = (
            bool(result.success)
            and math.isfinite(float(result.fun))
            and float(result.fun) < 1e99
            and np.all(np.isfinite(result.x))
        )
        if fallback and not success:
            fallback_result = minimizer(
                objective,
                start,
                method="Powell",
                bounds=bounds,
                options={"maxiter": maxiter * 2, "xtol": 1e-8, "ftol": 1e-10},
            )
            messages.append(f"Powell fallback: {fallback_result.message}")
            if (
                bool(fallback_result.success)
                and math.isfinite(float(fallback_result.fun))
                and float(fallback_result.fun) < 1e99
                and np.all(np.isfinite(fallback_result.x))
            ):
                result = fallback_result
                success = True
        if success:
            candidates.append((float(result.fun), result, "; ".join(messages)))
        else:
            failures.append("; ".join(messages))
    if not candidates:
        raise ValueError(
            "No optimizer start produced a finite fit: "
            + ("; ".join(dict.fromkeys(failures)) or "unknown optimizer failure")
        )
    _, selected, message = min(candidates, key=lambda item: item[0])
    return MultistartResult(
        x=np.asarray(selected.x, dtype=float),
        fun=float(selected.fun),
        success=bool(selected.success),
        message=message,
        starts=len(starts),
        converged_starts=len(candidates),
        failed_starts=len(failures),
    )
