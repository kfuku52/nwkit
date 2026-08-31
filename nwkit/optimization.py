"""Memory-bounded caching for scalar likelihood searches."""

import math
from functools import lru_cache


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
