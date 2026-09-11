"""Bounded global search for a single Gaussian variance component.

On independent contrasts, twice the negative log likelihood (up to a
constant) is D(r) + Q(r), where D=log|r C + S| increases and Q decreases.
For 0 <= shift <= lambda_min(C^-1/2 S C^-1/2), put u=log(r+shift).
Writing the whitened terms as log(exp(u)+lambda-shift) + w^2/(exp(u)+lambda-shift)
shows f''(u) <= d/8 + Q(left)/2 on an interval. A chord minus this curvature
bound is a lower bound on f throughout the interval. Combining it with
(D(left)+Q(right))/2 lets branch-and-bound check basins that a grid misses.

Bounds and convergence are evaluated to floating-point likelihood tolerance;
failure to resolve an interval or exhaustively bound the search is an error.
Only scalar evaluations are retained, not per-node pruning states.
"""

import heapq
import math
from dataclasses import dataclass
from itertools import count

_MAX_EVALUATIONS = 20000
_LIKELIHOOD_TOLERANCE = 1e-9


@dataclass(frozen=True)
class RatePoint:
    rate: float
    log_variance: float
    quadratic: float

    @property
    def score(self):
        return 0.5 * (self.log_variance + self.quadratic)


def _log_sum(first, second):
    larger, smaller = max(first, second), min(first, second)
    return math.log(larger) + math.log1p(smaller / larger)


def _log_width(left, right, shift):
    scale = max(left, shift)
    relative_width = ((right - left) / scale) / (left / scale + shift / scale)
    if math.isfinite(relative_width):
        return math.log1p(relative_width)
    return _log_sum(right, shift) - _log_sum(left, shift)


def _interpolate(left, right, shift, fraction):
    if fraction <= 0.0:
        return left
    if fraction >= 1.0:
        return right
    width = _log_width(left, right, shift)
    if width == 0.0:
        return left + (right - left) * fraction
    # Algebraically (right-left)*expm1(fraction*width)/expm1(width),
    # without overflowing exp(width) or underflowing a tiny ratio first.
    increment = math.exp(math.log(right - left) - width * (1.0 - fraction))
    increment *= -math.expm1(-width * fraction) / -math.expm1(-width)
    return min(right, max(left, left + increment))


def interval_lower_bound(left, right, shift, residual_df):
    width = _log_width(left.rate, right.rate, shift)
    curvature = residual_df / 8.0 + left.quadratic / 2.0
    gap = 0.5 * curvature * width * width
    chord_bound = -math.inf
    if math.isfinite(gap):
        weight = (
            min(1.0, max(0.0, 0.5 + 0.5 * (left.score - right.score) / gap))
            if gap > 0.0
            else 0.0
            if left.score <= right.score
            else 1.0
        )
        chord_bound = (
            (1.0 - weight) * left.score
            + weight * right.score
            - gap * weight * (1.0 - weight)
        )
    monotone_bound = 0.5 * (left.log_variance + right.quadratic)
    rounding = 64.0 * math.ulp(
        max(
            abs(left.log_variance),
            abs(right.log_variance),
            left.quadratic,
            right.quadratic,
            1.0,
        )
    )
    return max(chord_bound, monotone_bound) - rounding


def _initial_intervals(lower, upper, shift):
    pending = [(lower, upper)]
    intervals = []
    while pending:
        left, right = pending.pop()
        if _log_width(left, right, shift) <= 2.0:
            intervals.append((left, right))
            continue
        middle = _interpolate(left, right, shift, 0.5)
        if not left < middle < right:
            raise ValueError("Brownian REML search exceeds floating-point range.")
        pending.extend(((middle, right), (left, middle)))
    return intervals


class _RateSearch:
    def __init__(self, evaluate, shift, minimize_scalar):
        self.evaluate = evaluate
        self.shift = shift
        self.minimize_scalar = minimize_scalar
        self.points: dict[float, RatePoint] = {}

    def at(self, rate):
        rate = float(rate)
        if rate not in self.points:
            if len(self.points) >= _MAX_EVALUATIONS:
                raise ValueError("Brownian REML global search failed to converge.")
            log_variance, quadratic = self.evaluate(rate)
            point = RatePoint(rate, log_variance, quadratic)
            if not math.isfinite(point.score):
                raise ValueError("Brownian REML search exceeds floating-point range.")
            self.points[rate] = point
        return self.points[rate]

    def refine(self, left, right):
        def objective(fraction):
            return self.at(_interpolate(left, right, self.shift, fraction)).score

        result = self.minimize_scalar(
            objective,
            bounds=(0.0, 1.0),
            method="bounded",
            options={"xatol": 1e-12},
        )
        if not result.success or not math.isfinite(result.fun):
            raise ValueError("Brownian REML rate optimization failed to converge.")
        return self.at(_interpolate(left, right, self.shift, float(result.x)))


def minimize_brownian_rate(evaluate, lower, upper, shift, residual_df, minimize_scalar):
    if lower == upper:
        return lower
    if lower == 0.0 and shift <= 0.0:
        raise ValueError(
            "A Brownian REML search that includes zero requires a positive rate shift."
        )
    search = _RateSearch(evaluate, shift, minimize_scalar)
    intervals = _initial_intervals(lower, upper, shift)
    serial = count()
    pending: list[tuple[float, int, RatePoint, RatePoint]] = []

    def enqueue(left, right):
        bound = interval_lower_bound(left, right, shift, residual_df)
        heapq.heappush(pending, (bound, next(serial), left, right))

    for left, right in intervals:
        enqueue(search.at(left), search.at(right))
    best = min(search.points.values(), key=lambda point: (point.score, point.rate))
    while pending:
        bound, _, left, right = heapq.heappop(pending)
        # D and Q can cancel even after removing within-position constants.
        # Resolution follows those components, not their rounded sum alone.
        resolution = math.ulp(max(abs(best.log_variance), best.quadratic, 1.0))
        tolerance = _LIKELIHOOD_TOLERANCE + 128.0 * resolution
        if bound >= best.score - tolerance:
            break
        rate = _interpolate(left.rate, right.rate, shift, 0.5)
        if not left.rate < rate < right.rate:
            raise ValueError("Brownian REML optimum exceeds floating-point resolution.")
        middle = search.at(rate)
        best = min((best, middle), key=lambda point: (point.score, point.rate))
        enqueue(left, middle)
        enqueue(middle, right)
    neighbors = sorted(search.points)
    index = neighbors.index(best.rate)
    if 0 < index < len(neighbors) - 1:
        refined = search.refine(neighbors[index - 1], neighbors[index + 1])
        best = min((best, refined), key=lambda point: (point.score, point.rate))
    if lower == 0.0:
        zero = search.at(0.0)
        resolution = math.ulp(max(abs(best.log_variance), best.quadratic, 1.0))
        if zero.score <= best.score + 8.0 * resolution:
            return 0.0
    return best.rate
