"""Independent spectral checks of the bounded Brownian rate search."""

import math

import numpy as np
import pytest
from numpy.polynomial import Polynomial
from scipy.optimize import minimize_scalar

from nwkit import continuous_asr_optimize as optimize


def spectral_point(rate, noise, squares):
    variances = rate + noise
    return optimize.RatePoint(
        rate, float(np.log(variances).sum()), float((squares / variances).sum())
    )


@pytest.mark.parametrize("dimension", [1, 3, 8])
@pytest.mark.parametrize("exact", [False, True])
def test_interval_bounds_do_not_exclude_any_sampled_likelihood(dimension, exact):
    rng = np.random.default_rng(1762 + dimension)
    for _ in range(20):
        noise = 10 ** rng.uniform(-8, 8, size=dimension)
        squares = 10 ** rng.uniform(-10, 10, size=dimension)
        if exact:
            noise[0] = 0.0
        shift = float(noise.min()) / 2
        edges = sorted(10 ** rng.uniform(-10, 10, size=5))
        if not exact:
            edges[0] = 0.0
        for left, right in zip(edges[:-1], edges[1:], strict=True):
            bound = optimize.interval_lower_bound(
                spectral_point(left, noise, squares),
                spectral_point(right, noise, squares),
                shift,
                dimension,
            )
            # Build the sample grid independently of the production interpolator.
            rates = np.maximum(
                left,
                np.exp(
                    np.linspace(math.log(left + shift), math.log(right + shift), 101)
                )
                - shift,
            )
            variances = rates[:, None] + noise
            scores = 0.5 * (np.log(variances) + squares / variances).sum(axis=1)
            assert np.all(bound <= scores + 1e-12 * np.maximum(1.0, np.abs(scores)))


@pytest.mark.parametrize("exact", [False, True])
def test_global_search_matches_all_roots_of_the_spectral_score_polynomial(exact):
    rng = np.random.default_rng(581934)
    for _ in range(40):
        noise = 10 ** rng.uniform(-3, 3, size=3)
        squares = 10 ** rng.uniform(-3, 3, size=3)
        if exact:
            noise[0] = 0.0
        # Multiplying the score by prod(r+lambda_i)^2 gives a polynomial;
        # enumerating its roots checks all modes without any likelihood grid.
        polynomial = Polynomial([0.0])
        for index in range(3):
            term = Polynomial([noise[index] - squares[index], 1.0])
            for other in range(3):
                if other != index:
                    term *= Polynomial([noise[other], 1.0]) ** 2
            polynomial += term
        lower = float(squares[0] / 3) if exact else 0.0
        upper = float(squares.sum())
        candidates = [lower, upper] + [
            float(root.real)
            for root in polynomial.roots()
            if lower < root.real < upper and abs(root.imag) < 1e-8 * abs(root.real)
        ]
        expected = min(
            spectral_point(rate, noise, squares).score for rate in candidates
        )

        def evaluate(rate, noise=noise, squares=squares):
            point = spectral_point(rate, noise, squares)
            return point.log_variance, point.quadratic

        rate = optimize.minimize_brownian_rate(
            evaluate,
            lower,
            upper,
            min(upper, float(noise.min()) / 2),
            3,
            minimize_scalar,
        )
        assert spectral_point(rate, noise, squares).score <= expected + 2e-9


@pytest.mark.parametrize(
    "left, right, shift",
    [
        (math.ulp(0.0), 1e308, 0.0),
        (0.0, 1e308, math.ulp(0.0)),
        (1.0, math.nextafter(1.0, math.inf), 0.0),
        (0.0, math.ulp(0.0), 1e308),
    ],
)
def test_log_rate_interpolation_is_bounded_at_extreme_scales(left, right, shift):
    rates = [
        optimize._interpolate(left, right, shift, fraction)
        for fraction in [0, 0.1, 0.5, 0.9, 1]
    ]
    assert rates == sorted(rates)
    assert rates[0] == left
    assert rates[-1] == right
    assert all(math.isfinite(rate) and left <= rate <= right for rate in rates)


def test_global_search_budget_exhaustion_is_an_error(monkeypatch):
    monkeypatch.setattr(optimize, "_MAX_EVALUATIONS", 5)
    with pytest.raises(ValueError, match="global search failed to converge"):
        optimize.minimize_brownian_rate(
            lambda rate: (math.log(rate + 1), 4 / (rate + 1)),
            0.0,
            4.0,
            1.0,
            1,
            minimize_scalar,
        )


def test_convergence_accounts_for_cancellation_between_likelihood_components():
    dimension = 1_000_000

    def evaluate(rate):
        # At the optimum, the two million-sized components cancel to zero.
        # An ULP tolerance based only on their sum cannot bound roundoff.
        return dimension * (math.log(rate + 0.1) - 1), dimension / (rate + 0.1)

    rate = optimize.minimize_brownian_rate(
        evaluate, 0.0, float(dimension), 0.1, dimension, minimize_scalar
    )
    assert rate == pytest.approx(0.9, rel=1e-7)


def test_global_search_uses_only_one_local_refinement():
    calls = 0

    def counted_minimize(*args, **kwargs):
        nonlocal calls
        calls += 1
        return minimize_scalar(*args, **kwargs)

    rate = optimize.minimize_brownian_rate(
        lambda candidate: (math.log(candidate + 1.0), 4.0 / (candidate + 1.0)),
        0.0,
        4.0,
        1.0,
        1,
        counted_minimize,
    )
    assert rate == pytest.approx(3.0, rel=1e-7)
    assert calls == 1


def test_search_interval_that_includes_zero_requires_a_positive_shift():
    with pytest.raises(
        ValueError, match="includes zero requires a positive rate shift"
    ):
        optimize.minimize_brownian_rate(
            lambda candidate: (math.log1p(candidate), 1.0 / (candidate + 1.0)),
            0.0,
            4.0,
            0.0,
            1,
            minimize_scalar,
        )


def test_global_search_does_not_skip_numerically_unusable_candidates():
    def evaluate(rate):
        if rate > 0.0:
            raise ValueError("unresolvable observation scale")
        return 0.0, 4.0

    with pytest.raises(ValueError, match="unresolvable observation scale"):
        optimize.minimize_brownian_rate(evaluate, 0.0, 4.0, 1.0, 1, minimize_scalar)


@pytest.mark.parametrize("quadratic", [math.inf, math.nan])
def test_nonfinite_global_objective_is_an_error(quadratic):
    with pytest.raises(ValueError, match="floating-point range"):
        optimize.minimize_brownian_rate(
            lambda rate: (0.0, quadratic), 0.0, 4.0, 1.0, 1, minimize_scalar
        )
