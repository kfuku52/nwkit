from decimal import Decimal, localcontext

import numpy as np
import pytest

from nwkit.phylogenetic_glmm import _count_log_likelihood, _count_log_pmf


def _decimal_nb_log_probabilities(values, log_means, dispersion):
    """Independent finite-product reference for integer counts, in 80 digits."""
    probabilities = []
    zero_probabilities = []
    with localcontext() as context:
        context.prec = 80
        one = Decimal(1)
        size = one / Decimal(float(dispersion))
        for value, log_mean in zip(values, log_means, strict=True):
            count = int(value)
            eta = Decimal(float(log_mean))
            denominator = (one + eta.exp() / size).ln()
            factorial = Decimal(1)
            ratio = Decimal(0)
            for integer in range(count):
                ratio += (one + Decimal(integer) / size).ln()
                factorial *= integer + 1
            probabilities.append(
                ratio
                - factorial.ln()
                + Decimal(count) * eta
                - (size + count) * denominator
            )
            zero_probabilities.append(-size * denominator)
    return probabilities, zero_probabilities


@pytest.mark.parametrize(
    "dispersion", [0.4, 1 / 127, 1 / 128, np.exp(-12), 1e-8, 1e-12]
)
@pytest.mark.parametrize(
    "family",
    [
        "negative-binomial",
        "zero-inflated-negative-binomial",
        "hurdle-negative-binomial",
    ],
)
def test_count_likelihood_matches_high_precision_near_poisson_limit(family, dispersion):
    values = np.array([0.0, 1.0, 4.0, 17.0, 100.0])
    log_means = np.log([0.3, 1.0, 2.5, 12.0, 80.0])
    log_pmf, log_zero = _decimal_nb_log_probabilities(values, log_means, dispersion)
    actual_pmf, actual_zero = _count_log_pmf(
        values, np.exp(log_means), family, dispersion
    )
    np.testing.assert_allclose(
        actual_pmf, list(map(float, log_pmf)), rtol=0, atol=2e-12
    )
    np.testing.assert_allclose(
        actual_zero, list(map(float, log_zero)), rtol=0, atol=2e-12
    )

    expected = []
    with localcontext() as context:
        context.prec = 80
        pi = Decimal(0.23)
        for count, log_probability, zero in zip(values, log_pmf, log_zero, strict=True):
            if family == "negative-binomial":
                probability = log_probability
            elif count == 0:
                probability = (
                    (pi + (1 - pi) * zero.exp()).ln()
                    if family.startswith("zero-inflated")
                    else pi.ln()
                )
            else:
                probability = (1 - pi).ln() + log_probability
                if family.startswith("hurdle"):
                    probability -= (1 - zero.exp()).ln()
            expected.append(float(probability))
    actual = _count_log_likelihood(
        values, log_means, family, dispersion, 0.23, np.zeros(len(values))
    )
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-12)


def test_nb_dispersion_difference_is_not_corrupted_by_gamma_cancellation():
    values = np.array([0.0, 1.0, 4.0, 17.0, 100.0])
    log_means = np.log([0.3, 1.0, 2.5, 12.0, 80.0])
    step = 1e-4
    dispersions = np.exp(np.array([-12.0 + step, -12.0 - step]))
    actual = [
        _count_log_likelihood(
            values, log_means, "negative-binomial", dispersion, None, np.zeros(5)
        )
        for dispersion in dispersions
    ]
    expected = [
        _decimal_nb_log_probabilities(values, log_means, dispersion)[0]
        for dispersion in dispersions
    ]
    with localcontext() as context:
        context.prec = 80
        reference = [
            float((plus - minus) / Decimal(2 * step))
            for plus, minus in zip(*expected, strict=True)
        ]
    np.testing.assert_allclose(
        (actual[0] - actual[1]) / (2 * step), reference, rtol=0, atol=2e-9
    )
