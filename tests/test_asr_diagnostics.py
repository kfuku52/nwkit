import numpy as np
import pytest

from nwkit.asr_diagnostics import (
    gaussian_posterior_predictive,
    parametric_bootstrap,
    profile_likelihood_interval,
)
from nwkit.evolution import build_evolutionary_process
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def test_profile_likelihood_interval_recovers_normal_quadratic():
    interval = profile_likelihood_interval(
        lambda value: -0.5 * ((value - 2.0) / 0.5) ** 2,
        2.0,
        (-2.0, 6.0),
    )
    assert interval.lower == pytest.approx(2.0 - 1.9599639845 * 0.5, abs=1e-8)
    assert interval.upper == pytest.approx(2.0 + 1.9599639845 * 0.5, abs=1e-8)


@pytest.mark.parametrize("mle", [0.0, 4.0])
def test_profile_likelihood_interval_marks_boundary_mle(mle):
    interval = profile_likelihood_interval(
        lambda value: -0.5 * ((value - mle) / 0.5) ** 2,
        mle,
        (0.0, 4.0),
    )
    if mle == 0.0:
        assert interval.lower == 0.0
        assert interval.lower_boundary_limited
        assert not interval.upper_boundary_limited
    else:
        assert interval.upper == 4.0
        assert interval.upper_boundary_limited
        assert not interval.lower_boundary_limited


def test_profile_likelihood_marks_invalid_search_region_as_limited():
    def profile(value):
        if value > 3.0:
            raise ValueError("outside numerical support")
        return -0.5 * (value - 2.0) ** 2

    interval = profile_likelihood_interval(profile, 2.0, (-2.0, 6.0))
    assert interval.upper <= 3.0
    assert interval.upper_boundary_limited


def test_gaussian_posterior_predictive_is_reproducible():
    tree = tree_from("((A:1,B:1)I:1,(C:1,D:1)J:1)R;")
    process = build_evolutionary_process(
        tree, model="brownian", root_mode="flat"
    ).scaled_variance(0.7)
    values = {"A": -0.2, "B": 0.4, "C": 1.0, "D": 1.2}
    first = gaussian_posterior_predictive(process, values, num_simulations=200, seed=19)
    second = gaussian_posterior_predictive(
        process, values, num_simulations=200, seed=19
    )
    assert first.equals(second)
    assert first["p_two_sided"].between(0.0, 1.0).all()


def test_single_posterior_predictive_replicate_has_defined_sd():
    tree = tree_from("(A:1,B:1,C:1)R;")
    process = build_evolutionary_process(
        tree, model="brownian", root_mode="flat"
    ).scaled_variance(0.7)
    result = gaussian_posterior_predictive(
        process,
        {"A": -0.2, "B": 0.4, "C": 1.0},
        num_simulations=1,
        seed=19,
    )
    assert (result["replicate_sd"] == 0.0).all()


def test_callback_parametric_bootstrap_records_successes_and_failures():
    def simulator(seed):
        return int(np.random.default_rng(seed).integers(0, 10))

    def fitter(value):
        if value == 0:
            raise ValueError("synthetic failure")
        return value * 2

    successes, failures = parametric_bootstrap(
        simulator,
        fitter,
        lambda value: {"estimate": value},
        num_simulations=50,
        seed=8,
    )
    assert len(successes) + len(failures) == 50
    assert set(successes["fit_status"]) == {"ok"}
