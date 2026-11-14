"""Seeded scenario regressions, not an empirical false-positive/power study.

The opt-in negative calibration requests 59 real refitted bootstrap datasets: its
minimum plus-one p-value is 1/60, so rejection at alpha=.05 is possible. One
fixed dataset cannot establish composite-null error control; Monte Carlo SE is
at most about .065. Failed simulated refits must abort rather than be discarded.

All families originate at the root. Family rates are homogeneous here; these
scenarios do not establish robustness to correlated detection, de novo family
origins, gamma misspecification, or simultaneous genome events.
"""

import numpy as np
import pytest

from nwkit.wgd_count_fit import calibrate_scan, fit_counts, scan_counts
from nwkit.wgd_count_model import (
    CountLikelihood,
    CountTree,
    MultiplicationEvent,
    birth_death_transition,
)


def _tree():
    return CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 0.6, 0.4, 0.4, 1.0),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )


def _simulate(
    families,
    rates,
    root_mean,
    seed,
    *,
    detection=None,
    branch_groups=None,
    event=None,
    missing=False,
):
    masks = np.ones((families, 3))
    if missing:
        masks[::4, 0] = np.nan
        masks[1::4, 1] = np.nan
    template = CountLikelihood(
        _tree(), masks, detection=detection, branch_groups=branch_groups
    )
    counts = template.simulate(rates, root_mean, np.random.default_rng(seed), event)
    np.testing.assert_array_equal(np.isnan(counts), np.isnan(masks))
    assert np.all(np.nansum(counts, axis=1) > 0)
    return CountLikelihood(_tree(), counts, detection=detection)


def _scan(model):
    # Search every non-root branch; fraction is deliberately fixed to bound cost.
    return scan_counts(model, fractions=(0.5,), max_states=128)


def _assert_verified(scan):
    for fit in (
        scan.background,
        *(candidate.event_fit for candidate in scan.candidates),
        *(candidate.burst_fit for candidate in scan.candidates),
    ):
        assert fit.converged
        assert fit.state_error <= 1e-7
        assert np.isfinite(fit.log_likelihood)


@pytest.mark.slow
@pytest.mark.study
def test_no_wgd_search_bootstrap_does_not_support_genome_event():
    model = _simulate(80, [[0.12, 0.22]], 1.25, 2901)
    observed = _scan(model)
    _assert_verified(observed)
    assert {candidate.node for candidate in observed.candidates} == {1, 2, 3, 4}
    calibrated = calibrate_scan(
        model, observed, 59, 8101, fractions=(0.5,), max_states=128
    )
    statistics = np.asarray(calibrated.bootstrap_statistics)
    assert np.ptp(statistics) > 0.5
    assert all(candidate.p_value > 0.05 for candidate in calibrated.candidates)
    _assert_calibration_metadata(calibrated, 59)


def _assert_calibration_metadata(calibrated, draws):
    statistics = np.asarray(calibrated.bootstrap_statistics)
    assert len(statistics) == draws
    assert np.all(np.isfinite(statistics))
    assert calibrated.calibration == "plugin-parametric-bootstrap-search-maximum"
    for candidate in calibrated.candidates:
        exceedances = np.count_nonzero(statistics >= candidate.improvement - 1e-9)
        assert candidate.p_value == pytest.approx((1 + exceedances) / (draws + 1))
        assert candidate.p_value_mc_se == pytest.approx(
            np.sqrt(candidate.p_value * (1 - candidate.p_value) / (draws + 1))
        )


def test_search_bootstrap_refits_all_branches_and_reports_calibration():
    # Three draws exercise real null/candidate refits; the CLI test checks seed
    # replay. This verifies the calculation, not rejection at .05 or error rates.
    model = _simulate(80, [[0.12, 0.22]], 1.25, 2901)
    observed = _scan(model)
    _assert_verified(observed)
    assert {candidate.node for candidate in observed.candidates} == {1, 2, 3, 4}
    calibrated = calibrate_scan(model, observed, 3, 8101, max_states=128)
    _assert_calibration_metadata(calibrated, 3)
    assert np.ptp(calibrated.bootstrap_statistics) > 0


def test_no_wgd_null_fit_is_not_worse_than_known_feasible_parameters():
    model = _simulate(80, [[0.12, 0.22]], 2.0, 2905)
    fit = fit_counts(model, max_states=128)
    assert fit.converged
    # Compare the same truncated objective, not different state approximations.
    feasible = model.log_likelihood([[0.12, 0.22]], 2.0, fit.max_count)
    assert fit.log_likelihood >= feasible - 1e-6


def test_ssd_branch_burst_is_not_count_support_for_wgd():
    # Only AB has elevated continuous SSD; no retained multiplication is present.
    model = _simulate(
        100,
        [[0.08, 0.15], [0.9, 0.15]],
        1.2,
        2902,
        branch_groups=(0, 1, 0, 0, 0),
    )
    scan = _scan(model)
    _assert_verified(scan)
    candidate = scan.candidates[0]
    assert candidate.node == 1
    # The homogeneous null is misspecified: LR alone would be misleading here.
    assert candidate.improvement > 20
    assert candidate.burst_aic_difference < -10
    assert candidate.burst_fit.log_likelihood > candidate.event_fit.log_likelihood
    assert candidate.p_value is None


def test_strong_loss_missingness_scan_remains_finite():
    model = _simulate(
        80,
        [[0.12, 1.1]],
        1.4,
        2903,
        detection=np.array([0.65, 0.8, 0.9]),
        missing=True,
    )
    assert np.count_nonzero(np.isnan(model.counts)) == 40
    assert np.count_nonzero(model.counts == 0) > 40
    assert np.count_nonzero(model.counts > 1) > 0
    scan = _scan(model)
    _assert_verified(scan)


def test_retained_wgd_wins_against_continuous_ssd_with_incomplete_detection():
    model = _simulate(
        160,
        [[0.1, 0.2]],
        1.25,
        2904,
        detection=np.array([0.9, 0.85, 0.95]),
        event=MultiplicationEvent(1, 0.85),
    )
    scan = _scan(model)
    _assert_verified(scan)
    candidate = scan.candidates[0]
    assert candidate.node == 1
    assert 0.5 < candidate.event_fit.event.retention < 0.99
    assert candidate.improvement > 50
    assert candidate.burst_aic_difference > 20
    # Recovery and an AIC comparison are not calibrated significance or power.
    assert candidate.p_value is None


@pytest.mark.parametrize("max_count", [8, 64])
def test_pure_loss_survival_is_not_rounded_to_absence(max_count):
    transition = birth_death_transition(0.0, 40.0, 1.0, max_count)
    # Independent analytic reference; no truncation or simulation is involved.
    assert transition[1, 1] == pytest.approx(np.exp(-40.0), rel=1e-12, abs=0)
    tree = CountTree(
        (-1, 0, 0),
        (0.0, 1.0, 1.0),
        (1, 2),
        ("A", "B"),
        (0, 1, 2),
        ("AB", "A", "B"),
    )
    model = CountLikelihood(tree, np.array([[1.0, 0.0]]))
    survival = np.exp(-40.0)
    expected = np.log1p(-survival) - np.log(2 - survival)
    assert model.log_likelihood([[0.0, 40.0]], 1.0, max_count) == pytest.approx(
        expected, abs=1e-12
    )
