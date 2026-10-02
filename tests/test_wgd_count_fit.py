import numpy as np
import pytest

from nwkit.wgd_count_fit import calibrate_scan, fit_counts, scan_counts
from nwkit.wgd_count_model import CountLikelihood, CountTree, MultiplicationEvent


def model(counts):
    tree = CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 0.5, 0.5, 0.5, 1.0),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )
    return CountLikelihood(tree, np.asarray(counts))


def test_coordinated_multiplication_beats_stochastic_branch_burst():
    result = scan_counts(
        model(np.tile([2, 2, 1], (40, 1))), nodes=(1,), fractions=(0.5,), max_states=64
    )
    candidate = result.candidates[0]
    assert candidate.event_fit.event.retention > 0.99
    assert candidate.event_fit.root_mean == pytest.approx(1, abs=1e-4)
    assert candidate.improvement > 100
    assert candidate.burst_aic_difference > 50
    assert candidate.p_value is None
    assert result.calibration == "not-run"
    assert candidate.event_fit.state_error < 1e-7
    assert candidate.burst_fit.log_likelihood >= result.background.log_likelihood - 1e-7


def test_known_rate_wgd_simulation_is_located_without_input_event():
    template = model(np.ones((120, 3)))
    counts = template.simulate(
        [[0.08, 0.1]], 1, np.random.default_rng(150), MultiplicationEvent(1, 0.9)
    )
    result = scan_counts(model(counts), fractions=(0.5,), max_states=64)
    assert result.candidates[0].node == 1
    assert result.candidates[0].event_fit.event.retention > 0.6
    assert result.candidates[0].improvement > 20


def test_background_fit_time_units_and_truncation_diagnostic():
    import json

    from nwkit.wgd_count import _fit_json

    counts = np.array([[1, 1, 1], [1, 2, 1], [1, 1, 0], [2, 2, 2]] * 8)
    original = model(counts)
    fit = fit_counts(original, max_states=64)
    tree = original.tree
    scaled_tree = CountTree(
        tree.parents,
        tuple(t * 1000 for t in tree.lengths),
        tree.tip_nodes,
        tree.tip_names,
        tree.branch_ids,
        tree.clade_ids,
    )
    scaled = fit_counts(CountLikelihood(scaled_tree, counts), max_states=64)
    assert scaled.log_likelihood == pytest.approx(fit.log_likelihood, abs=1e-7)
    np.testing.assert_allclose(scaled.rates * 1000, fit.rates, rtol=1e-5)
    assert fit.converged
    assert fit.state_error < 1e-7
    assert type(fit.boundary) is bool
    assert (
        json.loads(json.dumps(_fit_json(fit)))["nuisance_bound_reached"] == fit.boundary
    )


def test_calibration_uses_maximum_and_repeats_null_fit(monkeypatch):
    # Algebraic calibration check; optimization/recovery has separate real tests.
    from dataclasses import replace

    import nwkit.wgd_count_fit as fitting

    data = model(np.tile([2, 2, 1], (10, 1)))
    observed = scan_counts(data, nodes=(1, 2), fractions=(0.5,), max_states=64)
    observed = replace(
        observed,
        candidates=tuple(
            replace(c, improvement=float(i + 1))
            for i, c in enumerate(observed.candidates)
        ),
    )
    calls = []

    def fake_scan(simulated, **kwargs):
        calls.append((simulated.counts.copy(), kwargs))
        return replace(
            observed,
            candidates=(
                replace(observed.candidates[0], improvement=0.5),
                replace(observed.candidates[1], improvement=1.5),
            ),
        )

    monkeypatch.setattr(fitting, "scan_counts", fake_scan)
    result = calibrate_scan(data, observed, 9, 123, fractions=(0.5,), max_states=64)
    assert len(calls) == 9
    assert all(call[1]["compare_bursts"] is False for call in calls)
    assert all(call[1]["nodes"] == (1, 2) for call in calls)
    assert result.bootstrap_statistics == (1.5,) * 9
    assert result.candidates[0].p_value == 1
    assert result.candidates[1].p_value == 0.1
    assert result.calibration == "plugin-parametric-bootstrap-search-maximum"


def test_invalid_scan_and_bound_fail_instead_of_returning_unverified_fit():
    data = model([[1, 1, 1]])
    with pytest.raises(ValueError, match="twice"):
        fit_counts(data, max_states=8)
    with pytest.raises(ValueError, match="fractions"):
        scan_counts(data, fractions=())
    with pytest.raises(ValueError, match="Candidate"):
        scan_counts(data, nodes=(0,))
    with pytest.raises(ValueError, match="State tolerance"):
        fit_counts(data, state_tolerance=0)


def test_optimizer_success_without_stationarity_is_rejected(monkeypatch):
    from scipy.optimize import OptimizeResult

    import nwkit.wgd_count_fit as fitting

    def false_success(objective, start, **kwargs):
        return OptimizeResult(
            x=start,
            fun=objective(start),
            jac=np.zeros_like(start),
            success=True,
            message="False relative-function convergence",
            nit=0,
        )

    monkeypatch.setattr(fitting, "minimize", false_success)
    data = model(np.tile([2, 2, 1], (20, 1)))
    with pytest.raises(ValueError, match="did not converge.*projected gradient"):
        fit_counts(data, max_states=64)


def test_nonfinite_likelihood_is_not_replaced_by_flat_objective(monkeypatch):
    data = model([[1, 1, 1]])
    monkeypatch.setattr(data, "log_likelihood", lambda *args: -np.inf)
    with pytest.raises(ValueError, match="nonfinite likelihood"):
        fit_counts(data, max_states=64)


def test_richardson_gradient_removes_verified_stiff_root_step_bias():
    from nwkit.wgd_count_fit import _objective_gradient

    root = 0.00146
    parameters, bounds = np.array([root]), [(0.0, 1.0)]

    def objective(point):
        return -np.log(point[0]) + point[0] / root

    step = np.cbrt(np.finfo(float).eps)
    legacy = (objective(parameters + step) - objective(parameters - step)) / (2 * step)
    assert abs(legacy) > 1e-3
    assert abs(_objective_gradient(objective, parameters, bounds)[0]) < 1e-5
    assert (
        abs(_objective_gradient(objective, parameters, bounds, step_factor=2.0)[0])
        < 1e-5
    )


def test_gradient_subtracts_per_family_values_before_summing():
    from nwkit.wgd_count_fit import _objective_gradient

    root = 0.00146
    parameters, bounds = np.array([root + 1e-6]), [(0.0, 1.0)]

    def family_objective(point):
        return np.array([1e16, -np.log(point[0]) + point[0] / root])

    expected = 1 / root - 1 / parameters[0]
    actual = _objective_gradient(family_objective, parameters, bounds)[0]
    assert actual == pytest.approx(expected, abs=1e-6)
    scalar = _objective_gradient(
        lambda p: family_objective(p).sum(), parameters, bounds
    )
    assert scalar[0] == 0.0


@pytest.mark.parametrize("point", [0.0, 1e-6, 0.4, 1 - 1e-6, 1.0])
def test_richardson_gradient_keeps_both_stencils_inside_bounds(point):
    from nwkit.wgd_count_fit import _objective_gradient

    def objective(parameters):
        x = parameters[0]
        assert 0 <= x <= 1
        return np.array([12.0, x + x**2 + x**3])

    actual = _objective_gradient(objective, np.array([point]), [(0.0, 1.0)])[0]
    assert actual == pytest.approx(1 + 2 * point + 3 * point**2, abs=1e-8)


def test_stiff_root_profile_is_polished_without_relaxing_stationarity(monkeypatch):
    from scipy.optimize import OptimizeResult, minimize

    import nwkit.wgd_count_fit as fitting

    data = model(np.ones((500, 3)))
    optimum = 0.00146
    target = np.array([-0.5, -0.7, optimum + 2e-9])
    phases = []

    def family_likelihood(
        rates, root_mean, max_count, event=None, *, _split_zero_event=False
    ):
        parameters = np.log(np.asarray(rates).ravel() * 0.5)
        root_parameter = np.log(root_mean)
        values = np.full(500, -3.0)
        values[0] -= (
            np.sum((parameters - target[:2]) ** 2)
            - np.log(root_parameter)
            + root_parameter / optimum
        )
        return values

    def stalled_initial(objective, start, **kwargs):
        phases.append(kwargs["options"]["maxiter"])
        if kwargs["jac"] == "3-point":
            return OptimizeResult(
                x=target.copy(),
                fun=objective(target),
                success=True,
                nit=3,
                message="RELATIVE REDUCTION OF F <= FACTR*EPSMCH",
            )
        return minimize(objective, start, **kwargs)

    monkeypatch.setattr(data, "family_log_likelihoods", family_likelihood)
    monkeypatch.setattr(fitting, "minimize", stalled_initial)
    fit = fit_counts(data, max_states=32, max_iterations=40)
    root_parameter = np.log(fit.root_mean)
    assert abs(1 / optimum - 1 / root_parameter) <= 1e-5
    assert fit.converged and fit.state_error == 0.0
    assert phases == [40, 40, 37, 37]


def test_bootstrap_inherits_exact_observed_search_and_rejects_changes(monkeypatch):
    import nwkit.wgd_count_fit as fitting

    data = model(np.tile([2, 2, 1], (10, 1)))
    observed = scan_counts(
        data, nodes=(1,), fractions=(0.1, 0.8), multiplicity=3, max_states=64
    )
    calls = []

    def fake_scan(simulated, **kwargs):
        calls.append(kwargs)
        return observed

    monkeypatch.setattr(fitting, "scan_counts", fake_scan)
    calibrate_scan(data, observed, 1, 91, max_states=64)
    assert calls[0]["fractions"] == (0.1, 0.8)
    assert calls[0]["multiplicity"] == 3
    assert calls[0]["nodes"] == (1,)
    for options in ({"fractions": (0.5,)}, {"multiplicity": 2}):
        with pytest.raises(ValueError, match="must match the observed search"):
            calibrate_scan(data, observed, 1, 91, **options)


def test_zero_retention_fit_refines_hidden_intermediate_count_bound(monkeypatch):
    import nwkit.wgd_count_fit as fitting

    design = CountTree(
        (-1, 0, 0), (0.0, 1.0, 1.0), (1, 2), ("A", "B"), (0, 1, 2), ("AB", "A", "B")
    )
    data = CountLikelihood(design, np.array([[1.0, 1.0]]))
    calls = []

    def fixed_reference_fit(model, max_count, **kwargs):
        calls.append(max_count)
        rates, event = np.array([[20.0, 20.0]]), MultiplicationEvent(1, 0.0)
        return fitting.CountFit(
            rates,
            1.0,
            model.log_likelihood(rates, 1.0, max_count, event),
            4,
            max_count,
            float("inf"),
            False,
            True,
            "Fixed reference parameters",
            event,
        )

    monkeypatch.setattr(fitting, "_fit_at_bound", fixed_reference_fit)
    fit = fit_counts(data, node=1, max_states=256)
    assert calls == [8, 16, 32, 64, 128]
    assert fit.max_count == 128
    assert fit.converged and fit.state_error <= 1e-7
