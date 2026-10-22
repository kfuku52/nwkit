"""Untimed call-path and retained post-factor finite-check audit."""

import json

import benchmark as bench


def count_conditioner():
    calls = 0

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return bench.CANDIDATE(*args, **kwargs)

    inputs = bench.fixture(64, 4)
    bench.call(counted, inputs)
    return calls


def count_fit():
    from nwkit.vector_fit import fit_pruning_mvbm

    process, observed, errors = bench.fixture(8, 2)
    standard_errors = {
        name: bench.np.sqrt(bench.np.diag(covariance)).tolist()
        for name, covariance in errors.items()
    }
    calls = 0

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return bench.CANDIDATE(*args, **kwargs)

    with bench.patched(counted):
        result = fit_pruning_mvbm(
            process.tree, observed, ("x", "y"), standard_errors=standard_errors
        )
    return {"propagate_calls": calls, "fit_status": result[1].fit_status}


def finite_factor_failure():
    original = bench.vg.cho_factor

    def broken_factor(*args, **kwargs):
        return bench.np.array([[bench.np.nan, 0.0], [0.0, 1.0]]), True

    def failure(function):
        try:
            function(bench.np.eye(2))
        except Exception as exc:
            return {"type": type(exc).__name__, "message": str(exc)}
        raise AssertionError("Nonfinite factor was not rejected")

    bench.vg.cho_factor = broken_factor
    try:
        before = failure(lambda matrix: bench.vg._inverse(matrix, "Branch covariance"))
        after = failure(bench.validate_branch_covariance)
    finally:
        bench.vg.cho_factor = original
    assert before == after
    return before


def main():
    with bench.measurement_lock():
        result = {
            "environment": bench.environment(),
            "conditioner_64tips_4traits_propagate_calls": count_conditioner(),
            "complete_BM_fit_8tips_2traits": count_fit(),
            "retained_finite_factor_failure": finite_factor_failure(),
        }
        (bench.HERE / "call-audit.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
