"""Predeclared high-dimensional and complete-fit follow-up for candidate 06."""

import json
from dataclasses import fields, is_dataclass

import benchmark as bench


def compare(left, right, path="result"):
    if is_dataclass(left):
        for field in fields(left):
            compare(
                getattr(left, field.name),
                getattr(right, field.name),
                path + "." + field.name,
            )
    elif isinstance(left, dict):
        assert set(left) == set(right), path
        for name in left:
            compare(left[name], right[name], path + "." + str(name))
    elif isinstance(left, (tuple, list)):
        assert len(left) == len(right), path
        for index, (one, two) in enumerate(zip(left, right, strict=True)):
            compare(one, two, path + "." + str(index))
    elif isinstance(left, (bench.np.ndarray, bench.np.number, float)):
        bench.np.testing.assert_array_equal(left, right, err_msg=path)
    else:
        assert left == right, (path, left, right)


def condition_case(tips, dimension):
    inputs = bench.fixture(tips, dimension, mode="noisy")

    def baseline():
        return bench.call(bench.ORIGINAL, inputs)

    def candidate():
        return bench.call(bench.CANDIDATE, inputs)

    return {
        "case": f"{tips}tips-{dimension}traits-noisy-proper",
        "scope": "complete condition_vector_tree; fixture excluded",
        "comparison": bench.compare_results(baseline(), candidate()),
        "timing": bench.paired(baseline, candidate, minimum_seconds=0.15),
        "peak_bytes": {
            "baseline": bench.peak_bytes(baseline),
            "candidate": bench.peak_bytes(candidate),
        },
    }


def complete_fit():
    from nwkit.vector_fit import fit_pruning_mvbm

    process, observed, errors = bench.fixture(8, 2)
    standard_errors = {
        name: bench.np.sqrt(bench.np.diag(covariance)).tolist()
        for name, covariance in errors.items()
    }

    def baseline():
        with bench.patched(bench.ORIGINAL):
            return fit_pruning_mvbm(
                process.tree, observed, ("x", "y"), standard_errors=standard_errors
            )

    def candidate():
        with bench.patched(bench.CANDIDATE):
            return fit_pruning_mvbm(
                process.tree, observed, ("x", "y"), standard_errors=standard_errors
            )

    first, second = baseline(), candidate()
    compare(first, second)
    return {
        "case": "8tips-2traits-partial-noisy-complete-BM-fit",
        "scope": "complete fit_pruning_mvbm including preparation, deterministic multistart optimization and all posterior fields",
        "comparison": {
            "every_fit_and_posterior_field_exact": True,
            "fit_fields": [field.name for field in fields(first[1])],
            "posterior_nodes": len(first[0]),
            "status": first[1].fit_status,
        },
        "timing": bench.paired(baseline, candidate, minimum_seconds=0.15),
        "peak_bytes": {
            "baseline": bench.peak_bytes(baseline),
            "candidate": bench.peak_bytes(candidate),
        },
    }


def main():
    with bench.measurement_lock():
        results = {
            "environment": bench.environment(),
            "reason": "Initial 2-64 trait effects overlapped host noise; test O(d^3) unused solve at 128/256 traits and an actual complete fit without weakening checks.",
            "cases": [condition_case(16, 128), condition_case(16, 256), complete_fit()],
        }
        (bench.HERE / "additional-results.json").write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n"
        )
        for case in results["cases"]:
            print(
                case["case"], case["timing"]["speedup"], case["peak_bytes"], flush=True
            )


if __name__ == "__main__":
    main()
