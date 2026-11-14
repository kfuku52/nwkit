"""Rank-threshold and failure controls for the unchanged full-matrix arithmetic."""

import json

import benchmark as study


def threshold_case(gap):
    leaves = ",".join(f"T{i}:{1e-6 + i * gap:.17g}" for i in range(8))
    tree = study.read_tree("(" + leaves + ")R;", "1", True, quiet=True, rooted="yes")
    values = {f"T{i}": [float(i), float(i % 2)] for i in range(8)}
    data = study._prepare_observations(tree, values, ("x", "y"), None)
    design, complete = study.module.observation_moment_design(data)
    return {
        "base_branch_length": 1e-6,
        "gap": gap,
        "design_entries": len(design),
        "existing_relative_singular_threshold": 1e-6,
        **study.compare(study.parameters_for(2, initial=True), 2, 1, design, complete),
    }


def outcome(function, args):
    try:
        function(*args)
    except Exception as error:
        return {"exception": type(error).__name__, "message": str(error)}
    raise AssertionError("Expected a failure")


def main():
    with study.measurement_lock():
        threshold = [threshold_case(gap) for gap in (1e-20, 5e-20, 6e-20, 1e-19)]
        failures = []
        cases = [
            ("trait_index", (study.parameters_for(2), 2, 1, ((0.1, 0.2, 2, 0),))),
            ("short_parameters", (study.np.zeros(1), 2, 1, ((0.1, 0.2, 0, 0),))),
        ]
        for name, args in cases:
            first = outcome(study.ORIGINAL, args)
            second = outcome(study.consecutive_moments, args)
            assert first == second, (first, second)
            failures.append({"case": name, "exact": True, **first})
        for function in (study.ORIGINAL, study.consecutive_moments):
            with study.patched(function):
                result = outcome(
                    study.module.full_ou_identifiability,
                    (study.parameters_for(2), 2, 1, (), True),
                )
            if function is study.ORIGINAL:
                baseline = result
            else:
                assert result == baseline
        failures.append({"case": "empty_identifiability", "exact": True, **baseline})
        output = {
            "environment": study.environment(),
            "source_path": study.module.__file__,
            "threshold_controls": threshold,
            "failures": failures,
        }
        (study.HERE / "additional-results.json").write_text(
            json.dumps(output, indent=2, allow_nan=False) + "\n"
        )
        for row in threshold:
            print(row, flush=True)
        print("failure_controls", failures, flush=True)


if __name__ == "__main__":
    main()
