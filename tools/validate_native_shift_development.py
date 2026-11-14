"""Frozen small development smoke cases; not a production-adoption study."""

import argparse
import hashlib
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from benchmark_native_shift import balanced_tree

from nwkit.shift_native_bootstrap import calibrate_native_search, simulate_native_data
from nwkit.shift_native_fit import NativeFitOptions
from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_model import ShiftData, ShiftLayout, evaluate_trait
from nwkit.shift_native_provenance import native_implementation
from nwkit.util import read_tree

CASES = [
    {
        "name": "null_complete",
        "shifts": [],
        "groups": [[0]],
        "known_error": 0.0,
        "missing": False,
    },
    {
        "name": "null_known_error_missing",
        "shifts": [],
        "groups": [[0]],
        "known_error": 0.05,
        "missing": True,
    },
    {
        "name": "single_terminal_shift",
        "shifts": [7],
        "groups": [[0], [7]],
        "known_error": 0.0,
        "missing": False,
    },
    {
        "name": "shared_terminal_shifts",
        "shifts": [7, 14],
        "groups": [[0], [7, 14]],
        "known_error": 0.05,
        "missing": False,
    },
]


def generate(case, seed):
    tree = read_tree(balanced_tree(8), "auto", True, quiet=True)
    values = np.random.default_rng(1).normal(size=(8, 2))
    if case["missing"]:
        values[2, 1] = np.nan
    errors = case["known_error"] * np.arange(1, 9)[:, None] * np.ones((1, 2))
    data = ShiftData.build(tree, values, ["x", "y"], errors)
    layout = ShiftLayout.build(data.tree, case["shifts"], case["groups"])
    fits = []
    for trait, alpha in enumerate([0.7, 1.4]):
        fit = evaluate_trait(
            data,
            layout,
            trait,
            alpha,
            1 / data.scales[trait] ** 2,
            0.1 / data.scales[trait] ** 2,
        )
        coefficients = np.full(len(layout.groups), 8.0 / data.scales[trait])
        coefficients[0] = -data.centers[trait] / data.scales[trait]
        fits.append(
            replace(
                fit,
                coefficients=coefficients,
                predicted=layout.design(data.tree, alpha) @ coefficients,
            )
        )
    return simulate_native_data(
        data, {"fits": fits}, np.random.default_rng(seed)
    ), layout


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--replicates", type=int, default=19)
    parser.add_argument("--cases", nargs="+", choices=[case["name"] for case in CASES])
    args = parser.parse_args(argv)
    destination = Path(args.output)
    destination.mkdir(exist_ok=False, parents=True)
    specification = {
        "scope": "development_smoke_only; not_independent_confirmation_or_power_noninferiority",
        "seed": 2026091002,
        "cases": [
            case for case in CASES if args.cases is None or case["name"] in args.cases
        ],
        "datasets_per_case": 1,
        "alpha_height_truth": [0.7, 1.4],
        "process_tip_variance_truth": 1.0,
        "additional_measurement_variance_truth": 0.1,
        "scaled_shift_offset_truth": 8.0,
        "replicates": args.replicates,
        "level": 0.05,
        "search": {
            "max_shifts": 2,
            "convergence": True,
            "candidate_pool": 14,
            "refit_budget": 12,
        },
        "fit": {"root_model": "OUfixedRoot", "estimate_measurement_error": True},
        "implementation": native_implementation(),
        "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (destination / "specification.json").write_text(
        json.dumps(specification, indent=2, allow_nan=False) + "\n"
    )
    results = []
    for index, case in enumerate(CASES):
        if args.cases is not None and case["name"] not in args.cases:
            continue
        data, truth = generate(case, specification["seed"] + index)
        options = NativeSearchOptions(**specification["search"])

        def search(sample, options=options):
            return heuristic_native_search(
                sample,
                options=options,
                fit_arguments={"options": NativeFitOptions(**specification["fit"])},
            )

        started = time.perf_counter()
        record = {
            "case": case["name"],
            "truth_shifts": list(truth.shifts),
            "truth_groups": [list(g) for g in truth.groups],
        }
        try:
            searched = search(data)
            selected, calibration = calibrate_native_search(
                data,
                searched,
                search,
                replicates=args.replicates,
                seed=specification["seed"] + 100 + index,
            )
            record.update(
                {
                    "status": "complete",
                    "selected_shifts": list(selected["layout"].shifts),
                    "selected_groups": [list(g) for g in selected["layout"].groups],
                    "any_shift_selected": bool(selected["layout"].shifts),
                    "exact_location_recovery": selected["layout"].shifts
                    == truth.shifts,
                    "exact_layout_recovery": selected["layout"] == truth,
                    "calibration": calibration,
                }
            )
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            record.update({"status": "failed", "error": str(exc)})
        record["wall_seconds"] = time.perf_counter() - started
        results.append(record)
        (destination / "results.json").write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n"
        )
        print(
            json.dumps(
                {key: value for key, value in record.items() if key != "calibration"}
            ),
            flush=True,
        )
    (destination / "status.json").write_text(
        json.dumps(
            {
                "completed_cases": sum(r["status"] == "complete" for r in results),
                "failed_cases": sum(r["status"] == "failed" for r in results),
                "production_adoption": "not_evaluated; development_only; no_baseline_noninferiority_or_confirmatory_sample",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
