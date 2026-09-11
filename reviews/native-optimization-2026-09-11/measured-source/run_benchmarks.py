"""Run identical harnesses sequentially against frozen before/after sources."""

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--before-results", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[2]
    sources = {"before": args.before.resolve(), "after": args.after.resolve()}
    env = dict(os.environ)
    for key in (
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        env[key] = "1"
    manifest = {
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "commands": [],
        "source_sha256": {
            label: {
                str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted((source / "nwkit").glob("*.py"))
            }
            for label, source in sources.items()
        },
    }
    if args.before_results:
        prior = json.loads((args.before_results / "manifest.json").read_text())
        assert prior["source_sha256"]["before"] == manifest["source_sha256"]["before"]
        manifest["baseline_measurements_from"] = str(args.before_results.resolve())

    def run(label, tool, options, output):
        command = [
            sys.executable,
            str(repo / "tools" / tool),
            *options,
            "--output",
            str(output.resolve()),
        ]
        manifest["commands"].append(
            {"label": label, "argv": command, "PYTHONPATH": str(sources[label])}
        )
        completed = subprocess.run(
            command,
            env={**env, "PYTHONPATH": str(sources[label])},
            capture_output=True,
            text=True,
            check=True,
        )
        print(label, output.name, completed.stdout.strip(), flush=True)
        return json.loads(output.read_text())

    comparisons = {}
    search_cases = [
        ("balanced128-2-est", "balanced", 128, 2, True),
        ("pectinate128-2-est", "pectinate", 128, 2, True),
        ("balanced512-4-est", "balanced", 512, 4, True),
        ("balanced1000-1-fixed", "balanced", 1000, 1, False),
    ]
    for name, shape, tips, traits, estimate in search_cases:
        options = [
            "--tips",
            str(tips),
            "--traits",
            str(traits),
            "--shape",
            shape,
            "--shifts",
            "3",
            "--candidate-pool",
            "24",
            "--refit-budget",
            "24",
            "--screening-budget",
            "1000",
        ]
        if estimate:
            options.append("--fit-covariance")
        results = {label: [] for label in sources}
        for repeat in range(3):
            for label in sources if repeat % 2 == 0 else reversed(sources):
                filename = f"{name}-{label}-{repeat}.json"
                if label == "before" and args.before_results:
                    reference = args.before_results / filename
                    shutil.copyfile(reference, args.output / filename)
                    results[label].append(json.loads(reference.read_text()))
                    continue
                results[label].append(
                    run(
                        label,
                        "benchmark_native_search.py",
                        options,
                        args.output / filename,
                    )
                )
        fields = (
            "input_sha256",
            "best_log_likelihood",
            "candidates",
            "metadata",
            "selected_shifts_without_penalty_or_calibration",
            "largest_fitted_shift_count",
        )
        for result in results["before"] + results["after"]:
            for field in fields:
                assert result[field] == results["before"][0][field], (name, field)
        comparisons[name] = {
            "all_candidate_records_and_search_metadata_exactly_equal": True,
            **{
                label: {
                    "median_wall_seconds": float(
                        np.median([r["wall_seconds"] for r in records])
                    ),
                    "peak_rss_bytes": [r["peak_rss_bytes"] for r in records],
                    "wall_seconds": [r["wall_seconds"] for r in records],
                }
                for label, records in results.items()
            },
        }
    with tempfile.TemporaryDirectory(prefix="nwkit-screen-matrices-") as temporary:
        for shape in ("balanced", "pectinate"):
            for increments in (False, True):
                name = (
                    f"screen-{shape}-{'increments' if increments else 'standardized'}"
                )
                results = {label: [] for label in sources}
                maximum = 0.0
                for repeat in range(3):
                    for label in sources if repeat % 2 == 0 else reversed(sources):
                        output = Path(temporary) / f"{name}-{label}-{repeat}.json"
                        options = ["--shape", shape, "--tips", "1000", "--traits", "4"]
                        if increments:
                            options.append("--optimum-increments")
                        results[label].append(
                            run(label, "benchmark_native_screen.py", options, output)
                        )
                        (args.output / output.name).write_text(output.read_text())
                    paths = [
                        Path(temporary) / f"{name}-{label}-{repeat}.npz"
                        for label in sources
                    ]
                    with np.load(paths[0]) as before, np.load(paths[1]) as after:
                        for key in before.files:
                            np.testing.assert_allclose(
                                before[key], after[key], rtol=1e-12, atol=1e-12
                            )
                            maximum = max(
                                maximum, float(np.max(np.abs(before[key] - after[key])))
                            )
                comparisons[name] = {"maximum_absolute_difference": maximum, **results}
    (args.output / "equivalence.json").write_text(
        json.dumps(comparisons, indent=2) + "\n"
    )
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
