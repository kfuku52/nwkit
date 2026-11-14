"""Compare scientific outputs, retaining structural and convergence decisions."""

import json
from pathlib import Path

import numpy as np


def compare(left, right, path="root"):
    if isinstance(left, dict):
        if left.keys() != right.keys():
            raise ValueError(f"Different keys at {path}")
        for key in left:
            if key in {"evaluations", "optimizer"}:
                continue
            compare(left[key], right[key], f"{path}.{key}")
    elif isinstance(left, list):
        if len(left) != len(right):
            raise ValueError(f"Different lengths at {path}")
        for i, (first, second) in enumerate(zip(left, right, strict=True)):
            compare(first, second, f"{path}[{i}]")
    elif isinstance(left, float):
        if not np.isfinite(left) or not np.isfinite(right):
            raise ValueError(f"Nonfinite value at {path}")
        np.testing.assert_allclose(left, right, rtol=1e-6, atol=1e-6, err_msg=path)
    elif left != right:
        raise ValueError(f"Mismatch at {path}: {left!r} != {right!r}")


def check(before, after):
    for key in ["status", "input_sha256", "truth", "fit_options", "search_options"]:
        compare(before[key], after[key], key)
    for key in ["fit", "candidates", "best_unpenalized"]:
        if key in before:
            compare(before[key], after[key], key)
    summary = "fit" if "fit" in before else "best_unpenalized"
    for first, second in zip(
        before[summary]["traits"], after[summary]["traits"], strict=True
    ):

        def modes(trait):
            return [
                [row["alpha_mode"], row["success"]]
                for row in trait["optimizer"]["alpha_candidates"]
            ]

        compare(modes(first), modes(second), "covariance mode decisions")
        for key in ["complete_alpha_modes", "nuisance_variance_at_numerical_bound"]:
            compare(first["optimizer"][key], second["optimizer"][key], key)
    for first, second in zip(
        before["search_calls"], after["search_calls"], strict=True
    ):
        compare(first["metadata"], second["metadata"], "search metadata")
        compare(
            first["largest_fitted_shift_count"],
            second["largest_fitted_shift_count"],
            "largest fitted shifts",
        )


if __name__ == "__main__":
    import sys

    root = Path(sys.argv[1])
    checks = []
    for path in sorted(root.glob("*-before-*.json")):
        before = json.loads(path.read_text())
        target = path.with_name(path.name.replace("-before-", "-after-"))
        if not target.exists():
            raise ValueError(f"Missing comparison: {target}")
        check(before, json.loads(target.read_text()))
        checks.append([path.name, target.name])
    if not checks:
        raise ValueError("No comparisons found")
    print(json.dumps({"passed": checks, "rtol": 1e-6, "atol": 1e-6}, indent=2))
