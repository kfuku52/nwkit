"""Post-hoc independent null-model ML profile check; never changes selections.

For zero measurement error, profile the intercept and scalar variance exactly,
then check alpha on a log grid and refine every grid-local maximum. This is a
diagnostic of the fitted null likelihood, not a new selection rule or certificate.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar
from shift_alpha_audit import read_table

from nwkit.util import read_tree


def profile(tree, observations, root_model, lower, upper):
    leaves = list(tree.leaves())
    heights = np.array([tree.get_distance(tree, node) for node in leaves])
    total = heights[:, None] + heights[None, :]
    distance = np.array([[tree.get_distance(a, b) for b in leaves] for a in leaves])
    shared = (total - distance) / 2
    ones = np.ones(len(leaves))
    observed = np.asarray(observations)

    def evaluate(log_alpha):
        alpha = math.exp(log_alpha)
        covariance = np.exp(-alpha * distance) / (2 * alpha)
        if root_model == "OUfixedRoot":
            covariance *= -np.expm1(-2 * alpha * shared)
        factor = cho_factor(covariance)
        precision_one = cho_solve(factor, ones)
        intercept = (precision_one @ observed) / (precision_one @ ones)
        residual = observed - intercept
        sigma2 = residual @ cho_solve(factor, residual) / len(leaves)
        objective = 0.5 * (
            len(leaves) * (math.log(2 * math.pi * sigma2) + 1)
            + 2 * np.log(np.diag(factor[0])).sum()
        )
        return objective, float(intercept), float(sigma2)

    def negative_ll(log_alpha):
        return evaluate(log_alpha)[0]

    grid = np.linspace(math.log(lower), math.log(upper), 161)
    values = [negative_ll(x) for x in grid]
    candidates = [(values[0], grid[0]), (values[-1], grid[-1])]
    for i in range(1, len(grid) - 1):
        if values[i] <= min(values[i - 1], values[i + 1]):
            fit = minimize_scalar(
                negative_ll, bounds=(grid[i - 1], grid[i + 1]), method="bounded"
            )
            candidates.append((float(fit.fun), float(fit.x)))
    objective, optimum = min(candidates)
    _, intercept, sigma2 = evaluate(optimum)
    return {
        "profile_log_likelihood": -objective,
        "profile_alpha": math.exp(optimum),
        "profile_intercept": intercept,
        "profile_sigma2": sigma2,
    }


def main(root, output):
    rows = []
    for folder in sorted(root.glob("c*-r*")):
        case = json.loads((folder / "case.json").read_text())
        if case["family"] != "primary" or case["scenario"] != "null":
            continue
        truth = json.loads((folder / "truth.json").read_text())
        tree = read_tree(str(folder / "tree.nwk"), "auto", True, quiet=True)
        settings = read_table(folder / "settings.tsv")
        ledger = read_table(folder / "candidates-results.tsv")
        for setting in settings:
            backend = next(
                r
                for r in ledger
                if r["floor_id"] == setting["floor_id"] and r["candidate_id"] == "0"
            )
            result = profile(
                tree,
                truth["observations"],
                case["root_model"],
                float(setting["lower"]),
                float(setting["upper"]),
            )
            rows.append(
                {
                    "case_id": case["case_id"],
                    "root_model": case["root_model"],
                    "floor_id": setting["floor_id"],
                    **result,
                    "backend_log_likelihood": float(backend["log_likelihood"]),
                    "profile_minus_backend": result["profile_log_likelihood"]
                    - float(backend["log_likelihood"]),
                }
            )
    output.write_text(
        json.dumps(
            {
                "post_hoc": True,
                "grid_points": 161,
                "records": rows,
                "maximum_log_likelihood_improvement": max(
                    r["profile_minus_backend"] for r in rows
                ),
            },
            indent=2,
        )
        + "\n"
    )
    print(output.read_text()[-120:])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.input, args.output)
