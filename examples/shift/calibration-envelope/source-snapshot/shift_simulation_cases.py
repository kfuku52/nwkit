"""Independent branchwise OU generator and truth-based selection metrics."""

import math

import numpy as np
from scipy.stats import binomtest

from nwkit.util import assign_branch_ids, read_tree

SCENARIOS = ("null", "single", "distinct", "convergent")


def partition(labels):
    groups = {}
    for name, label in labels.items():
        groups.setdefault(label, []).append(name)
    return sorted(sorted(group) for group in groups.values())


def balanced_newick(tips):
    if tips < 8 or tips & (tips - 1):
        raise ValueError("tips must be a power of two >= 8")
    nodes = [f"t{i}:1" for i in range(tips)]
    while len(nodes) > 1:
        nodes = [f"({nodes[i]},{nodes[i + 1]}):1" for i in range(0, len(nodes), 2)]
    return nodes[0][:-2] + ";"


def simulate(*, tips, scenario, root_model, se, effect, seed, alpha=0.7, sigma2=0.25):
    """Draw innovations directly; no fitted model or backend simulation is used."""
    if scenario not in SCENARIOS or root_model not in ("OUfixedRoot", "OUrandomRoot"):
        raise ValueError("Unknown simulation scenario or root treatment")
    if (
        not all(math.isfinite(x) for x in (se, effect, alpha, sigma2))
        or min(alpha, sigma2, effect) <= 0
        or se < 0
    ):
        raise ValueError("Invalid generating parameters")
    newick = balanced_newick(tips)
    tree = read_tree(newick, "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    quarter = tips // 4
    targets = [
        set(f"t{i}" for i in range(start, start + quarter)) for start in (0, tips // 2)
    ]
    shifts = {}
    for node in tree.traverse():
        clade = set(node.leaf_names())
        if scenario != "null" and clade == targets[0]:
            shifts[node] = effect
        if scenario in ("distinct", "convergent") and clade == targets[1]:
            shifts[node] = -effect if scenario == "distinct" else effect
    rng = np.random.default_rng(seed)
    root_var = sigma2 / (2 * alpha) if root_model == "OUrandomRoot" else 0.0
    states, means, optima, regimes = (
        {tree: rng.normal(0, math.sqrt(root_var))},
        {tree: 0.0},
        {tree: 0.0},
        {tree: 0},
    )
    # Innovation loadings provide an independent exact covariance for generator QA.
    loadings = {tree: np.eye(1, len(ids), ids[tree]).ravel() * math.sqrt(root_var)}
    for node in tree.traverse("preorder"):
        if node.is_root:
            continue
        optima[node] = shifts.get(node, optima[node.up])
        regimes[node] = ids[node] if node in shifts else regimes[node.up]
        decay = math.exp(-alpha * node.dist)
        sd = math.sqrt(sigma2 * (-math.expm1(-2 * alpha * node.dist)) / (2 * alpha))
        offset = (1 - decay) * optima[node]
        states[node] = decay * states[node.up] + offset + rng.normal(0, sd)
        means[node] = decay * means[node.up] + offset
        loadings[node] = decay * loadings[node.up]
        loadings[node][ids[node]] += sd
    leaves = list(tree.leaves())
    matrix = np.asarray([loadings[node] for node in leaves])
    observations = [states[node] + rng.normal(0, se) for node in leaves]
    truth = {
        "seed": seed,
        "scenario": scenario,
        "tips": tips,
        "root_model": root_model,
        "standard_error": se,
        "effect": effect,
        "alpha": alpha,
        "sigma2": sigma2,
        "shift_branch_ids": sorted(ids[node] for node in shifts),
        "shift_optima": {str(ids[node]): value for node, value in shifts.items()},
        "tip_names": [node.name for node in leaves],
        "observations": observations,
        "tip_mean": {node.name: means[node] for node in leaves},
        "tip_optimum": {node.name: optima[node] for node in leaves},
        "ancestry_partition": partition({node.name: regimes[node] for node in leaves}),
        "shared_partition": partition({node.name: optima[node] for node in leaves}),
        "covariance": (matrix @ matrix.T + se**2 * np.eye(tips)).tolist(),
    }
    return newick, truth


def selection_metrics(model, truth):
    tips = model["tip_predictions"]
    parents = {row["branch_id"]: row["parent"] for row in model["branches"]}
    selected = set(model["shift_branch_ids"])
    ancestry = {}
    for row in tips:
        branch = row["branch_id"]
        while branch != 0 and branch not in selected:
            branch = parents[branch]
        ancestry[row["leaf_name"]] = branch
    result = {
        "selected_shifts": len(selected),
        "any_shift": bool(selected),
        "exact_edges": sorted(selected) == truth["shift_branch_ids"],
        "ancestry_recovered": partition(ancestry) == truth["ancestry_partition"],
        "shared_recovered": model["parameters"]["alpha"] > 0
        and partition({row["leaf_name"]: row["regime"] for row in tips})
        == truth["shared_partition"],
        "tip_mean_rmse": float(
            np.sqrt(
                np.mean(
                    [
                        (row["predicted"] - truth["tip_mean"][row["leaf_name"]]) ** 2
                        for row in tips
                    ]
                )
            )
        ),
        "alpha": model["parameters"]["alpha"],
        "any_merge": bool(model["convergence"] and model["convergence"]["merges"]),
        "backend_version": model["backend_version"],
    }
    boot = model.get("bootstrap")
    if boot:
        result.update(
            {
                f"bootstrap_{key}": boot[key]
                for key in ("attempted", "successful", "failed")
            }
        )
        for field, target in (
            ("tip_partition_frequencies", "ancestry_partition"),
            ("shared_optimum_partition_frequencies", "shared_partition"),
        ):
            if field in boot:
                result[f"bootstrap_truth_{target}_frequency"] = sum(
                    row["frequency"]
                    for row in boot[field]
                    if row["groups"] == truth[target]
                )
    return result


def rate(count, denominator):
    if not denominator:
        return {"count": count, "denominator": 0, "rate": None, "wilson_95": None}
    ci = binomtest(count, denominator).proportion_ci(method="wilson")
    return {
        "count": count,
        "denominator": denominator,
        "rate": count / denominator,
        "wilson_95": [ci.low, ci.high],
    }


def summarize(rows):
    good = [row for row in rows if row["status"] == "completed"]
    result = {
        "attempted": len(rows),
        "completed": len(good),
        "failure": rate(len(rows) - len(good), len(rows)),
    }
    for key in (
        "any_shift",
        "any_merge",
        "exact_edges",
        "ancestry_recovered",
        "shared_recovered",
    ):
        count = sum(row[key] for row in good)
        result[key] = rate(count, len(good))
        if key.endswith("recovered") or key == "exact_edges":
            result[f"{key}_returned"] = rate(count, len(rows))
    result["mean_tip_mean_rmse"] = (
        float(np.mean([row["tip_mean_rmse"] for row in good])) if good else None
    )
    result["bootstrap_counts"] = {
        key: sum(row.get(f"bootstrap_{key}", 0) for row in good)
        for key in ("attempted", "successful", "failed")
    }
    for key in ("ancestry_partition", "shared_partition"):
        values = [
            row[f"bootstrap_truth_{key}_frequency"]
            for row in good
            if f"bootstrap_truth_{key}_frequency" in row
        ]
        result[f"bootstrap_truth_{key}"] = {
            "outer_fits_available": len(values),
            "mean_frequency": float(np.mean(values)) if values else None,
        }
    return result
