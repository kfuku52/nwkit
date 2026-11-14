"""Design invariants and paired-denominator checks for the OU sensitivity study."""

import importlib
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def alpha_modules(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "tools"))
    return (
        importlib.import_module("shift_alpha_design"),
        importlib.import_module("summarize_shift_alpha"),
    )


def test_protocol_has_independent_seeds_and_prespecified_denominators(alpha_modules):
    design, _ = alpha_modules
    specification = design.protocol()
    cases = list(design.cases(specification))
    assert len(cases) == len({r["seed"] for r in cases}) == 520
    assert sum(r["family"] == "primary" for r in cases) == 400
    assert sum(r["tips"] == 16 for r in cases) == 40
    assert cases == list(design.cases(specification))
    assert all(
        0.1 < r["alpha_height"] < specification["alpha_upper_height"] for r in cases
    )


def test_dimensionless_generator_mean_and_covariance(alpha_modules):
    from nwkit.shift_reference import evaluate_shift_model
    from nwkit.util import assign_branch_ids, read_tree

    design, _ = alpha_modules
    specification = design.protocol(1, 1)
    for case in design.cases(specification):
        newick, truth = design.generate_case(case)
        tree = read_tree(newick, "auto", True, quiet=True)
        branches = assign_branch_ids(tree)
        heights = {
            node: max(tree.get_distance(node, tip) for tip in node.leaves())
            for node in tree.traverse()
        }
        optima, effects = {tree: 0.0}, {}
        for node in tree.traverse("preorder"):
            if node.is_root:
                continue
            optima[node] = truth["shift_optima"].get(
                str(branches[node]), optima[node.up]
            )
            if branches[node] in truth["shift_branch_ids"]:
                effects[node] = (optima[node] - optima[node.up]) * -np.expm1(
                    -truth["alpha"] * heights[node.up]
                )
        mean, covariance, _ = evaluate_shift_model(
            tree,
            observations=truth["observations"],
            standard_errors=[case["standard_error"]] * case["tips"],
            alpha=truth["alpha"],
            sigma2=truth["sigma2"],
            intercept=0.0,
            mean_effects=effects,
            root_model=case["root_model"],
        )
        np.testing.assert_allclose(mean, list(truth["tip_mean"].values()), atol=1e-12)
        np.testing.assert_allclose(covariance, truth["covariance"], atol=1e-12)


def test_paired_counts_exclude_failed_pair_without_losing_attempt(alpha_modules):
    _, summary = alpha_modules
    rows = []
    for i in range(3):
        for floor in ("small", "raised"):
            rows.append(
                {
                    "case_id": str(i),
                    "floor_id": floor,
                    "method": "joint",
                    "criterion": "BIC",
                    "status": "failed" if i == 2 and floor == "raised" else "completed",
                    "shared_partition": [["a", "b"]]
                    if floor == "small"
                    else [["a"], ["b"]],
                    "any_shift": floor == "raised",
                    "shared_recovered": floor == "small",
                    "exact_edges": False,
                    "tip_mean_rmse": 1.0 if floor == "small" else 2.0,
                }
            )
    result = summary.paired(rows, "floor_id", "small", "raised")
    assert result["attempted"] == 3 and result["completed"] == 2
    assert (
        result["partition_changed"]["count"]
        == result["partition_changed"]["denominator"]
        == 2
    )
    assert result["any_shift"] == {"left_only": 0, "right_only": 2}
    assert result["mean_rmse_right_minus_left"] == 1.0


def test_search_comparison_excludes_incomplete_or_unmatched_fits(alpha_modules):
    _, summary = alpha_modules
    rows = []
    for i in range(4):
        for method in ("two_stage", "joint"):
            rows.append(
                {
                    "case_id": str(i),
                    "floor_id": "small",
                    "method": method,
                    "criterion": "pBIC",
                    "status": "completed",
                    "shared_partition": [["a"]],
                    "any_shift": False,
                    "shared_recovered": True,
                    "exact_edges": True,
                    "tip_mean_rmse": 0.0,
                    "audit_passed": True,
                    "candidate_failed": int(i == 1),
                    "baseline_in_candidate_set": i != 2,
                    "same_model_refit_score_gap": 0.01 if i == 3 else 0.0,
                    "score": 10.0 if method == "two_stage" else 9.0,
                }
            )
    result = summary.paired(rows, "method", "two_stage", "joint")
    assert result["completed"] == 4
    assert result["search_score_comparison"] == {
        "eligible": 1,
        "excluded": 3,
        "joint_lower": 1,
        "tied": 0,
        "joint_higher": 0,
    }


@pytest.mark.parametrize("root", ["OUfixedRoot", "OUrandomRoot"])
def test_null_profile_matches_independent_tree_likelihood(alpha_modules, root):
    from shift_alpha_null_diagnostic import profile

    from nwkit.shift_reference import evaluate_shift_model
    from nwkit.util import read_tree

    design, _ = alpha_modules
    case = next(
        c for c in design.cases(design.protocol(1, 0)) if c["root_model"] == root
    )
    newick, truth = design.generate_case(case)
    tree = read_tree(newick, "auto", True, quiet=True)
    result = profile(tree, truth["observations"], root, 1e-7 / 3, 10 / 3)
    _, _, ll = evaluate_shift_model(
        tree,
        observations=truth["observations"],
        alpha=result["profile_alpha"],
        sigma2=result["profile_sigma2"],
        intercept=result["profile_intercept"],
        mean_effects={},
        root_model=root,
    )
    assert ll == pytest.approx(result["profile_log_likelihood"], abs=1e-7)
