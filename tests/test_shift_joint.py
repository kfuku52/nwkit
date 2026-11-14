"""Candidate-space checks and public-backend reference integration."""

import csv
import importlib.util
import itertools
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nwkit.util import assign_branch_ids, read_tree

TOOLS = Path(__file__).resolve().parents[1] / "tools"
TREE = "((A:1,B:1):1,(C:1,D:1):1);"


@pytest.fixture
def joint(monkeypatch):
    monkeypatch.syspath_prepend(str(TOOLS))
    spec = importlib.util.spec_from_file_location(
        "validate_shift_joint", TOOLS / "validate_shift_joint.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("tips,expected", [(4, 35), (8, 195), (16, 899)])
def test_candidate_counts_and_unique_keys(joint, tips, expected):
    if tips == 4:
        text = TREE
    else:
        from shift_simulation_cases import balanced_newick

        text = balanced_newick(tips)
    candidates, counts = joint.enumerate_candidates(
        read_tree(text, "auto", True, quiet=True)
    )
    assert len(candidates) == expected
    assert counts["unidentifiable_configurations"] == 1
    keys = [
        (tuple(row["shift_branch_ids"]), joint.canonical_groups(row["groups"]))
        for row in candidates
    ]
    assert len(set(keys)) == expected


def test_independent_label_assignment_enumeration(joint):
    tree = read_tree(TREE, "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    all_branches = sorted(branch for branch in ids.values() if branch)
    expected = set()
    for size in range(3):
        for selected in itertools.combinations(all_branches, size):
            for colors in itertools.product(range(size + 1), repeat=size):
                aliases = dict(zip(selected, colors, strict=True))
                aliases[0] = 0
                labels, ancestral = {}, {}
                valid = True
                for node in tree.traverse("preorder"):
                    branch = ids[node]
                    labels[node] = aliases.get(branch, labels.get(node.up))
                    ancestral[node] = (
                        branch if branch in aliases else ancestral[node.up]
                    )
                    if (
                        not node.is_root
                        and branch in aliases
                        and labels[node] == labels[node.up]
                    ):
                        valid = False
                if len({ancestral[node] for node in tree.leaves()}) != size + 1:
                    valid = False
                if valid:
                    groups = {}
                    for branch, color in aliases.items():
                        groups.setdefault(color, []).append(branch)
                    expected.add((selected, joint.canonical_groups(groups.values())))
    candidates, _ = joint.enumerate_candidates(tree)
    actual = {
        (tuple(row["shift_branch_ids"]), joint.canonical_groups(row["groups"]))
        for row in candidates
    }
    assert actual == expected


def test_disconnected_convergence_and_nested_return_are_kept(joint):
    tree = read_tree(TREE, "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    leaves = {node.name: node for node in tree.leaves()}
    a, c = ids[leaves["A"]], ids[leaves["C"]]
    parent = ids[leaves["A"].up]
    candidates, _ = joint.enumerate_candidates(tree)
    keys = {
        (tuple(row["shift_branch_ids"]), joint.canonical_groups(row["groups"]))
        for row in candidates
    }
    assert ((a, c), joint.canonical_groups([[0], [a, c]])) in keys
    assert ((parent, a), joint.canonical_groups([[0, a], [parent]])) in keys
    assert joint.tip_groups(tree, [a, c], [[0], [a, c]]) == (("A", "C"), ("B", "D"))


@pytest.mark.parametrize(
    "tree",
    ["(A:1,B:1,C:1);", "((A:1,B:2):1,(C:1,D:1):1);", "((A:0,B:1):1,(C:1,D:1):1);"],
)
def test_invalid_trees_rejected(joint, tree):
    with pytest.raises(ValueError):
        joint.enumerate_candidates(read_tree(tree, "auto", True, quiet=True))


def test_limit_rejected_before_fitting(joint):
    with pytest.raises(ValueError, match="limit"):
        joint.enumerate_candidates(read_tree(TREE, "auto", True, quiet=True), limit=1)


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"),
    reason="Set NWKIT_TEST_RSCRIPT for real joint enumeration",
)
@pytest.mark.parametrize("root_model", ["OUfixedRoot", "OUrandomRoot"])
@pytest.mark.parametrize("se", [0.0, 0.2])
def test_real_joint_reference_covers_candidates_and_matches_density(
    joint, tmp_path, root_model, se
):
    case = tmp_path / "case"
    case.mkdir()
    (case / "tree.nwk").write_text(TREE)
    (case / "traits.tsv").write_text(
        "leaf_name\tvalue\tse\n"
        + "".join(
            f"{name}\t{value}\t{se}\n"
            for name, value in zip("ABCD", [1, 1.2, 3, 2.7], strict=True)
        )
    )
    output = tmp_path / "result"
    options = SimpleNamespace(
        case=case,
        output=output,
        rscript=os.environ["NWKIT_TEST_RSCRIPT"],
        criterion="BIC",
        root_model=root_model,
        max_shifts=1,
        candidate_limit=100,
    )
    joint.run(options)
    summary = json.loads((output / "summary.json").read_text())
    assert summary["score_audit"]["comparison_validated"]
    assert summary["attempted"] == 7
    assert summary["successful"] + summary["failed"] == 7
    assert summary["score_improvement"] >= -1e-5
    assert summary["same_candidate_refit_score"] is not None
    from summarize_shift_joint import verify_fit

    with pytest.raises(ValueError, match="comparison summary"):
        verify_fit(output, {**summary, "score_improvement": 123456})
    np.testing.assert_allclose(
        summary["joint"]["reference_log_likelihood"],
        summary["joint"]["log_likelihood"],
        rtol=1e-6,
    )


@pytest.mark.parametrize(
    "penalty_gap,anchor,expected", [(0, 12, True), (69, 12, False), (0, 13, False)]
)
def test_score_audit_rejects_incomparable_penalties_and_stale_baseline(
    joint, tmp_path, penalty_gap, anchor, expected
):
    (tmp_path / "baseline-refit.tsv").write_text(f"score\terror\n{anchor}\t\n")
    results = [
        {
            "status": "completed",
            "score": 12 + penalty_gap,
            "free_score": 14,
            "log_likelihood": -3,
            "free_log_likelihood": -4,
        }
    ]
    audit = joint._score_audit(tmp_path, results, {"parameters": {"score": 12}})
    assert audit["comparison_validated"] == expected
    assert audit["maximum_information_penalty_gap"] == penalty_gap


def test_candidate_ledger_preserves_quoted_warnings_and_failures(joint, tmp_path):
    fields = [
        "candidate_id",
        "status",
        "error",
        "warnings",
        "score",
        "log_likelihood",
        "alpha",
        "sigma2",
        "free_score",
        "free_log_likelihood",
    ]
    warning = 'bounds "upper.bound"\nsecond line'
    with (tmp_path / "candidate-results.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(fields)
        writer.writerow([0, "completed", "", warning, 12, -3, 1, 1, 12, -3])
        writer.writerow(
            [1, "failed", "singular covariance", "", "NA", "NA", "NA", "NA", "NA", "NA"]
        )
    results = joint._collect(tmp_path, [{"candidate_id": 0}, {"candidate_id": 1}])
    assert results[0]["warnings"] == warning
    assert results[1]["error"] == "singular covariance"
    assert "score" not in results[1]
    with pytest.raises(ValueError, match="missing, duplicated or reordered"):
        joint._collect(tmp_path, [{"candidate_id": 0}])


def test_joint_export_rejects_missing_and_duplicate_jobs(joint):
    from summarize_shift_joint import verify_job_coverage

    manifest = {"grid": [[8, "convergent", "OUfixedRoot", 0]]}
    rows = [
        {"case": "case-00", "criterion": criterion} for criterion in ("BIC", "pBIC")
    ]
    verify_job_coverage(rows, manifest)
    for changed in (rows[:1], rows + rows[:1]):
        with pytest.raises(ValueError, match="pilot jobs"):
            verify_job_coverage(changed, manifest)
