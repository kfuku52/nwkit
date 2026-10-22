import numpy as np
import pytest
from ete4 import Tree

from nwkit.time_tree import (
    _credible_interval,
    _hpd_interval,
    annotate_mcmctree_calibrations,
    infer_node_ages_from_branch_lengths,
    paml_node_mapping,
    parse_mcmctree_calibration,
    posterior_sample_tree,
    prepare_time_tree_annotations,
    read_mcmctree_posterior,
    summarize_mcmctree_posterior,
)


def test_rounded_short_branches_use_whole_tree_height_without_losing_intervals():
    tree = Tree("((A:0.1000,B:0.1003):252.2315,C:252.3318);", parser=1)
    clade = tree.common_ancestor(["A", "B"])
    clade.props.update(
        age_ci_low=0.08, age_ci_high=0.12, age_ci_kind="HPD", age_ci_level=0.95
    )
    before = {node: (node.dist, dict(node.props)) for node in tree.traverse()}

    assert prepare_time_tree_annotations(tree)["credible_interval_count"] == 1

    assert clade.props["age"] == pytest.approx(0.10015)
    for node, (distance, properties) in before.items():
        assert node.dist == distance
        assert all(node.props[key] == value for key, value in properties.items())


@pytest.mark.parametrize(
    "text",
    [
        "((A:0.1,B:0.2):252.2315,C:252.3315);",
        # Every old local comparison passes, but the full path spread is too large.
        "((A:100,B:100.0009):0.0009,C:100);",
    ],
)
def test_non_ultrametric_tree_is_rejected_before_any_age_is_attached(text):
    tree = Tree(text, parser=1)
    before = {node: dict(node.props) for node in tree.traverse()}
    with pytest.raises(ValueError, match="not ultrametric"):
        infer_node_ages_from_branch_lengths(tree)
    assert all(node.props == properties for node, properties in before.items())


@pytest.mark.parametrize("distance", [-1.0, float("nan"), float("inf")])
def test_dated_ages_reject_invalid_branch_lengths(distance):
    tree = Tree("(A:1,B:1);", parser=1)
    tree["A"].dist = distance
    with pytest.raises(ValueError, match="finite and nonnegative"):
        infer_node_ages_from_branch_lengths(tree)


@pytest.mark.parametrize("tolerance", [-1.0, float("nan"), float("inf")])
def test_dated_ages_reject_invalid_tolerance(tolerance):
    with pytest.raises(ValueError, match="tolerance"):
        infer_node_ages_from_branch_lengths(Tree("(A:1,B:1);", parser=1), tolerance)


def test_explicit_ages_and_zero_length_branches_are_preserved():
    tree = Tree("((A:1,B:1):0,C:1);", parser=1)
    tree.props["age"] = 1.000001
    infer_node_ages_from_branch_lengths(tree)
    assert tree.props["age"] == 1.000001
    assert tree.common_ancestor(["A", "B"]).dist == 0.0


def test_mcmctree_interval_indices_match_paml_order_statistics():
    values = np.arange(10, dtype=float)

    assert _hpd_interval(values, 0.8) == (0.0, 8.0)
    assert _credible_interval(values, 0.8, "equal-tail") == (1.0, 9.0)


@pytest.mark.parametrize("use_file", [False, True])
def test_posterior_rejects_duplicate_node_age_headers(tmp_path, use_file):
    source = "# posterior samples\n\nGen t_n3 t_n3\n1 2 99\n2 3 100\n"
    if use_file:
        path = tmp_path / "mcmc.txt"
        path.write_text(source)
        source = path
    with pytest.raises(ValueError, match="duplicate columns"):
        read_mcmctree_posterior(source, Tree("(A,B);", parser=1))


@pytest.mark.parametrize(
    ("text", "kind", "lower", "upper"),
    [
        ("@12.5", "point", 12.5, 12.5),
        ("B(10, 20, 0.025, 0.025)", "bounded", 10.0, 20.0),
        ("L(10, 0.1, 1, 0.025)", "lower", 10.0, None),
        ("U(20, 0.025)", "upper", None, 20.0),
        (">.06<.08", "bounded", 0.06, 0.08),
        (">.06", "lower", 0.06, None),
        ("<.08", "upper", None, 0.08),
    ],
)
def test_parse_mcmctree_calibration_variants(text, kind, lower, upper):
    record = parse_mcmctree_calibration(text)

    assert record["type"] == kind
    assert record.get("lower") == lower
    assert record.get("upper") == upper


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("B(20,10,0.025,0.025)", "must not exceed"),
        ("B(10,20,0,0.025)", "tail probabilities"),
        ("L(10,0,-1,0.025)", "scales must be positive"),
        ("@1e999", "must be finite"),
    ],
)
def test_parse_mcmctree_calibration_rejects_invalid_values(text, message):
    with pytest.raises(ValueError, match=message):
        parse_mcmctree_calibration(text)


def test_paml_node_mapping_is_tips_then_internal_preorder():
    tree = Tree("(((A,B),C),(D,E));", parser=1)
    mapping = paml_node_mapping(tree)

    assert [mapping[index].name for index in range(1, 6)] == ["A", "B", "C", "D", "E"]
    assert mapping[6].is_root
    assert set(mapping[7].leaf_names()) == {"A", "B", "C"}
    assert set(mapping[8].leaf_names()) == {"A", "B"}
    assert set(mapping[9].leaf_names()) == {"D", "E"}


def test_read_and_summarize_posterior_reconstructs_dated_tree(tmp_path):
    tree = Tree("((A,B)'B(2,6,0.025,0.025)',(C,D));", parser=1)
    posterior_path = tmp_path / "mcmc.txt"
    posterior_path.write_text(
        "Gen\tt_n5\tt_n6\tt_n7\tlnL\n"
        "1\t10\t4\t3\t-1\n"
        "2\t12\t5\t4\t-1\n"
        "3\t11\t6\t5\t-1\n"
        "4\t13\t7\t6\t-1\n"
    )

    annotate_mcmctree_calibrations(tree)
    posterior = read_mcmctree_posterior(posterior_path, tree, burnin=1, thin=2)
    summarize_mcmctree_posterior(
        tree, posterior, point="median", ci="equal-tail", level=0.8
    )

    assert posterior.sample_count == 2
    assert posterior.node_ids == (5, 6, 7)
    assert tree.props["age"] == pytest.approx(12.5)
    assert tree.common_ancestor(["A", "B"]).props["age"] == pytest.approx(6.0)
    assert tree.common_ancestor(["A", "B"]).props["calibration_type"] == "bounded"
    assert tree["A"].dist == pytest.approx(6.0)
    assert tree.common_ancestor(["A", "B"]).dist == pytest.approx(6.5)
    assert tree.props["age_ci_kind"] == "equal-tail"
    assert np.isfinite(float(tree.props["age_ci_low"]))

    sample = posterior_sample_tree(tree, posterior, 0)
    assert sample["A"].dist == pytest.approx(5.0)
    assert sample.common_ancestor(["A", "B"]).dist == pytest.approx(7.0)


def test_posterior_rejects_chronologically_invalid_sample(tmp_path):
    tree = Tree("((A,B),(C,D));", parser=1)
    path = tmp_path / "mcmc.txt"
    path.write_text("Gen t_n5 t_n6 t_n7\n1 4 5 2\n")

    with pytest.raises(ValueError, match="older than its parent"):
        read_mcmctree_posterior(path, tree)


def test_posterior_rejects_nonbinary_topology(tmp_path):
    tree = Tree("(A,B,C);", parser=1)
    path = tmp_path / "mcmc.txt"
    path.write_text("Gen t_n4\n1 4\n")

    with pytest.raises(ValueError, match="rooted and binary"):
        read_mcmctree_posterior(path, tree)
