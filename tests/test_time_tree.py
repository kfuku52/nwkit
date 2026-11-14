import numpy as np
import pytest
from ete4 import Tree

from nwkit.time_tree import (
    _age_order_statistics,
    annotate_mcmctree_calibrations,
    infer_node_ages_from_branch_lengths,
    paml_node_mapping,
    parse_mcmctree_calibration,
    posterior_sample_tree,
    read_mcmctree_posterior,
    summarize_mcmctree_posterior,
)


def test_mcmctree_interval_indices_match_paml_order_statistics():
    values = np.arange(9, -1, -1, dtype=float)

    assert _age_order_statistics(values, 0.8, "hpd") == (4.5, 0.0, 8.0)
    assert _age_order_statistics(values, 0.8, "equal-tail") == (4.5, 1.0, 9.0)


@pytest.mark.parametrize("point", ["mean", "median"])
@pytest.mark.parametrize("ci", ["hpd", "equal-tail"])
def test_mcmctree_dated_tree_preserves_finite_large_node_ages(point, ci):
    tree = Tree("((A,B),C);", parser=1)
    posterior = read_mcmctree_posterior(
        "Gen t_n4 t_n5\n1 1.6e308 9e307\n2 1.4e308 1e308\n", tree
    )
    summarize_mcmctree_posterior(tree, posterior, point=point, ci=ci, level=0.8)

    assert tree.props["age_mean"] == pytest.approx(1.5e308)
    assert tree.props["age_median"] == pytest.approx(1.5e308)
    assert tree.props["age_ci_low"] == pytest.approx(1.4e308)
    assert tree.props["age_ci_high"] == pytest.approx(1.6e308)
    internal = tree.common_ancestor(["A", "B"])
    assert internal.dist == pytest.approx(5.5e307)
    assert tree["A"].dist == pytest.approx(9.5e307)


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


@pytest.mark.parametrize("scale", [1e-180, 1e180])
def test_posterior_age_order_validation_does_not_depend_on_units(scale):
    tree = Tree("((A,B),C);", parser=1)
    with pytest.raises(ValueError, match="older than its parent"):
        read_mcmctree_posterior(f"Gen t_n4 t_n5\n1 {scale} {2 * scale}\n", tree)
    posterior = read_mcmctree_posterior(f"Gen t_n4 t_n5\n1 {2 * scale} {scale}\n", tree)
    summarize_mcmctree_posterior(tree, posterior)
    assert tree["A"].dist / scale == pytest.approx(1.0)


def test_inferred_dated_ages_avoid_overflow_in_child_path_average():
    tree = Tree("(A:1.7e308,B:1.7e308);", parser=1)
    infer_node_ages_from_branch_lengths(tree)
    assert tree.props["age"] == 1.7e308


def test_posterior_rejects_nonbinary_topology(tmp_path):
    tree = Tree("(A,B,C);", parser=1)
    path = tmp_path / "mcmc.txt"
    path.write_text("Gen t_n4\n1 4\n")

    with pytest.raises(ValueError, match="rooted and binary"):
        read_mcmctree_posterior(path, tree)
