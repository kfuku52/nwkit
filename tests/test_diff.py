import pandas as pd
import pytest

from nwkit.diff import _rows_have_differences, compare_trees, diff_main
from nwkit.util import read_tree
from tests.helpers import make_args


def test_rooted_diff_reports_common_and_tree_specific_clades():
    target = read_tree("(((A:1,B:1):1,C:1):1,D:1);", "1", True, quiet=True)
    source = read_tree("(((A:1,C:1):1,B:1):1,D:1);", "1", True, quiet=True)
    rows = compare_trees(target, source, comparison="rooted", target_class="intnode")
    statuses = {row["status"] for row in rows if row["record_type"] == "node"}
    assert "exact_match" in statuses
    assert "unmatched" in statuses
    assert "source_only" in statuses
    root_row = next(row for row in rows if row["comparison"] == "root_split")
    assert root_row["reason"] == "projected_root_bipartition"


def test_unrooted_diff_with_partial_taxon_overlap(tmp_nwk, tmp_path):
    target_path = tmp_nwk("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", "target.nwk")
    source_path = tmp_nwk("((A:1,B:1):1,(C:1,D:1):1);", "source.nwk")
    outfile = tmp_path / "diff.tsv"
    args = make_args(
        infile=target_path,
        infile2=source_path,
        outfile=str(outfile),
        format="1",
        format2="1",
        taxon_mode="intersection",
        comparison="unrooted",
        target="intnode",
        property=[],
        fail_on_difference=False,
    )
    diff_main(args)
    rows = pd.read_csv(outfile, sep="\t")
    leaf_summary = rows[rows["comparison"] == "leaf_set"].iloc[0]
    assert leaf_summary["status"] == "different"
    assert "E" in leaf_summary["target_taxa"]
    assert set(rows["comparison"]) >= {"leaf_set", "root_split", "unrooted_split"}


def test_unrooted_diff_does_not_report_topology_change_for_rerooted_tree():
    target = read_tree(
        "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);",
        "1",
        True,
        quiet=True,
    )
    source = read_tree(
        "(A:0.5,(B:1,(C:1,(D:1,E:1):1):2):0.5);",
        "1",
        True,
        quiet=True,
    )

    rows = compare_trees(
        target,
        source,
        comparison="unrooted",
        target_class="intnode",
    )

    root_row = next(row for row in rows if row["comparison"] == "root_split")
    assert root_row["status"] == "different"
    split_rows = [row for row in rows if row["comparison"] == "unrooted_split"]
    assert len(split_rows) == 2
    assert {row["status"] for row in split_rows} == {"exact_match"}


def test_unrooted_fail_ignores_arbitrary_newick_root_position(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        "(A:1,B:1,(C:1,(D:1,E:1):1):1);",
        "target.nwk",
    )
    source_path = tmp_nwk(
        "(D:1,E:1,(C:1,(A:1,B:1):1):1);",
        "source.nwk",
    )
    outfile = tmp_path / "diff.tsv"
    args = _diff_args(
        target_path,
        source_path,
        outfile,
        comparison="unrooted",
    )

    diff_main(args)

    rows = pd.read_csv(outfile, sep="\t")
    root_row = rows.loc[rows["comparison"] == "root_split"].iloc[0]
    assert root_row["status"] == "unresolved"
    split_rows = rows.loc[rows["comparison"] == "unrooted_split"]
    assert set(split_rows["status"]) == {"exact_match"}


def _diff_args(target_path, source_path, outfile, **overrides):
    values = {
        "infile": target_path,
        "infile2": source_path,
        "outfile": str(outfile),
        "format": "1",
        "format2": "1",
        "taxon_mode": "exact",
        "comparison": "rooted",
        "target": "intnode",
        "property": [],
        "fail_on_difference": True,
    }
    values.update(overrides)
    return make_args(**values)


@pytest.mark.parametrize(
    ("target_tree", "source_tree", "format", "properties", "target_class"),
    [
        (
            "((A:1,B:1)AB:2,C:1)target_root;",
            "((A:1,B:1)AB:2,C:1)source_root;",
            "1",
            [],
            "root",
        ),
        (
            "((A:1,B:1)90:2,C:1);",
            "((A:1,B:1)91:2,C:1);",
            "0",
            [],
            "intnode",
        ),
        (
            "((A:1,B:1):2,C:1);",
            "((A:1,B:1):3,C:1);",
            "1",
            [],
            "intnode",
        ),
        (
            "((A:1,B:1)[&&NHX:state=present]:2,C:1);",
            "((A:1,B:1)[&&NHX:state=absent]:2,C:1);",
            "1",
            ["state"],
            "intnode",
        ),
    ],
    ids=("root-name", "support", "length", "custom-property"),
)
def test_fail_on_difference_includes_matched_node_values(
    tmp_nwk, tmp_path, target_tree, source_tree, format, properties, target_class
):
    target_path = tmp_nwk(target_tree, "target.nwk")
    source_path = tmp_nwk(source_tree, "source.nwk")
    args = _diff_args(
        target_path,
        source_path,
        tmp_path / "diff.tsv",
        format=format,
        format2=format,
        property=properties,
        target=target_class,
    )

    with pytest.raises(ValueError, match="Tree differences were detected"):
        diff_main(args)


def test_fail_on_difference_uses_numeric_tolerance(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        "((A:1,B:1)90:2,C:1):1;",
        "target.nwk",
    )
    source_path = tmp_nwk(
        "((A:1.0000000000005,B:1)90.00000000001:2.0000000000005,C:1):1.0000000000005;",
        "source.nwk",
    )
    args = _diff_args(
        target_path,
        source_path,
        tmp_path / "diff.tsv",
        format="0",
        format2="0",
        target="all",
    )

    diff_main(args)


def test_missing_support_is_normalized_before_reporting_and_comparison():
    target = read_tree(
        "((A:1,B:1):2,C:1);",
        "0",
        True,
        quiet=True,
    )
    source = read_tree(
        "((A:1,B:1):2,C:1);",
        "1",
        True,
        quiet=True,
    )
    target.common_ancestor(["A", "B"]).support = None

    rows = compare_trees(
        target,
        source,
        comparison="rooted",
        target_class="intnode",
    )

    ab_row = next(
        row
        for row in rows
        if row["record_type"] == "node" and row["target_taxa"] == "A,B"
    )
    assert target.common_ancestor(["A", "B"]).support is None
    assert source.common_ancestor(["A", "B"]).support == -999999.0
    assert ab_row["target_support"] == ""
    assert ab_row["source_support"] == ""
    assert ab_row["support_delta"] == ""
    assert not _rows_have_differences(rows)


def test_unrooted_fail_sums_root_edge_halves_across_rerooting(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        "((A:1,B:1):2,(C:1,(D:1,E:1):7):3);",
        "target.nwk",
    )
    source_path = tmp_nwk(
        "((D:1,E:1):4,(C:1,(A:1,B:1):5):3);",
        "source.nwk",
    )
    outfile = tmp_path / "diff.tsv"
    args = _diff_args(
        target_path,
        source_path,
        outfile,
        comparison="unrooted",
    )

    diff_main(args)

    rows = pd.read_csv(outfile, sep="\t").to_dict("records")
    split_rows = {
        row["shared_taxa_or_split"]: row
        for row in rows
        if row["comparison"] == "unrooted_split"
    }
    assert split_rows["A,B|C,D,E"]["target_length"] == 5
    assert split_rows["A,B|C,D,E"]["source_length"] == 5
    assert split_rows["D,E|A,B,C"]["target_length"] == 7
    assert split_rows["D,E|A,B,C"]["source_length"] == 7


def test_unrooted_fail_detects_physical_root_edge_length_difference(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        "((A:1,B:1):2,(C:1,(D:1,E:1):7):3);",
        "target.nwk",
    )
    source_path = tmp_nwk(
        "((A:1,B:1):4,(C:1,(D:1,E:1):7):4);",
        "source.nwk",
    )
    outfile = tmp_path / "diff.tsv"
    args = _diff_args(
        target_path,
        source_path,
        outfile,
        comparison="unrooted",
    )

    with pytest.raises(ValueError, match="Tree differences were detected"):
        diff_main(args)

    rows = pd.read_csv(outfile, sep="\t")
    root_edge = rows.loc[rows["shared_taxa_or_split"] == "A,B|C,D,E"].iloc[0]
    assert root_edge["target_length"] == 5
    assert root_edge["source_length"] == 8


def test_unrooted_root_edge_retains_unique_or_equal_annotations():
    target = read_tree(
        "((A:1,B:1):2,(C:1,D:1):3);",
        "1",
        True,
        quiet=True,
    )
    source = read_tree(
        "(A:0.5,(B:1,(C:1,D:1):5):0.5);",
        "1",
        True,
        quiet=True,
    )
    for node in target.children:
        node.name = "edge"
        node.support = 95
        node.add_prop("state", "present")
    source_edge = source.common_ancestor(["C", "D"])
    source_edge.name = "edge"
    source_edge.support = 95
    source_edge.add_prop("state", "present")

    rows = compare_trees(
        target,
        source,
        comparison="unrooted",
        target_class="intnode",
        properties=[("state", "state")],
    )

    edge_row = next(row for row in rows if row["comparison"] == "unrooted_split")
    assert edge_row["status"] == "exact_match"
    assert edge_row["target_name"] == edge_row["source_name"] == "edge"
    assert edge_row["target_support"] == edge_row["source_support"] == 95
    assert edge_row["target_properties"] == '{"state": "present"}'
    assert edge_row["source_properties"] == '{"state": "present"}'


@pytest.mark.parametrize("annotation", ("name", "support", "state"))
def test_unrooted_root_edge_annotation_conflicts_are_ambiguous(annotation):
    target = read_tree(
        "((A:1,B:1):2,(C:1,D:1):3);",
        "1",
        True,
        quiet=True,
    )
    source = read_tree(
        "((A:1,B:1):2,(C:1,D:1):3);",
        "1",
        True,
        quiet=True,
    )
    properties = []
    if annotation == "name":
        target.children[0].name = "left"
        target.children[1].name = "right"
    elif annotation == "support":
        target.children[0].support = 90
        target.children[1].support = 91
    else:
        target.children[0].add_prop("state", "present")
        target.children[1].add_prop("state", "absent")
        properties = [("state", "state")]

    rows = compare_trees(
        target,
        source,
        comparison="unrooted",
        target_class="intnode",
        properties=properties,
    )

    edge_row = next(row for row in rows if row["comparison"] == "unrooted_split")
    assert edge_row["status"] == "ambiguous"
    assert "target={}".format(annotation) in edge_row["reason"]
    assert _rows_have_differences(rows, comparison="unrooted")


def test_unrooted_root_edge_length_is_missing_if_either_half_is_missing():
    target = read_tree(
        "((A:1,B:1),(C:1,D:1):3);",
        "1",
        True,
        quiet=True,
    )
    source = read_tree(
        "((A:1,B:1),(C:1,D:1):3);",
        "1",
        True,
        quiet=True,
    )

    rows = compare_trees(
        target,
        source,
        comparison="unrooted",
        target_class="intnode",
    )

    edge_row = next(row for row in rows if row["comparison"] == "unrooted_split")
    assert edge_row["target_length"] == ""
    assert edge_row["source_length"] == ""


def test_unrooted_diff_rejects_unrepresentable_physical_root_edge_total():
    target = read_tree(
        "((A:1,B:1):1e308,(C:1,D:1):1e308);",
        "1",
        True,
        quiet=True,
    )
    source = read_tree(
        "((A:1,B:1):1e308,(C:1,D:1):1e308);",
        "1",
        True,
        quiet=True,
    )

    with pytest.raises(
        ValueError,
        match="physical root-edge length total is too large",
    ):
        compare_trees(
            target,
            source,
            comparison="unrooted",
            target_class="intnode",
        )


def test_diff_rejects_an_unrepresentable_numeric_delta():
    target = read_tree(
        "(A:1e308,B:1);",
        "1",
        True,
        quiet=True,
    )
    source = read_tree(
        "(A:-1e308,B:1);",
        "1",
        True,
        quiet=True,
    )

    with pytest.raises(
        ValueError,
        match="numeric delta is too large",
    ):
        compare_trees(
            target,
            source,
            comparison="rooted",
            target_class="leaf",
        )
