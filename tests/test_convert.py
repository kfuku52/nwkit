"""Real conversion contracts, including annotation loss and malformed inputs."""

import io
from decimal import Decimal

import pytest

from nwkit.cli import main
from nwkit.convert import convert_tree_text
from nwkit.util import read_tree


def parsed(text):
    return read_tree(text, "auto", True, quiet=True)


def test_property_labels_are_internal_only_and_precede_removal():
    source = "((A:1[&&NHX:S=tip],B:2)old:3[&&NHX:S=AB:D=Y],C:4)root[&&NHX:S=ancestor];"
    result = convert_tree_text(
        source, target="newick", node_label="S", properties="drop"
    )
    tree = read_tree(result, 1, True, quiet=True)
    assert set(tree.leaf_names()) == {"A", "B", "C"}
    assert {n.name for n in tree.traverse() if not n.is_leaf} == {"AB", "ancestor"}
    assert tree["A"].dist == 1 and tree["AB"].dist == 3
    assert "NHX" not in result


def test_label_copy_retains_properties_and_missing_labels():
    source = "((A:1,B:1)unchanged:1,C:2)root[&&NHX:S=two words];"
    result = convert_tree_text(source, node_label="S")
    tree = read_tree(result, 1, True, quiet=True)
    assert tree.name == "two words"
    assert tree.props["S"] == "two words"
    assert tree["unchanged"].name == "unchanged"


def test_drop_properties_preserves_rooting_comments_and_precision():
    source = "[&R](A:0.123456789012345,B:2)[note][&&NHX:D=Y:age=2];"
    result = convert_tree_text(source, target="newick", properties="drop")
    assert "[&R]" in result and "[note]" in result
    assert "0.123456789012345" in result and "NHX" not in result


@pytest.mark.parametrize("label", ["95", "O'Brien", "a,b", "two words"])
def test_property_labels_are_quoted_and_replace_support(label):
    result = convert_tree_text(
        f"(A:1,B:1)90[&&NHX:S={label}];",
        tree_format=0,
        node_label="S",
        properties="drop",
        target="newick",
    )
    assert read_tree(result, 1, True, quiet=True).name == label


def test_age_property_label_is_copied_before_scaling():
    result = convert_tree_text(
        "(A:1,B:1)[&&NHX:age=2];", node_label="age", time_factor=1000
    )
    tree = read_tree(result, 1, True, quiet=True)
    assert tree.name == "2"
    assert float(tree.props["age"]) == 2000


@pytest.mark.parametrize("target", ["nhx", "newick", "figtree"])
def test_generated_quotes_are_independent_of_input_quote_policy(target):
    result = convert_tree_text(
        "(A:1,B:1)[&&NHX:S=two words];",
        target=target,
        quoted_node_names=False,
        node_label="S",
        properties="drop",
    )
    assert read_tree(result, 1, True, quiet=True).name == "two words"
    with pytest.raises(ValueError, match="Quoted node names"):
        convert_tree_text(
            "('A':1,B:1)[&&NHX:S=ancestor];",
            node_label="S",
            quoted_node_names=False,
        )


@pytest.mark.parametrize("target", ["nhx", "newick", "figtree"])
@pytest.mark.parametrize("rooted", ["auto", "no"])
def test_rooting_property_is_copied_before_canonicalization(target, rooted):
    result = convert_tree_text(
        "(A:1,B:1)[&&NHX:nwkit_rooted=yes];",
        target=target,
        rooted=rooted,
        node_label="nwkit_rooted",
        properties="drop",
    )
    assert read_tree(result, 1, True, quiet=True).name == "yes"
    assert ("[&U]" if rooted == "no" else "[&R]") in result


@pytest.mark.parametrize("property_value", ["name=old", "dist=9", "support=95"])
@pytest.mark.parametrize("properties", ["keep", "drop"])
@pytest.mark.parametrize("target", ["nhx", "newick", "figtree"])
def test_reserved_properties_cannot_silently_override_converted_fields(
    property_value, properties, target
):
    source = f"((A:1,B:1)original:2[&&NHX:S=ancestor:{property_value}],C:3);"
    with pytest.raises(ValueError, match="Reserved NHX properties"):
        convert_tree_text(
            source, target=target, properties=properties, node_label="S", time_factor=10
        )


def test_reserved_tip_property_failure_preserves_output(tmp_path):
    source = tmp_path / "input.nhx"
    output = tmp_path / "output.nwk"
    source.write_text("(A:1[&&NHX:name=X],B:1);")
    output.write_text("previous result")
    with pytest.raises(ValueError, match="Reserved NHX properties"):
        main(["convert", "-i", str(source), "-o", str(output), "--properties", "drop"])
    assert output.read_text() == "previous result"


def test_property_drop_includes_intervals_but_validates_input():
    source = "(A:1,B:1)[&95%HPD={0.5,1.5}];"
    assert (
        convert_tree_text(source, target="newick", properties="drop") == "(A:1,B:1);\n"
    )
    with pytest.raises(ValueError):
        convert_tree_text("(A:1,B:1)[&95%HPD={2,1}];", properties="drop")
    with pytest.raises(ValueError, match="properties"):
        convert_tree_text("(A:1,B:1);", properties="invalid")


def test_convert_cli_replaces_nhx2nwk_label_workflow(monkeypatch, capsys):
    monkeypatch.setattr("sys.stdin", io.StringIO("(A:1,B:1)[&&NHX:S=ancestor];"))
    main(["convert", "--to", "newick", "--node-label", "S", "--properties", "drop"])
    assert read_tree(capsys.readouterr().out, 1, True, quiet=True).name == "ancestor"


def test_multiline_nexus_roundtrip_scales_lengths_ages_and_intervals():
    source = "#NEXUS\nBEGIN TREES;\nTREE dated = [&R] (A:0.1,\nB:0.1)[&95%HPD={0.05,0.15},age=0.1];\nEND;\n"
    nhx = convert_tree_text(source, time_factor=1000)
    tree = parsed(nhx)
    assert tree["A"].dist == 100
    assert float(tree.props["age"]) == 100
    assert float(tree.props["age_ci_low"]) == 50
    assert float(tree.props["age_ci_high"]) == 150
    assert tree.props["age_ci_kind"] == "HPD"
    restored = parsed(
        convert_tree_text(nhx, target="figtree", time_factor=Decimal("0.001"))
    )
    assert restored["A"].dist == pytest.approx(0.1)
    assert float(restored.props["age_ci_low"]) == pytest.approx(0.05)


def test_names_support_and_non_time_comments_are_not_scaled():
    text = "(('A:1,B[&95%={2,3}]':1,C:1)95:1,D:2)[note='x:2,y'][&rate=0.3];"
    result = convert_tree_text(text, time_factor=1000)
    assert "'A:1,B[&95%={2,3}]':1000" in result or "'A:1,B[&95%={2,3}]':1E+3" in result
    assert ")95:" in result
    assert "[note='x:2,y']" in result and "[&rate=0.3]" in result


def test_plain_newick_requires_explicit_interval_removal():
    text = "(A:1,B:1)[&95%={0.5,1.5}];"
    with pytest.raises(ValueError, match="cannot retain"):
        convert_tree_text(text, target="newick")
    assert convert_tree_text(text, target="newick", age_ci="drop") == "(A:1,B:1);\n"
    with pytest.raises(ValueError, match="cannot retain"):
        convert_tree_text("(A:1,B:1)[&&NHX:D=Y];", target="newick", age_ci="drop")


def test_credible_interval_kind_and_level_are_preserved():
    text = "(A:1,B:1)[&&NHX:age_ci_low=0.5:age_ci_high=1.5:age_ci_kind=equal-tail:age_ci_level=0.9];"
    result = convert_tree_text(text, target="figtree")
    assert "[&90%={0.5,1.5}]" in result or "[&9E+1%={0.5,1.5}]" in result
    assert float(parsed(result).props["age_ci_level"]) == 0.9
    assert parsed(result).props["age_ci_kind"] == "equal-tail"


def test_mcmctree_known_views_select_annotated_tree_but_distinct_trees_are_ambiguous():
    header = "Species tree for FigTree.\n"
    source = header + "(1_A,2_B)3;\n(A:1,B:1);\n(A:1,B:1)[&95%HPD={0.5,1.5}];\n"
    assert "age_ci_low" in convert_tree_text(source)
    ambiguous = header + "(A:1,B:1);\n(A:2,B:2);\n"
    with pytest.raises(ValueError, match="tree-index"):
        convert_tree_text(ambiguous)
    assert parsed(convert_tree_text(ambiguous, tree_index=2))["A"].dist == 2
    with pytest.raises(ValueError, match="topology"):
        convert_tree_text(source, tree_index=1)


def test_generic_multi_tree_selection_counts_input_statements():
    with pytest.raises(ValueError, match="tree-index"):
        convert_tree_text("(A:1,B:1);(C:2,D:2);")
    assert set(
        parsed(convert_tree_text("(A:1,B:1);(C:2,D:2);", tree_index=2)).leaf_names()
    ) == {"C", "D"}


@pytest.mark.parametrize(
    "text",
    [
        "(A:1,);",
        "(A:NaN,B:1);",
        "(A:-1,B:1);",
        "(A:1,A:1);",
        "(A B:1,C:1);",
        "(A:1,B:1)[&95%HPD={nan,1}];",
        "(A:1,B:1)[&95%={2,1}];",
        "(A:1,B:1)[&95%={oops}];",
        "(A:1,B:1);garbage",
        "#NEXUS\nBEGIN TREES;\nTRANSLATE 1 A,2 B;\nTREE t=(1:1,2:1);\nEND;",
    ],
)
def test_invalid_inputs_are_rejected(text):
    with pytest.raises(Exception):  # noqa: B017 - ETE and NWKIT expose different parse errors
        convert_tree_text(text)


@pytest.mark.parametrize("factor", ["NaN", "sNaN", "Infinity", "-Infinity", "0", "-1"])
def test_invalid_time_factor(factor):
    with pytest.raises(ValueError, match="finite and positive"):
        convert_tree_text("(A:1,B:1);", time_factor=factor)


def test_cli_stdio_and_failed_conversion_preserves_output(
    tmp_path, monkeypatch, capsys
):
    monkeypatch.setattr("sys.stdin", io.StringIO("(A:1,B:1)[&95%HPD={0.5,1.5}];"))
    main(["convert", "--to", "nhx"])
    assert "age_ci_low" in capsys.readouterr().out
    source, target = tmp_path / "source.tre", tmp_path / "result.tre"
    source.write_text("(A:1,B:1)[&95%HPD={0.5,1.5}];")
    target.write_text("original")
    with pytest.raises(ValueError):
        main(["convert", "-i", str(source), "-o", str(target), "--to", "newick"])
    assert target.read_text() == "original"
    main(
        [
            "convert",
            "-i",
            str(source),
            "-o",
            str(target),
            "--to",
            "newick",
            "--age-ci",
            "drop",
        ]
    )
    assert target.read_text() == "(A:1,B:1);\n"


@pytest.mark.parametrize(
    "source",
    [
        "#NEXUS\nBEGIN TREES; TRANSLATE 1 A,2 B; TREE t=(1:1,2:1); END;",
        "#NEXUS\nBEGIN TREES; [comment] TRANSLATE 1 A,2 B; TREE t=(1:1,2:1); END;",
    ],
)
def test_inline_nexus_translation_is_rejected(source):
    with pytest.raises(ValueError, match="TRANSLATE"):
        convert_tree_text(source)


def test_nexus_leading_comments_and_numeric_names_are_preserved():
    result = convert_tree_text(
        "#NEXUS\nBEGIN TREES; [comment] TREE t=(A:1,B:1)'0.2,0.3'; END;",
        time_factor=1000,
    )
    assert "'0.2,0.3'" in result
    assert "age_ci_low" not in result


def test_explicit_rooting_and_override_survive_conversion():
    from nwkit.rooting_state import get_rooting_info

    source = "[&U](A:1,B:1,C:1);"
    assert (
        get_rooting_info(parsed(convert_tree_text(source, target="figtree"))).state
        == "unrooted"
    )
    assert (
        get_rooting_info(
            parsed(convert_tree_text(source, target="newick", rooted="yes"))
        ).state
        == "rooted"
    )
    with pytest.raises(ValueError, match="Conflicting"):
        convert_tree_text("[&R](A:1,B:1)[&&NHX:nwkit_rooted=no];")


def test_scaling_preserves_missing_lengths_and_scales_explicit_ages():
    result = parsed(convert_tree_text("((A:1,B):2,C)[&age=3];", time_factor=1000))
    assert result["A"].dist == 1000
    assert result["B"].dist is None and result["C"].dist is None
    assert float(result.props["age"]) == 3000
