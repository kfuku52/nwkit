"""Shared and direct writers preserve rooting while tokens default to off."""

import io
from types import SimpleNamespace

import pytest

from nwkit.cli import main, subparsers
from nwkit.convert import convert_tree_text
from nwkit.gene_tree_search_generax import tree_text
from nwkit.mul_msc import dated_text
from nwkit.mul_reconcile import topology_text
from nwkit.rooting_state import get_rooting_info, rooting_output_policy
from nwkit.util import read_tree, write_tree


def read(text, format=1):
    return read_tree(text, format, True, quiet=True)


def signature(tree):
    return [
        (
            tuple(sorted(node.leaf_names())),
            node.name,
            node.dist,
            node.support,
            node.props.get("tag"),
        )
        for node in tree.traverse()
    ]


@pytest.mark.parametrize(
    "source",
    [
        "[&R]((A:1.25,B:2.5)inner:0.5,'C[&U]':3)root:0;",
        "((A:1.25,B:2.5)inner:0.5,'C[&U]':3)root:0;",
        "(A:1,B:1,C:1)root;",
    ],
)
@pytest.mark.parametrize("token", [False, True])
def test_writer_no_unconditional_nhx_and_known_state_opt_in(source, token):
    tree = read(source)
    props = [dict(node.props) for node in tree.traverse()]
    output = io.StringIO()
    write_tree(
        tree, SimpleNamespace(outfile=output, rooting_token=token), 1, quiet=True
    )
    text = output.getvalue()
    assert "NHX" not in text
    assert text.startswith("[&R]") is (token and get_rooting_info(tree).rooted is True)
    assert get_rooting_info(read(text)).rooted is get_rooting_info(tree).rooted
    assert signature(read(text)) == signature(tree)
    assert [dict(node.props) for node in tree.traverse()] == props


@pytest.mark.parametrize(
    "source,expected",
    [
        ("[&U](A:1,B:2)root;", "no"),
        ("[&R](A:1,B:2,C:3)root;", "yes"),
        ("(A:1,B:2)[&&NHX:nwkit_rooted=unknown];", "unknown"),
    ],
)
@pytest.mark.parametrize("token", [False, True])
def test_non_inferable_states_and_existing_nhx_are_preserved(source, expected, token):
    tree = read(source)
    output = io.StringIO()
    write_tree(
        tree, SimpleNamespace(outfile=output, rooting_token=token), 1, quiet=True
    )
    text = output.getvalue()
    if not token or expected == "unknown":
        assert not text.startswith("[&")
        assert f"nwkit_rooted={expected}" in text
    else:
        assert text.startswith("[&R]" if expected == "yes" else "[&U]")
    assert get_rooting_info(read(text)).rooted is get_rooting_info(tree).rooted
    assert signature(read(text)) == signature(tree)


def test_nhx_switch_precedes_token_and_preserves_selected_annotations():
    tree = read("[&R]((A:1,B:2)inner:3,C:4)root[&&NHX:tag=kept];")
    output = io.StringIO()
    write_tree(
        tree,
        SimpleNamespace(outfile=output, rooting_token=True, rooting_nhx=True),
        1,
        quiet=True,
        props=["tag"],
    )
    assert not output.getvalue().startswith("[&")
    assert "nwkit_rooted=yes" in output.getvalue()
    assert signature(read(output.getvalue())) == signature(tree)


@pytest.mark.parametrize("target", ["newick", "nhx", "figtree"])
@pytest.mark.parametrize("token", [False, True])
def test_conversion_roundtrip_preserves_fields_and_format_declarations(target, token):
    source = "[&R](('A[&U]':0.123456789012345,B:2)95:3,C:5)100:0;"
    result = convert_tree_text(
        source, target=target, tree_format=0, rooting_token=token
    )
    restored = read(result, 0)
    assert signature(restored) == signature(read(source, 0))
    assert get_rooting_info(restored).rooted is True
    if target == "figtree":
        assert "TREE 1 = [&R]" in result
    else:
        assert result.startswith("[&R]") is token
        assert "NHX" not in result
    again = convert_tree_text(
        result, target="newick", tree_format=0, rooting_token=token
    )
    assert signature(read(again, 0)) == signature(read(source, 0))
    assert again.startswith("[&R]") is token


@pytest.mark.parametrize(
    "source",
    ["[&U](A:1,B:1);", "[&R](A:1,B:1,C:1);", "(A:1,B:1)[&&NHX:nwkit_rooted=unknown];"],
)
def test_plain_newick_rejects_rooting_loss_before_replacing_output(tmp_path, source):
    destination = tmp_path / "result.nwk"
    destination.write_text("previous output\n")
    with pytest.raises(ValueError, match="cannot preserve this rooting"):
        main(
            [
                "convert",
                "-i",
                source,
                "--to",
                "newick",
                "--properties",
                "drop",
                "-o",
                str(destination),
            ]
        )
    assert destination.read_text() == "previous output\n"
    restored = convert_tree_text(source, target="nhx", properties="drop")
    assert (
        get_rooting_info(read(restored)).rooted is get_rooting_info(read(source)).rooted
    )


@pytest.mark.parametrize("token", [False, True])
def test_unrooted_nexus_and_existing_annotations_survive_conversion(token):
    source = "[&U]((A:1,B:2)inner:3,C:4)root[&&NHX:tag=kept];"
    nhx = convert_tree_text(source, target="nhx", rooting_token=token)
    nexus = convert_tree_text(nhx, target="figtree", rooting_token=token)
    assert "UTREE 1 = [&U]" in nexus
    assert signature(read(nexus)) == signature(read(source))
    assert get_rooting_info(read(nexus)).rooted is False
    again = convert_tree_text(nexus, target="nhx", rooting_token=token)
    assert again.startswith("[&U]") is token
    assert signature(read(again)) == signature(read(source))


@pytest.mark.parametrize("writer", [dated_text, topology_text, tree_text])
@pytest.mark.parametrize("token", [False, True])
def test_direct_serializers_share_policy_without_mutating_tree(writer, token):
    tree = read("[&R]((A:0.123456789012345,B:2)inner:3,C:4)root;")
    original = [dict(node.props) for node in tree.traverse()]
    text = writer(tree, args=SimpleNamespace(rooting_token=token))
    assert text.startswith("[&R]") is token
    assert "NHX" not in text
    assert get_rooting_info(read(text)).rooted is True
    assert list(read(text).leaf_names()) == list(tree.leaf_names())
    assert [dict(node.props) for node in tree.traverse()] == original
    if writer is dated_text:
        assert [n.dist for n in read(text).traverse()] == [
            n.dist for n in tree.traverse()
        ]


@pytest.mark.parametrize("token", [False, True])
@pytest.mark.parametrize(
    "command", ["label", "sample", "shuffle", "prune", "rescale", "nhx2nwk"]
)
@pytest.mark.parametrize("stdout", [False, True])
def test_real_cli_writers_file_and_stdout(tmp_path, capsys, command, token, stdout):
    source = "[&R]((A:1,B:2):1,C:3);"
    destination = tmp_path / "result.nwk"
    options = {
        "label": [],
        "sample": ["--n", "3"],
        "shuffle": ["--seed", "1"],
        "prune": ["--pattern", "never-matches"],
        "rescale": ["--factor", "1"],
        "nhx2nwk": [],
    }
    main(
        [
            command,
            "-i",
            source,
            "-o",
            "-" if stdout else str(destination),
            "--rooting-token",
            "yes" if token else "no",
            *options[command],
        ]
    )
    text = capsys.readouterr().out if stdout else destination.read_text()
    lines = [line for line in text.splitlines() if line]
    assert len(lines) == 1
    for line in lines:
        assert line.startswith("[&R]") is token
        assert "NHX" not in line
        assert get_rooting_info(read(line)).rooted is True


def test_all_public_tree_output_parsers_expose_policy_defaults_and_aliases():
    commands = []
    for command, parser in subparsers.choices.items():
        destinations = {action.dest for action in parser._actions}
        if destinations & {"outformat", "tree_outformat", "tree_out"} or command in {
            "convert",
            "radte",
            "mcmctree",
        }:
            commands.append(command)
            action = next(
                action for action in parser._actions if action.dest == "rooting_token"
            )
            assert action.default is False
            assert action.option_strings == ["--rooting-token", "--rooting_token"]
            assert action.type("yes") is True and action.type("no") is False
    assert {
        "label",
        "sample",
        "convert",
        "radte",
        "mul-reconcile",
        "gene-tree-search",
        "asr",
    } <= set(commands)
    assert rooting_output_policy(read("(A,B);")) == ("", False)


def test_snake_case_rooting_token_alias_is_deprecated_but_works(capsys):
    main(["label", "-i", "(A:1,B:1);", "--rooting_token", "yes"])
    result = capsys.readouterr()
    assert result.out.startswith("[&R]")
    assert "--rooting-token" in result.err


@pytest.mark.parametrize("token", [False, True])
@pytest.mark.parametrize("stdout", [False, True])
def test_root_candidate_multiple_tree_bundle_preserves_output_policy(
    tmp_path, capsys, token, stdout
):
    from nwkit.root import _write_reconciliation_candidates

    tree = read("[&R](A:1,(B:1,C:1)inner:1)root;")
    candidate_path = tmp_path / "candidates.nwk"
    primary = tmp_path / "selected.nwk"
    args = SimpleNamespace(
        outfile="-" if stdout else str(primary),
        candidates_out=str(candidate_path),
        outformat=1,
        rooting_token=token,
    )
    evaluation = SimpleNamespace(
        candidates=[
            SimpleNamespace(split=(("A",), ("B", "C"))),
            SimpleNamespace(split=(("B",), ("A", "C"))),
        ]
    )
    _write_reconciliation_candidates(tree, evaluation, args, [])
    primary_text = capsys.readouterr().out if stdout else primary.read_text()
    candidates = candidate_path.read_text().splitlines()
    assert len(candidates) == 2
    for text in [primary_text, *candidates]:
        assert text.startswith("[&R]") is token
        assert "NHX" not in text
        assert set(read(text).leaf_names()) == {"A", "B", "C"}
        assert get_rooting_info(read(text)).rooted is True


@pytest.mark.parametrize("token", [False, True])
def test_fixed_msc_cli_tree_columns_and_model_scores_share_policy(tmp_path, token):
    import json

    import pandas as pd

    from tests.test_mul_msc import invoke

    result = invoke(
        tmp_path,
        [
            "--tree-out",
            str(tmp_path / "best.nwk"),
            "--model-out",
            str(tmp_path / "model.json"),
            "--rooting-token",
            "yes" if token else "no",
        ],
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "best.nwk").read_text().startswith("[&R]") is token
    scores = pd.read_csv(tmp_path / "scores.tsv", sep="\t")
    for text in scores["dated.tree"].dropna():
        assert text.startswith("[&R]") is token
    model = json.loads((tmp_path / "model.json").read_text())
    for row in model["scores"]:
        if row["dated.tree"] is not None:
            assert row["dated.tree"].startswith("[&R]") is token


@pytest.mark.parametrize("token", [False, True])
def test_locus_json_population_trees_share_policy_and_remain_readable(tmp_path, token):
    import json

    from tests.test_mul_locus import arguments

    args, paths = arguments(tmp_path)
    main([*args, "--rooting-token", "yes" if token else "no"])
    model = json.loads(paths["model"].read_text())
    texts = [
        model["hypothesis_scope"]["species_tree"],
        *(bank["population_tree"] for bank in model["banks"]),
    ]
    assert texts and all(text.startswith("[&R]") is token for text in texts)
    assert all(get_rooting_info(read(text)).rooted is True for text in texts)


@pytest.mark.parametrize("token", [False, True])
def test_dl_mapping_and_topology_columns_share_policy(tmp_path, token):
    import pandas as pd

    from tests.test_mul_reconcile import invoke

    result = invoke(
        tmp_path,
        [
            "--tree-out",
            str(tmp_path / "best.nwk"),
            "--report",
            str(tmp_path / "detail.tsv"),
            "--rooting-token",
            "yes" if token else "no",
        ],
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "best.nwk").read_text().startswith("[&R]") is token
    scores = pd.read_csv(tmp_path / "scores.tsv", sep="\t")
    mappings = pd.read_csv(tmp_path / "detail.tsv", sep="\t")
    for text in [*scores["labeled.tree"], *mappings["maps"]]:
        assert text.startswith("[&R]") is token
