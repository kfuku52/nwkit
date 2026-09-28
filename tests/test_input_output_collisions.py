"""CLI outputs must not replace files still used as inputs."""

import argparse

import pytest

from nwkit.cli import main, parser
from nwkit.provenance import INPUT_PATH_ARGUMENTS, OUTPUT_ARGUMENTS


@pytest.mark.parametrize(
    ("command", "other_option", "output_option", "output_suffix"),
    [
        ("sample", "--trait", "--report", ""),
        ("skim", "--trait", "--group-table-prefix", ".all.tsv"),
        ("transfer", "--infile2", "--report", ""),
        ("compose", "--name-source", "--report", ""),
        ("root", "--infile2", "--outfile", ""),
        ("mcmctree", "--posterior", "--outfile", ""),
        ("consensus", "--weight-tsv", "--outfile", ""),
        ("rename", "--name-tsv", "--outfile", ""),
        ("intersection", "--seqin", "--outfile", ""),
    ],
)
def test_secondary_inputs_are_preserved(
    tmp_path, command, other_option, output_option, output_suffix
):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,C:1);")
    other = tmp_path / ("other" + output_suffix)
    other.write_text("original input\n")
    output = other if not output_suffix else tmp_path / "other"
    options = [command, "--infile", str(tree), other_option, str(other)]
    if command == "root":
        options.extend(["--method", "transfer"])
    options.extend([output_option, str(output)])
    with pytest.raises(ValueError, match="must not overwrite input"):
        main(options)
    assert other.read_text() == "original input\n"


@pytest.mark.parametrize(
    "command",
    [
        "info",
        "printlabel",
        "nwk2table",
        "validate",
        "monophyly",
        "cladefreq",
        "diff",
        "dist",
    ],
)
def test_non_tree_output_cannot_replace_primary_tree(tmp_path, command):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,C:1);")
    options = [command, "--infile", str(tree), "--outfile", str(tree)]
    if command in {"diff", "dist"}:
        options.extend(["--infile2", str(tmp_path / "other.nwk")])
    with pytest.raises(ValueError, match="must not overwrite input"):
        main(options)
    assert tree.read_text() == "((A:1,B:1):1,C:1);"


@pytest.mark.parametrize("command", ["table2nwk", "image"])
def test_cross_format_output_cannot_replace_primary_input(tmp_path, command):
    source = tmp_path / "input.tsv"
    source.write_text("original input\n")
    if command == "image":
        options = [
            "image",
            "--infile",
            str(source),
            "--out-dir",
            str(tmp_path / "images"),
            "--manifest-out",
            str(source),
        ]
    else:
        options = ["table2nwk", "--infile", str(source), "--outfile", str(source)]
    with pytest.raises(ValueError, match="must not overwrite input"):
        main(options)
    assert source.read_text() == "original input\n"


def test_tree_editor_can_update_its_primary_tree(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A:1,B:1);")
    main(["rescale", "--infile", str(tree), "--outfile", str(tree), "--factor", "2"])
    assert "A:2" in tree.read_text()


@pytest.mark.parametrize("alias_kind", ["symlink", "hardlink"])
def test_output_alias_cannot_replace_secondary_input(tmp_path, alias_kind):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A:1,B:1);")
    trait = tmp_path / "trait.tsv"
    trait.write_text("leaf_name\tx\nA\t1\nB\t2\n")
    alias = tmp_path / "report.tsv"
    try:
        if alias_kind == "symlink":
            alias.symlink_to(trait)
        else:
            alias.hardlink_to(trait)
    except OSError:
        pytest.skip("Filesystem links are unavailable")
    with pytest.raises(ValueError, match="must not overwrite input"):
        main(
            [
                "sample",
                "--infile",
                str(tree),
                "--trait",
                str(trait),
                "--report",
                str(alias),
            ]
        )
    assert trait.read_text() == "leaf_name\tx\nA\t1\nB\t2\n"


def test_all_parser_path_arguments_have_provenance_roles():
    path_arguments = set()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for subparser in action.choices.values():
                path_arguments.update(
                    option.dest
                    for option in subparser._actions
                    if option.metavar == "PATH"
                )
    # Cache roots and a directory used only to resolve manifest-relative assets
    # are not direct input/output files in an analysis result.
    path_arguments -= {"download_dir", "tip_image_root"}
    assert path_arguments <= INPUT_PATH_ARGUMENTS | OUTPUT_ARGUMENTS


@pytest.mark.parametrize(
    ("command", "input_option"),
    [("asr", "--tip-likelihoods"), ("mcmctree", "--posterior")],
)
def test_audit_path_cannot_replace_newly_registered_input(
    tmp_path, command, input_option
):
    source = tmp_path / "source.tsv"
    source.write_text("original input\n")
    options = [command, "--infile", "(A:1,B:1);", input_option, str(source)]
    options.extend(["--audit", str(source)])
    with pytest.raises(ValueError, match="must not overwrite input|distinct"):
        main(options)
    assert source.read_text() == "original input\n"


@pytest.mark.parametrize(
    ("command", "output_option"),
    [("root", "--candidates-out"), ("asr", "--posterior-samples-out")],
)
def test_audit_path_cannot_replace_newly_registered_output(
    tmp_path, command, output_option
):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,C:1);")
    target = tmp_path / "output.tsv"
    target.write_text("previous output\n")
    options = [command, "--infile", str(tree), output_option, str(target)]
    if command == "root":
        options.extend(["--method", "reconciliation", "--species-tree", str(tree)])
    options.extend(["--audit", str(target)])
    with pytest.raises(ValueError, match="Output paths must be distinct"):
        main(options)
    assert target.read_text() == "previous output\n"
