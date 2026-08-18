import argparse
import subprocess
import sys
from pathlib import Path

import pytest

from nwkit import __version__
from nwkit.cli import main, parser


def _subcommand_parser(command):
    subparsers_action = next(
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers_action.choices[command]


def test_version_option(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"nwkit {__version__}"


def test_regress_replaces_pgls_without_compatibility_alias(capsys):
    assert _subcommand_parser("regress").prog == "nwkit regress"

    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["pgls"])

    assert exc_info.value.code == 2
    assert "invalid choice: 'pgls'" in capsys.readouterr().err


@pytest.mark.parametrize(
    "command",
    ["annotate", "asr", "draw", "monophyly", "sample", "skim"],
)
def test_shared_tip_table_commands_default_to_unmatched_warn(command):
    assert _subcommand_parser(command).get_default("unmatched") == "warn"


def test_contrast_has_independent_unmatched_default_and_allows_override():
    assert _subcommand_parser("contrast").get_default("unmatched") == "error"
    assert (
        parser.parse_args(
            [
                "contrast",
                "--trait",
                "traits.tsv",
                "--columns",
                "expression",
                "--unmatched",
                "ignore",
            ]
        ).unmatched
        == "ignore"
    )
    assert parser.parse_args(["draw", "--unmatched", "error"]).unmatched == "error"


@pytest.mark.integration
def test_python_module_entry_point():
    project_root = Path(__file__).parents[1]
    result = subprocess.run(
        [sys.executable, "-m", "nwkit", "--version"],
        cwd=project_root,
        capture_output=True,
        check=False,
        text=True,
    )

    assert result.returncode == 0
    assert result.stdout.strip() == f"nwkit {__version__}"


@pytest.mark.parametrize(
    "arguments",
    [
        ["rescale", "--factor", "nan"],
        ["draw", "--figure-width", "inf"],
        ["consensus", "--min-freq", "-inf"],
        ["consensus", "--min-freq", "1.01"],
        ["mcmctree", "--min-clade-prop", "-0.01"],
        ["subtree", "--dup-conf-score-threshold", "1.01"],
    ],
)
def test_float_options_reject_non_finite_or_out_of_range_values(arguments, capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(arguments)

    assert exc_info.value.code == 2
    assert "error:" in capsys.readouterr().err


def test_console_handler_error_is_concise_without_traceback(
    tmp_path,
    monkeypatch,
    capsys,
):
    infile = tmp_path / "broken.nwk"
    infile.write_text("(A:1,B:1)")
    monkeypatch.setattr(
        sys,
        "argv",
        ["nwkit", "info", "-i", str(infile)],
    )

    assert main() == 2
    stderr = capsys.readouterr().err
    assert "nwkit: error:" in stderr
    assert "Traceback" not in stderr


def test_console_debug_option_reraises_for_traceback(tmp_path, monkeypatch):
    infile = tmp_path / "broken.nwk"
    infile.write_text("(A:1,B:1)")
    monkeypatch.setattr(
        sys,
        "argv",
        ["nwkit", "info", "-i", str(infile), "--debug"],
    )

    with pytest.raises(Exception, match="parse"):
        main()
