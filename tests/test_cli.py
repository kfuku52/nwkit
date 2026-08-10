import subprocess
import sys
from pathlib import Path

import pytest

from nwkit import __version__
from nwkit.cli import main


def test_version_option(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"nwkit {__version__}"


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
