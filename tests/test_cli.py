import subprocess
import sys
from pathlib import Path

import pytest

from nwkit import __version__
from nwkit.cli import main


def test_version_option(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(['--version'])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out.strip() == f'nwkit {__version__}'


def test_python_module_entry_point():
    project_root = Path(__file__).parents[1]
    result = subprocess.run(
        [sys.executable, '-m', 'nwkit', '--version'],
        cwd=project_root,
        capture_output=True,
        check=False,
        text=True,
    )

    assert result.returncode == 0
    assert result.stdout.strip() == f'nwkit {__version__}'
