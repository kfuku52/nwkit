#!/usr/bin/env python3
"""Generate version-locked inventory data for the NWKIT manuscript.

The report is written to standard output so it can be inspected directly or
captured by a higher-level reproducibility workflow.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import nwkit  # noqa: E402
from nwkit.cli import parser as nwkit_parser  # noqa: E402


def git_value(*args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def command_inventory() -> list[dict[str, object]]:
    subparser_action = next(
        action
        for action in nwkit_parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    commands = []
    for name, command_parser in sorted(subparser_action.choices.items()):
        if name == "help":
            continue
        options = {
            action.dest: action.default
            for action in command_parser._actions
            if action.dest != "help"
        }
        commands.append(
            {
                "name": name,
                "reads_stdin_by_default": options.get("infile") == "-",
                "writes_stdout_by_default": options.get("outfile") == "-",
                "has_seed_option": "seed" in options,
            }
        )
    return commands


def collected_test_count() -> int:
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    match = re.search(r"(\d+) tests? collected", result.stdout)
    if match is None:
        raise RuntimeError("Could not parse the pytest collection count")
    return int(match.group(1))


def build_report() -> dict[str, object]:
    commands = command_inventory()
    return {
        "nwkit_version": nwkit.__version__,
        "git_commit": git_value("rev-parse", "HEAD"),
        "git_commit_short": git_value("rev-parse", "--short=12", "HEAD"),
        "python_version": sys.version.split()[0],
        "functional_subcommand_count": len(commands),
        "functional_subcommands": commands,
        "reads_stdin_by_default_count": sum(
            bool(command["reads_stdin_by_default"]) for command in commands
        ),
        "writes_stdout_by_default_count": sum(
            bool(command["writes_stdout_by_default"]) for command in commands
        ),
        "seeded_subcommands": [
            command["name"] for command in commands if command["has_seed_option"]
        ],
        "pytest_test_count": collected_test_count(),
        "pytest_module_count": len(list((PROJECT_ROOT / "tests").glob("test_*.py"))),
    }


def main() -> None:
    print(json.dumps(build_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
