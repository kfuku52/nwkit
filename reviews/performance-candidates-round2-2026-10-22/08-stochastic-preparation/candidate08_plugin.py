"""Inject candidate08 before collecting existing focused regressions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import benchmark  # noqa: E402


def pytest_configure(config):
    for (module, name), function in benchmark.CANDIDATE.items():
        setattr(module, name, function)


def pytest_unconfigure(config):
    for (module, name), function in benchmark.ORIGINAL.items():
        setattr(module, name, function)
