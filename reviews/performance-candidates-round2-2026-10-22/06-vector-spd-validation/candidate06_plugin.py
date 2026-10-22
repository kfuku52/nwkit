"""Install the candidate only within a dedicated pytest process."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark  # noqa: E402


def pytest_configure(config):
    benchmark.vg._propagate = benchmark.CANDIDATE


def pytest_unconfigure(config):
    benchmark.vg._propagate = benchmark.ORIGINAL
