"""Install the fit-local moment-product prototype for full-OU regressions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import benchmark


def pytest_configure(config):
    benchmark.module._moment_vector = benchmark.consecutive_moments


def pytest_unconfigure(config):
    benchmark.module._moment_vector = benchmark.ORIGINAL
