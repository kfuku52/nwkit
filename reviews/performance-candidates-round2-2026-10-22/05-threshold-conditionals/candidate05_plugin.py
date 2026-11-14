"""Install only the per-fit threshold conditional-preparation prototype."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark


def pytest_configure(config):
    benchmark.module.compute_threshold_marginals = benchmark.CANDIDATE


def pytest_unconfigure(config):
    benchmark.module.compute_threshold_marginals = benchmark.ORIGINAL
