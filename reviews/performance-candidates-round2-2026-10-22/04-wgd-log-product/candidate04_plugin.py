"""Use candidate 04 for existing snapshot WGD regressions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import benchmark


def pytest_configure(config):
    benchmark.model._log_matrix_product = benchmark.batched_product


def pytest_unconfigure(config):
    benchmark.model._log_matrix_product = benchmark.ORIGINAL
