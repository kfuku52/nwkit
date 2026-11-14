"""Force the ER algorithm for tiny existing references as well as large cases."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import benchmark


def pytest_configure(config):
    benchmark.asr._compute_log_inside_likelihoods = benchmark.build_candidate(2)


def pytest_unconfigure(config):
    benchmark.asr._compute_log_inside_likelihoods = benchmark.ORIGINAL
