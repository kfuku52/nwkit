"""Install the covariance prototype for existing focused regressions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import benchmark


def pytest_configure(config):
    # Force the exact same block algorithm even for existing tiny reference
    # fixtures, so these checks exercise the prototype rather than its guard.
    benchmark.gt.GaussianTreeProcess.covariance = benchmark._candidate_function(0)


def pytest_unconfigure(config):
    benchmark.gt.GaussianTreeProcess.covariance = benchmark.ORIGINAL
