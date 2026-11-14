"""Apply the final-whitening PCA prototype before tests import its alias."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark  # noqa: E402


def pytest_configure(config):
    from nwkit import pca

    benchmark.module.fit_pca = benchmark.CANDIDATE
    pca.fit_pca = benchmark.CANDIDATE


def pytest_unconfigure(config):
    from nwkit import pca

    benchmark.module.fit_pca = benchmark.ORIGINAL
    pca.fit_pca = benchmark.ORIGINAL
