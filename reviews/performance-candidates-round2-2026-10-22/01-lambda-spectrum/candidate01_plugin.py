"""Patch only this pytest process, before test modules import lambda_fit."""

import benchmark

from nwkit import signal_stats


def pytest_configure(config):
    assert "/nwkit-perf-round2-20261022-rkzoGe/baseline/" in signal_stats.__file__
    signal_stats.lambda_fit = benchmark.spectral_lambda_fit
    from nwkit import signal

    signal.lambda_fit = benchmark.spectral_lambda_fit
