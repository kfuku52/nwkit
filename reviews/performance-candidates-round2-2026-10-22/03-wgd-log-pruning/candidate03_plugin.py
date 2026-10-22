"""Scope the experimental log-pruning method to this pytest process."""

from benchmark import build_candidate


def pytest_configure(config):
    from nwkit.wgd_count_model import CountLikelihood

    CountLikelihood._prune_log = build_candidate()
