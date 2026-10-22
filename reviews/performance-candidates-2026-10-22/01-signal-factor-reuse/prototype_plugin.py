"""Load the scoped prototype before pytest collects the existing suite."""

from benchmark import candidate_permutation_test


def pytest_configure(config):
    import nwkit.signal
    import nwkit.signal_stats

    nwkit.signal_stats.permutation_test = candidate_permutation_test
    nwkit.signal.permutation_test = candidate_permutation_test
