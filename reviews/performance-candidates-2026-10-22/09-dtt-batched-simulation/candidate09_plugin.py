"""Install only candidate09's scoped DTT chunk implementation for pytest."""

import benchmark


def pytest_configure(config):
    import nwkit.disparity as disparity

    disparity._simulate_chunk = benchmark.candidate_chunk
