"""Install exactly the scoped observed-covariance prototype during pytest."""

from benchmark import prototype

PATCH = None


def pytest_configure(config):
    global PATCH
    PATCH = prototype()
    PATCH.__enter__()


def pytest_unconfigure(config):
    if PATCH is not None:
        PATCH.__exit__(None, None, None)
