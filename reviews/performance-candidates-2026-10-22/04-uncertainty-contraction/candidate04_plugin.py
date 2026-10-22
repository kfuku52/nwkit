"""Install the isolated projection prototype for existing regression tests."""

from benchmark import prototype

_context = None


def pytest_configure(config):
    global _context
    _context = prototype()
    _context.__enter__()


def pytest_unconfigure(config):
    if _context is not None:
        _context.__exit__(None, None, None)
