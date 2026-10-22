"""Load the exact isolated prototype into existing pytest regressions."""

from benchmark import prototype

_context = None


def pytest_configure(config):
    global _context
    _context = prototype()
    _context.__enter__()


def pytest_unconfigure(config):
    if _context is not None:
        _context.__exit__(None, None, None)
