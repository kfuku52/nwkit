"""Install candidate 06 before existing test modules are collected."""

import importlib.util
from pathlib import Path

_context = None


def pytest_configure(config):
    global _context
    path = Path(__file__).with_name("benchmark.py")
    spec = importlib.util.spec_from_file_location("candidate06_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _context = module.prototype()
    _context.__enter__()


def pytest_unconfigure(config):
    if _context is not None:
        _context.__exit__(None, None, None)
