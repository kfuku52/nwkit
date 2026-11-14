"""Apply the guarded local experiment to existing tests in this process only."""

import importlib

from benchmark import guarded_is_diagonal


def pytest_configure(config):
    for name in (
        "nwkit.gaussian",
        "nwkit.regress",
        "nwkit.measurement_error",
        "nwkit.ordinary_regression",
    ):
        module = importlib.import_module(name)
        module.is_diagonal = guarded_is_diagonal
