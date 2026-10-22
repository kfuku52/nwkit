"""Install only candidate 03's experimental helper during pytest startup."""


def pytest_configure(config):
    from benchmark import candidate_function

    from nwkit import regress

    regress._fit_model = candidate_function()
