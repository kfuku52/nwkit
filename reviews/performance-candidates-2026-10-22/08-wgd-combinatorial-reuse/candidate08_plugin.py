"""Enable only the candidate 08 transition kernel during existing regressions."""


def pytest_configure(config):
    from benchmark import build_candidate

    from nwkit import wgd_count_model

    wgd_count_model._birth_death_log_transition = build_candidate()
