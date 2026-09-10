"""Signal formulas, independent R references and real CLI contracts."""

from io import StringIO

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize
from scipy.stats import multivariate_normal

from nwkit.cli import main
from nwkit.evolution import build_evolutionary_covariance
from nwkit.signal_stats import bh_adjust, k_fit, k_statistic, lambda_fit, rate_fit
from nwkit.util import read_tree

TREE = "(((A:0.7,B:0.7):0.8,(C:0.9,D:0.9):0.6):0.5,((E:0.8,F:0.8):0.7,(G:1.1,H:1.1):0.4):0.5);"
VALUES = np.array([1, 1.5, 2, 1.8, 4, 4.2, 3.8, 4.5])
ERRORS = np.array([0.1, 0.2] * 4)


def covariance():
    return build_evolutionary_covariance(read_tree(TREE, 0, False), list("ABCDEFGH"))


def test_phytools_reference_k():
    # R phytools 2.3.0, using examples/signal; no R dependency at runtime.
    # Error-aware rate reoptimized with optimize(..., tol=1e-12).
    c = covariance()
    assert k_statistic(c, VALUES) == pytest.approx(2.16451564383, rel=1e-10)
    k, _ = k_fit(c, VALUES - VALUES[0], ERRORS)
    assert k == pytest.approx(2.141105106703, rel=1e-8)
    assert k_statistic(c * 17, VALUES * -3 + 200) == pytest.approx(
        k_statistic(c, VALUES)
    )


@pytest.mark.parametrize("errors", [ERRORS, np.zeros(8), np.array([0.0, 0.2] * 4)])
def test_likelihood_against_independent_direct_optimization(errors):
    c = covariance()
    y = VALUES - VALUES[0]
    likelihood, rate, mean = rate_fit(c, y, errors)

    def objective(params):
        return -multivariate_normal.logpdf(
            y,
            mean=np.full(8, params[0]),
            cov=np.exp(params[1]) * c + np.diag(errors**2),
        )

    independent = minimize(objective, [2, -1], method="BFGS", tol=1e-8)
    assert likelihood == pytest.approx(-independent.fun, abs=1e-7)
    assert rate == pytest.approx(np.exp(independent.x[1]), rel=1e-5)
    assert mean == pytest.approx(independent.x[0], abs=1e-5)


def test_lambda_boundary_and_profile():
    fit = lambda_fit(covariance(), VALUES, np.zeros(8))
    assert fit["estimate"] == 1
    assert fit["status"] == "boundary"
    assert fit["log_likelihood"] == pytest.approx(-10.2136362027)
    assert fit["null_log_likelihood"] == pytest.approx(-13.5439938997)
    assert fit["ci_lower"] == pytest.approx(0.395258805358)
    assert fit["ci_upper"] == 1
    assert fit["p_value"] == pytest.approx(0.00985613398069)


def test_star_has_no_identifiable_lambda():
    assert (
        lambda_fit(np.eye(4), np.arange(4), np.zeros(4))["status"]
        == "unidentifiable_lambda"
    )
    assert k_statistic(np.eye(4), np.arange(4)) == pytest.approx(1)


def test_bh():
    assert bh_adjust([0.04, 0.001, 0.03]) == pytest.approx([0.04, 0.003, 0.04])
    assert len(bh_adjust([])) == 0


def run_cli(tmp_path, capsys, extra=(), frame=None, tree=TREE):
    path = tmp_path / "traits.tsv"
    if frame is None:
        frame = pd.DataFrame({"leaf_name": list("ABCDEFGH"), "x": VALUES, "se": ERRORS})
    frame.to_csv(path, sep="\t", index=False)
    main(
        [
            "signal",
            "-i",
            tree,
            "--trait",
            str(path),
            "--columns",
            "x",
            "--n-sim",
            "9",
            *extra,
        ]
    )
    return pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")


def test_cli_se_and_reproducible_permutations(tmp_path, capsys):
    result = run_cli(tmp_path, capsys, ["--standard-error-column", "se"])
    again = run_cli(tmp_path, capsys, ["--standard-error-column", "se"])
    pd.testing.assert_frame_equal(result, again)
    assert result.iloc[0].estimate == pytest.approx(2.141105106703, rel=1e-8)
    assert result.p_value.between(0.1, 1).iloc[0]
    assert result.measurement_error.tolist() == ["yes", "yes"]


def test_missing_constant_and_trait_order(tmp_path, capsys):
    frame = pd.DataFrame(
        {"leaf_name": list("ABCDEFGH"), "x": VALUES, "constant": 1, "missing": np.nan}
    )
    first = run_cli(tmp_path, capsys, ["--columns", "x,constant,missing"], frame)
    assert (
        first.status.tolist()[2:] == ["constant_trait"] * 2 + ["insufficient_taxa"] * 2
    )
    second = run_cli(tmp_path, capsys, ["--columns", "missing,x"], frame)
    assert first.iloc[0].p_value == second.iloc[2].p_value
    assert first.iloc[0].num_tests == 1


def test_missing_uses_retained_mrca(tmp_path, capsys):
    frame = pd.DataFrame(
        {"leaf_name": list("ABCDEFGH"), "x": [1, 2, 3, 4, None, None, None, None]}
    )
    result = run_cli(tmp_path, capsys, ["--test", "no"], frame)
    reduced = run_cli(
        tmp_path,
        capsys,
        ["--test", "no"],
        frame.iloc[:4],
        "((A:0.7,B:0.7):0.8,(C:0.9,D:0.9):0.6);",
    )
    assert result.estimate.tolist() == pytest.approx(reduced.estimate.tolist())
    assert result.num_taxa.tolist() == [4, 4]
    assert result.p_value.isna().all()


def test_singular_covariance_and_root_polytomy(tmp_path, capsys):
    frame = pd.DataFrame({"leaf_name": list("ABC"), "x": [1, 2, 3]})
    result = run_cli(
        tmp_path, capsys, ["--input-rooted", "yes"], frame, "((A:0,B:0):1,C:1);"
    )
    assert set(result.status) == {"singular_covariance"}
    result = run_cli(
        tmp_path, capsys, ["--input-rooted", "yes"], frame, "(A:1,B:1,C:1);"
    )
    assert result.iloc[0].estimate == pytest.approx(1)
    assert result.iloc[0].p_value == 1
    assert result.iloc[1].status == "unidentifiable_lambda"


@pytest.mark.parametrize(
    "extra",
    [
        ["--ci-level", "nan"],
        ["--n-sim", "0"],
        ["--seed", "-1"],
        ["--columns", "x,x"],
        ["--input-rooted", "no"],
    ],
)
def test_invalid_options(tmp_path, capsys, extra):
    with pytest.raises(ValueError):
        run_cli(tmp_path, capsys, extra)


def test_output_cannot_replace_input(tmp_path):
    path = tmp_path / "traits.tsv"
    text = "leaf_name\tx\nA\t1\nB\t2\nC\t3\n"
    path.write_text(text)
    with pytest.raises(ValueError):
        main(
            [
                "signal",
                "-i",
                "((A:1,B:1):1,C:2);",
                "--trait",
                str(path),
                "--columns",
                "x",
                "-o",
                str(path),
            ]
        )
    assert path.read_text() == text


def test_interior_lambda_with_errors_independent_optimization():
    c = covariance()
    y = (VALUES + np.array([1, 4, 2, 3, 1.5, 4.2, 1.8, 3.8])) / 2
    y -= y[0]
    diagonal = np.diag(np.diag(c))
    fit = lambda_fit(c, y, ERRORS)

    def objective(params):
        matrix = np.exp(params[1]) * (diagonal + params[2] * (c - diagonal)) + np.diag(
            ERRORS**2
        )
        return -multivariate_normal.logpdf(y, mean=np.full(8, params[0]), cov=matrix)

    oracle = minimize(
        objective,
        [1, -1, 0.5],
        method="L-BFGS-B",
        bounds=[(None, None), (None, None), (0, 1)],
        options={"ftol": 1e-14, "gtol": 1e-8},
    )
    assert fit["status"] == "ok"
    assert fit["estimate"] == pytest.approx(oracle.x[2], abs=1e-5)
    assert fit["log_likelihood"] == pytest.approx(-oracle.fun, abs=1e-7)
    assert fit["ci_lower"] <= fit["estimate"] <= fit["ci_upper"]


def test_zero_rate_and_singular_zero_rate():
    c = covariance()
    fit = lambda_fit(c, VALUES, np.full(8, 100.0))
    assert fit["status"] == "unidentifiable_lambda"
    assert fit["sigma2"] == 0
    assert "p_value" not in fit
    with pytest.raises(ValueError, match="Unbounded"):
        rate_fit(c, VALUES, np.array([0.0] + [1.0] * 7))


def test_signal_stdin_and_audit(tmp_path, capsys, monkeypatch):
    import json

    audit = tmp_path / "audit.json"
    monkeypatch.setattr("sys.stdin", StringIO("leaf_name\tx\nA\t1\nB\t2\nC\t3\n"))
    main(
        [
            "signal",
            "-i",
            "((A:1,B:1):1,C:2);",
            "--trait",
            "-",
            "--columns",
            "x",
            "--n-sim",
            "3",
            "--audit",
            str(audit),
        ]
    )
    result = pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")
    assert len(result) == 2
    assert json.loads(audit.read_text())["status"] == "ok"
    with pytest.raises(ValueError, match="only one input"):
        main(["signal", "-i", "-", "--trait", "-", "--columns", "x"])


@pytest.mark.parametrize(
    "values,errors",
    [
        ([1, 2, 3], [0.1, -1, 0.1]),
        ([1, 2, 3], [0.1, None, 0.1]),
        ([1, "oops", 3], [0.1] * 3),
        ([1, np.inf, 3], [0.1] * 3),
    ],
)
def test_invalid_observations(tmp_path, capsys, values, errors):
    frame = pd.DataFrame({"leaf_name": list("ABC"), "x": values, "se": errors})
    with pytest.raises(ValueError):
        run_cli(
            tmp_path,
            capsys,
            ["--standard-error-column", "se"],
            frame,
            "((A:1,B:1):1,C:2);",
        )


def test_rate_far_below_noisy_observation_scale():
    y = np.array([1e-10, -1e-10, 1, -1, 2, -2, 3, -3])
    se = np.array([0, 0, 10, 10, 10, 10, 10, 10])
    likelihood, rate, mean = rate_fit(np.eye(8), y, se)
    assert rate == pytest.approx(1e-20, rel=1e-5, abs=0)
    assert likelihood == pytest.approx(23.744683036279252, abs=1e-8)
    assert mean == pytest.approx(0, abs=1e-20)


def test_extreme_finite_trait_ranges_do_not_produce_false_significance(
    tmp_path, capsys
):
    values = np.array([-1, -0.5, 0.5, 1, -0.8, -0.4, 0.4, 0.8]) * 1e308
    frame = pd.DataFrame({"leaf_name": list("ABCDEFGH"), "x": values})
    with np.errstate(all="raise"):
        result = run_cli(tmp_path, capsys, ["--method", "K"], frame)
    assert result.iloc[0].estimate == pytest.approx(
        k_statistic(covariance(), values / 1e308)
    )
    assert result.iloc[0].status == "ok"


def test_rate_units_avoid_intermediate_square_overflow():
    from nwkit.signal import _restore_units

    result = _restore_units({"sigma2": 0.5}, 1e200, 1e200, 0, 8)
    assert result["sigma2"] == pytest.approx(0.5e200)
    result = _restore_units({"sigma2": 0.5}, 1e-200, 1e-200, 0, 8)
    assert result["sigma2"] == pytest.approx(0.5e-200, rel=1e-15, abs=0)


@pytest.mark.parametrize("seed", range(4))
def test_error_rate_stress_against_direct_covariance(seed):
    from scipy.linalg import cho_solve
    from scipy.optimize import minimize_scalar

    rng = np.random.default_rng(seed)
    c = covariance()
    values = rng.normal(size=8)
    se = 10.0 ** rng.uniform(-3, 1, 8)
    if seed % 2:
        se[:2] = 0
    actual, rate, _ = rate_fit(c, values, se)

    def independent(lograte):
        matrix = np.exp(lograte) * c + np.diag(se**2)
        lower = np.linalg.cholesky(matrix)
        ones = cho_solve((lower, True), np.ones(8))
        mean = ones @ values / ones.sum()
        delta = values - mean
        return 0.5 * (
            8 * np.log(2 * np.pi)
            + 2 * np.log(np.diag(lower)).sum()
            + delta @ cho_solve((lower, True), delta)
        )

    grid = np.linspace(-25, 5, 80)
    best = int(np.argmin([independent(x) for x in grid]))
    reference = minimize_scalar(
        independent,
        bounds=(grid[max(0, best - 1)], grid[min(79, best + 1)]),
        method="bounded",
    )
    assert actual == pytest.approx(-reference.fun, abs=1e-7)
    assert rate == pytest.approx(np.exp(reference.x), rel=1e-4)


def test_k_permutations_invariant_to_input_order(tmp_path, capsys):
    frame = pd.DataFrame({"leaf_name": list("ABCDEFGH"), "x": VALUES})
    first = run_cli(tmp_path, capsys, ["--method", "K"], frame)
    reverse_tree = "(((H:1.1,G:1.1):0.4,(F:0.8,E:0.8):0.7):0.5,((D:0.9,C:0.9):0.6,(B:0.7,A:0.7):0.8):0.5);"
    second = run_cli(
        tmp_path, capsys, ["--method", "K"], frame.iloc[::-1], reverse_tree
    )
    pd.testing.assert_frame_equal(first, second)
