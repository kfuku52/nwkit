"""Constrained means, shared regime maps and full-selection bootstrap checks."""

import json
import math
import os

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.shift_convergence import decode_groups
from nwkit.shift_reference import evaluate_shift_model
from nwkit.util import assign_branch_ids, read_tree
from tests import test_shift as shift_support

shift_inputs = shift_support.shift_inputs


def nested_convergence_backend(directory, args):
    result = shift_support.fake_backend(directory, args)
    (directory / "unconstrained-model.tsv").write_bytes(
        (directory / "model.tsv").read_bytes()
    )
    (directory / "convergence.tsv").write_text("clades\n@root;t1\nt0/t1\n")
    shifts = pd.read_csv(directory / "shifts.tsv", sep="\t")
    delta = float(shifts.iloc[0].optimum_effect)
    mean = -delta * -math.expm1(-0.5)
    shifts.loc[1, ["mean_effect", "optimum_effect"]] = [mean, -delta]
    shifts.to_csv(directory / "shifts.tsv", sep="\t", index=False)
    tips = pd.read_csv(directory / "tips.tsv", sep="\t")
    tips.loc[tips.token == "t1", ["predicted", "residual", "optimum"]] = [
        2.4 + mean,
        2 - (2.4 + mean),
        2,
    ]
    tips.to_csv(directory / "tips.tsv", sep="\t", index=False)
    return result


def test_nested_shift_returns_to_background(shift_inputs, tmp_path, monkeypatch):
    monkeypatch.setattr("nwkit.shift.run_backend", nested_convergence_backend)
    output = tmp_path / "regimes.tsv"
    main([*shift_inputs, "--convergence", "-o", str(output)])
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["convergence"]["groups"] == [
        {"regime": "baseline", "branch_ids": [0, 4]},
        {"regime": "shift_1", "branch_ids": [1]},
    ]
    assert len(model["regime_parameters"]) == 2
    assert (
        pd.read_csv(output, sep="\t").set_index("branch_id").loc[4, "regime"]
        == "baseline"
    )
    assert model["shift_effects"][1]["parent_regime"] == "shift_1"


@pytest.mark.parametrize(
    "values",
    [
        ["@root;t0;t0"],
        ["t0"],
        ["@root;unknown"],
        ["@root", "rootclade;t0"],
        ["@root", ""],
    ],
)
def test_invalid_convergence_groups(values):
    with pytest.raises(ValueError, match="[Cc]onvergence"):
        decode_groups(values, {"t0": 1, "rootclade": 0}, [1])


@pytest.mark.parametrize("criterion", ["mBIC", "pBICess"])
def test_unsupported_convergence_criterion_preflight(
    shift_inputs, monkeypatch, criterion
):
    def backend(*args):
        pytest.fail("Must reject before R runs")

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    with pytest.raises(ValueError, match="supports only"):
        main([*shift_inputs, "--convergence", "--criterion", criterion])


def test_unequal_optima_cannot_be_merged(shift_inputs, tmp_path, monkeypatch):
    def backend(directory, args):
        result = nested_convergence_backend(directory, args)
        (directory / "convergence.tsv").write_text("clades\n@root;t0/t1;t1\n")
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    model = tmp_path / "model.json"
    model.write_text("original")
    with pytest.raises(ValueError, match="optima disagree"):
        main([*shift_inputs, "--convergence"])
    assert model.read_text() == "original"


def convergent_command(tmp_path, root_model, with_error):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);")
    trait = tmp_path / "traits.tsv"
    values = [6, 6.2, 1, 1.2, 5.8, 6.1, 0.9, 1.1]
    pd.DataFrame(
        {
            "leaf_name": list("ABCDEFGH"),
            "value": values,
            "se": [0, 0.02, 0.03, 0.01, 0.02, 0.01, 0.03, 0.02],
        }
    ).iloc[::-1].to_csv(trait, sep="\t", index=False)
    return [
        "shift",
        "--selection",
        "ic",
        "-i",
        str(tree),
        "--trait",
        str(trait),
        "--state-column",
        "value",
        "--model-out",
        str(tmp_path / "model.json"),
        "--rscript",
        os.environ["NWKIT_TEST_RSCRIPT"],
        "--max-shifts",
        "2",
        "--criterion",
        "BIC",
        "--search-strategy",
        "exhaustive",
        "--root-model",
        root_model,
        "--convergence",
        "--fit-out",
        str(tmp_path / "fit.rds"),
        "-o",
        str(tmp_path / "regimes.tsv"),
        *(["--standard-error-column", "se"] if with_error else []),
    ]


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
@pytest.mark.parametrize("root_model", ["OUfixedRoot", "OUrandomRoot"])
@pytest.mark.parametrize("with_error", [False, True])
def test_constrained_reference_and_bootstrap(tmp_path, root_model, with_error):
    command = convergent_command(tmp_path, root_model, with_error)
    main([*command, "--bootstrap", "2", "--bootstrap-seed", "71"])
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["convergence_searched"]
    assert model["convergence"]["merges"] >= 1
    assert len(model["regime_parameters"]) < len(model["shift_branch_ids"]) + 1
    assert model["parameters"]["score"] < model["unconstrained_parameters"]["score"]
    tree = read_tree(str(tmp_path / "tree.nwk"), "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    by_id = {branch: node for node, branch in ids.items()}
    tips = {row["leaf_name"]: row for row in model["tip_predictions"]}
    ordered = [tips[node.name] for node in tree.leaves()]
    params = model["parameters"]
    mean, covariance, loglik = evaluate_shift_model(
        tree,
        observations=[row["observed"] for row in ordered],
        alpha=params["alpha"],
        sigma2=params["sigma2"],
        intercept=params["intercept"],
        root_model=root_model,
        mean_effects={
            by_id[row["branch_id"]]: row["mean_effect"]
            for row in model["shift_effects"]
        },
        standard_errors=[row["standard_error"] for row in ordered],
    )
    np.testing.assert_allclose(
        mean, [row["predicted"] for row in ordered], rtol=1e-8, atol=1e-8
    )
    assert loglik == pytest.approx(params["log_likelihood"], abs=1e-6)
    assert np.all(np.linalg.eigvalsh(covariance) > 0)
    bootstrap = model["bootstrap"]
    assert bootstrap["selection"] == "shift_and_convergence"
    assert bootstrap["successful"] == 2
    assert len(bootstrap["successful_convergence_groups"]) == 2
    assert (
        sum(row["count"] for row in bootstrap["shared_optimum_partition_frequencies"])
        == 2
    )
    if root_model == "OUfixedRoot" and not with_error:
        main([*command, "--bootstrap", "2", "--bootstrap-seed", "71"])
        assert (
            json.loads((tmp_path / "model.json").read_text())["bootstrap"] == bootstrap
        )
        main(
            [
                "asr",
                "-i",
                str(tmp_path / "tree.nwk"),
                "--trait",
                str(tmp_path / "traits.tsv"),
                "--state-column",
                "value",
                "--trait-type",
                "continuous",
                "--model",
                "OUM",
                "--regime-map",
                str(tmp_path / "regimes.tsv"),
                "-o",
                str(tmp_path / "asr.tsv"),
            ]
        )
        assert (tmp_path / "asr.tsv").is_file()


def test_convergence_can_leave_all_regimes_separate(
    shift_inputs, tmp_path, monkeypatch
):
    def backend(directory, args):
        result = shift_support.fake_backend(directory, args)
        (directory / "unconstrained-model.tsv").write_bytes(
            (directory / "model.tsv").read_bytes()
        )
        (directory / "convergence.tsv").write_text("clades\n@root\nt0/t1\nt1\n")
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    main([*shift_inputs, "--convergence"])
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["convergence"]["merges"] == 0
    assert len(model["regime_parameters"]) == 3


def test_unidentifiable_convergence_rejected(shift_inputs, tmp_path, monkeypatch):
    def backend(directory, args):
        result = nested_convergence_backend(directory, args)
        model = pd.read_csv(directory / "model.tsv", sep="\t")
        model["alpha"] = 0
        model.to_csv(directory / "model.tsv", sep="\t", index=False)
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    with pytest.raises(ValueError, match="identifiable OU optima"):
        main([*shift_inputs, "--convergence"])
    assert not (tmp_path / "model.json").exists()


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
def test_no_shifts_convergence(tmp_path):
    main([*convergent_command(tmp_path, "OUfixedRoot", False), "--max-shifts", "0"])
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["shift_branch_ids"] == []
    assert model["convergence"]["groups"] == [{"regime": "baseline", "branch_ids": [0]}]
    assert model["convergence"]["merges"] == 0


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
@pytest.mark.parametrize("criterion", ["pBIC", "AICc"])
def test_other_supported_convergence_criteria(tmp_path, criterion):
    main(
        [*convergent_command(tmp_path, "OUfixedRoot", False), "--criterion", criterion]
    )
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["criterion"] == criterion
    assert model["parameters"]["alpha"] > 0


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
def test_partial_failure_in_convergence_stage(tmp_path, monkeypatch, capsys):
    from nwkit import shift_backend

    script = shift_backend.R_SCRIPT.replace(
        "converge <- function(model) {",
        "convergence.calls <- 0L\nconverge <- function(model) {\n"
        "convergence.calls <<- convergence.calls + 1L\n"
        "if (convergence.calls == 2L) stop('a \"quoted\" failure\\nwith newline')\n",
    )
    monkeypatch.setattr(shift_backend, "R_SCRIPT", script)
    main(
        [
            *convergent_command(tmp_path, "OUfixedRoot", False),
            "--bootstrap",
            "2",
            "--bootstrap-seed",
            "71",
        ]
    )
    boot = json.loads((tmp_path / "model.json").read_text())["bootstrap"]
    assert (boot["successful"], boot["failed"]) == (1, 1)
    assert boot["failure_messages"] == [
        {"message": 'a "quoted" failure\nwith newline', "count": 1}
    ]
    assert len(boot["successful_convergence_groups"]) == 1
    assert boot["shared_optimum_partition_frequencies"][0]["frequency"] == 1
    assert "1/2 bootstrap refits failed" in capsys.readouterr().err


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
def test_rds_retains_original_fit(tmp_path):
    import subprocess

    main(convergent_command(tmp_path, "OUfixedRoot", False))
    script = tmp_path / "inspect.R"
    script.write_text("""
fit <- readRDS("fit.rds")
original <- attr(fit, "nwkit.unconstrained.fit")
stopifnot(inherits(fit, "l1ou"), inherits(original, "l1ou"),
          isTRUE(fit$convergent), !isTRUE(original$convergent),
          original$score > fit$score,
          length(original$shift.configuration) == length(fit$shift.configuration))
""")
    subprocess.run(
        [os.environ["NWKIT_TEST_RSCRIPT"], "--vanilla", str(script)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
