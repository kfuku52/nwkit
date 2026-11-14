"""Cancellation and locality checks for shift result validation."""

import os

import pandas as pd
import pytest

from nwkit.cli import main
from tests import test_shift as shift_support

shift_inputs = shift_support.shift_inputs


def cancelling_backend(directory, args, discrepancy):
    result = shift_support.fake_backend(directory, args)
    model = pd.read_csv(directory / "model.tsv", sep="\t")
    model[["alpha", "intercept"]] = 0
    model.to_csv(directory / "model.tsv", sep="\t", index=False)
    (directory / "shifts.tsv").write_text(
        "clade\tmean_effect\toptimum_effect\nt0/t1\t1e16\tNA\nt1\t-1e16\tNA\n"
    )
    tips = pd.read_csv(directory / "tips.tsv", sep="\t")
    tips["predicted"] = [1e16, discrepancy, 0, 0]
    tips["residual"] = tips.observed - tips.predicted
    tips["optimum"] = "NA"
    tips.to_csv(directory / "tips.tsv", sep="\t", index=False)
    return result


def test_cancellation_roundoff_is_not_invalid_fit(shift_inputs, monkeypatch):
    monkeypatch.setattr(
        "nwkit.shift.run_backend", lambda d, a: cancelling_backend(d, a, 1.0)
    )
    main(shift_inputs)


def test_large_cancellation_error_still_rejected(shift_inputs, monkeypatch):
    monkeypatch.setattr(
        "nwkit.shift.run_backend", lambda d, a: cancelling_backend(d, a, 1000.0)
    )
    with pytest.raises(ValueError, match="inconsistent"):
        main(shift_inputs)


def test_unrelated_large_effect_does_not_hide_constraint_error():
    from nwkit.shift_convergence import shared_regime_rows
    from nwkit.shift_math import ancestral_effect_scales
    from nwkit.util import assign_branch_ids, read_tree

    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1);", "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    effects = {1: {"optimum_effect": 1.0}, 2: {"optimum_effect": 1e100}}
    scales = ancestral_effect_scales(tree, ids, 0.0, effects, "optimum_effect")
    by_id = {ids[node]: scale for node, scale in scales.items()}
    rows = [
        {"regime": "baseline", "branch_id": 0, "optimum": 0.0},
        {"regime": "baseline", "branch_id": 1, "optimum": 1e-10},
        {"regime": "other", "branch_id": 2, "optimum": 1e100},
    ]
    with pytest.raises(ValueError, match="optima disagree"):
        shared_regime_rows(rows, by_id)


def test_shared_optimum_cancellation_still_allowed():
    from nwkit.shift_convergence import shared_regime_rows

    rows = [
        {"regime": "baseline", "branch_id": 0, "optimum": 0.0},
        {"regime": "baseline", "branch_id": 4, "optimum": 1.0},
    ]
    assert len(shared_regime_rows(rows, {0: 0.0, 4: 1e16})) == 1


def test_shared_scale_includes_intermediate_ancestral_values():
    from nwkit.shift_math import ancestral_effect_scales
    from nwkit.util import assign_branch_ids, read_tree

    tree = read_tree("((A:1,B:1):1,C:2);", "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    effects = {
        ids[node]: {"mean_effect": 1.0} for node in tree.traverse() if not node.is_root
    }
    scales = ancestral_effect_scales(tree, ids, 1.0, effects, "mean_effect")
    assert scales[next(node for node in tree.leaves() if node.name == "A")] == 3.0


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
def test_convergence_unit_changes_preserve_fit(tmp_path):
    import json
    import math

    import numpy as np

    from tests.test_shift_convergence import convergent_command

    baseline = None
    for scale, time in [(1.0, 1.0), (1e-4, 1e3), (1e4, 1e-3)]:
        directory = tmp_path / str(scale)
        directory.mkdir()
        command = convergent_command(directory, "OUfixedRoot", True)
        data = pd.read_csv(directory / "traits.tsv", sep="\t")
        data[["value", "se"]] *= scale
        data.to_csv(directory / "traits.tsv", sep="\t", index=False)
        tree = directory / "tree.nwk"
        tree.write_text(tree.read_text().replace(":1", f":{time}"))
        main(command)
        model = json.loads((directory / "model.json").read_text())
        if baseline is None:
            baseline = model
        expected, actual = baseline["parameters"], model["parameters"]
        assert actual["alpha"] * time == pytest.approx(expected["alpha"], rel=1e-5)
        assert actual["sigma2"] * time / scale**2 == pytest.approx(
            expected["sigma2"], rel=1e-5
        )
        assert actual["log_likelihood"] + 8 * math.log(scale) == pytest.approx(
            expected["log_likelihood"], abs=1e-5
        )
        np.testing.assert_allclose(
            [row["predicted"] / scale for row in model["tip_predictions"]],
            [row["predicted"] for row in baseline["tip_predictions"]],
            rtol=1e-5,
        )
        assert model["convergence"]["merges"] == baseline["convergence"]["merges"]


def test_direct_shared_optimum_survives_cancelling_deltas(
    shift_inputs, tmp_path, monkeypatch
):
    import json
    import math

    from tests.test_shift_convergence import nested_convergence_backend

    def backend(directory, args):
        result = nested_convergence_backend(directory, args)
        model = pd.read_csv(directory / "model.tsv", sep="\t")
        model["intercept"] = 1.0
        model.to_csv(directory / "model.tsv", sep="\t", index=False)
        (directory / "unconstrained-model.tsv").write_bytes(
            (directory / "model.tsv").read_bytes()
        )
        delta = 1e16 / -math.expm1(-1.0)
        back = -delta * -math.expm1(-0.5)
        pd.DataFrame(
            {
                "clade": ["t0/t1", "t1"],
                "mean_effect": [1e16, back],
                "optimum_effect": [delta, -delta],
            }
        ).to_csv(directory / "shifts.tsv", sep="\t", index=False)
        tips = pd.read_csv(directory / "tips.tsv", sep="\t")
        tips["predicted"] = [1 + 1e16, 1 + 1e16 + back, 1, 1]
        tips["residual"] = tips.observed - tips.predicted
        tips["optimum"] = [1 + delta, 1, 1, 1]
        tips.to_csv(directory / "tips.tsv", sep="\t", index=False)
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    main([*shift_inputs, "--convergence"])
    model = json.loads((tmp_path / "model.json").read_text())
    assert (
        next(row for row in model["tip_predictions"] if row["leaf_name"] == "B")[
            "optimum"
        ]
        == 1.0
    )
    assert model["regime_parameters"][0]["optimum"] == 1.0


def test_direct_optima_cannot_disagree_within_regime():
    from nwkit.shift_results import align_regime_optima

    regimes = [{"regime": "baseline", "branch_id": 0, "optimum": 1.0}]
    tips = [{"regime": "baseline", "branch_id": 4, "optimum": 1.01}]
    with pytest.raises(ValueError, match="different optima"):
        align_regime_optima(regimes, tips)


def test_regime_without_extant_tips_retains_its_optimum():
    from nwkit.shift_results import align_regime_optima

    regimes = [
        {"regime": "baseline", "branch_id": 0, "optimum": 1.0},
        {"regime": "historical", "branch_id": 1, "optimum": 3.0},
        {"regime": "current", "branch_id": 2, "optimum": 2.0},
    ]
    tips = [{"regime": "current", "branch_id": 4, "optimum": 2.0000000000000004}]
    align_regime_optima(regimes, tips)
    assert regimes[1]["optimum"] == 3.0
    assert regimes[2]["optimum"] == tips[0]["optimum"]
