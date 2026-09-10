"""The evidence auditor must reject altered results, not merely valid JSON."""

import copy
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from shift_alpha_design import cases, generate_case, protocol  # noqa: E402
from shift_calibration_audit import audit_record, audit_records  # noqa: E402
from shift_simulation_cases import rate  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.shift_candidates import tip_groups  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


@pytest.fixture(scope="module")
def evidence():
    spec = protocol(1, 0, 20260917)
    spec.update(calibration_replicates=19, calibration_level=0.05)
    case = next(cases(spec))
    tree, truth = generate_case(case)
    search = CalibratedSearch(read_tree(tree, "auto", True, quiet=True))
    fit = search.fit(
        truth["observations"], seed=case["seed"] + 9000000000, replicates=19
    )
    row = dict(
        case=case,
        tree=tree,
        truth=truth,
        fit=fit,
        status="completed",
        any_shift=bool(fit["model"]["shift_branch_ids"]),
        partition_recovered=sorted(
            map(
                list,
                tip_groups(
                    search.tree,
                    fit["model"]["shift_branch_ids"],
                    fit["model"]["groups"],
                ),
            )
        )
        == truth["shared_partition"],
        mean_rmse=float(
            np.sqrt(
                np.mean(
                    (np.asarray(fit["predicted"]) - list(truth["tip_mean"].values()))
                    ** 2
                )
            )
        ),
    )
    summary = [
        dict(
            family=case["family"],
            scenario=case["scenario"],
            root_model=case["root_model"],
            attempted=1,
            completed=1,
            any_shift=rate(int(row["any_shift"]), 1),
            partition_recovered=rate(int(row["partition_recovered"]), 1),
            alpha_limit_supported=rate(int(fit["alpha_limit_supported"]), 1),
        )
    ]
    return spec, row, summary, search


def test_valid_evidence_reconstructs(evidence):
    spec, row, summary, _ = evidence
    assert audit_records([row], spec, summary) == summary


@pytest.mark.parametrize(
    "change",
    [
        "any_shift",
        "partition_recovered",
        "mean_rmse",
        "summary",
        "statistic",
        "alpha_profile",
        "seed",
    ],
)
def test_modified_evidence_is_rejected(evidence, change):
    spec, original, original_summary, search = evidence
    row, summary = copy.deepcopy(original), copy.deepcopy(original_summary)
    if change in ("any_shift", "partition_recovered"):
        row[change] = not row[change]
    elif change == "mean_rmse":
        row[change] += 1
    elif change == "summary":
        summary[0]["any_shift"]["count"] += 1
    elif change == "statistic":
        row["fit"]["tests"][0]["statistic"] += 2
    elif change == "alpha_profile":
        row["fit"]["alpha_profile"][3]["log_likelihood"] -= 1
    else:
        row["fit"]["seed"] += 1
    with pytest.raises(ValueError):
        if change == "summary":
            audit_records([row], spec, summary)
        else:
            audit_record(row, spec, search)


def test_changed_alpha_status_rejected(evidence):
    spec, original, _, search = evidence
    row = copy.deepcopy(original)
    row["fit"]["alpha_status"] = "invalid"
    with pytest.raises(ValueError, match="Alpha estimate"):
        audit_record(row, spec, search)


def test_seeded_bootstrap_replay_rejects_plausible_altered_p(evidence):
    spec, original, _, search = evidence
    row = copy.deepcopy(original)
    p = row["fit"]["tests"][-1]["p_value"]
    assert p > spec["calibration_level"]
    row["fit"]["tests"][-1]["p_value"] = 0.9 if p != 0.9 else 0.8
    with pytest.raises(ValueError):
        audit_record(row, spec, search, replay_bootstrap=True)


def test_replay_checks_nuisance_probabilities_even_when_upper_bound_is_one(evidence):
    spec, original, _, search = evidence
    row = copy.deepcopy(original)
    test = row["fit"]["tests"][0]
    assert test["p_value_kind"] == "conservative_upper_bound"
    value = 0.9 if test["p_value_lower_bound"] != 0.9 else 0.8
    for entry in test["null_alpha_evaluations"]:
        entry["p_value"] = value
    test["p_value_lower_bound"] = value
    with pytest.raises(ValueError, match="seeded replay"):
        audit_record(row, spec, search, replay_bootstrap=True)
