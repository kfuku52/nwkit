"""Verify simulation lineage before publishing aggregate evidence."""

import csv
import itertools
import json
import math
from collections import defaultdict

import numpy as np
from shift_simulation_cases import SCENARIOS, selection_metrics, simulate, summarize

from nwkit.util import assign_branch_ids, read_tree

CELL_KEYS = ("tips", "scenario", "root_model", "standard_error", "effect", "mode")
TRUTH_KEYS = ("seed", "scenario", "tips", "root_model", "standard_error", "effect")


def equivalent(actual, expected):
    """Allow floating-point rounding, never missing keys or changed labels."""
    if isinstance(expected, dict):
        return (
            isinstance(actual, dict)
            and actual.keys() == expected.keys()
            and all(equivalent(actual[k], v) for k, v in expected.items())
        )
    if isinstance(expected, list):
        return (
            isinstance(actual, list)
            and len(actual) == len(expected)
            and all(equivalent(a, b) for a, b in zip(actual, expected, strict=True))
        )
    if isinstance(expected, float):
        return type(actual) in (int, float) and math.isclose(
            actual, expected, rel_tol=1e-12, abs_tol=1e-14
        )
    return type(actual) is type(expected) and actual == expected


def expected_cases(options):
    cells = itertools.product(
        options["tips"],
        SCENARIOS,
        ("OUfixedRoot", "OUrandomRoot"),
        options["standard_errors"],
        options["effects"],
    )
    for cell_index, (tips, scenario, root_model, se, effect) in enumerate(cells):
        for replicate in range(options["replicates"]):
            seed = int(
                np.random.SeedSequence(
                    [options["seed"], cell_index, replicate]
                ).generate_state(1)[0]
            )
            yield (
                f"c{cell_index:03d}-r{replicate:04d}",
                dict(
                    tips=tips,
                    scenario=scenario,
                    root_model=root_model,
                    se=se,
                    effect=effect,
                    seed=seed,
                ),
            )


def audit_inputs(directory, parameters):
    newick, expected = simulate(**parameters)
    truth = json.loads((directory / "truth.json").read_text())
    if (
        not equivalent(truth, expected)
        or (directory / "tree.nwk").read_text().strip() != newick
    ):
        raise ValueError(
            f"Generating truth or tree disagrees with manifest: {directory.name}"
        )
    with (directory / "traits.tsv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if reader.fieldnames != ["leaf_name", "value", "se"]:
            raise ValueError("Unexpected simulation trait columns")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("Malformed simulation trait row")
    if [row["leaf_name"] for row in rows] != truth["tip_names"]:
        raise ValueError("Simulation trait tips disagree with truth")
    for row, observed in zip(rows, truth["observations"], strict=True):
        if (
            float(row["value"]) != observed
            or float(row["se"]) != truth["standard_error"]
        ):
            raise ValueError("Simulation observations or SEs disagree with truth")
    return {"tree": newick, "truth": truth}


def audit_model(model, item, options, mode):
    truth = item["truth"]
    if (
        model["criterion"] != options["criterion"]
        or model["root_model"] != truth["root_model"]
        or model["convergence_searched"] != (mode != "shift")
    ):
        raise ValueError("Fitted model settings disagree with manifest")
    if (
        model["max_shifts"] != 2
        or model["search_strategy"] != "exhaustive"
        or model["trait"] != "value"
        or model["standard_error_column"] != "se"
    ):
        raise ValueError(
            "Fitted search or observation settings disagree with simulation protocol"
        )
    tree = read_tree(item["tree"], "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    parents = {ids[node]: -1 if node.is_root else ids[node.up] for node in ids}
    branches = model["branches"]
    if (
        len(branches) != len(parents)
        or {r["branch_id"]: r["parent"] for r in branches} != parents
    ):
        raise ValueError("Fitted branch topology disagrees with generating tree")
    tips = model["tip_predictions"]
    if len(tips) != truth["tips"] or {r["leaf_name"] for r in tips} != set(
        truth["tip_names"]
    ):
        raise ValueError("Fitted tip coverage disagrees with generating tree")
    tip_ids = {node.name: ids[node] for node in tree.leaves()}
    observed = dict(zip(truth["tip_names"], truth["observations"], strict=True))
    for row in tips:
        if (
            row["branch_id"] != tip_ids[row["leaf_name"]]
            or not math.isclose(
                row["observed"],
                observed[row["leaf_name"]],
                rel_tol=1e-14,
                abs_tol=1e-15,
            )
            or not equivalent(row["standard_error"], truth["standard_error"])
        ):
            raise ValueError(
                "Fitted observations, SEs or tip coordinates disagree with inputs"
            )
    selected = model["shift_branch_ids"]
    if len(selected) != len(set(selected)) or not set(selected) <= set(parents) - {0}:
        raise ValueError("Invalid fitted shift configuration")
    boot = model.get("bootstrap")
    if (boot is not None) != (mode == "bootstrap"):
        raise ValueError("Bootstrap availability disagrees with run mode")
    if boot is not None and boot["attempted"] != options["bootstrap"]:
        raise ValueError("Bootstrap count disagrees with manifest")


def audit_run(run, manifest, rows):
    options = manifest["options"]
    inputs = {
        case: audit_inputs(run / case, parameters)
        for case, parameters in expected_cases(options)
    }
    groups = defaultdict(list)
    for row in rows:
        item = inputs[row["case"]]
        if any(not equivalent(row.get(key), item["truth"][key]) for key in TRUTH_KEYS):
            raise ValueError("Trial conditions or seed disagree with generating truth")
        directory = run / row["case"] / row["mode"]
        if json.loads((directory / "result.json").read_text()) != row:
            raise ValueError("Trial record disagrees with its saved result")
        if row["status"] == "completed":
            model = json.loads((directory / "model.json").read_text())
            audit_model(model, item, options, row["mode"])
            metrics = selection_metrics(model, item["truth"])
            if any(
                not equivalent(row.get(key), value) for key, value in metrics.items()
            ):
                raise ValueError("Trial metrics disagree with fitted model and truth")
        elif not isinstance(row.get("error"), str) or not row["error"]:
            raise ValueError("Failed trial is missing its error")
        groups[tuple(row[key] for key in CELL_KEYS)].append(row)
    cells = [
        {**dict(zip(CELL_KEYS, key, strict=True)), **summarize(group)}
        for key, group in groups.items()
    ]
    saved = json.loads((run / "summary.json").read_text())

    def key(row):
        return tuple(row[field] for field in CELL_KEYS)

    if not equivalent(sorted(saved, key=key), sorted(cells, key=key)):
        raise ValueError("Saved cell summary disagrees with audited trial records")
    return inputs, cells
