"""Read-only independent audit and full seeded replay of the bounded experiment."""

import argparse
import gzip
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shift_alpha_design import generate_case  # noqa: E402
from shift_confidence_calibration import restricted_envelope  # noqa: E402
from shift_continuous_reference import DenseOUReference  # noqa: E402
from try_shift_confidence_calibration import (  # noqa: E402
    SCENARIOS,
    execute,
    generating_case,
    summarize,
)

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


def check_condition(tree, search, condition, data_seed, B):
    _, truth = generate_case(generating_case(data_seed, condition["scenario"]))
    if condition["truth"] != truth:
        raise ValueError("Independent generating observations/truth disagree")
    values = truth["observations"]
    result = condition["result"]
    model = search.models[result["full_winner"]]
    if model != result["full_model"]:
        raise ValueError("Full winner candidate identity disagrees")
    reference = DenseOUReference(tree, model)
    full = reference.at(values, search.grid[result["full_alpha_index"]])[
        "log_likelihood"
    ]
    null = DenseOUReference(tree, search.models[0])
    points = np.array([null.at(values, a)["log_likelihood"] for a in search.grid])
    error = max(
        abs(full - result["full_log_likelihood"]),
        float(np.max(abs(points - result["null_grid_log_likelihood"]))),
    )
    if error > 1e-7 or not np.isclose(
        2 * (full - points.max()), result["shift_statistic"], atol=1e-7
    ):
        raise ValueError("Independent likelihood or shift statistic disagrees")
    if not np.allclose(
        2 * (points.max() - points), result["confidence_statistics"], atol=1e-7
    ):
        raise ValueError("Confidence-set likelihood statistic disagrees")
    for prefix in ("shift", "confidence"):
        counts = np.asarray(result[prefix + "_exceedances"])
        if (
            counts.shape != (len(search.grid),)
            or not np.issubdtype(counts.dtype, np.integer)
            or np.any((counts < 0) | (counts > B))
        ):
            raise ValueError("Invalid Monte Carlo exceedance counts")
        if not np.array_equal(
            (1 + counts) / (B + 1), result[prefix + "_probabilities"]
        ):
            raise ValueError("Monte Carlo count/probability disagreement")
    if result["full_envelope_p"] != max(result["shift_probabilities"]) or result[
        "full_envelope_reject"
    ] != (result["full_envelope_p"] <= 0.05):
        raise ValueError("Full-grid envelope decision disagrees")
    for candidate in result["candidates"]:
        expected = restricted_envelope(
            result["shift_probabilities"],
            result["confidence_probabilities"],
            candidate["beta"],
        )
        if expected != candidate:
            raise ValueError(
                "Confidence-set restriction, correction or decision disagrees"
            )
    return error


def replay_case(payload):
    index, specification, expected = payload
    replay = execute((index, specification))
    replay.pop("seconds")
    expected = {k: v for k, v in expected.items() if k != "seconds"}
    if replay != expected:
        raise ValueError(f"Seeded full-bank replay disagrees in block {index}")
    return 1


def audit(directory, workers=4):
    specification = json.loads((directory / "protocol.json").read_text())
    for name, sha in specification["source_sha256"].items():
        snapshot = directory / "source-snapshot" / Path(name).name
        if hashlib.sha256(snapshot.read_bytes()).hexdigest() != sha:
            raise ValueError(f"Frozen source altered: {name}")
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != sha:
            raise ValueError(
                f"Active source differs from the frozen experiment: {name}"
            )
    with gzip.open(directory / "records.jsonl.gz", "rt") as stream:
        rows = list(map(json.loads, stream))
    if [r["block"] for r in rows] != list(range(specification["blocks"])):
        raise ValueError("Missing, duplicated or reordered blocks")
    search, tree = None, None
    maximum = 0.0
    for row in rows:
        seeds = np.random.SeedSequence(
            [specification["master_seed"], row["block"]]
        ).spawn(2)
        expected_seeds = [int(s.generate_state(1)[0]) for s in seeds]
        if [row["data_seed"], row["bootstrap_seed"]] != expected_seeds:
            raise ValueError("Paired generating/Monte Carlo seeds disagree")
        if row["status"] == "failed":
            if not row.get("error"):
                raise ValueError("Failed block has no reason")
            continue
        expected_tree, _ = generate_case(generating_case(row["data_seed"], "null"))
        if row["tree"] != expected_tree or [
            c["scenario"] for c in row["conditions"]
        ] != list(SCENARIOS):
            raise ValueError("Tree or paired scenarios disagree")
        if search is None:
            tree = read_tree(row["tree"], "auto", True, quiet=True)
            search = CalibratedSearch(tree, convergence=False)
        for condition in row["conditions"]:
            maximum = max(
                maximum,
                check_condition(
                    tree,
                    search,
                    condition,
                    row["data_seed"],
                    specification["bootstrap_replicates"],
                ),
            )
    if summarize(rows, specification["phase"]) != json.loads(
        (directory / "summary.json").read_text()
    ):
        raise ValueError("Reconstructed screening/futility summary disagrees")
    with ProcessPoolExecutor(max_workers=workers) as executor:
        replays = sum(
            executor.map(
                replay_case, [(i, specification, row) for i, row in enumerate(rows)]
            )
        )
    return dict(
        status="passed",
        independent_blocks=len(rows),
        paired_observations=4 * len(rows),
        generating_truth="exactly regenerated",
        selected_likelihood_and_null_profile="independently checked",
        maximum_likelihood_error=maximum,
        all_block_banks_replayed=replays,
        confidence_sets_and_corrections="reconstructed",
        screening_summary="matched",
        record_sha256=hashlib.sha256(
            (directory / "records.jsonl.gz").read_bytes()
        ).hexdigest(),
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Positive workers required")
    print(json.dumps(audit(args.directory, args.workers), indent=2))
