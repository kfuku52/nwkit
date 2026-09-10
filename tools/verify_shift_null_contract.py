"""Read-only audit of null-contract data, independent fits and cellwise bounds."""

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shift_calibration_audit import check_probability_metadata  # noqa: E402
from shift_continuous_reference import DenseOUReference  # noqa: E402
from validate_shift_null_contract import cells, summarize  # noqa: E402

from nwkit.shift_calibration import ALPHA_HEIGHT_GRID, CalibratedSearch  # noqa: E402
from nwkit.shift_candidates import enumerate_candidates  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


def check_winner(tree, y, variances, fit, models):
    winner = fit["winner"]
    if not 0 <= winner < len(models) or fit["model"] != models[winner]:
        raise ValueError("Selected model disagrees with its candidate id")
    if fit["candidate_count"] != len(models):
        raise ValueError("Incorrect candidate count")
    alpha = np.inf if fit["alpha_height"] is None else fit["alpha_height"]
    if alpha not in ALPHA_HEIGHT_GRID:
        raise ValueError("Selected alpha is outside the declared grid")
    reference = DenseOUReference(tree, fit["model"], variances)
    point = reference.at(
        y, alpha, fit["process_tip_variance"] if np.any(variances) else None
    )
    error = abs(point["log_likelihood"] - fit["contrast_log_likelihood"])
    if error > 1e-6 or not np.isclose(
        point["process_tip_variance"],
        fit["process_tip_variance"],
        rtol=1e-8,
        atol=1e-10,
    ):
        raise ValueError("Independent winner likelihood or process variance disagrees")
    design, K = reference.branch_moments(alpha)
    mean = design @ point["coefficients"]
    covariance = point["process_tip_variance"] * K + np.diag(variances)
    weights = (
        (variances == 0).astype(float)
        if point["process_tip_variance"] == 0 and np.any(variances == 0)
        else np.linalg.solve(covariance, np.ones(len(y)))
    )
    intercept = float(weights @ (y - mean) / weights.sum())
    if not np.allclose(mean + intercept, fit["predicted"], rtol=1e-8, atol=1e-8):
        raise ValueError("Independent mean predictions disagree")
    return error


def check_stages(fit, known_error, B, level):
    tests = fit["tests"]
    if not tests:
        raise ValueError("Missing no-shift test")
    search = SimpleNamespace(grid=ALPHA_HEIGHT_GRID, known_error=known_error)
    for index, test in enumerate(tests):
        p = test["p_value"]
        if (
            test["stage"] != index
            or not 0 < p <= 1
            or abs(p * (B + 1) - round(p * (B + 1))) > 1e-10
        ):
            raise ValueError("Invalid stage index or Monte Carlo probability")
        check_probability_metadata(test, index, search, B, level)
        if index < len(tests) - 1 and p > level:
            raise ValueError("Search continued after nonrejection")
    shifted = bool(fit["model"]["shift_branch_ids"])
    if shifted != (tests[0]["p_value"] <= level):
        raise ValueError("Any-shift decision disagrees with the null test")
    return shifted


def audit(directory, replay_stride=0, *, frozen_engine=False):
    specification = json.loads((directory / "protocol.json").read_text())
    if (
        frozen_engine
        and "nwkit/shift_calibration.py" not in specification["source_sha256"]
    ):
        raise ValueError("Archived engine lacks a declared source hash")
    for name, sha in specification["source_sha256"].items():
        archived = directory / "source-snapshot" / Path(name).name
        if hashlib.sha256(archived.read_bytes()).hexdigest() != sha:
            raise ValueError(f"Snapshot hash mismatch: {name}")
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != sha:
            if not frozen_engine or name != "nwkit/shift_calibration.py":
                raise ValueError(f"Active generator/fitting hash mismatch: {name}")
    engine_class = CalibratedSearch
    if frozen_engine:
        frozen = directory / "source-snapshot" / "shift_calibration.py"
        module = ModuleType("nwkit._archived_shift_calibration")
        module.__file__ = str(frozen)
        # Compile the verified source directly: import loaders would create a
        # __pycache__ inside the input evidence, violating read-only auditing.
        exec(compile(frozen.read_bytes(), str(frozen), "exec"), module.__dict__)
        if not np.array_equal(module.ALPHA_HEIGHT_GRID, ALPHA_HEIGHT_GRID):
            raise ValueError("Archived alpha grid differs from the audit grid")
        engine_class = module.CalibratedSearch
    designs = {c["cell_id"]: c for c in cells(specification["suite"])}
    if any(c != designs[c["cell_id"]] for c in specification["cells"]):
        raise ValueError("Protocol cells disagree with the frozen design")
    expected = [
        (c["cell_id"], r)
        for c in specification["cells"]
        for r in range(specification["replicates"])
    ]
    with gzip.open(directory / "records.jsonl.gz", "rt") as stream:
        rows = [json.loads(line) for line in stream]
    if [(r["cell_id"], r["replicate"]) for r in rows] != expected:
        raise ValueError("Missing, duplicated or reordered datasets")
    geometries, models, replay_engines = {}, {}, {}
    maximum_error = 0.0
    replays = 0
    for index, row in enumerate(rows):
        cell = designs[row["cell_id"]]
        key = cell["tree"], cell["known_error"], cell["alpha_height"]
        if key not in geometries:
            tree = read_tree(cell["tree"], "auto", True, quiet=True)
            variances = (
                np.geomspace(0.01, 0.25, cell["tips"])
                if cell["known_error"]
                else np.zeros(cell["tips"])
            )
            null = {"groups": [[0]], "shift_branch_ids": []}
            reference = DenseOUReference(tree, null, variances)
            a = np.inf if cell["alpha_height"] is None else cell["alpha_height"]
            _, K = reference.branch_moments(a)
            covariance = cell["process_tip_variance"] * K + np.diag(variances)
            geometries[key] = tree, variances, covariance
        tree, variances, covariance = geometries[key]
        seed = np.random.SeedSequence(
            [specification["seed"], cell["cell_id"], row["replicate"]]
        )
        data_seed, bootstrap_seed = seed.spawn(2)
        boot = int(bootstrap_seed.generate_state(1)[0] % (2**31))
        y = np.linalg.cholesky(covariance) @ np.random.default_rng(data_seed).normal(
            size=cell["tips"]
        )
        if (
            row["bootstrap_seed"] != boot
            or not np.array_equal(y, row["observations"])
            or not np.array_equal(covariance, row["generating_covariance"])
            or not np.array_equal(variances, row["known_variances"])
        ):
            raise ValueError("Generated observations/covariance/seeds disagree")
        if [lane["convergence"] for lane in row["lanes"]] != specification[
            "convergence_lanes"
        ]:
            raise ValueError("Missing or reordered search lanes")
        for lane in row["lanes"]:
            if lane["status"] == "failed":
                if not lane.get("error"):
                    raise ValueError("Failure lacks an explanation")
                continue
            if lane["status"] != "completed":
                raise ValueError("Unknown fit status")
            model_key = cell["tree"], lane["convergence"]
            if model_key not in models:
                models[model_key], _ = enumerate_candidates(
                    tree, convergence=lane["convergence"]
                )
            fit = lane["fit"]
            B = specification["bootstrap_replicates"]
            if (
                fit["seed"] != boot
                or fit["calibration_replicates"] != B
                or fit["calibration_level"] != 0.05
            ):
                raise ValueError("Fit settings disagree with protocol")
            if lane["any_shift"] != check_stages(fit, bool(np.any(variances)), B, 0.05):
                raise ValueError("Saved any-shift flag disagrees with the model")
            maximum_error = max(
                maximum_error, check_winner(tree, y, variances, fit, models[model_key])
            )
            if replay_stride and index % replay_stride == 0:
                engine_key = model_key, cell["known_error"]
                if engine_key not in replay_engines:
                    replay_engines[engine_key] = engine_class(
                        tree, convergence=lane["convergence"], variances=variances
                    )
                replay = replay_engines[engine_key].fit(y, seed=boot, replicates=B)
                if replay["model"] != fit["model"] or replay["tests"] != fit["tests"]:
                    raise ValueError("Seeded complete-search replay disagrees")
                replays += 1
    summary = summarize(rows, specification)
    if summary != json.loads((directory / "summary.json").read_text()):
        raise ValueError("Summary differs from reconstructed counts")
    for result in summary["cell_results"]:
        upper = result["simultaneous_one_sided_95_upper"]
        worst = result["any_shift"] + result["failed"]
        if worst < result["attempted"] and not np.isclose(
            binom.cdf(worst, result["attempted"], upper),
            0.05 / summary["comparisons"],
            rtol=1e-8,
            atol=1e-12,
        ):
            raise ValueError("Cellwise binomial upper bound is inconsistent")
    return {
        "status": "passed",
        "datasets": len(rows),
        "fits": summary["fits"],
        "source_snapshots": "matched",
        "engine_scope": "archived engine; not current CLI validation"
        if frozen_engine
        else "active engine",
        "inputs": "exactly_regenerated",
        "independent_winner_likelihood_max_error": maximum_error,
        "independent_predictions": "matched",
        "complete_search_replays": replays,
        "cellwise_counts_and_bounds": "matched",
        "all_cells_meet_acceptance": summary["all_cells_pass"],
        "record_sha256": hashlib.sha256(
            (directory / "records.jsonl.gz").read_bytes()
        ).hexdigest(),
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "--replay-stride",
        type=int,
        default=0,
        help="Also replay every Nth dataset; zero disables bootstrap replay",
    )
    parser.add_argument(
        "--frozen-engine",
        action="store_true",
        help="Audit the hash-verified archived fitting engine; all other sources must match active files",
    )
    args = parser.parse_args()
    if args.replay_stride < 0:
        parser.error("--replay-stride must be nonnegative")
    print(
        json.dumps(
            audit(args.directory, args.replay_stride, frozen_engine=args.frozen_engine),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
