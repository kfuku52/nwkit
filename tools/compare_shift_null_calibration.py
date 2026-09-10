"""Paired replay of the archived plug-in and current grid-envelope procedures."""

import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from nwkit.shift_candidates import tip_groups  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

BASELINE = Path(
    "examples/shift/calibration-review/source-snapshot/shift_calibration.py"
)
BASELINE_SHA = "e30e9aab8c4886811c315817e0b730e7d3966bfddb743899ec18c6cded6809cf"
ENGINES = {}


def replay(row):
    spec = importlib.util.spec_from_file_location(
        "archived_shift_calibration", BASELINE
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if row["tree"] not in ENGINES:
        ENGINES[row["tree"]] = module.CalibratedSearch(
            read_tree(row["tree"], "auto", True, quiet=True)
        )
    search = ENGINES[row["tree"]]
    result = {"case": row["case"], "envelope": row}
    try:
        current = row["fit"]
        fit = search.fit(
            row["truth"]["observations"],
            seed=current["seed"],
            replicates=current["calibration_replicates"],
            level=current["calibration_level"],
        )
        model = fit["model"]
        groups = sorted(
            list(g)
            for g in tip_groups(search.tree, model["shift_branch_ids"], model["groups"])
        )
        result["plugin"] = dict(
            status="completed",
            fit=fit,
            any_shift=bool(model["shift_branch_ids"]),
            partition_recovered=groups == row["truth"]["shared_partition"],
            mean_rmse=float(
                np.sqrt(
                    np.mean(
                        (
                            np.array(fit["predicted"])
                            - np.array(list(row["truth"]["tip_mean"].values()))
                        )
                        ** 2
                    )
                )
            ),
        )
    except (ValueError, np.linalg.LinAlgError) as exc:
        result["plugin"] = dict(status="failed", error=str(exc))
    return result


def summarize(rows):
    cells = defaultdict(list)
    for row in rows:
        c = row["case"]
        cells[c["family"], c["scenario"], c["root_model"]].append(row)
    result = []
    for key, members in cells.items():
        good = [r for r in members if r["plugin"]["status"] == "completed"]
        lost = sum(
            r["plugin"]["any_shift"] and not r["envelope"]["any_shift"] for r in good
        )
        gained = sum(
            r["envelope"]["any_shift"] and not r["plugin"]["any_shift"] for r in good
        )
        failures = len(members) - len(good)
        # An upper bound for loss-minus-gain is the exact upper bound for loss.
        # Include every failed pair as a loss; this is conservative and paired.
        worst_loss = lost + failures
        upper = (
            1.0
            if worst_loss == len(members)
            else float(beta.ppf(0.95, worst_loss + 1, len(members) - worst_loss))
        )
        result.append(
            dict(
                family=key[0],
                scenario=key[1],
                root_model=key[2],
                attempted=len(members),
                failed=failures,
                plugin_detected=sum(r["plugin"]["any_shift"] for r in good),
                envelope_detected=sum(r["envelope"]["any_shift"] for r in good),
                paired_loss=lost,
                paired_gain=gained,
                power_loss_upper_one_sided_95=upper,
                power_noninferiority_5pp_demonstrated=upper <= 0.05,
                plugin_partition=sum(r["plugin"]["partition_recovered"] for r in good),
                envelope_partition=sum(
                    r["envelope"]["partition_recovered"] for r in good
                ),
                plugin_mean_rmse=float(
                    np.mean([r["plugin"]["mean_rmse"] for r in good])
                )
                if good
                else None,
                envelope_mean_rmse=float(
                    np.mean([r["envelope"]["mean_rmse"] for r in good])
                )
                if good
                else None,
            )
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if hashlib.sha256(BASELINE.read_bytes()).hexdigest() != BASELINE_SHA:
        raise ValueError("Unexpected archived baseline")
    source = Path("examples/shift/calibration-envelope/records.jsonl.gz")
    with gzip.open(source, "rt") as stream:
        rows = [json.loads(line) for line in stream]
    rows = [r for r in rows if r["case"]["standard_error"] == 0]
    if any(r["status"] != "completed" for r in rows):
        raise ValueError("Incomplete envelope evidence")
    args.output.mkdir(parents=True, exist_ok=False)
    sources = [BASELINE, Path(__file__), Path("nwkit/shift_candidates.py")]
    protocol = dict(
        description="Paired software-effect replay of existing evidence, not fresh validation data",
        records_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_sha256={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
        },
        convergence=True,
        B=199,
        attempted=len(rows),
        known_error="Excluded: the first-stage method did not change",
        power_interval="Per-cell conservative one-sided exact 95% bound for paired losses; no multiplicity adjustment",
    )
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    snapshot = args.output / "source-snapshot"
    snapshot.mkdir()
    for source_file in sources:
        (snapshot / source_file.name).write_bytes(source_file.read_bytes())
    results = []
    with (
        ProcessPoolExecutor(max_workers=args.workers) as executor,
        gzip.open(args.output / "records.jsonl.gz", "wt") as stream,
    ):
        for index, row in enumerate(executor.map(replay, rows, chunksize=4)):
            results.append(row)
            stream.write(json.dumps(row, allow_nan=False) + "\n")
            stream.flush()
            if index % 25 == 0:
                print(f"Completed {index + 1}/{len(rows)} paired datasets", flush=True)
    (args.output / "summary.json").write_text(
        json.dumps(summarize(results), indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
