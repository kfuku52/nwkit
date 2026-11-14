"""Read-only independent likelihood and pairing audit of calibration comparisons."""

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_shift_null_calibration import summarize  # noqa: E402
from verify_shift_null_contract import check_winner  # noqa: E402

from nwkit.shift_candidates import enumerate_candidates  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


def read_rows(path):
    with gzip.open(path, "rt") as stream:
        return [json.loads(line) for line in stream]


def audit(directory):
    protocol = json.loads((directory / "protocol.json").read_text())
    source = Path("examples/shift/calibration-envelope/records.jsonl.gz")
    if hashlib.sha256(source.read_bytes()).hexdigest() != protocol["records_sha256"]:
        raise ValueError("Envelope evidence changed")
    for name, sha in protocol["source_sha256"].items():
        archived = directory / "source-snapshot" / Path(name).name
        source_file = archived if archived.exists() else Path(name)
        if hashlib.sha256(source_file.read_bytes()).hexdigest() != sha:
            raise ValueError(f"Paired replay source changed: {name}")
    expected = {
        row["case"]["case_id"]: row
        for row in read_rows(source)
        if row["case"]["standard_error"] == 0
    }
    rows = read_rows(directory / "records.jsonl.gz")
    if [row["case"]["case_id"] for row in rows] != list(expected):
        raise ValueError("Paired cases missing, duplicated or reordered")
    geometries = {}
    maximum = 0.0
    completed = 0
    for row in rows:
        original = expected[row["case"]["case_id"]]
        if row["envelope"] != original or row["case"] != original["case"]:
            raise ValueError("Paired envelope record differs from original")
        if original["tree"] not in geometries:
            tree = read_tree(original["tree"], "auto", True, quiet=True)
            models, _ = enumerate_candidates(tree, convergence=True)
            geometries[original["tree"]] = tree, models
        tree, models = geometries[original["tree"]]
        for lane in (row["plugin"], row["envelope"]):
            if lane["status"] != "completed":
                continue
            fit = lane["fit"]
            if (
                fit["seed"] != original["fit"]["seed"]
                or fit["calibration_replicates"] != 199
                or fit["calibration_level"] != 0.05
                or lane["any_shift"] != bool(fit["model"]["shift_branch_ids"])
            ):
                raise ValueError("Paired settings or shift flag disagree")
            maximum = max(
                maximum,
                check_winner(
                    tree,
                    np.array(original["truth"]["observations"]),
                    np.zeros(original["case"]["tips"]),
                    fit,
                    models,
                ),
            )
            completed += 1
    if summarize(rows) != json.loads((directory / "summary.json").read_text()):
        raise ValueError("Paired summary disagrees with raw records")
    return dict(
        status="passed",
        pairs=len(rows),
        independent_winner_fits=completed,
        maximum_likelihood_error=maximum,
        record_sha256=hashlib.sha256(
            (directory / "records.jsonl.gz").read_bytes()
        ).hexdigest(),
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.directory), indent=2))
