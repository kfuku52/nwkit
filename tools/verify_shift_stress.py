"""Replay every random-tree stress dataset and bootstrap fit from its saved seed."""

import argparse
import gzip
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from stress_shift_calibration import execute


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    root = args.directory
    protocol = json.loads((root / "protocol.json").read_text())
    with gzip.open(root / "records.jsonl.gz", "rt") as stream:
        rows = [json.loads(line) for line in stream]
    if [row["case"] for row in rows] != protocol["cases"]:
        raise ValueError("Stress case identities do not match the frozen design")
    changes = {}
    for source, expected in protocol["source_sha256"].items():
        snapshot = root / "source-snapshot" / Path(source).name
        if hashlib.sha256(snapshot.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Frozen source snapshot mismatch: {source}")
        current = hashlib.sha256(Path(source).read_bytes()).hexdigest()
        if current != expected:
            changes[source] = {"frozen": expected, "current": current}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for start in range(0, len(rows), 8):
            originals = rows[start : start + 8]
            replayed = executor.map(execute, [row["case"] for row in originals])
            for original, replay in zip(originals, replayed, strict=True):
                if replay != original:
                    raise ValueError(
                        f"Stress replay differs: case {original['case']['case_id']}"
                    )
            print(f"Stress replay: {min(start + 8, len(rows))}/{len(rows)}", flush=True)
    summary = json.loads((root / "summary.json").read_text())
    if (
        summary["attempted"] != len(rows)
        or summary["completed"] != sum(r["status"] == "completed" for r in rows)
        or summary["false_selection"]["count"]
        != sum(r.get("any_shift", False) for r in rows)
    ):
        raise ValueError("Stress summary differs from replayed records")
    audit = {
        "cases_replayed": len(rows),
        "inputs_covariances_and_fits": "exact_match",
        "source_changes": changes,
        "records_sha256": hashlib.sha256(
            (root / "records.jsonl.gz").read_bytes()
        ).hexdigest(),
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (root / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")


if __name__ == "__main__":
    main()
