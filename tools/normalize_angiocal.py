"""Reproduce or verify the bundled records from the pinned AngioCal XLS."""

import argparse
import csv
import hashlib
import io
import sys
from dataclasses import asdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
from nwkit.angiocal import (  # noqa: E402
    ANGIOCAL_MAX_BYTES,
    ANGIOCAL_NORMALIZED_SOURCE_SHA256,
    FossilCalibration,
    _read_xls,
)


def normalized_records(path):
    with Path(path).open("rb") as handle:
        data = handle.read(ANGIOCAL_MAX_BYTES + 1)
    if hashlib.sha256(data).hexdigest() != ANGIOCAL_NORMALIZED_SOURCE_SHA256:
        raise ValueError("Expected the exact pinned AngioCal v1.0 XLS workbook.")
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=list(FossilCalibration.__dataclass_fields__),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    # Read the workbook independently of the bundled-record shortcut.
    writer.writerows(asdict(record) for record in _read_xls(data))
    return output.getvalue().encode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = normalized_records(args.workbook)
    target = PROJECT_ROOT / "nwkit" / "data_angiocal" / "v1.0.tsv"
    if args.check:
        if target.read_bytes() != expected:
            raise ValueError("Bundled AngioCal records differ from the official XLS.")
        print("All bundled AngioCal v1.0 records match the official XLS.")
    else:
        target.write_bytes(expected)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
