"""Diagnostic call counts (instrumented; not benchmark timing evidence)."""

import argparse
import contextlib
import cProfile
import json
import sys
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.checkout.resolve()))
    from nwkit.cli import main as cli

    case = next(
        c
        for c in json.loads(args.manifest.read_text())["cases"]
        if c["name"] == args.case
    )
    profile = cProfile.Profile(timer=time.process_time)
    with open("stdout.txt", "w") as stream:
        with contextlib.redirect_stdout(stream):
            profile.runcall(cli, case["args"])
    profile.create_stats()
    rows = [
        dict(
            file=file,
            line=line,
            function=name,
            primitive_calls=values[0],
            calls=values[1],
            self_cpu=values[2],
            cumulative_cpu=values[3],
        )
        for (file, line, name), values in profile.stats.items()
    ]
    rows.sort(key=lambda row: row["cumulative_cpu"], reverse=True)
    args.output.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
