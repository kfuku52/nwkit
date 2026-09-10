"""Check every measured output against its first same-revision and paired result."""

import argparse
import importlib.util
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location(
        "benchmark", Path(__file__).with_name("benchmark.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    results = json.loads(args.results.read_text())
    verified = {}
    specs = {case["name"]: case for case in results["manifest"]["cases"]}
    for name, case in results["cases"].items():
        rows = [
            row
            for row in case["records"]
            if row["repeat"] >= 0 and row["exit_code"] == 0
        ]
        first = {
            label: next((r for r in rows if r["label"] == label), None)
            for label in ["before", "after"]
        }
        details = []
        for row in rows:
            for baseline_label in [row["label"], "before"]:
                baseline = first[baseline_label]
                if baseline is None:
                    continue
                for table in specs[name]["tables"]:
                    check = module.compare_tables(
                        Path(baseline["directory"]) / table,
                        Path(row["directory"]) / table,
                    )
                    details.append(
                        dict(
                            label=row["label"],
                            repeat=row["repeat"],
                            baseline=baseline_label,
                            table=table,
                            check=check,
                        )
                    )
        verified[name] = details
    args.output.write_text(json.dumps(verified, indent=2) + "\n")


if __name__ == "__main__":
    main()
