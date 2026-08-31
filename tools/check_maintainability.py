"""Enforce per-function cyclomatic-complexity ratchets without an average gate."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

MAX_NEW_FUNCTION_COMPLEXITY = 40
PROJECT_ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = Path(__file__).with_name("complexity_baseline.json")


def function_complexities(results):
    """Include methods and nested functions, each with a stable qualified name."""
    functions = {}

    def visit(path, block, prefix=""):
        name = prefix + block["name"]
        if block["type"] == "class":
            for method in block.get("methods", []):
                visit(path, method, name + ".")
            for child in block.get("inner_classes", []):
                visit(path, child, name + ".")
            return
        functions[f"{path}:{name}"] = int(block["complexity"])
        for closure in block.get("closures", []):
            visit(path, closure, name + ".<locals>.")

    for path, blocks in results.items():
        for block in blocks:
            visit(path.replace("\\", "/"), block)
    return functions


def complexity_violations(current, baseline):
    violations = []
    for name, complexity in sorted(current.items()):
        limit = baseline.get(name, MAX_NEW_FUNCTION_COMPLEXITY)
        if complexity > limit:
            label = (
                "existing function ceiling"
                if name in baseline
                else "new function limit"
            )
            violations.append(
                f"{name}: complexity {complexity} exceeds {label} {limit}"
            )
    return violations


def collect_complexities():
    completed = subprocess.run(
        [sys.executable, "-m", "radon", "cc", "nwkit", "--json"],
        cwd=PROJECT_ROOT,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    current = function_complexities(json.loads(completed.stdout))
    if not current:
        raise RuntimeError("Radon did not find any functions to analyze.")
    return current


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--update-baseline",
        action="store_true",
        help="Record reductions/new functions and remove deleted functions; never raise an existing ceiling.",
    )
    args = parser.parse_args(argv)
    baseline = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
    current = collect_complexities()
    violations = complexity_violations(current, baseline)
    if violations:
        raise RuntimeError(
            "Cyclomatic-complexity budget failed:\n- " + "\n- ".join(violations)
        )
    if args.update_baseline:
        BASELINE_PATH.write_text(
            json.dumps(current, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    average = sum(current.values()) / len(current)
    print(
        f"Radon: {len(current)} functions, {average:.2f} average (informational), {max(current.values())} maximum; all function ceilings respected."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
