"""Report complexity growth and enforce reviewed per-function upper limits."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

MAX_FUNCTION_COMPLEXITY = 40
PROJECT_ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = Path(__file__).with_name("complexity_baseline.json")
EXCEPTIONS_PATH = Path(__file__).with_name("complexity_exceptions.json")


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


def complexity_violations(current, exceptions=None):
    # The baseline is a comparison point, never implicit permission to exceed
    # the common limit. Only a documented exception changes the hard ceiling.
    exceptions = {} if exceptions is None else exceptions
    violations = []
    for name, complexity in sorted(current.items()):
        limit = exceptions.get(name, {}).get("limit", MAX_FUNCTION_COMPLEXITY)
        if complexity > limit:
            violations.append(f"{name}: complexity {complexity} exceeds limit {limit}")
    return violations


def complexity_increases(current, baseline):
    return [
        f"{name}: complexity increased from {baseline[name]} to {complexity}; review the added branching."
        for name, complexity in sorted(current.items())
        if name in baseline and complexity > baseline[name]
    ]


def validate_exceptions(exceptions, current):
    if not isinstance(exceptions, dict):
        raise ValueError("Complexity exceptions must be an object keyed by function.")
    for name, record in exceptions.items():
        if name not in current:
            raise ValueError(f"Stale complexity exception: {name}")
        if (
            not isinstance(record, dict)
            or set(record) != {"limit", "reason", "tests"}
            or type(record["limit"]) is not int
            or record["limit"] <= MAX_FUNCTION_COMPLEXITY
            or not isinstance(record["reason"], str)
            or not record["reason"].strip()
            or not isinstance(record["tests"], list)
            or not record["tests"]
        ):
            raise ValueError(
                f"Exception requires limit, reason, and test paths: {name}"
            )
        for test in record["tests"]:
            if (
                not isinstance(test, str)
                or not test.startswith("tests/")
                or ".." in Path(test).parts
                or not test.endswith(".py")
                or not (PROJECT_ROOT / test).is_file()
            ):
                raise ValueError(f"Exception has an invalid test path: {name}: {test}")


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
        help="Record reviewed current values for future growth warnings; hard limits and exceptions stay unchanged.",
    )
    args = parser.parse_args(argv)
    baseline = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
    current = collect_complexities()
    exceptions = json.loads(EXCEPTIONS_PATH.read_text(encoding="utf-8"))
    validate_exceptions(exceptions, current)
    for warning in complexity_increases(current, baseline):
        print("Warning: " + warning, file=sys.stderr)
    violations = complexity_violations(current, exceptions)
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
        f"Radon: {len(current)} functions, {average:.2f} average (informational), {max(current.values())} maximum; all hard limits respected."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
