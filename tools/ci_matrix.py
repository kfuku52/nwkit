"""Select CI coverage from changed paths; scheduled/manual runs cover every platform."""

import ast
import json
import os
import subprocess
from pathlib import Path

DEPENDENCY_PATHS = {"pyproject.toml", "setup.py", "constraints-dev.txt", "MANIFEST.in"}
NUMERICAL_MODULES = {
    "clade_index",
    "clade_mapping",
    "contrast",
    "evolution",
    "gaussian",
    "measurement_error",
    "model_matrix",
    "model_specs",
    "multivariate_pgls",
    "optimization",
    "ordinary_regression",
    "phylogenetic_glmm",
    "reconcile",
    "reconciliation_properties",
    "regress",
    "regression_pipeline",
    "replicates",
    "root_algorithms",
    "root_evaluation",
    "rsc_diagnostics",
    "sparse_laplace",
    "numerical_invariance",
}


def documentation_path(path):
    return path.endswith((".md", ".rst")) or path in {"LICENSE", ".gitignore"}


def numerical_path(path):
    parent, stem = Path(path).parent.as_posix(), Path(path).stem
    if parent == "tests" and stem.startswith("test_"):
        stem = stem[5:]
    return parent in {"nwkit", "tests"} and stem in NUMERICAL_MODULES


def select_coverage(paths, *, full=False):
    changed = [path for path in paths if not documentation_path(path)]
    source_checks = bool(full or changed)
    all_versions = full or any(
        path in DEPENDENCY_PATHS
        or path.startswith(".github/")
        or path == "tools/ci_matrix.py"
        for path in changed
    )
    platform_checks = bool(full or any(not numerical_path(path) for path in changed))
    matrix = []
    if source_checks:
        for version in ["3.10", "3.11", "3.12", "3.13"] if all_versions else ["3.10"]:
            matrix.append(
                {
                    "os": "ubuntu-latest",
                    "python-version": version,
                    "extras": "test,image",
                }
            )
    if platform_checks:
        for system in ("macos-latest", "windows-latest"):
            matrix.append({"os": system, "python-version": "3.14", "extras": "test"})
    return {
        "source_checks": source_checks,
        "macos_clean": platform_checks,
        "matrix": {"include": matrix},
    }


def version_only_change(before, after):
    def normalized(source):
        tree = ast.parse(source)
        versions = []
        for node in tree.body:
            if (
                isinstance(node, ast.Assign)
                and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == "__version__"
            ):
                if not isinstance(node.value, ast.Constant) or not isinstance(
                    node.value.value, str
                ):
                    return None, None
                versions.append(node.value.value)
                node.value = ast.Constant(value="VERSION")
        return ast.dump(tree), versions[0] if len(versions) == 1 else None

    old, old_version = normalized(before)
    new, new_version = normalized(after)
    return (
        old_version is not None and new_version is not None and old == new,
        new_version,
    )


def changed_paths(event_name, event):
    if event_name in {"schedule", "workflow_dispatch"}:
        return [], True
    base = event.get("pull_request", {}).get("base", {}).get("sha") or event.get(
        "before"
    )
    if not base or set(base) == {"0"}:
        return [], True
    # Diff against the checked-out merge result for PRs, never trust paths
    # interpolated into a shell command. Renames report both sides of the move.
    result = subprocess.run(
        ["git", "diff", "--name-only", "--no-renames", "-z", base, "HEAD"],
        check=True,
        stdout=subprocess.PIPE,
    )
    paths = [os.fsdecode(path) for path in result.stdout.split(b"\0") if path]
    full = False
    if "nwkit/__init__.py" in paths:
        previous = subprocess.run(
            ["git", "show", f"{base}:nwkit/__init__.py"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        )
        current = Path("nwkit/__init__.py").read_text(encoding="utf-8")
        version_only, version = version_only_change(previous.stdout, current)
        if version_only:
            paths.remove("nwkit/__init__.py")
            # A major/minor release receives complete coverage before tagging.
            full = version.split(".")[-1] == "0"
    return paths, full


def main():
    event_name = os.environ["GITHUB_EVENT_NAME"]
    event = json.loads(
        Path(os.environ["GITHUB_EVENT_PATH"]).read_text(encoding="utf-8")
    )
    paths, full = changed_paths(event_name, event)
    coverage = select_coverage(paths, full=full)
    with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as handle:
        for name, value in coverage.items():
            handle.write(f"{name}={json.dumps(value, separators=(',', ':'))}\n")
    print(json.dumps({"changed_paths": paths, **coverage}, indent=2))


if __name__ == "__main__":
    main()
