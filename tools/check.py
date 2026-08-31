"""Single entry point for local and CI validation."""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
# Distribution checks must work before installing runtime dependencies.
from nwkit import __version__  # noqa: E402

PYTHON = sys.executable


def run(*command: str, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=PROJECT_ROOT, check=True, env=env)


def lint_and_typecheck(*, incremental: bool = False) -> None:
    run(PYTHON, "-m", "ruff", "check", "nwkit", "tests", "setup.py", "tools")
    run(
        PYTHON,
        "-m",
        "ruff",
        "format",
        "--check",
        "nwkit",
        "tests",
        "setup.py",
        "tools",
    )
    run(PYTHON, "-m", "mypy", *([] if incremental else ["--no-incremental"]))


def run_tests(pytest_args: tuple[str, ...] = (), *, quick: bool = False) -> None:
    marker_selected = any(
        arg == "--markexpr" or arg.startswith("--markexpr=") or arg.startswith("-m")
        for arg in pytest_args
    )
    markers = ["-m", "not slow"] if quick and not marker_selected else []
    run(
        PYTHON,
        "-m",
        "pytest",
        "-q",
        "--durations=20",
        "--durations-min=0.05",
        *markers,
        *pytest_args,
    )


def run_full_checks() -> None:
    lint_and_typecheck()
    run(PYTHON, "-m", "pip", "check")
    run(PYTHON, "-m", "bandit", "-r", "nwkit", "-ll", "-ii", "-q")
    run(PYTHON, "-m", "pip_audit", ".")
    run(PYTHON, "-m", "coverage", "erase")
    run(PYTHON, "-m", "coverage", "run", "-m", "pytest", "tests/", "-q")
    run(PYTHON, "-m", "coverage", "report")
    run(PYTHON, "tools/check_maintainability.py")


def build_and_check_distributions() -> None:
    for directory_name in ("build", "dist", "direct-dist"):
        directory = PROJECT_ROOT / directory_name
        if directory.parent != PROJECT_ROOT:
            raise RuntimeError(f"Refusing to clean artifact path: {directory}")
        if directory.exists():
            shutil.rmtree(directory)
    stale_module = PROJECT_ROOT / "build" / "lib" / "nwkit" / "_mad.py"
    stale_module.parent.mkdir(parents=True, exist_ok=True)
    stale_module.write_text('raise RuntimeError("stale build artifact")\n')
    environment = os.environ.copy()
    if "SOURCE_DATE_EPOCH" not in environment:
        completed = subprocess.run(
            ["git", "log", "-1", "--pretty=%ct"],
            cwd=PROJECT_ROOT,
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        environment["SOURCE_DATE_EPOCH"] = completed.stdout.strip()
    run(
        PYTHON,
        "-m",
        "build",
        "--wheel",
        "--outdir",
        "direct-dist",
        env=environment,
    )
    run(
        PYTHON,
        "-m",
        "build",
        "--sdist",
        "--outdir",
        "direct-dist",
        env=environment,
    )
    direct_sdist = "direct-dist/nwkit-{}.tar.gz".format(__version__)
    run(PYTHON, "tools/normalize_sdist.py", direct_sdist, env=environment)
    run(PYTHON, "-m", "build", env=environment)
    sdist = "dist/nwkit-{}.tar.gz".format(__version__)
    run(PYTHON, "tools/normalize_sdist.py", sdist, env=environment)
    run(PYTHON, "tools/check_dist.py")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        choices=("test", "quick", "full", "dist", "release"),
        help="Validation depth to run.",
    )
    parser.add_argument(
        "pytest_args",
        nargs=argparse.REMAINDER,
        help="Optional pytest targets/options for test or quick (after --).",
    )
    args = parser.parse_args(argv)
    pytest_args = tuple(args.pytest_args)
    if pytest_args[:1] == ("--",):
        pytest_args = pytest_args[1:]
    if pytest_args and args.mode not in {"test", "quick"}:
        parser.error(
            "pytest selection is only supported by test and quick; full/release always run the entire suite"
        )
    if args.mode == "test":
        run_tests(pytest_args)
    elif args.mode == "quick":
        lint_and_typecheck(incremental=True)
        run_tests(pytest_args, quick=True)
    elif args.mode == "full":
        run_full_checks()
    elif args.mode == "dist":
        build_and_check_distributions()
    else:
        run_full_checks()
        build_and_check_distributions()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
