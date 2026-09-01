"""Validate reproducibility and contents of built NWKIT distributions."""

import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
from nwkit import __version__  # noqa: E402

WHEEL_NAME = f"nwkit-{__version__}-py3-none-any.whl"


def _wheel_members(path: Path) -> set[str]:
    with zipfile.ZipFile(path) as wheel:
        return set(wheel.namelist())


def main() -> int:
    direct_wheel = PROJECT_ROOT / "direct-dist" / WHEEL_NAME
    wheel = PROJECT_ROOT / "dist" / WHEEL_NAME
    direct_sdist = PROJECT_ROOT / "direct-dist" / f"nwkit-{__version__}.tar.gz"
    sdist = PROJECT_ROOT / "dist" / f"nwkit-{__version__}.tar.gz"
    for artifact in (direct_wheel, wheel, direct_sdist, sdist):
        if not artifact.is_file():
            raise FileNotFoundError(f"Expected distribution was not built: {artifact}")
    if direct_wheel.read_bytes() != wheel.read_bytes():
        raise RuntimeError("Direct and sdist-built wheels are not byte-for-byte equal.")
    if direct_sdist.read_bytes() != sdist.read_bytes():
        raise RuntimeError(
            "Independent source distributions are not byte-for-byte equal."
        )

    required_wheel = {
        "nwkit/__init__.py",
        "nwkit/contrast.py",
        "nwkit/evolution.py",
        "nwkit/gaussian.py",
        "nwkit/image_metadata.py",
        "nwkit/measurement_error.py",
        "nwkit/model_matrix.py",
        "nwkit/multivariate_pgls.py",
        "nwkit/ordinary_regression.py",
        "nwkit/regress.py",
        "nwkit/regression_pipeline.py",
        "nwkit/phylogenetic_glmm.py",
        "nwkit/reconcile.py",
        "nwkit/replicates.py",
        "nwkit/root.py",
        "nwkit/root_compare.py",
        "nwkit/root_evaluation.py",
        "nwkit/sparse_laplace.py",
        "nwkit/data_tree/apgiv.nwk",
    }
    required_wheel.update(
        path.relative_to(PROJECT_ROOT).as_posix()
        for path in (PROJECT_ROOT / "nwkit").glob("*.py")
    )
    forbidden_names = {
        "nwkit/_mad.py",
        "nwkit/ordinary_pgls.py",
        "nwkit/pgls.py",
        "nwkit/pgls_pipeline.py",
        "THIRD_PARTY_NOTICES",
    }
    for candidate in (direct_wheel, wheel):
        members = _wheel_members(candidate)
        missing = required_wheel - members
        if missing:
            raise RuntimeError(f"Wheel is missing required members: {sorted(missing)}")
        if any(
            member in forbidden_names or member.endswith("/THIRD_PARTY_NOTICES")
            for member in members
        ):
            raise RuntimeError(f"Wheel contains a forbidden member: {candidate}")

    with tarfile.open(sdist, "r:gz") as archive:
        sdist_members = set(archive.getnames())
    required_sdist = {"/" + member for member in required_wheel} | {
        "/ASR.md",
        "/CHANGELOG.md",
        "/CLI_TSV_CONVENTIONS.md",
        "/DEVELOPMENT.md",
        "/PHYLOGENETIC_REGRESSION.md",
        "/RECONCILED_SPECIATION_CONTRAST_MATH.md",
        "/RELEASING.md",
        "/constraints-dev.txt",
        "/setup.py",
        "/tests/test_properties.py",
        "/tests/test_measurement_error.py",
        "/tests/test_distribution_reproducibility.py",
        "/tests/test_cli_contracts.py",
        "/tests/test_numerical_invariance.py",
        "/tests/test_output_transaction.py",
        "/tools/check_dist.py",
        "/tools/normalize_sdist.py",
        "/tools/complexity_baseline.json",
        "/tools/benchmark.py",
    }
    for suffix in sorted(required_sdist):
        if not any(member.endswith(suffix) for member in sdist_members):
            raise RuntimeError(f"Source distribution is missing {suffix}.")
    forbidden_sdist_suffixes = tuple(
        "/{}".format(name) for name in sorted(forbidden_names)
    )
    if any(member.endswith(forbidden_sdist_suffixes) for member in sdist_members):
        raise RuntimeError("Source distribution contains a forbidden member.")

    subprocess.run(
        [sys.executable, "-m", "twine", "check", str(wheel), str(sdist)],
        check=True,
    )
    subprocess.run(["check-wheel-contents", str(wheel)], check=True)
    print("Distribution reproducibility and contents are valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
