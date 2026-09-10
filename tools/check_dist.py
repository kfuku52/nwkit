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
        "nwkit/data_model/lg.txt",
        "nwkit/data_iqtree/worker.cpp",
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
            Path(member).name in {"iqtree3", "iqtree3.exe", "nwkit-iqtree-worker"}
            or Path(member).name.startswith("libiqtree.")
            for member in members
        ):
            raise RuntimeError(
                "NWKIT wheels must not bundle IQ-TREE binaries or libraries."
            )
        if any(
            member in forbidden_names or member.endswith("/THIRD_PARTY_NOTICES")
            for member in members
        ):
            raise RuntimeError(f"Wheel contains a forbidden member: {candidate}")

    with tarfile.open(sdist, "r:gz") as archive:
        sdist_members = set(archive.getnames())
    required_sdist = {"/" + member for member in required_wheel} | {
        "/BRANCH_GAUSSIAN.md",
        "/examples/branch_gaussian/mixed_process.py",
        "/DTT.md",
        "/examples/dtt/README.md",
        "/examples/dtt/dtt.png",
        "/examples/dtt/heatmap.png",
        "/examples/dtt/reference.R",
        "/examples/dtt/traits.tsv",
        "/examples/dtt/tree.nwk",
        "/ASR.md",
        "/SHIFT.md",
        "/SHIFT_VALIDATION.md",
        "/SHIFT_JOINT.md",
        "/SHIFT_PBIC.md",
        "/SHIFT_ALPHA.md",
        "/SHIFT_CALIBRATION.md",
        "/reviews/shift-calibration-review.md",
        "/tools/shift_calibration_audit.py",
        "/tools/stress_shift_calibration.py",
        "/tools/verify_shift_stress.py",
        "/examples/shift/calibration-envelope-stress/audit.json",
        "/examples/shift/calibration-envelope/audit.json",
        "/examples/shift/calibration-envelope/records.jsonl.gz",
        "/examples/shift/calibration-envelope-stress/summary.json",
        "/examples/shift/calibration-weak-null/audit.json",
        "/examples/shift/calibration-weak-null/records.jsonl.gz",
        "/examples/shift/calibration-validation/protocol.json",
        "/examples/shift/calibration-validation/summary.json",
        "/examples/shift/calibration-validation/records.jsonl.gz",
        "/tools/validate_shift_calibration.py",
        "/tools/verify_shift_calibration.py",
        "/examples/shift/calibration-validation/audit.json",
        "/examples/shift/alpha-validation/protocol.json",
        "/examples/shift/alpha-validation/audit.json",
        "/examples/shift/alpha-validation/summary.json",
        "/examples/shift/alpha-validation/paired.json",
        "/examples/shift/alpha-validation/records.jsonl.gz",
        "/examples/shift/alpha-validation/candidate-ledger.jsonl.gz",
        "/examples/shift/alpha-validation/report.html",
        "/tools/validate_shift_alpha.py",
        "/tools/shift_alpha_backend.R",
        "/tools/shift_alpha_audit.py",
        "/tools/summarize_shift_alpha.py",
        "/examples/shift/pbic-correction/before.csv",
        "/examples/shift/pbic-correction/after.csv",
        "/examples/shift/pbic-correction/backend-provenance.json",
        "/examples/shift/joint-validation-pbic-fixed/records.json",
        "/examples/shift/joint-validation-pbic-fixed/candidates.json",
        "/examples/shift/joint-validation-pbic-fixed/backend-provenance.json",
        "/tools/shift_joint_candidates.py",
        "/tools/shift_joint_backend.py",
        "/tools/validate_shift_joint.py",
        "/tools/validate_shift_joint_pilot.py",
        "/tools/summarize_shift_joint.py",
        "/examples/shift/joint-validation/records.json",
        "/examples/shift/joint-validation/candidates.json",
        "/examples/shift/joint-validation/inputs.json",
        "/examples/shift/joint-validation/models.json",
        "/examples/shift/joint-validation/source-snapshot.json",
        "/tools/shift_simulation_audit.py",
        "/tools/shift_simulation_cases.py",
        "/tools/validate_shift_simulations.py",
        "/tools/summarize_shift_validation.py",
        "/examples/shift/validation/records.json",
        "/examples/shift/validation/inputs.json",
        "/examples/shift/validation/cell-summary.json",
        "/examples/shift/validation/source-snapshot.json",
        "/examples/shift/validation-16tip/records.json",
        "/examples/shift/validation-16tip/inputs.json",
        "/examples/shift/validation-16tip/cell-summary.json",
        "/examples/shift/tree.nwk",
        "/examples/shift/traits.tsv",
        "/examples/shift/traits-with-se.tsv",

        "/RADTE.md",
        "/IQTREE_LIBRARY.md",
        "/RADTE_MATH.md",
        "/RADTE_VALIDATION.md",
        "/examples/radte/gene.nwk",
        "/examples/radte/species.nwk",
        "/examples/radte/species-map.tsv",
        "/examples/radte/README.md",
        "/ASR_PERFORMANCE.md",
        "/CHANGELOG.md",
        "/CLI_TSV_CONVENTIONS.md",
        "/DEVELOPMENT.md",
        "/PHYLOGENETIC_REGRESSION.md",
        "/REGRESSION_SELECTION.md",
        "/examples/regression_selection/data.tsv",
        "/examples/regression_selection/tree.nwk",
        "/examples/regression_selection/folds.tsv",
        "/examples/regression_selection/predictors.txt",
        "/RECONCILED_SPECIATION_CONTRAST_MATH.md",
        "/RELEASING.md",
        "/constraints-dev.txt",
        "/examples/asr_figure/tree.nwk",
        "/examples/asr_figure/traits.tsv",
        "/examples/asr_figure/regimes.tsv",
        "/examples/asr_figure/parameters.tsv",
        "/examples/asr_figure/README.md",
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
