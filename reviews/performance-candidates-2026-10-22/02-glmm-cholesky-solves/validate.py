"""Run affected existing regressions with the isolated candidate installed."""

import fcntl
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
targets = [
    "tests/test_ordinary_regression.py::test_conventional_phylogenetic_glmm_response_families",
    "tests/test_ordinary_regression.py::test_conventional_scalar_non_gaussian_response_families",
    "tests/test_ordinary_regression.py::test_sparse_poisson_glmm_matches_dense_fit",
    "tests/test_ordinary_regression.py::test_sparse_binomial_glmm_matches_dense_fit",
    "tests/test_ordinary_regression.py::test_sparse_multilevel_glmm_matches_dense_fit",
    "tests/test_ordinary_regression.py::test_sparse_glmm_with_group_effect_matches_dense_fit",
    "tests/test_ordinary_regression.py::test_sparse_glmm_with_grouped_predictor_uncertainty_matches_dense",
    "tests/test_regression_selection.py::test_scalar_mode_verifies_score_instead_of_solver_flag",
    "tests/test_regression_selection.py::test_small_penalty_matches_existing_glmm",
    "tests/test_regression_inference_contract.py::test_penalized_wald_does_not_report_frequentist_p_or_intervals",
    "tests/test_regression_inference_contract.py::test_penalized_poisson_null_bootstrap_refits_both_models",
]
env = os.environ.copy()
env["PYTHONPATH"] = str(HERE) + os.pathsep + str(ROOT)
for name in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    env[name] = "1"
env["NWKIT_CI_DIAGNOSTICS"] = "0"
command = [
    sys.executable,
    "tools/check.py",
    "test",
    "--",
    "-p",
    "candidate02_plugin",
    "-rs",
    *targets,
]
with Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock").open("a") as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    with (HERE / "validation.txt").open("w") as output:
        output.write("Prototype validation command: " + repr(command) + "\n")
        output.flush()
        completed = subprocess.run(
            command, cwd=ROOT, env=env, stdout=output, stderr=subprocess.STDOUT
        )
    print((HERE / "validation.txt").read_text())
    raise SystemExit(completed.returncode)
