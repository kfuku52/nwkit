"""Serialize the affected snapshot regressions with the scoped prototype."""

import os
import subprocess
import sys

import benchmark as candidate

with candidate.measurement_lock():
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        (str(candidate.HERE), str(candidate.SOURCE), str(candidate.HERE.parent))
    )
    targets = sys.argv[1:] or [
        "tests/test_wgd_count_model.py",
        "-k",
        "pulse or pruning or event or convergence or fully_missing or mixture or root_clade",
        "-rs",
    ]
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate04_plugin",
        *targets,
    ]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=candidate.SOURCE, env=environment, check=True)
