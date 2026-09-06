"""Reproduce the 0.43.3 multivariate audit workloads (warmup + 3 timed runs).

Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1.
Reported timings are workload-specific, not general speedup guarantees.
"""

import argparse
import json
import platform
import statistics
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from nwkit import __version__  # noqa: E402
from nwkit.asr_compare import (  # noqa: E402
    ComparisonCandidate,
    ComparisonContext,
    _fit_continuous,
)
from nwkit.multivariate_pgls import fit_multivariate_pgls  # noqa: E402
from nwkit.replicates import _expand_trait_estimates  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


def measure(function):
    function()
    elapsed = []
    result = None
    for _ in range(3):
        start = time.perf_counter()
        result = function()
        elapsed.append(time.perf_counter() - start)
    return statistics.median(elapsed), result


def pgls_workload(rng):
    responses = rng.normal(size=(80, 2))
    design = np.ones((80, 1))
    covariance = 0.2 * np.ones((80, 80)) + 0.8 * np.eye(80)
    seconds, fit = measure(
        lambda: fit_multivariate_pgls(responses, design, {"tree": covariance})
    )
    return {
        "seconds": seconds,
        "logL": fit.log_likelihood,
        "coefficients": fit.coefficients.tolist(),
        "covariance": fit.component_trait_covariances["tree"].tolist(),
    }


def replicate_workload():
    names = [f"t{i}" for i in range(8000)]
    observed_names = names[::2]
    values = np.ones(len(observed_names))
    seconds, expanded = measure(
        lambda: _expand_trait_estimates(
            names, observed_names, values, values, values, values
        )
    )
    np.testing.assert_array_equal(expanded[0][::2], values)
    assert np.isnan(expanded[0][1::2]).all()
    return {"seconds": seconds}


def comparison_workload(rng):
    tree = read_tree(
        "(" + ",".join(f"t{i}:{1 + (i % 7) * 0.1}" for i in range(400)) + ")R;",
        "1",
        True,
        quiet=True,
        rooted="yes",
    )
    values = rng.normal(size=(20, 2))
    data = pd.DataFrame(
        {
            "leaf_name": [f"t{i}" for i in range(20)],
            "x": values[:, 0],
            "y": values[:, 1],
        }
    )

    def compare():
        context = ComparisonContext(
            tree, data, "continuous", ("x", "y"), None, SimpleNamespace(alpha=0.7)
        )
        return _fit_continuous(
            context, ComparisonCandidate("MV-OU", "stationary", "model-default")
        )

    seconds, fit = measure(compare)
    return {
        "seconds": seconds,
        "logL": fit.log_likelihood,
        "sigma": fit.sigma.tolist(),
        "theta": fit.theta.tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional JSON result file.")
    args = parser.parse_args()
    rng = np.random.default_rng(241)
    result = {
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "nwkit": __version__,
        },
        "pgls": pgls_workload(rng),
        "replicates": replicate_workload(),
        "asrcompare": comparison_workload(rng),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    print(encoded, end="")
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()
