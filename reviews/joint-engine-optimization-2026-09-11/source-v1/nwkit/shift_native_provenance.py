"""Content-based replay guards for completed native output artifacts."""

import hashlib
import json
import os
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np


def _digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def native_implementation():
    import scipy

    sources = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(Path(__file__).parent.glob("*.py"))
    }
    return {
        "sources": sources,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "ete4": version("ete4"),
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "thread_environment": {
            name: os.environ.get(name)
            for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        },
    }


def native_implementation_sha256():
    return _digest(native_implementation())


def native_provenance(data, args, layout=None, names=None):
    tree = data.tree
    inputs = {
        "branch_ids": tree.branch_ids,
        "parents": [int(i) for i in tree.compiled.parents],
        "names": [str(node.name or "") for node in tree.compiled.nodes],
        "times": tree.times.tolist(),
        "height": tree.height,
        "traits": data.trait_names,
        "values": [
            [None if np.isnan(x) else float(x) for x in row] for row in data.values
        ],
        "variances": data.variances.tolist(),
        "centers": data.centers.tolist(),
        "scales": data.scales.tolist(),
        "layout": None if layout is None else [layout.shifts, layout.groups, names],
    }
    configuration = {
        name: getattr(args, name)
        for name in (
            "max_shifts",
            "convergence",
            "search_strategy",
            "exhaustive_max_configurations",
            "candidate_pool",
            "refit_budget",
            "screening_budget",
            "beam_width",
            "lasso_iterations",
            "search_memory_mb",
            "alpha",
            "process_tip_variance",
            "measurement_variance",
            "root_model",
            "estimate_measurement_error",
            "optimizer_starts",
            "max_iterations",
            "seed",
            "calibration_replicates",
            "calibration_level",
            "global_null_gate",
            "criterion",
            "bootstrap",
            "bootstrap_seed",
        )
    }
    configuration.update(
        trait_covariance=getattr(args, "trait_covariance", "diagonal"),
        alpha_model=getattr(args, "alpha_model", "trait-specific"),
        covariance_engine=getattr(args, "covariance_engine", "auto"),
    )
    implementation = native_implementation()
    return {
        "analysis_input_sha256": _digest(inputs),
        "configuration_sha256": _digest(configuration),
        "implementation_sha256": _digest(implementation),
        "configuration": configuration,
        "implementation": implementation,
    }


def read_native_resume(path, provenance):
    model = json.loads(Path(path).read_text(encoding="utf-8"))
    if model.get("schema_version") != 8 or model.get("selection") != "native":
        raise ValueError("Resume requires a completed schema-8 native model.")
    for field in (
        "analysis_input_sha256",
        "configuration_sha256",
        "implementation_sha256",
    ):
        if model.get(field) != provenance[field]:
            raise ValueError(f"Native resume rejected: {field} differs.")
    if model.get("completion_status") != "complete":
        raise ValueError("Native resume requires a completed analysis.")
    return model
