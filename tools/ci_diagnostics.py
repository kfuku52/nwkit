"""Print allowlisted runtime metadata, without host identity or environment dumps."""

import importlib.metadata
import json
import os
import platform
import re

import numpy as np


def runtime_metadata():
    packages = {}
    for name in ("numpy", "scipy", "ete4", "pytest", "coverage"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    thread_settings = {
        name: value
        for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        if re.fullmatch(r"[1-9][0-9]{0,3}", value := os.environ.get(name, ""))
    }
    blas = (
        getattr(np.__config__, "CONFIG", {})
        .get("Build Dependencies", {})
        .get("blas", {})
    )
    build = {
        key: str(blas[key])[:256]
        for key in ("name", "version", "openblas configuration")
        if key in blas
    }
    pools = []
    try:
        from threadpoolctl import threadpool_info
    except ImportError:
        pass
    else:
        pools = [
            {
                key: pool[key]
                for key in ("internal_api", "version", "architecture", "num_threads")
                if key in pool
            }
            for pool in threadpool_info()
        ]
    return {
        "python": platform.python_version(),
        "system": platform.system(),
        "machine": platform.machine(),
        "packages": packages,
        "thread_settings": thread_settings,
        "numpy_blas_build": build,
        "threadpools": pools,
    }


if __name__ == "__main__":
    print(json.dumps(runtime_metadata(), sort_keys=True), flush=True)
