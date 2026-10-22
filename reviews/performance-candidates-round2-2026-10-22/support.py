"""Small shared helpers for exploratory, serialized candidate measurements."""

import fcntl
import gc
import platform
import statistics
import sys
import time
import tracemalloc
from contextlib import contextmanager


@contextmanager
def measurement_lock():
    with open("/tmp/nwkit-perf-round2-20261022-rkzoGe/measurement.lock", "a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def paired(baseline, candidate, *, repeats=7, minimum_seconds=0.03):
    baseline()
    candidate()
    start = time.perf_counter()
    baseline()
    elapsed = time.perf_counter() - start
    loops = max(1, min(10000, int(minimum_seconds / max(elapsed, 1e-9)) + 1))
    times = {"baseline": [], "candidate": []}
    for repeat in range(repeats):
        order = (
            ("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")
        )
        for name in order:
            function = baseline if name == "baseline" else candidate
            start = time.perf_counter()
            for _ in range(loops):
                function()
            times[name].append((time.perf_counter() - start) / loops)
    summaries = {
        name: {
            "raw_seconds": values,
            "median_seconds": statistics.median(values),
            "min_seconds": min(values),
            "max_seconds": max(values),
        }
        for name, values in times.items()
    }
    return {
        "loops_per_sample": loops,
        "repeats": repeats,
        **summaries,
        "speedup": summaries["baseline"]["median_seconds"]
        / summaries["candidate"]["median_seconds"],
    }


def peak_bytes(function):
    gc.collect()
    tracemalloc.start()
    result = function()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    del result
    return peak


def environment():
    import numpy
    import pandas
    import scipy

    import nwkit

    return {
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "pandas": pandas.__version__,
        "nwkit": nwkit.__version__,
        "nwkit_source": nwkit.__file__,
        "baseline_commit": "86c8f2c537769d852c17e2c789340f6c9de4ef54",
        "timing_limit": "Exploratory host-local measurements; unrelated CPU jobs are active.",
        "memory_limit": "tracemalloc peak is not complete RSS or all native allocations.",
    }
