"""Compare full versus overflow detection integration with equivalent target mass."""

import json
import statistics
import time
import tracemalloc

from nwkit.mul_locus import signature_size
from nwkit.mul_locus_integral import OVERFLOW, detection_distribution


def gene_tree(indices):
    if len(indices) == 1:
        index = indices[0]
        return ("tip", ("A", "X", "B")[index % 3], index)
    half = len(indices) // 2
    return ("node", gene_tree(indices[:half]), gene_tree(indices[half:]))


def measure(function):
    function()
    timings = []
    for _ in range(3):
        started = time.perf_counter()
        for _ in range(20):
            function()
        timings.append((time.perf_counter() - started) / 20)
    tracemalloc.start()
    function()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return {
        "seconds_per_call": timings,
        "median": statistics.median(timings),
        "peak_python_bytes": peak,
    }


def benchmark(tips):
    gene = gene_tree(tuple(range(tips)))
    detection = {"A": 0.6, "X": 0.6, "B": 0.6}

    def full():
        return detection_distribution(gene, detection, max_states=100000)

    def overflow():
        return detection_distribution(
            gene, detection, max_states=100000, max_observed_tips=4
        )

    all_sizes = full()
    selected = {
        key: mass for key, mass in all_sizes.items() if signature_size(key) <= 4
    }
    selected[OVERFLOW] = sum(
        mass for key, mass in all_sizes.items() if signature_size(key) > 4
    )
    bounded = overflow()
    if selected.keys() != bounded.keys():
        raise ArithmeticError("Detection integration changed the observation universe.")
    error = max(abs(selected[key] - bounded[key]) for key in selected)
    if error > 3e-13:
        raise ArithmeticError(
            "Detection overflow did not preserve target probabilities."
        )
    return {
        "tips": tips,
        "max_abs_probability_difference": error,
        "warmups": 1,
        "repetitions": 3,
        "calls_per_repetition": 20,
        "full": measure(full),
        "overflow": measure(overflow),
    }


if __name__ == "__main__":
    for tips in (8, 12):
        print(json.dumps(benchmark(tips)))
