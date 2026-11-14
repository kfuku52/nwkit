"""Compare complete ill-conditioned lambda fits with Cholesky profiles."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from nwkit import signal_stats


def cholesky_only(*args):
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    provider = signal_stats._zero_error_lambda_profile
    rows = []
    try:
        for size in (8, 32, 128):
            for gap in (1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10):
                for seed in (314, 412):
                    half = size // 2
                    covariance = np.zeros((size, size))
                    covariance[:half, :half] = covariance[half:, half:] = 1.0
                    covariance += gap * np.eye(size)
                    values = np.r_[np.zeros(half), np.ones(half)]
                    values += np.random.default_rng(seed).normal(size=size) * np.sqrt(
                        gap
                    )
                    errors = np.zeros(size)
                    signal_stats._zero_error_lambda_profile = provider
                    spectral = provider(covariance, values, errors) is not None
                    actual = signal_stats.lambda_fit(covariance, values, errors)
                    signal_stats._zero_error_lambda_profile = cholesky_only
                    reference = signal_stats.lambda_fit(covariance, values, errors)
                    differences = {
                        key: abs(actual[key] - value)
                        for key, value in reference.items()
                        if isinstance(value, float)
                    }
                    assert actual["status"] == reference["status"]
                    assert all(
                        value <= (1e-5 if key == "estimate" else 1e-7)
                        for key, value in differences.items()
                    ), (size, gap, seed, differences)
                    rows.append(
                        {
                            "tips": size,
                            "gap": gap,
                            "seed": seed,
                            "spectral": spectral,
                            "exact": actual == reference,
                            "differences": differences,
                        }
                    )
    finally:
        signal_stats._zero_error_lambda_profile = provider
    source = Path(signal_stats.__file__)
    args.output.write_text(
        json.dumps(
            {
                "source": str(source),
                "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "controls": rows,
                "result": "passed",
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    print(f"{len(rows)} complete lambda fits match Cholesky reference allowances")


if __name__ == "__main__":
    main()
