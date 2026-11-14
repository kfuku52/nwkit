"""Check small/few-event controls against the archived cluster implementation."""

import argparse
import ast
import fcntl
import json
import os
import tempfile
from pathlib import Path

for key in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"

from benchmark import fixture_module, flatten  # noqa: E402


def main():
    import numpy as np

    from nwkit import regress

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    args = parser.parse_args()
    source = (args.baseline / "nwkit/regress.py").read_text()
    node = next(
        x
        for x in ast.parse(source).body
        if isinstance(x, ast.FunctionDef) and x.name == "_fit_model"
    )
    namespace = dict(regress.__dict__)
    exec(
        compile(ast.get_source_segment(source, node), "archived_cluster", "exec"),
        namespace,
    )
    baseline = namespace["_fit_model"]
    fixture = fixture_module("03-cluster-score-aggregation")
    results = []
    with Path(tempfile.gettempdir(), "nwkit-performance-adoption.lock").open(
        "a"
    ) as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for n, g, p in ((200, 50, 3), (6000, 12, 5), (6000, 1500, 5)):
            frame, terms, options, _ = fixture.fixture(n, g, p, "event")

            def old(frame=frame, terms=terms, options=options):
                return baseline(frame, terms, **options)

            def new(frame=frame, terms=terms, options=options):
                return regress._fit_model(frame, terms, **options)

            left, right = flatten(old()), flatten(new())
            exact = True
            for key, value in left.items():
                other = right[key]
                if value.dtype.kind in "fc":
                    exact &= np.array_equal(value, other, equal_nan=True)
                    np.testing.assert_allclose(
                        value, other, rtol=1e-12, atol=1e-12, equal_nan=True
                    )
                else:
                    np.testing.assert_array_equal(value, other)
            measured = fixture.measure({"baseline": old, "candidate": new})
            results.append(
                {
                    "observations": n,
                    "events": g,
                    "predictors": p,
                    "exact": bool(exact),
                    "measurements": measured,
                }
            )
    Path(__file__).with_name("cluster-controls.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(
        [
            (
                x["observations"],
                x["events"],
                x["measurements"]["candidate"]["speedup_over_baseline"],
            )
            for x in results
        ]
    )


if __name__ == "__main__":
    main()
