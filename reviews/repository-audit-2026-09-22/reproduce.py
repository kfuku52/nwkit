"""Read-only source audit; all destructive scenarios use temporary files.

Run from the repository root: python reviews/repository-audit-2026-09-22/reproduce.py
This records observed behavior, rather than asserting that bugs must persist.
"""

import json
import os
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from ete4 import Tree

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from nwkit.cli import main
from nwkit.optimization import global_bounded_scalar_minimize
from nwkit.rename import read_name_tsv
from nwkit.util import write_tree


def observations():
    results = {}
    centers = np.array([1.0, 3.0, 5.0, 7.0, 9.25])
    offsets = np.array([0.0, 0.0, 0.0, 0.0, -10.0])

    def objective(x):
        return float(np.min(100 * (x - centers) ** 2 + offsets))

    fit = global_bounded_scalar_minimize(objective, (0, 16))
    results["optimizer"] = {
        "selected_x": fit.x,
        "selected_objective": fit.fun,
        "success": fit.success,
        "known_optimal_x": 9.25,
        "known_optimal_objective": objective(9.25),
    }
    with TemporaryDirectory(prefix="nwkit-review-") as directory:
        root = Path(directory)
        tree = root / "tree.nwk"
        tree.write_text("((A:1,B:1):1,(C:1,D:1):1);")
        mapping = root / "mapping.tsv"
        mapping.write_text("old_name\tnew_name\tnew_name\nA\tB\tC\n")
        try:
            results["duplicate_headers"] = read_name_tsv(mapping)
        except ValueError:
            results["duplicate_headers"] = "rejected"
        traits = root / "traits.tsv"
        original = "leaf_name\tx\nA\t1\nB\t2\nC\t3\nD\t4\n"
        traits.write_text(original)
        annotate_command = [
            "annotate",
            "--infile",
            str(tree),
            "--table",
            str(traits),
            "--columns",
            "x",
            "--report",
            str(traits),
            "--outfile",
            str(root / "annotated.nwk"),
        ]
        try:
            main(annotate_command)
        except ValueError:
            pass
        results["annotation_input_replaced"] = traits.read_text() != original

        target = root / "existing.nwk"
        target.write_text("ORIGINAL")

        class FailingWriter:
            def __init__(self, handle):
                self.handle = handle

            def __enter__(self):
                return self

            def write(self, data):
                self.handle.write(data[:5])
                self.handle.flush()
                raise OSError("simulated disk failure")

            def __exit__(self, *args):
                self.handle.close()

        original_fdopen = os.fdopen

        def fdopen(descriptor, mode, **kwargs):
            handle = original_fdopen(descriptor, mode, **kwargs)
            return FailingWriter(handle) if mode == "w" else handle

        with patch("nwkit.output_transaction.os.fdopen", side_effect=fdopen):
            try:
                write_tree(
                    Tree("(A:1,B:1);", parser=1),
                    SimpleNamespace(outfile=str(target)),
                    format=1,
                    quiet=True,
                )
            except OSError:
                pass
        results["existing_output_after_write_failure"] = target.read_text()

        source = root / "parameters.json"
        spec = {
            "trait_names": ["x"],
            "alpha": 1.0,
            "regime_optima": [[0.0]],
            "process_tip_covariance": [[1.0]],
        }
        truth = root / "truth.json"
        command = [
            "shift-simulate",
            "--infile",
            str(tree),
            "--input-rooted",
            "yes",
            "--parameters",
            str(source),
            "--truth-out",
            str(truth),
            "--outfile",
            str(root / "simulation.tsv"),
        ]
        audit = root / "audit.jsonl"
        source.write_text(json.dumps(spec))
        main(command + ["--audit", str(audit)])
        record = json.loads(audit.read_text())
        results["simulation_audit"] = {
            "inputs": [item["argument"] for item in record["inputs"]],
            "outputs": [item["argument"] for item in record["outputs"]],
        }
        for label, path in [("parameters", source), ("truth", truth)]:
            source.write_text(json.dumps(spec))
            try:
                main(command + ["--audit", str(path)])
            except ValueError:
                pass
            try:
                json.loads(path.read_text())
                valid = True
            except json.JSONDecodeError:
                valid = False
            results[f"audit_collision_{label}_remains_valid_json"] = valid
    return results


if __name__ == "__main__":
    print(json.dumps(observations(), indent=2))
