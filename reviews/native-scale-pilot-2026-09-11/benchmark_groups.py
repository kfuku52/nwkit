"""Check regime propagation on terminal and nested 1,000-tip layouts."""

import argparse
import hashlib
import json
import time
from pathlib import Path

from benchmark_native_shift import balanced_tree, pectinate_tree

from nwkit.shift_native_model import ShiftLayout, ShiftTree
from nwkit.util import read_tree

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", required=True)
args = parser.parse_args()
path = Path(args.output)
if path.exists():
    parser.error("Use a fresh output.")
records = []
for shape, newick in [
    ("balanced", balanced_tree(1000)),
    ("pectinate", pectinate_tree(1000)),
]:
    tree = ShiftTree.build(read_tree(newick, "auto", True, quiet=True))
    for kind in ("terminal", "internal"):
        indices = (
            tree.compiled.leaf_indices
            if kind == "terminal"
            else [
                i for i, children in enumerate(tree.compiled.children) if i and children
            ]
        )
        for count in (10, 100):
            layout = ShiftLayout.build(
                tree, [tree.branch_ids[i] for i in indices[:count]]
            )
            layout.node_groups(tree)
            times = []
            for _ in range(3):
                start = time.perf_counter()
                for _ in range(100):
                    labels = layout.node_groups(tree)
                times.append((time.perf_counter() - start) / 100)
            records.append(
                {
                    "shape": shape,
                    "kind": kind,
                    "shifts": count,
                    "seconds_per_assignment": times,
                    "labels_sha256": hashlib.sha256(labels.tobytes()).hexdigest(),
                }
            )
path.write_text(json.dumps(records, indent=2) + "\n")
