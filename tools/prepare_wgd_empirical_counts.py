"""Adapt a literal species-count matrix for an audited root-spanning study.

This prepares study inputs only; it does not fit models or label known events.
"""

import argparse
import hashlib
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit.gene_family_input import read_gene_family_table
from nwkit.util import read_tree
from nwkit.wgd_count import _count_tree


def prepare_counts(counts_path, tree_path, output, max_families, seed):
    if max_families < 0 or seed < 0:
        raise ValueError("Family limit and seed must be nonnegative integers.")
    tree = read_tree(str(tree_path), "auto", True)
    design, nodes, _ = _count_tree(tree)
    header, rows = read_gene_family_table(str(counts_path), design.tip_names)
    if set(header) != set(design.tip_names):
        raise ValueError("Raw count columns must exactly match species-tree tips.")
    try:
        numbers = [[Decimal(row[name]) for name in design.tip_names] for row in rows]
    except InvalidOperation as exc:
        raise ValueError("The complete-count study requires numeric counts.") from exc
    count_limit = int(np.iinfo(np.int64).max)
    if not numbers or any(
        not value.is_finite()
        or not 0 <= value <= count_limit
        or value != value.to_integral_value()
        for row in numbers
        for value in row
    ):
        raise ValueError(
            "This complete-count study needs finite nonnegative int64 counts."
        )
    counts = np.array(
        [[int(value) for value in row] for row in numbers], dtype=np.int64
    )
    groups = [
        set(leaf.name for leaf in nodes[child].leaves()) for child in design.children[0]
    ]
    eligible = np.ones(len(counts), dtype=bool)
    for names in groups:
        columns = [i for i, name in enumerate(design.tip_names) if name in names]
        eligible &= np.any(counts[:, columns] > 0, axis=1)
    indices = np.flatnonzero(eligible)
    if not len(indices):
        raise ValueError("No family spans all immediate species-root child clades.")
    if max_families and len(indices) > max_families:
        indices = np.sort(
            np.random.default_rng(seed).choice(indices, max_families, replace=False)
        )
    selected = np.zeros(len(counts), dtype=bool)
    selected[indices] = True
    family_ids = [f"source_row_{i + 1:05d}" for i in range(len(counts))]
    table = pd.DataFrame(counts[indices].astype(np.int64), columns=design.tip_names)
    table.insert(0, "family_id", [family_ids[i] for i in indices])
    audit = pd.DataFrame(
        {
            "family_id": family_ids,
            "source_data_row": np.arange(1, len(counts) + 1),
            "root_spanning": eligible.astype(int),
            "selected": selected.astype(int),
        }
    )
    provenance = {
        "schema_version": 1,
        "source_counts": str(Path(counts_path).resolve()),
        "source_tree": str(Path(tree_path).resolve()),
        "source_counts_sha256": hashlib.sha256(
            Path(counts_path).read_bytes()
        ).hexdigest(),
        "source_tree_sha256": hashlib.sha256(Path(tree_path).read_bytes()).hexdigest(),
        "source_family_ids": "one-based source data row; original family labels are unavailable",
        "raw_families": len(counts),
        "root_spanning_families": int(eligible.sum()),
        "selected_families": len(indices),
        "max_families": max_families,
        "seed": seed,
        "sampling": "uniform without replacement after root-clade selection; no copy-number cap/filter",
        "selected_max_copy_number": int(counts[indices].max()),
        "root_child_clades": [sorted(names) for names in groups],
        "known_event_locations_used": False,
    }
    output = Path(output)
    if output.exists():
        raise ValueError(
            "Study output directory must be new; existing evidence is preserved."
        )
    output.mkdir(parents=True)
    table.to_csv(output / "counts.tsv", sep="\t", index=False)
    audit.to_csv(output / "family_selection.tsv", sep="\t", index=False)
    (output / "provenance.json").write_text(
        json.dumps(provenance, indent=2, allow_nan=False) + "\n"
    )
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--species-tree", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--max-families",
        type=int,
        default=500,
        help="Uniform family limit; 0 includes all eligible families.",
    )
    parser.add_argument("--seed", type=int, default=620)
    args = parser.parse_args()
    print(
        json.dumps(
            prepare_counts(
                args.counts,
                args.species_tree,
                args.output,
                args.max_families,
                args.seed,
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
