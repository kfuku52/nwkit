"""Select among eligible root edges with explicit method priority."""

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any


def _topology_splits(tree, leaf_bits, leaf_mask):
    """Encode each unrooted edge using a shared, unique tip-to-bit mapping."""
    descendants: dict[Any, int] = {}
    splits = set()
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            mask = leaf_bits[node.name]
        else:
            mask = 0
            for child in node.children:
                mask |= descendants[child]
        descendants[node] = mask
        if not node.is_root:
            # Either side identifies the same unordered edge bipartition.
            splits.add(min(mask, leaf_mask ^ mask))
    return splits


def select_eligible_root(
    args, before_publish=None, *, priority=("mad", "midpoint", "candidate")
):
    if len(priority) != 3 or set(priority) != {"mad", "midpoint", "candidate"}:
        raise ValueError(
            "Root priority must order mad, midpoint and candidate exactly once"
        )
    from nwkit.clade_mapping import projected_root_split
    from nwkit.cli import parser as nwkit_parser
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction, validate_output_targets
    from nwkit.root import root_main
    from nwkit.root_compare import root_compare_main
    from nwkit.util import (
        get_tree_property_names,
        read_tree,
        read_tree_strings,
        write_tree,
    )

    candidate_strings = read_tree_strings(str(args.candidate_trees))
    if not candidate_strings:
        raise ValueError("Reconciliation candidate tree collection is empty.")
    outputs = [args.out_tree, args.comparison_table, args.comparison_plot]
    validate_outputs_do_not_replace_inputs(
        [
            ("input tree", args.in_tree),
            ("reconciliation candidates", args.candidate_trees),
        ],
        [("output", path) for path in outputs],
    )
    validate_output_targets(outputs)
    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="nwkit-root-", dir=args.out_tree.parent
    ) as temporary:
        work = Path(temporary)
        candidates = []
        for index, text in enumerate(candidate_strings):
            path = work / f"candidate.{index}.nwk"
            path.write_text(text + "\n", encoding="utf-8")
            candidates.append(path)
        native_trees: dict[str, Any] = {}
        comparison_args = nwkit_parser.parse_args(
            [
                "rootcompare",
                "--infile",
                str(args.in_tree),
                "--methods",
                "midpoint,mad",
                "--species-parser",
                args.species_parser,
                "--outfile",
                str(work / "comparison.tsv"),
                "--figure-out",
                str(work / "comparison.pdf"),
            ]
        )
        root_compare_main(comparison_args, rooted_trees=native_trees)
        original_tree = read_tree(str(args.in_tree), "auto", True, quiet=True)
        preserve_unary = any(
            len(node.children) == 1 and not node.is_root
            for node in original_tree.traverse()
        )
        roots = {}
        for method, tree in native_trees.items():
            roots[method] = work / f"{method}.nwk"
            output_args = nwkit_parser.parse_args(
                [
                    "root",
                    "--method",
                    method,
                    "--infile",
                    str(args.in_tree),
                    "--outfile",
                    str(roots[method]),
                ]
            )
            if preserve_unary:
                # rootcompare contracts internal singleton nodes for comparison.
                # Native selection must retain the original CLI tree shape.
                root_main(output_args)
                native_trees[method] = read_tree(
                    str(roots[method]), "auto", True, quiet=True
                )
            else:
                write_tree(
                    tree,
                    output_args,
                    format=output_args.outformat,
                    props=get_tree_property_names(tree),
                )
        leaves = frozenset(native_trees["mad"].leaf_names())
        if any(
            len(tree) != len(leaves) or frozenset(tree.leaf_names()) != leaves
            for tree in native_trees.values()
        ):
            raise ValueError("NWKIT root trees must have matching, unique tip labels.")
        native_splits = {
            method: projected_root_split(tree, leaves)
            for method, tree in native_trees.items()
        }
        if any(split is None for split in native_splits.values()):
            raise ValueError("NWKIT did not produce a resolved root edge.")
        leaf_bits = {name: 1 << index for index, name in enumerate(sorted(leaves))}
        leaf_mask = (1 << len(leaves)) - 1
        expected_topology = _topology_splits(native_trees["mad"], leaf_bits, leaf_mask)
        reconciliation_splits = []
        for path in candidates:
            candidate = read_tree(str(path), "auto", True)
            if (
                len(candidate) != len(leaves)
                or frozenset(candidate.leaf_names()) != leaves
            ):
                raise ValueError(
                    f"Reconciliation candidate tips are duplicated or differ from the input tree: {path}"
                )
            if _topology_splits(candidate, leaf_bits, leaf_mask) != expected_topology:
                raise ValueError(
                    f"Reconciliation candidate has a different unrooted topology: {path}"
                )
            split = projected_root_split(candidate, leaves)
            if split is None:
                raise ValueError(
                    f"Reconciliation candidate does not have a resolved root edge: {path}"
                )
            reconciliation_splits.append(split)
        compatible = {
            method: split in reconciliation_splits
            for method, split in native_splits.items()
        }
        method = next(
            method for method in priority if method == "candidate" or compatible[method]
        )
        if method == "candidate":
            selected, reason = candidates[0], "first_reconciliation_candidate"
        else:
            selected, reason = roots[method], method + "_compatible_with_reconciliation"
        print(f"Number of optimal reconciliation roots: {len(candidates)}", flush=True)
        for method in ("mad", "midpoint"):
            print(
                f"is_{method}_compatible_with_reconciliation: {compatible[method]}",
                flush=True,
            )
        print(f"Selected root: {reason}", flush=True)
        subprocess.run(
            [
                "nwkit",
                "drop",
                "--infile",
                str(selected),
                "--outfile",
                str(work / "selected.nwk"),
                "--target",
                "intnode",
                "--name",
                "yes",
                "--support",
                "yes",
            ],
            check=True,
        )
        with output_transaction(outputs, create_parents=True) as staged:
            for source, target in zip(
                ("selected.nwk", "comparison.tsv", "comparison.pdf"),
                outputs,
                strict=True,
            ):
                shutil.copyfile(work / source, staged[target])
            if before_publish is not None:
                before_publish()
