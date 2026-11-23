"""Check tip retention and unrooted partial constraints in linear working space."""

import json
import sys
from typing import Any

from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction
from nwkit.util import read_item_per_line_file, read_tree, validate_unique_named_leaves


def _from_anchor(anchor):
    """Orient edges away from the same leaf, independent of input Newick roots."""
    ordered = []
    stack = [(anchor, None)]
    while stack:
        node, parent = stack.pop()
        ordered.append((node, parent))
        if node.up is not None and node.up is not parent:
            stack.append((node.up, node))
        stack.extend(
            (child, node) for child in reversed(node.children) if child is not parent
        )
    return ordered


def _intervals(ordered, indices):
    # In constraint DFS order every clade opposite the anchor is contiguous.
    # The leaf count detects holes without storing each clade's descendant set.
    summaries: dict[Any, tuple[int, int, int]] = {}
    intervals = set()
    n = len(indices)
    for node, parent in reversed(ordered):
        low, high, size = summaries.pop(node, (n, -1, 0))
        if node.is_leaf and node.name in indices:
            i = indices[node.name]
            low, high, size = min(i, low), max(i, high), size + 1
        if 1 < size <= n - 2 and high - low + 1 == size:
            intervals.add((low, high))
        if parent is not None and size:
            a, b, count = summaries.get(parent, (n, -1, 0))
            summaries[parent] = (min(a, low), max(b, high), count + size)
    return intervals


def check_constraint(inferred, constraint, expected):
    for tree, name in ((inferred, "--infile"), (constraint, "--constraint")):
        validate_unique_named_leaves(tree, option_name=name)
    if not expected or len(expected) != len(set(expected)):
        raise ValueError(
            "--species-list requires a nonempty list of unique input tip labels."
        )
    actual = {leaf.name: leaf for leaf in inferred.leaves()}
    fixed = {leaf.name: leaf for leaf in constraint.leaves()}
    wanted = set(expected)
    absent = sorted(fixed.keys() - actual.keys())
    result: dict[str, Any] = dict(
        missing_tips=sorted(wanted - actual.keys()),
        extra_tips=sorted(actual.keys() - wanted),
        free_tips=sorted(wanted - fixed.keys()),
        unexpected_constraint_tips=sorted(fixed.keys() - wanted),
        missing_constraint_tips=absent,
    )
    anchor = min(fixed)
    ordered = _from_anchor(fixed[anchor])
    order = [node.name for node, _ in ordered if node.is_leaf]
    indices = {name: i for i, name in enumerate(order)}
    required = _intervals(ordered, indices)
    incompatible = None
    if not absent:
        reproduced = _intervals(_from_anchor(actual[anchor]), indices)
        incompatible = sorted(required - reproduced)
    result.update(
        constraint_tip_order=order,
        required_splits=len(required),
        displayed_splits=len(required) - len(incompatible)
        if incompatible is not None
        else None,
        violated_intervals=incompatible,
        ok=not any(
            (
                absent,
                result["missing_tips"],
                result["extra_tips"],
                result["unexpected_constraint_tips"],
                incompatible,
            )
        ),
    )
    return result


def check_main(args):
    validate_outputs_do_not_replace_inputs(
        [
            ("--infile", args.infile),
            ("--constraint", args.constraint),
            ("--species-list", args.species_list),
        ],
        [("--outfile", args.outfile)],
    )
    result = check_constraint(
        read_tree(args.infile, args.format, args.quoted_node_names),
        read_tree(args.constraint, "auto", args.quoted_node_names),
        read_item_per_line_file(args.species_list),
    )
    if args.outfile == "-":
        json.dump(result, sys.stdout, indent=2)
        print()
    else:
        with output_transaction([args.outfile]) as staged:
            staged.write_text(
                args.outfile, lambda handle: json.dump(result, handle, indent=2)
            )
    if not result["ok"]:
        raise ValueError(
            "Inferred tree failed the tip or constraint split check. See JSON diagnostics."
        )
    return 0
