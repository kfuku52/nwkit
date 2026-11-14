from bisect import bisect_left

import pandas as pd

from nwkit.clade_index import LcaIndex
from nwkit.util import (
    get_species_group_records,
    read_tip_table,
    read_tree,
    validate_unique_named_leaves,
)

MONOPHYLY_COLUMNS = (
    "group",
    "status",
    "is_monophyletic",
    "num_target_taxa",
    "num_mrca_taxa",
    "num_intruder_taxa",
    "target_taxa",
    "intruder_taxa",
)


class _GroupClades:
    def __init__(self, tree, unrooted):
        self._index = LcaIndex(tree)
        self._leaves = tuple(node for node in self._index.nodes if node.is_leaf)
        self._ranks = {node.name: rank for rank, node in enumerate(self._leaves)}
        self._intervals = {}
        for node in tree.traverse("postorder"):
            if node.is_leaf:
                start = self._ranks[node.name]
                self._intervals[node] = (start, start + 1)
            else:
                self._intervals[node] = (
                    self._intervals[node.children[0]][0],
                    self._intervals[node.children[-1]][1],
                )
        self._edge_order: dict[tuple[int, int], int] = {}
        if unrooted:
            for rank, node in enumerate(self._index.nodes):
                self._edge_order.setdefault(self._intervals[node], rank)
        self._complement_sides = sorted(
            (len(self._leaves) - (end - start), rank, start, end)
            for (start, end), rank in self._edge_order.items()
        )

    def _ancestor_interval(self, positions):
        # The ancestor of the first and last DFS leaves contains all between them.
        node = self._index.common_ancestor(
            self._leaves[positions[0]], self._leaves[positions[-1]]
        )
        return self._intervals[node]

    def classify(self, names):
        positions = sorted(self._ranks[name] for name in names)
        start, end = self._ancestor_interval(positions)
        if end - start == len(positions):
            return True, "monophyletic", []
        ranges: tuple[range, ...] = (range(start, end),)
        if self._complement_sides:
            rooted_rank = (end - start, self._edge_order[start, end])
            # The smallest side wins, with the original preorder breaking ties.
            for size, rank, left, right in self._complement_sides:
                if (size, rank) >= rooted_rank:
                    break
                first = bisect_left(positions, left)
                if first == len(positions) or positions[first] >= right:
                    ranges = (range(left), range(right, len(self._leaves)))
                    break
        members = set(positions)
        intruders = [
            position for span in ranges for position in span if position not in members
        ]
        if not intruders:
            return True, "monophyletic", []
        left, right = self._ancestor_interval(intruders)
        status = "paraphyletic" if right - left == len(intruders) else "polyphyletic"
        return False, status, sorted(self._leaves[i].name for i in intruders)


def _read_trait_groups(
    path, tree_leaf_name_set, group_by, unmatched="warn", missing_values=None
):
    trait_df, _, _ = read_tip_table(
        path,
        option_name="--trait",
        tree_leaf_names=tree_leaf_name_set,
        required_columns=(group_by,),
        unmatched=unmatched,
        missing_values=missing_values,
    )
    trait_df = trait_df[trait_df["leaf_name"].isin(set(tree_leaf_name_set))].copy()
    trait_df = trait_df[~trait_df[group_by].isna()].copy()
    trait_df[group_by] = trait_df[group_by].astype(str)
    group_to_leaf_names = dict()
    for group_name, group_df in trait_df.groupby(group_by, sort=True):
        leaf_names = sorted(group_df["leaf_name"].astype(str).tolist())
        if len(leaf_names) > 0:
            group_to_leaf_names[str(group_name)] = leaf_names
    return group_to_leaf_names


def _get_groups(tree, args):
    if args.trait not in ["", None]:
        if args.group_by in ["", None]:
            raise ValueError("'--group-by' is required when '--trait' is specified.")
        return _read_trait_groups(
            path=args.trait,
            tree_leaf_name_set=set(tree.leaf_names()),
            group_by=args.group_by,
            unmatched=getattr(args, "unmatched", "warn"),
            missing_values=getattr(args, "missing_values", None),
        )
    _, species_to_leaf_names, _ = get_species_group_records(
        tree,
        option_name="--infile",
        context=" for 'monophyly'",
        args=args,
        require_monophyly=False,
    )
    return species_to_leaf_names


def monophyly_main(args):
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    validate_unique_named_leaves(
        tree, option_name="--infile", context=" for 'monophyly'"
    )
    group_to_leaf_names = _get_groups(tree, args)
    clades = _GroupClades(tree, args.unrooted)
    rows = list()
    non_monophyletic_groups = list()
    for group_name in sorted(group_to_leaf_names.keys()):
        leaf_names = sorted(str(name) for name in group_to_leaf_names[group_name])
        if len(leaf_names) == 0:
            continue
        is_monophyletic, status, intruder_leaf_names = clades.classify(leaf_names)
        row = {
            "group": group_name,
            "status": status,
            "is_monophyletic": bool(is_monophyletic),
            "num_target_taxa": len(leaf_names),
            "num_mrca_taxa": len(leaf_names) + len(intruder_leaf_names),
            "num_intruder_taxa": len(intruder_leaf_names),
            "target_taxa": ",".join(leaf_names),
            "intruder_taxa": ",".join(intruder_leaf_names),
        }
        rows.append(row)
        if not is_monophyletic:
            non_monophyletic_groups.append(group_name)
    out = pd.DataFrame(rows, columns=MONOPHYLY_COLUMNS)
    if args.outfile == "-":
        print(out.to_csv(sep="\t", index=False), end="")
    else:
        out.to_csv(args.outfile, sep="\t", index=False)
    if getattr(args, "fail_on_non_monophyly", False) and non_monophyletic_groups:
        raise ValueError(
            "Non-monophyletic group(s): {}".format(
                ", ".join(sorted(non_monophyletic_groups))
            )
        )
