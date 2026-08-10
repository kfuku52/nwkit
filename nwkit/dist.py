import csv
import math
import sys
from itertools import combinations
from typing import Any

from nwkit.util import (
    get_subtree_leaf_name_sets,
    is_rooted,
    read_tree,
    validate_unique_named_leaves,
)

SUPPORTED_METRICS = (
    "rf",
    "normalized-rf",
    "weighted-rf",
    "branch-score",
    "path-topological",
    "path-length",
)
DISTANCE_COLUMNS = (
    "metric",
    "comparison",
    "num_taxa",
    "distance",
    "max_distance",
)
SPLIT_METRICS = frozenset(
    {
        "rf",
        "normalized-rf",
        "weighted-rf",
        "branch-score",
    }
)
BRANCH_LENGTH_METRICS = frozenset(
    {
        "weighted-rf",
        "branch-score",
        "path-length",
    }
)
PATH_MATRIX_MAX_TAXA = 2000
PATH_MATRIX_MAX_DYNAMIC_RANGE = 10**12
PATH_MATRIX_ROUNDOFF_FACTOR = 16


def _parse_metric_values(values):
    metrics = list()
    has_all = False
    for value in values:
        for item in str(value).split(","):
            metric = item.strip().lower().replace("_", "-")
            if metric == "":
                raise ValueError("Empty metric in `--metric`.")
            if metric == "all":
                has_all = True
                continue
            if metric not in SUPPORTED_METRICS:
                raise ValueError(
                    "Unsupported metric '{}'. Choose from: all, {}.".format(
                        item.strip(),
                        ", ".join(SUPPORTED_METRICS),
                    )
                )
            if metric not in metrics:
                metrics.append(metric)
    if has_all:
        return SUPPORTED_METRICS
    if not metrics:
        return SUPPORTED_METRICS
    return tuple(metrics)


def _resolve_metrics(args):
    canonical_values = getattr(args, "metric", None)
    legacy_value = getattr(args, "dist", None)
    if canonical_values is not None and legacy_value is not None:
        raise ValueError("Use either '--metric' or deprecated '--dist', not both.")
    if legacy_value is not None:
        metrics = _parse_metric_values([legacy_value])
        legacy_output = metrics == ("rf",)
        return metrics, legacy_output
    if canonical_values is None:
        return SUPPORTED_METRICS, False
    if isinstance(canonical_values, str):
        canonical_values = [canonical_values]
    return _parse_metric_values(canonical_values), False


def _validate_branch_lengths(tree, option_name, metrics):
    required_by = sorted(BRANCH_LENGTH_METRICS.intersection(metrics))
    if not required_by:
        return
    for node in tree.traverse():
        if node.is_root:
            continue
        if node.dist is None:
            raise ValueError(
                "Metric(s) {} require a branch length on every non-root edge; "
                "a length is missing in '{}'.".format(
                    ", ".join(required_by), option_name
                )
            )
        try:
            branch_length = float(node.dist)
        except (TypeError, ValueError):
            branch_length = math.nan
        if not math.isfinite(branch_length) or branch_length < 0:
            raise ValueError(
                "Metric(s) {} require finite, nonnegative branch lengths; "
                "found {!r} in '{}'.".format(
                    ", ".join(required_by),
                    node.dist,
                    option_name,
                )
            )


def _canonical_unrooted_split(descendant_taxa, all_taxa):
    other_taxa = all_taxa.difference(descendant_taxa)
    if len(descendant_taxa) < len(other_taxa):
        return frozenset(descendant_taxa)
    if len(other_taxa) < len(descendant_taxa):
        return frozenset(other_taxa)
    descendant_key = tuple(sorted(descendant_taxa))
    other_key = tuple(sorted(other_taxa))
    return frozenset(descendant_taxa if descendant_key <= other_key else other_taxa)


def _edge_length_vector(tree, comparison):
    subtree_taxa = get_subtree_leaf_name_sets(tree)
    all_taxa = frozenset(tree.leaf_names())
    vector: dict[Any, list[float]] = {}
    for node in tree.traverse():
        if node.is_root:
            continue
        if comparison == "rooted":
            key = frozenset(subtree_taxa[node])
        else:
            key = _canonical_unrooted_split(subtree_taxa[node], all_taxa)
        vector.setdefault(key, []).append(float(node.dist))
    return {key: tuple(lengths) for key, lengths in vector.items()}


def _stable_component_difference(positive, negative):
    scale = max(
        (abs(value) for value in (*positive, *negative)),
        default=0.0,
    )
    if scale == 0.0:
        return 0.0
    difference = scale * math.fsum(
        [
            *(value / scale for value in positive),
            *(-value / scale for value in negative),
        ]
    )
    if not math.isfinite(difference):
        raise ValueError("A branch-length difference is too large to represent.")
    return difference


def _vector_distances(vector1, vector2):
    keys = vector1.keys() | vector2.keys()
    differences = [
        _stable_component_difference(
            vector1.get(key, ()),
            vector2.get(key, ()),
        )
        for key in keys
    ]
    try:
        weighted_rf = math.fsum(abs(difference) for difference in differences)
    except OverflowError as exc:
        raise ValueError("The weighted RF distance is too large to represent.") from exc
    branch_score = 0.0
    for difference in differences:
        branch_score = math.hypot(branch_score, difference)
    if not math.isfinite(weighted_rf) or not math.isfinite(branch_score):
        raise ValueError("The branch-length distance is too large to represent.")
    return weighted_rf, branch_score


def _binary_root_sides(tree):
    root_children = tree.get_children()
    if len(root_children) != 2:
        return None
    sides = dict()
    for side, child in enumerate(root_children):
        for leaf in child.leaves():
            sides[leaf.name] = side
    return sides


def _path_matrix_is_safe(tree1, tree2, num_taxa, topological):
    if num_taxa > PATH_MATRIX_MAX_TAXA:
        return False
    if topological:
        return True
    positive_lengths = list()
    max_length = 0.0
    edge_count = 0
    for tree in (tree1, tree2):
        for node in tree.traverse():
            if node.is_root:
                continue
            edge_count += 1
            length = float(node.dist)
            max_length = max(max_length, length)
            if length > 0.0:
                positive_lengths.append(length)
    if max_length == 0.0:
        return True
    if max_length > (sys.float_info.max / max(4, edge_count * 2)):
        return False
    min_length = min(positive_lengths)
    return (max_length / min_length) <= PATH_MATRIX_MAX_DYNAMIC_RANGE


def _ordered_distance_matrix(tree, leaf_names, topological):
    import numpy as np

    leaves = list(tree.leaves())
    matrix = np.asarray(
        tree.distance_matrix(topological=topological, squared=True),
        dtype=float,
    )
    if topological:
        root_sides = _binary_root_sides(tree)
        if root_sides is not None:
            sides = np.asarray([root_sides[leaf.name] for leaf in leaves])
            matrix -= sides[:, None] != sides[None, :]
    index_by_name = {str(leaf.name): index for index, leaf in enumerate(leaves)}
    order = np.asarray([index_by_name[str(name)] for name in leaf_names])
    return matrix[np.ix_(order, order)]


def _path_distance_matrix(tree1, tree2, leaf_names, topological):
    import numpy as np

    matrix1 = _ordered_distance_matrix(tree1, leaf_names, topological)
    matrix2 = _ordered_distance_matrix(tree2, leaf_names, topological)
    upper = (matrix1 - matrix2)[np.triu_indices(len(leaf_names), 1)]
    if upper.size == 0:
        return 0.0
    scale = float(np.max(np.abs(upper)))
    if scale == 0.0:
        return 0.0
    matrix_scale = max(
        float(np.max(matrix1)),
        float(np.max(matrix2)),
    )
    # Subtracting independently accumulated path sums is ill-conditioned when
    # the trees are nearly identical. Fall back before roundoff can dominate
    # the actual path differences.
    roundoff_bound = (
        PATH_MATRIX_ROUNDOFF_FACTOR
        * sys.float_info.epsilon
        * max(4, 2 * len(leaf_names))
        * matrix_scale
    )
    if scale <= roundoff_bound:
        return None
    distance = scale * math.sqrt(float(np.sum((upper / scale) ** 2)))
    if not math.isfinite(distance):
        raise ValueError("The path distance is too large to represent.")
    return distance


def _path_distance(tree1, tree2, leaf_names, topological):
    if _path_matrix_is_safe(
        tree1,
        tree2,
        num_taxa=len(leaf_names),
        topological=topological,
    ):
        matrix_distance = _path_distance_matrix(
            tree1,
            tree2,
            leaf_names,
            topological,
        )
        if matrix_distance is not None:
            return matrix_distance
    leaves1 = {leaf.name: leaf for leaf in tree1.leaves()}
    leaves2 = {leaf.name: leaf for leaf in tree2.leaves()}
    root_sides1 = _binary_root_sides(tree1) if topological else None
    root_sides2 = _binary_root_sides(tree2) if topological else None

    def path_edge_lengths(first_leaf, second_leaf):
        first_ancestors = set()
        cursor = first_leaf
        while cursor is not None:
            first_ancestors.add(cursor)
            cursor = cursor.up
        common_ancestor = second_leaf
        while common_ancestor not in first_ancestors:
            common_ancestor = common_ancestor.up
        lengths = []
        for leaf in (first_leaf, second_leaf):
            cursor = leaf
            while cursor is not common_ancestor:
                lengths.append(float(cursor.dist))
                cursor = cursor.up
        return lengths

    def differences():
        for leaf_name1, leaf_name2 in combinations(leaf_names, 2):
            leaf1_tree1 = leaves1[leaf_name1]
            leaf2_tree1 = leaves1[leaf_name2]
            leaf1_tree2 = leaves2[leaf_name1]
            leaf2_tree2 = leaves2[leaf_name2]
            if topological:
                distance1 = tree1.get_distance(
                    leaf1_tree1,
                    leaf2_tree1,
                    topological=True,
                )
                distance2 = tree2.get_distance(
                    leaf1_tree2,
                    leaf2_tree2,
                    topological=True,
                )
                if (
                    root_sides1 is not None
                    and root_sides1[leaf_name1] != root_sides1[leaf_name2]
                ):
                    distance1 -= 1
                if (
                    root_sides2 is not None
                    and root_sides2[leaf_name1] != root_sides2[leaf_name2]
                ):
                    distance2 -= 1
                yield float(distance1) - float(distance2)
                continue
            yield _stable_component_difference(
                path_edge_lengths(leaf1_tree1, leaf2_tree1),
                path_edge_lengths(leaf1_tree2, leaf2_tree2),
            )

    distance = 0.0
    for difference in differences():
        distance = math.hypot(distance, difference)
    if not math.isfinite(distance):
        raise ValueError("The path distance is too large to represent.")
    return distance


def _calculate_distances(tree1, tree2, leaf_names, metrics, comparison):
    results = dict()
    if {"rf", "normalized-rf"}.intersection(metrics):
        rf_result = tree1.robinson_foulds(
            t2=tree2,
            unrooted_trees=(comparison == "unrooted"),
        )
        rf_distance, max_rf_distance = rf_result[:2]
        # ETE reports -2 as the rooted maximum for two identical one-tip trees.
        # Mathematically, no split can differ, so the maximum is zero.
        max_rf_distance = max(0, max_rf_distance)
        results["rf"] = (rf_distance, max_rf_distance)
        normalized_rf = 0.0 if max_rf_distance == 0 else rf_distance / max_rf_distance
        results["normalized-rf"] = (normalized_rf, 1.0)
    if {"weighted-rf", "branch-score"}.intersection(metrics):
        vector1 = _edge_length_vector(tree1, comparison)
        vector2 = _edge_length_vector(tree2, comparison)
        weighted_rf, branch_score = _vector_distances(vector1, vector2)
        results["weighted-rf"] = (weighted_rf, "")
        results["branch-score"] = (branch_score, "")
    if "path-topological" in metrics:
        results["path-topological"] = (
            _path_distance(tree1, tree2, leaf_names, topological=True),
            "",
        )
    if "path-length" in metrics:
        results["path-length"] = (
            _path_distance(tree1, tree2, leaf_names, topological=False),
            "",
        )
    return results


def _write_legacy_rf(outfile, result):
    rf_distance, max_rf_distance = result
    text = "rf_dist\tmax_rf_dist\n{}\t{}\n".format(rf_distance, max_rf_distance)
    if outfile == "-":
        sys.stdout.write(text)
        return
    with open(outfile, mode="w") as handle:
        handle.write(text)


def _write_distance_rows(outfile, rows):
    if outfile == "-":
        writer = csv.DictWriter(
            sys.stdout,
            fieldnames=DISTANCE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        return
    with open(outfile, mode="w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=DISTANCE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def dist_main(args):
    metrics, legacy_output = _resolve_metrics(args)
    if getattr(args, "infile2", None) in ["", None]:
        raise ValueError("'--infile2' is required for 'dist'.")
    comparison = getattr(args, "comparison", "rooted")
    if comparison not in {"rooted", "unrooted"}:
        raise ValueError("'--comparison' must be 'rooted' or 'unrooted'.")

    tree1 = read_tree(args.infile, args.format, args.quoted_node_names)
    tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
    validate_unique_named_leaves(tree1, option_name="--infile", context=" for 'dist'")
    validate_unique_named_leaves(tree2, option_name="--infile2", context=" for 'dist'")
    leaf_names1 = sorted(tree1.leaf_names())
    leaf_names2 = sorted(tree2.leaf_names())
    if set(leaf_names1) != set(leaf_names2):
        raise ValueError("Leaf name(s) did not match.")

    if comparison == "rooted" and SPLIT_METRICS.intersection(metrics):
        if (not is_rooted(tree1)) or (not is_rooted(tree2)):
            raise ValueError(
                "Rooted comparison requires rooted trees for split metrics. Root the input trees "
                "first or use '--comparison unrooted'."
            )
    _validate_branch_lengths(tree1, "--infile", metrics)
    _validate_branch_lengths(tree2, "--infile2", metrics)

    results = _calculate_distances(
        tree1,
        tree2,
        leaf_names1,
        metrics,
        comparison,
    )
    if legacy_output:
        _write_legacy_rf(args.outfile, results["rf"])
        return

    rows = list()
    for metric in metrics:
        distance, max_distance = results[metric]
        rows.append(
            {
                "metric": metric,
                "comparison": (
                    "root-independent"
                    if metric in {"path-topological", "path-length"}
                    else comparison
                ),
                "num_taxa": len(leaf_names1),
                "distance": distance,
                "max_distance": max_distance,
            }
        )
    _write_distance_rows(args.outfile, rows)
