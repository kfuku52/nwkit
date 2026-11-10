import csv
import math
import sys
from collections import defaultdict
from itertools import chain
from typing import Any

from nwkit.rf_distance import robinson_foulds
from nwkit.rooting_state import require_rooted
from nwkit.util import (
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


def _edge_postorder(tree, anchor):
    """Yield each node with the original edge and children in the chosen direction."""
    if anchor is None:
        for node in tree.traverse(strategy="postorder"):
            if node is not tree:
                yield node, node, node.children
        return
    if anchor is tree:
        return
    # Orient away from the anchor tip without modifying the input tree.
    order = []
    pending = [(anchor.up, anchor)]
    while pending:
        node, parent = pending.pop()
        if node.up is parent:
            children = node.children
            edge = node
        else:
            children = [child for child in node.children if child is not parent]
            if node is not tree:
                children.append(node.up)
            edge = parent
        order.append((node, edge, children))
        pending.extend((child, node) for child in children)
    yield from reversed(order)


def _edge_interval_vector(tree, positions, anchor):
    extents: dict[Any, tuple[int, int, int, Any]] = {}
    keys = {}
    for node, edge, children in _edge_postorder(tree, anchor):
        key: Any
        if node.is_leaf:
            lo = hi = positions[node.name]
            size = 1
            key = (lo, hi)
        else:
            lo, hi, size = len(positions), -1, 0
            components = 0
            key = ()
            for child in children:
                first, last, count, child_key = extents.pop(child)
                if count:
                    lo = min(lo, first)
                    hi = max(hi, last)
                    size += count
                    components += 1
                    key = child_key
            if components > 1:
                # Noncontiguous clades cannot occur in the reference. Within
                # this tree, only unary paths can repeat them, so reuse the
                # sole child's key above and use node identity for new clades.
                key = (lo, hi) if hi - lo + 1 == size else node
        extents[node] = (lo, hi, size, key)
        keys[edge] = key

    # Keep the input traversal order when a split has multiple edge lengths.
    vector: dict[Any, list[float]] = defaultdict(list)
    for node in tree.traverse():
        if node in keys:
            vector[keys[node]].append(float(node.dist))
    return vector


def _matched_edge_vectors(tree1, tree2, comparison):
    tips = list(tree1.leaves())
    positions = {node.name: index for index, node in enumerate(tips)}
    first_anchor = second_anchor = None
    if comparison == "unrooted":
        # The side excluding the first reference tip is a contiguous range
        # in the reference DFS order, including complements along its path.
        first_anchor = tips[0]
        second_anchor = next(
            node for node in tree2.leaves() if node.name == first_anchor.name
        )
    return (
        _edge_interval_vector(tree1, positions, first_anchor),
        _edge_interval_vector(tree2, positions, second_anchor),
    )


def _stable_component_difference(positive, negative):
    if len(positive) <= 1 and len(negative) <= 1:
        # Preserve the scaling but sum the two signed terms directly.
        left = positive[0] if positive else 0.0
        right = negative[0] if negative else 0.0
        scale = max(abs(left), abs(right)) if positive else abs(right)
        difference = scale * (left / scale - right / scale) if scale else 0.0
    else:
        scale = max(abs(value) for value in chain(positive, negative))
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
    keys = chain(vector1, (key for key in vector2 if key not in vector1))
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
    branch_score = math.hypot(*differences)
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
        # Rounded path sums can agree even when a branch length changed.
        return 0.0 if topological else None
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


def _ancestor_rows(tip_depths, adjacent_lca_depths, order):
    """Find the shallowest adjacent ancestor between requested DFS tips."""
    import numpy as np

    ancestors = np.empty_like(tip_depths)
    for row, left in enumerate(order[:-1]):
        # The shallowest adjacent-tip ancestor between two DFS tips is their LCA.
        np.minimum.accumulate(
            adjacent_lca_depths[:left][::-1], out=ancestors[:left][::-1]
        )
        np.minimum.accumulate(adjacent_lca_depths[left:], out=ancestors[left + 1 :])
        ancestors[left] = tip_depths[left]
        right = order[row + 1 :]
        yield left, right, ancestors[right]


def _topological_path_rows(tree, leaf_names):
    """Yield upper-triangle rows in the requested tip order."""
    import numpy as np

    depths: list[int] = []
    adjacent_lca_depths: list[int] = []
    tip_positions = {}
    depth = shallowest = 0
    for is_post, node in tree.iter_prepostorder():
        if is_post:
            depth -= 1
            # The next subtree is reached via this node's parent.
            shallowest = min(shallowest, depth - 1)
        elif node.children:
            depth += 1
        else:
            tip_positions[node.name] = len(depths)
            if depths:
                adjacent_lca_depths.append(shallowest)
            depths.append(depth)
            shallowest = depth - 1

    tip_depths = np.asarray(depths, dtype=np.int64)
    gaps = np.asarray(adjacent_lca_depths, dtype=np.int64)
    order = np.asarray([tip_positions[name] for name in leaf_names], dtype=np.intp)
    suppress_binary_root = len(tree.children) == 2
    for left, right, lca_depths in _ancestor_rows(tip_depths, gaps, order):
        distances = tip_depths[left] + tip_depths[right] - 2 * lca_depths
        if suppress_binary_root:
            # Crossing a degree-two root uses one edge after suppression.
            distances -= lca_depths == 0
        yield distances


def _topological_path_distance(tree1, tree2, leaf_names):
    import numpy as np

    # Each path uses fewer edges than the number of nodes in its tree.
    max_nodes = max(sum(1 for _ in tree.traverse()) for tree in (tree1, tree2))
    int64_safe = len(leaf_names) * max_nodes**2 <= np.iinfo(np.int64).max
    rows1 = _topological_path_rows(tree1, leaf_names)
    rows2 = _topological_path_rows(tree2, leaf_names)
    squared = 0
    for first, second in zip(rows1, rows2, strict=True):
        differences = first - second
        if int64_safe:
            squared += int(differences @ differences)
        else:
            squared += sum(int(value) ** 2 for value in differences)
    return _rounded_dyadic_norm(squared, 0)


def _weighted_path_coordinates(tree, leaf_names, units, exponent):
    import numpy as np

    tip_depths: list[int] = []
    adjacent_lca_depths: list[int] = []
    tip_positions = {}
    depth: dict[Any, int] = {tree: 0}
    shallowest = 0
    for is_post, node in tree.iter_prepostorder():
        if is_post:
            if node is not tree:
                shallowest = min(shallowest, depth[node.up])
            continue
        if node is not tree:
            numerator, edge_exponent = units[node]
            depth[node] = depth[node.up] + (numerator << (exponent - edge_exponent))
        if node.is_leaf:
            tip_positions[node.name] = len(tip_depths)
            if tip_depths:
                adjacent_lca_depths.append(shallowest)
            tip_depths.append(depth[node])
            if node is not tree:
                shallowest = depth[node.up]

    # Even when depths need Python integers, their ranks fit in NumPy arrays.
    coordinates = sorted(set(tip_depths).union(adjacent_lca_depths))
    rank = {value: index for index, value in enumerate(coordinates)}
    tip_ranks = np.asarray([rank[value] for value in tip_depths], dtype=np.intp)
    gap_ranks = np.asarray(
        [rank[value] for value in adjacent_lca_depths], dtype=np.intp
    )
    order = np.asarray([tip_positions[name] for name in leaf_names], dtype=np.intp)
    return coordinates, tip_ranks, gap_ranks, order


def _weighted_path_rows(plan, dtype):
    import numpy as np

    coordinates, tip_ranks, gap_ranks, order = plan
    levels = np.asarray(coordinates, dtype=dtype)
    tips = levels[tip_ranks]
    for left, right, ancestors in _ancestor_rows(tip_ranks, gap_ranks, order):
        yield tips[left] + tips[right] - 2 * levels[ancestors]


def _rounded_dyadic_norm(squared, unit_exponent):
    """Round sqrt(squared) / 2**unit_exponent directly to a float."""
    if squared == 0:
        return 0.0
    # Round sqrt(squared) in units of the result's floating-point ulp.
    root_exponent = (squared.bit_length() - 1) // 2
    mantissa_bits = sys.float_info.mant_dig
    ulp_exponent = max(
        sys.float_info.min_exp - mantissa_bits,
        root_exponent - unit_exponent - mantissa_bits + 1,
    )
    shift = unit_exponent + ulp_exponent
    if shift < 0:
        squared <<= -2 * shift
        shift = 0
    lower = math.isqrt(squared) >> shift
    midpoint_squared = ((2 * lower + 1) << shift) ** 2
    scaled_squared = squared << 2
    if scaled_squared > midpoint_squared or (
        scaled_squared == midpoint_squared and lower % 2
    ):
        lower += 1
    return math.ldexp(lower, ulp_exponent)


def _weighted_path_distance(tree1, tree2, leaf_names):
    import numpy as np

    first_vector, second_vector = _matched_edge_vectors(tree1, tree2, "unrooted")
    if first_vector == second_vector:
        return 0.0
    del first_vector, second_vector

    # Finite float lengths are integer multiples of a common power of two.
    units: dict[Any, tuple[int, int]] = {}
    exponent = 0
    for tree in (tree1, tree2):
        for node in tree.traverse("preorder"):
            if node is tree:
                continue
            numerator, denominator = float(node.dist).as_integer_ratio()
            edge_exponent = denominator.bit_length() - 1
            exponent = max(exponent, edge_exponent)
            units[node] = (numerator, edge_exponent)
    first = _weighted_path_coordinates(tree1, leaf_names, units, exponent)
    second = _weighted_path_coordinates(tree2, leaf_names, units, exponent)
    bound = 2 * max(first[0][-1], second[0][-1])
    int64_safe = len(leaf_names) * bound**2 <= np.iinfo(np.int64).max
    dtype = np.int64 if int64_safe else object
    # Preserve the distinction between an overflowing pair and the total norm.
    half_ulp = 1 << (sys.float_info.max_exp - sys.float_info.mant_dig - 1)
    pair_limit = (int(sys.float_info.max) + half_ulp) << exponent
    check_pairs = bound >= pair_limit
    squared = 0
    for a, b in zip(
        _weighted_path_rows(first, dtype),
        _weighted_path_rows(second, dtype),
        strict=True,
    ):
        differences = a - b
        if check_pairs and any(abs(value) >= pair_limit for value in differences):
            raise ValueError("A branch-length difference is too large to represent.")
        squared += int(differences @ differences)
    try:
        return _rounded_dyadic_norm(squared, exponent)
    except OverflowError:
        raise ValueError("The path distance is too large to represent.") from None


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
    if topological:
        return _topological_path_distance(tree1, tree2, leaf_names)
    return _weighted_path_distance(tree1, tree2, leaf_names)


def _calculate_distances(tree1, tree2, leaf_names, metrics, comparison):
    results = dict()
    if {"rf", "normalized-rf"}.intersection(metrics):
        rf_distance, max_rf_distance = robinson_foulds(
            tree1,
            tree2,
            rooted=(comparison == "rooted"),
        )
        results["rf"] = (rf_distance, max_rf_distance)
        normalized_rf = 0.0 if max_rf_distance == 0 else rf_distance / max_rf_distance
        results["normalized-rf"] = (normalized_rf, 1.0)
    if {"weighted-rf", "branch-score"}.intersection(metrics):
        vector1, vector2 = _matched_edge_vectors(tree1, tree2, comparison)
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

    tree1 = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    tree2 = read_tree(
        args.infile2,
        args.format2,
        args.quoted_node_names,
        rooted=getattr(args, "infile2_rooted", "auto"),
    )
    validate_unique_named_leaves(tree1, option_name="--infile", context=" for 'dist'")
    validate_unique_named_leaves(tree2, option_name="--infile2", context=" for 'dist'")
    leaf_names1 = sorted(tree1.leaf_names())
    leaf_names2 = sorted(tree2.leaf_names())
    if set(leaf_names1) != set(leaf_names2):
        raise ValueError("Leaf name(s) did not match.")

    if comparison == "rooted" and SPLIT_METRICS.intersection(metrics):
        for tree, option in ((tree1, "--input-rooted"), (tree2, "--infile2-rooted")):
            require_rooted(
                tree,
                "Rooted comparison requires rooted trees for split metrics; "
                "use '--comparison unrooted' for root-independent splits.",
                option,
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
