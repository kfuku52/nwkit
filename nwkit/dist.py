import csv
import math
import sys
from itertools import combinations

from nwkit.util import (
    get_subtree_leaf_name_sets,
    is_rooted,
    read_tree,
    validate_unique_named_leaves,
)


SUPPORTED_METRICS = (
    'rf',
    'normalized-rf',
    'weighted-rf',
    'branch-score',
    'path-topological',
    'path-length',
)
DISTANCE_COLUMNS = (
    'metric',
    'comparison',
    'num_taxa',
    'distance',
    'max_distance',
)
SPLIT_METRICS = frozenset({
    'rf',
    'normalized-rf',
    'weighted-rf',
    'branch-score',
})
BRANCH_LENGTH_METRICS = frozenset({
    'weighted-rf',
    'branch-score',
    'path-length',
})


def _parse_metric_values(values):
    metrics = list()
    has_all = False
    for value in values:
        for item in str(value).split(','):
            metric = item.strip().lower().replace('_', '-')
            if metric == '':
                raise ValueError('Empty metric in `--metric`.')
            if metric == 'all':
                has_all = True
                continue
            if metric not in SUPPORTED_METRICS:
                raise ValueError(
                    "Unsupported metric '{}'. Choose from: all, {}.".format(
                        item.strip(),
                        ', '.join(SUPPORTED_METRICS),
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
    canonical_values = getattr(args, 'metric', None)
    legacy_value = getattr(args, 'dist', None)
    if canonical_values is not None and legacy_value is not None:
        raise ValueError("Use either '--metric' or deprecated '--dist', not both.")
    if legacy_value is not None:
        metrics = _parse_metric_values([legacy_value])
        legacy_output = metrics == ('rf',)
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
                "a length is missing in '{}'.".format(', '.join(required_by), option_name)
            )
        try:
            branch_length = float(node.dist)
        except (TypeError, ValueError):
            branch_length = math.nan
        if not math.isfinite(branch_length) or branch_length < 0:
            raise ValueError(
                "Metric(s) {} require finite, nonnegative branch lengths; "
                "found {!r} in '{}'.".format(
                    ', '.join(required_by),
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
    vector = dict()
    for node in tree.traverse():
        if node.is_root:
            continue
        if comparison == 'rooted':
            key = frozenset(subtree_taxa[node])
        else:
            key = _canonical_unrooted_split(subtree_taxa[node], all_taxa)
        vector[key] = vector.get(key, 0.0) + float(node.dist)
    return vector


def _vector_distances(vector1, vector2):
    keys = vector1.keys() | vector2.keys()
    differences = [vector1.get(key, 0.0) - vector2.get(key, 0.0) for key in keys]
    weighted_rf = math.fsum(abs(difference) for difference in differences)
    branch_score = math.sqrt(math.fsum(difference ** 2 for difference in differences))
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


def _path_distance(tree1, tree2, leaf_names, topological):
    leaves1 = {leaf.name: leaf for leaf in tree1.leaves()}
    leaves2 = {leaf.name: leaf for leaf in tree2.leaves()}
    root_sides1 = _binary_root_sides(tree1) if topological else None
    root_sides2 = _binary_root_sides(tree2) if topological else None

    def differences():
        for leaf_name1, leaf_name2 in combinations(leaf_names, 2):
            leaf1_tree1 = leaves1[leaf_name1]
            leaf2_tree1 = leaves1[leaf_name2]
            leaf1_tree2 = leaves2[leaf_name1]
            leaf2_tree2 = leaves2[leaf_name2]
            distance1 = tree1.get_distance(
                leaf1_tree1,
                leaf2_tree1,
                topological=topological,
            )
            distance2 = tree2.get_distance(
                leaf1_tree2,
                leaf2_tree2,
                topological=topological,
            )
            if topological:
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

    return math.sqrt(math.fsum(difference ** 2 for difference in differences()))


def _calculate_distances(tree1, tree2, leaf_names, metrics, comparison):
    results = dict()
    if {'rf', 'normalized-rf'}.intersection(metrics):
        rf_result = tree1.robinson_foulds(
            t2=tree2,
            unrooted_trees=(comparison == 'unrooted'),
        )
        rf_distance, max_rf_distance = rf_result[:2]
        # ETE reports -2 as the rooted maximum for two identical one-tip trees.
        # Mathematically, no split can differ, so the maximum is zero.
        max_rf_distance = max(0, max_rf_distance)
        results['rf'] = (rf_distance, max_rf_distance)
        normalized_rf = 0.0 if max_rf_distance == 0 else rf_distance / max_rf_distance
        results['normalized-rf'] = (normalized_rf, 1.0)
    if {'weighted-rf', 'branch-score'}.intersection(metrics):
        vector1 = _edge_length_vector(tree1, comparison)
        vector2 = _edge_length_vector(tree2, comparison)
        weighted_rf, branch_score = _vector_distances(vector1, vector2)
        results['weighted-rf'] = (weighted_rf, '')
        results['branch-score'] = (branch_score, '')
    if 'path-topological' in metrics:
        results['path-topological'] = (
            _path_distance(tree1, tree2, leaf_names, topological=True),
            '',
        )
    if 'path-length' in metrics:
        results['path-length'] = (
            _path_distance(tree1, tree2, leaf_names, topological=False),
            '',
        )
    return results


def _write_legacy_rf(outfile, result):
    rf_distance, max_rf_distance = result
    text = 'rf_dist\tmax_rf_dist\n{}\t{}\n'.format(rf_distance, max_rf_distance)
    if outfile == '-':
        sys.stdout.write(text)
        return
    with open(outfile, mode='w') as handle:
        handle.write(text)


def _write_distance_rows(outfile, rows):
    if outfile == '-':
        writer = csv.DictWriter(
            sys.stdout,
            fieldnames=DISTANCE_COLUMNS,
            delimiter='\t',
            lineterminator='\n',
        )
        writer.writeheader()
        writer.writerows(rows)
        return
    with open(outfile, mode='w', newline='') as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=DISTANCE_COLUMNS,
            delimiter='\t',
            lineterminator='\n',
        )
        writer.writeheader()
        writer.writerows(rows)


def dist_main(args):
    metrics, legacy_output = _resolve_metrics(args)
    if getattr(args, 'infile2', None) in ['', None]:
        raise ValueError("'--infile2' is required for 'dist'.")
    comparison = getattr(args, 'comparison', 'rooted')
    if comparison not in {'rooted', 'unrooted'}:
        raise ValueError("'--comparison' must be 'rooted' or 'unrooted'.")

    tree1 = read_tree(args.infile, args.format, args.quoted_node_names)
    tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
    validate_unique_named_leaves(tree1, option_name='--infile', context=" for 'dist'")
    validate_unique_named_leaves(tree2, option_name='--infile2', context=" for 'dist'")
    leaf_names1 = sorted(tree1.leaf_names())
    leaf_names2 = sorted(tree2.leaf_names())
    if set(leaf_names1) != set(leaf_names2):
        raise ValueError('Leaf name(s) did not match.')

    if comparison == 'rooted' and SPLIT_METRICS.intersection(metrics):
        if (not is_rooted(tree1)) or (not is_rooted(tree2)):
            raise ValueError(
                'Rooted comparison requires rooted trees for split metrics. Root the input trees '
                "first or use '--comparison unrooted'."
            )
    _validate_branch_lengths(tree1, '--infile', metrics)
    _validate_branch_lengths(tree2, '--infile2', metrics)

    results = _calculate_distances(
        tree1,
        tree2,
        leaf_names1,
        metrics,
        comparison,
    )
    if legacy_output:
        _write_legacy_rf(args.outfile, results['rf'])
        return

    rows = list()
    for metric in metrics:
        distance, max_distance = results[metric]
        rows.append({
            'metric': metric,
            'comparison': (
                'root-independent'
                if metric in {'path-topological', 'path-length'}
                else comparison
            ),
            'num_taxa': len(leaf_names1),
            'distance': distance,
            'max_distance': max_distance,
        })
    _write_distance_rows(args.outfile, rows)
