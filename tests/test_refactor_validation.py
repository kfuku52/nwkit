"""Independent numerical checks for refactored tree-distance statistics."""

import random

import pytest
from ete4 import Tree

import nwkit.root as root


def make_random_binary_tree(num_leaves, rng, leaf_prefix="L"):
    nodes = [f"{leaf_prefix}{i}:{rng.uniform(0.1, 2.0):.6f}" for i in range(num_leaves)]
    while len(nodes) > 1:
        rng.shuffle(nodes)
        left = nodes.pop()
        right = nodes.pop()
        nodes.append(f"({left},{right}):{rng.uniform(0.1, 2.0):.6f}")
    return Tree(nodes[0] + ";", parser=1)


def test_collect_leaf_distance_stats_matches_bruteforce():
    """The optimized accumulator must agree with an independent oracle."""
    rng = random.Random(53)
    for _ in range(30):
        tree = make_random_binary_tree(
            num_leaves=rng.randint(5, 10), rng=rng, leaf_prefix="R"
        )
        subtree_stats, all_stats = root._collect_leaf_distance_stats(tree)
        all_leaves = list(tree.leaves())
        for node in tree.traverse():
            sub_leaves = list(node.leaves())
            sub_dists = [tree.get_distance(node, leaf) for leaf in sub_leaves]
            s_count, s_sum, s_sumsq = subtree_stats[node]
            assert s_count == len(sub_dists)
            assert s_sum == pytest.approx(sum(sub_dists), abs=1e-9)
            assert s_sumsq == pytest.approx(sum(d * d for d in sub_dists), abs=1e-9)
            all_dists = [tree.get_distance(node, leaf) for leaf in all_leaves]
            a_count, a_sum, a_sumsq = all_stats[node]
            assert a_count == len(all_dists)
            assert a_sum == pytest.approx(sum(all_dists), abs=1e-9)
            assert a_sumsq == pytest.approx(sum(d * d for d in all_dists), abs=1e-9)
