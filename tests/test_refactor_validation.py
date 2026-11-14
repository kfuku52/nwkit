"""Independent numerical checks for refactored tree-distance statistics."""

import random
from statistics import mean, pvariance

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


def test_collect_leaf_distance_moments_matches_bruteforce():
    """The optimized accumulator must agree with an independent oracle."""
    rng = random.Random(53)
    for _ in range(30):
        tree = make_random_binary_tree(
            num_leaves=rng.randint(5, 10), rng=rng, leaf_prefix="R"
        )
        inside, outside = root._collect_leaf_distance_moments(tree)
        all_leaves = list(tree.leaves())
        for node in tree.traverse():
            sub_leaves = set(node.leaves())
            sub_dists = [tree.get_distance(node, leaf) for leaf in sub_leaves]
            count, average, ss = inside[node]
            assert count == len(sub_dists)
            assert average == pytest.approx(mean(sub_dists), abs=1e-9)
            assert ss == pytest.approx(count * pvariance(sub_dists), abs=1e-9)
            if node.is_root:
                assert outside[node][0] == 0
                continue
            other_dists = [
                tree.get_distance(node.up, leaf)
                for leaf in all_leaves
                if leaf not in sub_leaves
            ]
            count, average, ss = outside[node]
            assert count == len(other_dists)
            assert average == pytest.approx(mean(other_dists), abs=1e-9)
            assert ss == pytest.approx(count * pvariance(other_dists), abs=1e-9)
