import random

import pytest
from ete4 import Tree

from nwkit.rf_distance import robinson_foulds


def _random_binary_tree(names, rng):
    if len(names) == 1:
        return names[0] + ":1"
    names = rng.sample(names, len(names))
    split = rng.randrange(1, len(names))
    return "({},{})1:1".format(
        _random_binary_tree(names[:split], rng),
        _random_binary_tree(names[split:], rng),
    )


@pytest.mark.parametrize("rooted", [True, False])
def test_binary_rf_matches_ete(rooted):
    rng = random.Random(20)
    names = [f"T{i}" for i in range(16)]
    trees = [Tree(_random_binary_tree(names, rng) + ";", parser=0) for _ in range(8)]
    for first in trees:
        for second in trees:
            expected = first.robinson_foulds(second, unrooted_trees=not rooted)
            assert robinson_foulds(first, second, rooted=rooted) == tuple(expected[:2])


@pytest.mark.parametrize("rooted,expected", [(True, (2, 2)), (False, (1, 1))])
def test_rf_between_star_and_binary_uses_only_displayed_clades(rooted, expected):
    star = Tree("(A:1,B:1,C:1,D:1);", parser=1)
    resolved = Tree("(A:1,(B:1,(C:1,D:1):1):1);", parser=1)
    assert robinson_foulds(star, star, rooted=rooted) == (0, 0)
    assert robinson_foulds(star, resolved, rooted=rooted) == expected


@pytest.mark.parametrize("rooted", [True, False])
def test_rf_does_not_count_root_stems_or_duplicate_unary_clades(rooted):
    first = Tree("((A:1,B:1):1,(C:1,D:1):1):100;", parser=1)
    second = Tree("((((A:1,B:1):1):1,(C:1,D:1):1):1):200;", parser=1)
    assert robinson_foulds(first, second, rooted=rooted)[0] == 0
