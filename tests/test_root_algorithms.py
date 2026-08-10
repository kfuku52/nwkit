import math
import time

import pytest
from ete4 import Tree

import nwkit.root as root_mod
from nwkit.clade_mapping import canonical_split
from nwkit.root import (
    mad_rooting,
    midpoint_rooting,
    mv_rooting,
    outgroup_rooting,
    transfer_root,
)
from nwkit.util import is_rooted
from tests.helpers import make_deep_ladder_tree, safe_get_distance
from tests.root_test_support import (
    annotated_reroot_tree as _annotated_reroot_tree,
)
from tests.root_test_support import (
    assert_reroot_annotations as _assert_reroot_annotations,
)
from tests.root_test_support import (
    scaled_branch_profile as _scaled_branch_profile,
)


class TestMidpointRooting:
    @pytest.mark.slow
    def test_deep_ladder_does_not_exceed_python_recursion_limit(self):
        tree = make_deep_ladder_tree(1200)

        rooted = midpoint_rooting(tree)

        assert len(list(rooted.leaves())) == 1200
        assert set(rooted.leaf_names()) == set(tree.leaf_names())

    def test_basic(self):
        tree = Tree("(A:1,B:3,(C:2,D:4):2);", parser=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)

    def test_already_rooted(self):
        tree = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {"A", "B", "C", "D"}

    def test_symmetric_clock_exact_distances(self):
        """On a symmetric clock tree, all root-to-tip distances are equal."""
        tree = Tree("((A:2,B:2):1,(C:2,D:2):1);", parser=1)
        rooted = midpoint_rooting(tree)
        dists = [safe_get_distance(rooted, rooted, l) for l in rooted.leaves()]
        assert all(abs(d - 3.0) < 1e-6 for d in dists)

    def test_pairwise_distance_preservation(self):
        """Midpoint rooting must preserve all pairwise distances."""
        nwk = "((A:3,B:1):2,(C:4,D:2):1);"
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = midpoint_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert (
                        abs(
                            original.get_distance(l1.name, l2.name)
                            - rooted.get_distance(l1, l2)
                        )
                        < 1e-6
                    )

    def test_nonzero_root_distance_does_not_crash(self):
        tree = Tree("((A:1,B:3):1,(C:1,D:1):1):2;", parser=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {"A", "B", "C", "D"}

    def test_path_overflow_scale_invariance(self):
        subtree = "(T0:1,T1:1):1"
        for index in range(2, 23):
            subtree = "({},T{}:1):1".format(subtree, index)
        tree = Tree("({},T23:1);".format(subtree), parser=1)
        huge_tree = tree.copy(method="deepcopy")
        scale = 10**307
        for node in huge_tree.traverse():
            if not node.is_root:
                node.dist *= scale

        rooted = midpoint_rooting(tree)
        huge_rooted = midpoint_rooting(huge_tree)

        expected_root_sides = {
            frozenset(child.leaf_names()) for child in rooted.get_children()
        }
        observed_root_sides = {
            frozenset(child.leaf_names()) for child in huge_rooted.get_children()
        }
        assert observed_root_sides == expected_root_sides
        expected_profile = _scaled_branch_profile(rooted)
        observed_profile = _scaled_branch_profile(huge_rooted, scale=scale)
        assert observed_profile.keys() == expected_profile.keys()
        for split in expected_profile:
            assert observed_profile[split] == pytest.approx(
                expected_profile[split],
            )
        assert all(
            math.isfinite(node.dist)
            for node in huge_rooted.traverse()
            if not node.is_root
        )


class TestOutgroupRooting:
    def test_single_outgroup(self):
        tree = Tree("(A:1,(B:1,(C:1,D:1):1):1);", parser=1)
        rooted = outgroup_rooting(tree, "A")
        assert is_rooted(rooted)
        # A should be one of the subroot children
        subroot_children = rooted.get_children()
        outgroup_leaves = set()
        for child in subroot_children:
            outgroup_leaves.update(child.leaf_names())
        assert "A" in outgroup_leaves

    def test_multiple_outgroups(self):
        tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", parser=1)
        rooted = outgroup_rooting(tree, "D,E")
        assert is_rooted(rooted)
        leaf_names = set(rooted.leaf_names())
        assert leaf_names == {"A", "B", "C", "D", "E"}

    def test_multiple_outgroups_with_spaces(self):
        tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", parser=1)
        rooted = outgroup_rooting(tree, "D, E")
        children = rooted.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {"D", "E"} in child_leaf_sets

    def test_outgroup_not_found(self):
        tree = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)
        with pytest.raises(ValueError, match="Outgroup label"):
            outgroup_rooting(tree, "Z")

    def test_outgroup_partially_missing_raises(self):
        tree = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)
        with pytest.raises(ValueError, match="Outgroup label"):
            outgroup_rooting(tree, "A,Z")

    def test_outgroup_whole_tree_raises(self):
        tree = Tree("(A:1,B:1);", parser=1)
        with pytest.raises(ValueError, match="whole tree"):
            outgroup_rooting(tree, "A,B")

    def test_outgroup_single_leaf_tree_raises(self):
        tree = Tree("A;", parser=1)
        with pytest.raises(ValueError, match="whole tree"):
            outgroup_rooting(tree, "A")

    def test_nonzero_root_distance_does_not_crash(self):
        tree = Tree("((A:1,B:1):1,(C:1,D:1):1):2;", parser=1)
        rooted = outgroup_rooting(tree, "A")
        children = rooted.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {"A"} in child_leaf_sets

    def test_pairwise_distance_preservation(self):
        """Outgroup rooting must preserve all pairwise distances."""
        nwk = "((A:2,B:3):1,(C:4,D:5):1);"
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = outgroup_rooting(tree, "A")
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert (
                        abs(
                            original.get_distance(l1.name, l2.name)
                            - rooted.get_distance(l1, l2)
                        )
                        < 1e-6
                    )

    def test_multiple_outgroup_exact_bipartition(self):
        """Outgroup rooting with multiple tips creates correct bipartition."""
        tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", parser=1)
        rooted = outgroup_rooting(tree, "D,E")
        children = rooted.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {"D", "E"} in child_leaf_sets
        assert {"A", "B", "C"} in child_leaf_sets


class TestMadRooting:
    def test_basic(self):
        tree = Tree("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);", parser=1)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)

    def test_preserves_leaves(self):
        tree = Tree("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);", parser=1)
        rooted = mad_rooting(tree)
        assert set(rooted.leaf_names()) == {"A", "B", "C", "D", "E"}

    def test_asymmetric_branch_lengths(self):
        tree = Tree("((A:10,B:1):1,(C:1,D:1):1);", parser=1)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {"A", "B", "C", "D"}

    def test_preserves_quoted_leaf_names_with_punctuation(self):
        tree = Tree()
        ingroup = tree.add_child(dist=1.0)
        ingroup.add_child(name="A,B", dist=1.0)
        ingroup.add_child(name="C:D", dist=2.0)
        tree.add_child(name="E(1)'", dist=3.0)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {"A,B", "C:D", "E(1)'"}

    def test_requires_at_least_three_leaves(self):
        tree = Tree("(A:1,B:1);", parser=1)
        with pytest.raises(ValueError, match="at least 3 leaves"):
            mad_rooting(tree)

    def test_requires_at_least_three_positive_branches(self):
        tree = Tree("(A:0,B:1,C:1,D:0);", parser=1)
        with pytest.raises(ValueError, match="at least 3 positive branch lengths"):
            mad_rooting(tree)

    def test_requires_every_branch_length(self):
        tree = Tree("(A,B:1,C:2,D:3);", parser=1)

        with pytest.raises(ValueError, match="every branch length"):
            mad_rooting(tree)

    @pytest.mark.parametrize(
        "bad_length",
        [-1.0, float("inf"), float("nan")],
        ids=("negative", "infinite", "nan"),
    )
    def test_rejects_invalid_branch_lengths(self, bad_length):
        tree = Tree("(A:1,B:2,C:3,D:4);", parser=1)
        next(leaf for leaf in tree.leaves() if leaf.name == "A").dist = bad_length

        with pytest.raises(ValueError, match="finite, non-negative"):
            mad_rooting(tree)

    def test_zero_distance_tips_are_retained_as_one_effective_tip(self):
        tree = Tree("(E:3,D:2,C:1,B:0,A:0);", parser=1)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {"A", "B", "C", "D", "E"}

    def test_tied_root_edges_use_a_deterministic_split(self):
        expected_split = canonical_split(
            frozenset({"A"}),
            frozenset({"B", "C", "D"}),
        )
        observed_splits = list()
        for newick in (
            "(D:1,C:1,B:1,A:1);",
            "(B:1,D:1,A:1,C:1);",
        ):
            rooted = mad_rooting(Tree(newick, parser=1))
            children = rooted.get_children()
            observed_splits.append(
                canonical_split(
                    frozenset(children[0].leaf_names()),
                    frozenset(children[1].leaf_names()),
                )
            )
        assert observed_splits == [expected_split, expected_split]

    def test_root_position_is_invariant_to_branch_length_scale(self):
        base_tree = Tree(
            "((A:1,B:2):1,(C:3,(D:1,E:2):1):1);",
            parser=1,
        )
        observations = list()
        for scale in (10**-160, 1.0, 10**160):
            scaled_tree = base_tree.copy(method="deepcopy")
            for node in scaled_tree.traverse():
                if not node.is_root:
                    node.dist = float(node.dist) * scale
            rooted = mad_rooting(scaled_tree)
            child_dist_by_taxa = {
                frozenset(child.leaf_names()): child.dist
                for child in rooted.get_children()
            }
            root_total = sum(child_dist_by_taxa.values())
            observations.append(
                (
                    frozenset(child_dist_by_taxa),
                    child_dist_by_taxa[frozenset({"A", "B"})] / root_total,
                )
            )

        assert all(
            split
            == frozenset(
                (
                    frozenset({"A", "B"}),
                    frozenset({"C", "D", "E"}),
                )
            )
            for split, _ in observations
        )
        assert [position for _, position in observations] == pytest.approx(
            [observations[1][1]] * 3,
            rel=10**-12,
            abs=10**-12,
        )

    def test_extreme_within_tree_length_range_does_not_overflow(self):
        rooted = mad_rooting(
            Tree(
                "(A:1e-160,B:1e-160,(C:1,D:1):1);",
                parser=1,
            )
        )

        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {"A", "B", "C", "D"}
        assert all(
            math.isfinite(float(node.dist))
            for node in rooted.traverse()
            if not node.is_root
        )

    def test_heterogeneous_lengths_do_not_corrupt_path_delta_cancellation(self):
        tree = Tree(
            "(X1:1086.5437592639944,"
            "(X3:0.00085007102927072819,"
            "(X2:1.0785237796461825e-06,"
            "X0:0.00057934556794108716):75531.807545509451)"
            ":0.00763325784066344);",
            parser=1,
        )

        rooted = mad_rooting(tree)

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist for child in rooted.get_children()
        }
        assert set(child_dist_by_taxa) == {
            frozenset({"X0", "X2"}),
            frozenset({"X1", "X3"}),
        }
        assert child_dist_by_taxa[frozenset({"X0", "X2"})] == pytest.approx(
            38033.662192492,
        )

    def test_deep_short_clade_does_not_lose_root_distance_precision(self):
        tree = Tree(
            "(X4:0.0063973377995121773,"
            "X3:9.9523327390797709e-06,"
            "(X0:87377874218.288742,"
            "X1:3.9277797271534321e-10):3.7230847581009794e-12,"
            "X2:1540148.2696119624);",
            parser=1,
        )

        rooted = mad_rooting(tree)

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist for child in rooted.get_children()
        }
        assert frozenset({"X0"}) in child_dist_by_taxa
        assert child_dist_by_taxa[frozenset({"X0"})] == pytest.approx(
            43689129622.58888,
        )

    def test_moderate_dynamic_range_uses_exact_root_position(self):
        tree = Tree(
            "((((T6:5.2438826413098062e-08,T5:0.012135336082789326)"
            ":2.6598233474091526e-05,T0:0.06529318198647692)"
            ":3.4741992558598983e-05,(T1:0.0083846338251795938,"
            "T8:0.69976274685535922):1.3509927071941638e-07)"
            ":0.0017859630432441663,((T7:4.5932841447158775e-08,"
            "T3:1.0848429873742806e-08):0.84197554627230875,"
            "(T2:0.001082125534087372,(T9:0.0044256953801056686,"
            "T4:0.00035731219490173113):3.930847287325036e-06)"
            ":0.22936708849868917):0.6857711026599238);",
            parser=1,
        )

        rooted = mad_rooting(tree)

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): float(child.dist)
            for child in rooted.get_children()
        }
        assert child_dist_by_taxa[
            frozenset({"T2", "T3", "T4", "T7", "T9"})
        ] == pytest.approx(0.2050488039025591, rel=10**-12)

    def test_finite_root_halves_are_joined_without_intermediate_overflow(self):
        maximum_float = float.fromhex("0x1.fffffffffffffp+1023")
        tree = Tree(
            "((A:1,B:2):{0},(C:3,D:4):{0});".format(maximum_float),
            parser=1,
        )

        rooted = mad_rooting(tree)

        assert all(
            math.isfinite(float(node.dist))
            for node in rooted.traverse()
            if not node.is_root
        )
        assert [child.dist for child in rooted.get_children()] == pytest.approx(
            [maximum_float, maximum_float],
        )

    def test_unrepresentable_normalized_dynamic_range_is_rejected(self):
        tree = Tree(
            "(A:5e-324,B:1,C:2,D:1e308);",
            parser=1,
        )

        with pytest.raises(ValueError, match="branch-length dynamic range"):
            mad_rooting(tree)

    def test_scoring_does_not_walk_tree_paths_for_every_pair(self):
        tree = Tree(
            "((A:1,B:2):1,(C:3,(D:1,E:2):1):1);",
            parser=1,
        )

        assert "get_distance" not in mad_rooting.__code__.co_names
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)

    @pytest.mark.slow
    def test_balanced_320_tip_tree_completes_within_generous_budget(self):
        def balanced_newick(labels, depth=0):
            if len(labels) == 1:
                index = int(labels[0][1:])
                length = 0.7 + ((index * 17 + depth * 3) % 19) / 10
                return "{}:{}".format(labels[0], length)
            midpoint = len(labels) // 2
            length = 0.5 + ((len(labels) * 11 + depth) % 13) / 10
            return "({},{}):{}".format(
                balanced_newick(labels[:midpoint], depth + 1),
                balanced_newick(labels[midpoint:], depth + 1),
                length,
            )

        labels = ["T{}".format(index) for index in range(320)]
        tree = Tree(balanced_newick(labels) + ";", parser=1)

        started = time.perf_counter()
        rooted = mad_rooting(tree)
        elapsed = time.perf_counter() - started

        assert set(rooted.leaf_names()) == set(labels)
        assert elapsed < 5.0

    def test_wiki_exact_root_split(self):
        """MAD rooting on wiki tree: verify root split and pairwise distances.

        Input: ((A:1,B:2):1,(C:3,(D:1,E:2):1):1);
        Wiki output: ((A:1,B:2):1.58461,(C:3,(D:1,E:2):1):0.415389);
        Total root branch ≈ 2.0
        """
        tree = Tree("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);", parser=1)
        original = Tree("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);", parser=1)
        rooted = mad_rooting(tree)
        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist for child in rooted.get_children()
        }
        assert set(child_dist_by_taxa) == {
            frozenset({"A", "B"}),
            frozenset({"C", "D", "E"}),
        }
        assert child_dist_by_taxa[frozenset({"A", "B"})] == pytest.approx(
            1.584611134,
            abs=10**-8,
        )
        assert child_dist_by_taxa[frozenset({"C", "D", "E"})] == pytest.approx(
            0.415388866,
            abs=10**-8,
        )
        # Total root branch should sum to ≈ 2.0
        children = rooted.get_children()
        total = sum(c.dist for c in children)
        assert abs(total - 2.0) < 0.01
        # Pairwise distances must be preserved
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    orig_d = original.get_distance(l1.name, l2.name)
                    new_d = rooted.get_distance(l1, l2)
                    assert abs(orig_d - new_d) < 0.01


class TestMvRooting:
    def test_basic(self):
        tree = Tree("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);", parser=1)
        rooted = mv_rooting(tree)
        assert is_rooted(rooted)

    def test_preserves_leaves(self):
        tree = Tree("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);", parser=1)
        rooted = mv_rooting(tree)
        assert set(rooted.leaf_names()) == {"A", "B", "C", "D", "E"}

    def test_clock_like_tree(self):
        """On a clock-like tree, MV root should achieve near-zero variance."""
        import numpy as np

        tree = Tree("((A:2,B:2):1,(C:2,D:2):1);", parser=1)
        rooted = mv_rooting(tree)
        dists = [safe_get_distance(rooted, rooted, leaf) for leaf in rooted.leaves()]
        assert np.var(dists) < 1e-10

    def test_asymmetric_branch_lengths(self):
        tree = Tree("((A:10,B:1):1,(C:1,D:1):1);", parser=1)
        rooted = mv_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {"A", "B", "C", "D"}

    def test_wiki_exact_root_position(self):
        """MV rooting on wiki tree: root at x=5/12 from (C,(D,E)) side.

        Input: ((A:1,B:2):1,(C:3,(D:1,E:2):1):1);
        Expected root-to-tip distances:
          C: 3 + 5/12 = 41/12, D: 1+1+5/12 = 29/12, E: 2+1+5/12 = 41/12,
          A: 1 + 19/12 = 31/12, B: 2 + 19/12 = 43/12
        """
        tree = Tree("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);", parser=1)
        rooted = mv_rooting(tree)
        dists = {l.name: safe_get_distance(rooted, rooted, l) for l in rooted.leaves()}
        assert abs(dists["C"] - 41 / 12) < 1e-4
        assert abs(dists["D"] - 29 / 12) < 1e-4
        assert abs(dists["E"] - 41 / 12) < 1e-4
        assert abs(dists["A"] - 31 / 12) < 1e-4
        assert abs(dists["B"] - 43 / 12) < 1e-4
        # Root children: (C,(D,E)) side=5/12, (A,B) side=19/12
        children = rooted.get_children()
        for c in children:
            if "C" in c.leaf_names():
                assert abs(c.dist - 5 / 12) < 1e-4
            else:
                assert abs(c.dist - 19 / 12) < 1e-4

    def test_pairwise_distance_preservation(self):
        """MV rooting must preserve all pairwise distances."""
        nwk = "((A:1,B:2):1,(C:3,(D:1,E:2):1):1);"
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = mv_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert (
                        abs(
                            original.get_distance(l1.name, l2.name)
                            - rooted.get_distance(l1, l2)
                        )
                        < 1e-4
                    )

    def test_three_taxa_exact(self):
        """MV rooting on 3-taxon tree: verify pairwise distances preserved."""
        tree = Tree("(A:1,B:3,C:5);", parser=1)
        rooted = mv_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    d = rooted.get_distance(l1, l2)
                    if {l1.name, l2.name} == {"A", "B"}:
                        assert abs(d - 4) < 1e-6
                    elif {l1.name, l2.name} == {"A", "C"}:
                        assert abs(d - 6) < 1e-6
                    elif {l1.name, l2.name} == {"B", "C"}:
                        assert abs(d - 8) < 1e-6

    def test_nonzero_root_distance_does_not_crash(self):
        tree = Tree("((A:1,B:2):1,C:3):2;", parser=1)
        rooted = mv_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {"A", "B", "C"}

    def test_squared_length_overflow_scale_invariance(self):
        tree = Tree(
            "((A:1,B:2):1,(C:3,(D:1,E:2):1):1);",
            parser=1,
        )
        huge_tree = tree.copy(method="deepcopy")
        scale = 10**154
        for node in huge_tree.traverse():
            if not node.is_root:
                node.dist *= scale

        rooted = mv_rooting(tree)
        huge_rooted = mv_rooting(huge_tree)

        assert is_rooted(huge_rooted)
        expected_root_sides = {
            frozenset(child.leaf_names()) for child in rooted.get_children()
        }
        observed_root_sides = {
            frozenset(child.leaf_names()) for child in huge_rooted.get_children()
        }
        assert observed_root_sides == expected_root_sides
        expected_profile = _scaled_branch_profile(rooted)
        observed_profile = _scaled_branch_profile(huge_rooted, scale=scale)
        assert observed_profile.keys() == expected_profile.keys()
        for split in expected_profile:
            assert observed_profile[split] == pytest.approx(
                expected_profile[split],
            )
        assert all(
            math.isfinite(node.dist)
            for node in huge_rooted.traverse()
            if not node.is_root
        )


class TestRerootAnnotationSafety:
    @pytest.mark.parametrize(
        "rooter",
        [
            midpoint_rooting,
            mv_rooting,
            lambda tree: outgroup_rooting(tree, "A"),
        ],
    )
    def test_split_key_rooters_reject_duplicate_leaf_names(self, rooter):
        tree = Tree("(A:1,A:2,(B:3,C:4):5);", parser=1)

        with pytest.raises(ValueError, match="Duplicated leaf labels"):
            rooter(tree)

    @pytest.mark.parametrize(
        ("method_name", "rooter"),
        [
            ("Midpoint", midpoint_rooting),
            ("MV", mv_rooting),
        ],
    )
    @pytest.mark.parametrize(
        "bad_length",
        [-1.0, float("inf"), float("nan")],
        ids=("negative", "infinite", "nan"),
    )
    def test_numeric_rooters_reject_invalid_nonroot_lengths(
        self, method_name, rooter, bad_length
    ):
        tree = Tree("(A:1,B:2,(C:3,D:4):5);", parser=1)
        next(leaf for leaf in tree.leaves() if leaf.name == "A").dist = bad_length

        with pytest.raises(
            ValueError,
            match="{} rooting requires finite, non-negative branch lengths".format(
                method_name
            ),
        ):
            rooter(tree)

    @pytest.mark.parametrize("rooter", [midpoint_rooting, mv_rooting])
    def test_numeric_rooters_preserve_root_stem(self, rooter):
        tree = Tree("((A:1,B:2):3,(C:4,D:5):6);", parser=1)
        tree.dist = 7.5 * (10**250)

        rooted = rooter(tree)

        assert rooted.dist == tree.dist

    @pytest.mark.parametrize(
        ("method_name", "rooter"),
        [
            ("midpoint", midpoint_rooting),
            ("outgroup", lambda tree: outgroup_rooting(tree, "A")),
            ("mad", mad_rooting),
            ("mv", mv_rooting),
            (
                "taxonomy-edge",
                lambda tree: root_mod._root_by_outgroup_set(tree, {"E", "F"}),
            ),
        ],
    )
    def test_preserves_branch_root_and_tip_annotations(self, method_name, rooter):
        tree, expected_by_split = _annotated_reroot_tree()
        original_root_children = [
            set(child.leaf_names()) for child in tree.get_children()
        ]

        rooted = rooter(tree)

        assert rooted is not tree
        _assert_reroot_annotations(rooted, expected_by_split)
        _assert_reroot_annotations(tree, expected_by_split)
        assert [
            set(child.leaf_names()) for child in tree.get_children()
        ] == original_root_children

    def test_singleton_root_collapse_preserves_child_metadata(self):
        tree = Tree("(((A:1,B:1):1,C:1)INNER:2):3;", parser=1)
        tree.get_children()[0].props["root_tag"] = "kept"
        source = Tree("((A:1,B:1):1,C:1);", parser=1)

        rooted = transfer_root(tree, source)

        assert rooted.name == "INNER"
        assert rooted.props["root_tag"] == "kept"
        assert rooted.dist == pytest.approx(5.0)

    def test_singleton_root_collapse_rejects_nonfinite_combined_length(self):
        tree = Tree(
            "(((A:1,B:1):1,C:1):1e308):1e308;",
            parser=1,
        )
        source = Tree("((A:1,B:1):1,C:1);", parser=1)

        with pytest.raises(ValueError, match="non-finite root branch length"):
            transfer_root(tree, source)

    @pytest.mark.parametrize(
        "rooter",
        [
            midpoint_rooting,
            lambda tree: outgroup_rooting(tree, "A"),
            mv_rooting,
        ],
    )
    def test_rerooting_preserves_missing_branch_lengths(self, rooter):
        tree = Tree("(A,(B,(C,D)));", parser=1)
        rooted = rooter(tree)
        assert all(node.dist is None for node in rooted.traverse() if not node.is_root)

    def test_unchanged_root_preserves_one_sided_missing_root_length(self):
        tree = Tree("((A:1,B:1),(C:1,D:1):0);", parser=1)

        rooted = midpoint_rooting(tree)

        child_by_taxa = {
            frozenset(child.leaf_names()): child for child in rooted.get_children()
        }
        assert child_by_taxa[frozenset({"A", "B"})].dist is None
        assert child_by_taxa[frozenset({"C", "D"})].dist == pytest.approx(0.0)


class TestTransferRoot:
    def test_transfer(self):
        # tree_from is rooted with (A,B) | (C,D)
        tree_from = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)
        # tree_to is unrooted
        tree_to = Tree("(A:1,B:1,(C:1,D:1):1);", parser=1)
        result = transfer_root(tree_to, tree_from)
        assert is_rooted(result)
        assert set(result.leaf_names()) == {"A", "B", "C", "D"}

    def test_pairwise_distance_preservation(self):
        """Transfer root must preserve pairwise distances from tree_to."""
        tree_from = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)
        tree_to = Tree("(A:2,B:3,(C:4,D:5):1);", parser=1)
        original = Tree("(A:2,B:3,(C:4,D:5):1);", parser=1)
        result = transfer_root(tree_to, tree_from)
        for l1 in result.leaves():
            for l2 in result.leaves():
                if l1.name != l2.name:
                    orig_d = original.get_distance(l1.name, l2.name)
                    new_d = result.get_distance(l1, l2)
                    assert abs(orig_d - new_d) < 1e-6, (
                        f"{l1.name}-{l2.name}: {orig_d} vs {new_d}"
                    )

    def test_preserves_internal_support_on_unrooted_splits(self):
        tree_from = Tree(
            "((A:1,B:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);",
            parser=1,
        )
        tree_to = Tree(
            "((E:1,F:1)40:0.5,(D:1,(C:1,(A:1,B:1)20:1)30:1)40:0.5);",
            parser=0,
        )

        result = transfer_root(tree_to, tree_from)
        all_taxa = frozenset(result.leaf_names())
        support_by_split = dict()
        for node in result.traverse():
            if node.is_root or node.is_leaf:
                continue
            side = frozenset(node.leaf_names())
            split = canonical_split(side, all_taxa - side)
            support_by_split.setdefault(split, set()).add(float(node.support))

        assert support_by_split[
            canonical_split(frozenset({"A", "B"}), frozenset({"C", "D", "E", "F"}))
        ] == {20.0}
        assert support_by_split[
            canonical_split(frozenset({"A", "B", "C"}), frozenset({"D", "E", "F"}))
        ] == {30.0}
        assert support_by_split[
            canonical_split(frozenset({"E", "F"}), frozenset({"A", "B", "C", "D"}))
        ] == {40.0}

    def test_zero_subroot_length_in_source_does_not_crash(self):
        tree_from = Tree("((A:0,B:0):0,(C:0,D:0):0);", parser=1)
        tree_to = Tree("((A:3,B:3):3,(C:1,D:1):1);", parser=1)
        result = transfer_root(tree_to, tree_from)
        subroot_dists = sorted([child.dist for child in result.get_children()])
        assert subroot_dists == [1.0, 3.0]

    def test_nonzero_root_distance_in_target_does_not_crash(self):
        tree_from = Tree("((A:1,B:1):1,(C:1,D:1):1):0;", parser=1)
        tree_to = Tree("(A:1,B:1,(C:1,D:1):1):2;", parser=1)
        result = transfer_root(tree_to, tree_from)
        assert is_rooted(result)
        assert set(result.leaf_names()) == {"A", "B", "C", "D"}

    def test_source_subroot_none_distance_does_not_redistribute(self):
        tree_from = Tree("((A:1,B:1),(C:1,D:1):1);", parser=1)
        tree_to = Tree("((A:1,B:1):2,(C:1,D:1):2);", parser=1)
        result = transfer_root(tree_to, tree_from)
        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist for child in result.get_children()
        }
        assert child_dist_by_taxa == {
            frozenset({"A", "B"}): 2.0,
            frozenset({"C", "D"}): 2.0,
        }

    @pytest.mark.parametrize("invalid_length", [float("nan"), float("inf"), -1.0])
    def test_invalid_source_root_length_does_not_redistribute(self, invalid_length):
        tree_from = Tree("((A:1,B:1):1,(C:1,D:1):3);", parser=1)
        tree_from.get_children()[0].dist = invalid_length
        tree_to = Tree("((A:1,B:1):2,(C:1,D:1):4);", parser=1)

        result = transfer_root(tree_to, tree_from)

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist for child in result.get_children()
        }
        assert child_dist_by_taxa == {
            frozenset({"A", "B"}): 2.0,
            frozenset({"C", "D"}): 4.0,
        }

    def test_singleton_root_in_target_does_not_create_unnamed_leaf(self):
        tree_from = Tree("((A:1,B:1):1,C:1);", parser=1)
        tree_to = Tree("(((A:1,B:1):1,C:1):1);", parser=1)
        result = transfer_root(tree_to, tree_from)
        leaf_names = list(result.leaf_names())
        assert set(leaf_names) == {"A", "B", "C"}
        assert None not in leaf_names

    def test_exact_transfer_rejects_different_leaf_sets(self):
        source = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)
        target = Tree("((A:1,B:1):1,C:1,D:1,X:1);", parser=1)
        with pytest.raises(ValueError, match="identical leaf labels"):
            transfer_root(target, source)

    def test_intersection_transfer_honors_root_length_redistribution(self):
        target = Tree("((A:1,B:1):2,C:1,D:1,X:1);", parser=1)
        source = Tree("((A:1,B:1):1,(C:1,D:1):3);", parser=1)
        redistributed = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode="intersection",
            redistribute_root_length=True,
        )
        unchanged = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode="intersection",
            redistribute_root_length=False,
        )
        assert sorted(
            child.dist for child in redistributed.get_children()
        ) == pytest.approx([0.5, 1.5])
        assert sorted(
            child.dist for child in unchanged.get_children()
        ) == pytest.approx([1.0, 1.0])

    def test_matching_intersection_root_still_honors_length_redistribution(self):
        target = Tree(
            "((A:1,B:1):2,(C:1,D:1,X:1):2);",
            parser=1,
        )
        source = Tree("((A:1,B:1):1,(C:1,D:1):3);", parser=1)

        redistributed = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode="intersection",
            redistribute_root_length=True,
        )
        unchanged = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode="intersection",
            redistribute_root_length=False,
        )

        redistributed_by_taxa = {
            frozenset(child.leaf_names()): child.dist
            for child in redistributed.get_children()
        }
        assert redistributed_by_taxa[frozenset({"A", "B"})] == pytest.approx(1.0)
        assert redistributed_by_taxa[frozenset({"C", "D", "X"})] == pytest.approx(3.0)
        assert sorted(child.dist for child in unchanged.get_children()) == [2.0, 2.0]

    def test_intersection_transfer_skips_missing_source_root_length(self):
        target = Tree(
            "((A:1,B:1):2,(C:1,D:1,X:1):4);",
            parser=1,
        )
        source = Tree("((A:1,B:1),(C:1,D:1):3);", parser=1)

        rooted = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode="intersection",
            redistribute_root_length=True,
        )

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist for child in rooted.get_children()
        }
        assert child_dist_by_taxa == {
            frozenset({"A", "B"}): 2.0,
            frozenset({"C", "D", "X"}): 4.0,
        }

    def test_matching_exact_root_does_not_fill_missing_target_root_length(self):
        target = Tree("((A:1,B:1),(C:1,D:1):4);", parser=1)
        source = Tree("((A:1,B:1):1,(C:1,D:1):3);", parser=1)

        rooted = transfer_root(target, source)

        child_by_taxa = {
            frozenset(child.leaf_names()): child for child in rooted.get_children()
        }
        assert child_by_taxa[frozenset({"A", "B"})].dist is None
        assert child_by_taxa[frozenset({"C", "D"})].dist == pytest.approx(4.0)
