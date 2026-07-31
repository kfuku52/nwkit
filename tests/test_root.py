import json as json_lib
import math
import time

import pandas as pd
import pytest
from ete4 import Tree

import nwkit.constrain as constrain_mod
import nwkit.root as root_mod
from nwkit.root import (
    DEFAULT_TAXONOMY_SOURCE_CHAIN,
    midpoint_rooting,
    outgroup_rooting,
    transfer_root,
    mad_rooting,
    mv_rooting,
    taxonomy_rooting,
    root_main,
)
from nwkit.clade_mapping import canonical_split
from nwkit.util import read_tree, is_rooted
from tests.helpers import make_args, safe_get_distance


def _annotated_reroot_tree():
    tree = Tree(
        '(((A:2,B:3):5,C:7):11,(D:13,(E:17,F:19):23):29);',
        parser=1,
    )
    all_taxa = frozenset(tree.leaf_names())
    expected_by_split = dict()
    for node in tree.traverse():
        if node.is_root:
            node.name = 'ORIGINAL_ROOT'
            node.support = 0.99
            node.props['root_tag'] = 'root_value'
        elif node.is_leaf:
            node.support = 0.5
            node.props['tip_tag'] = 'tip_{}'.format(node.name)
        else:
            side = frozenset(node.leaf_names())
            split = canonical_split(side, all_taxa - side)
            if split not in expected_by_split:
                index = len(expected_by_split) + 1
                expected_by_split[split] = (
                    'EDGE_{}'.format(index),
                    float(10 + index),
                    'edge_value_{}'.format(index),
                )
            node.name, node.support, node.props['edge_tag'] = expected_by_split[split]
    return tree, expected_by_split


def _assert_reroot_annotations(tree, expected_by_split):
    all_taxa = frozenset(tree.leaf_names())
    observed_by_split = dict()
    for node in tree.traverse():
        if node.is_root or node.is_leaf:
            continue
        side = frozenset(node.leaf_names())
        split = canonical_split(side, all_taxa - side)
        observed_by_split.setdefault(split, set()).add(
            (node.name, node.support, node.props.get('edge_tag'))
        )
    for split, expected in expected_by_split.items():
        assert observed_by_split[split] == {expected}
    assert tree.name == 'ORIGINAL_ROOT'
    assert tree.support == 0.99
    assert tree.props['root_tag'] == 'root_value'
    for leaf in tree.leaves():
        assert leaf.support == 0.5
        assert leaf.props['tip_tag'] == 'tip_{}'.format(leaf.name)


def _scaled_branch_profile(tree, scale=1.0):
    all_taxa = frozenset(tree.leaf_names())
    profile = dict()
    for node in tree.traverse():
        if node.is_root:
            continue
        side = frozenset(node.leaf_names())
        split = canonical_split(side, all_taxa - side)
        profile.setdefault(split, list()).append(float(node.dist) / scale)
    return {
        split: sorted(lengths)
        for split, lengths in profile.items()
    }


def install_fake_ncbi(monkeypatch, name_to_taxid, lineage_by_taxid, taxid_to_name=None, rank_by_taxid=None):
    if taxid_to_name is None:
        taxid_to_name = {int(taxid): str(taxid) for taxid in lineage_by_taxid.keys()}
    if rank_by_taxid is None:
        rank_by_taxid = dict()

    class FakeNCBI:
        def __init__(self, *args, **kwargs):
            self.db = None

        def get_name_translator(self, names):
            translated = dict()
            for name in names:
                if name in name_to_taxid:
                    translated[name] = [name_to_taxid[name]]
            return translated

        def get_lineage(self, taxid):
            return list(lineage_by_taxid[int(taxid)])

        def get_taxid_translator(self, taxids):
            return {
                int(taxid): taxid_to_name[int(taxid)]
                for taxid in taxids
                if int(taxid) in taxid_to_name
            }

        def get_rank(self, taxids):
            return {
                int(taxid): rank_by_taxid.get(int(taxid), 'no rank')
                for taxid in taxids
            }

    monkeypatch.setattr(constrain_mod.ete4, 'NCBITaxa', FakeNCBI)

def install_fake_timetree(monkeypatch, upload_html, newick_text, upload_status=200, newick_status=200):
    calls = list()

    class FakeResponse:
        def __init__(self, text, status_code):
            self.text = text
            self.status_code = status_code
            self.headers = {}
            self.encoding = 'utf-8'

        def iter_content(self, chunk_size):
            payload = self.text.encode(self.encoding)
            for offset in range(0, len(payload), chunk_size):
                yield payload[offset:offset + chunk_size]

        def close(self):
            pass

    class FakeSession:
        def post(self, url, files=None, data=None, timeout=None, stream=False):
            calls.append((url, stream))
            if url.endswith('/ajax/prune/load_names/'):
                return FakeResponse(upload_html, upload_status)
            if url.endswith('/ajax/newick/prunetree/download'):
                return FakeResponse(newick_text, newick_status)
            raise AssertionError('Unexpected TimeTree URL: {}'.format(url))

    monkeypatch.setattr(root_mod.requests, 'Session', FakeSession)
    return calls

def install_fake_opentree(monkeypatch, tnrs_json, induced_subtree_json, tnrs_status=200, induced_subtree_status=200):
    calls = list()

    class FakeResponse:
        def __init__(self, text, status_code, json_data=None):
            self.text = text
            self.status_code = status_code
            self._json_data = json_data
            self.headers = {}
            self.encoding = 'utf-8'

        def json(self):
            if self._json_data is None:
                raise AssertionError('JSON payload was not configured')
            return self._json_data

        def iter_content(self, chunk_size):
            payload = self.text.encode(self.encoding)
            for offset in range(0, len(payload), chunk_size):
                yield payload[offset:offset + chunk_size]

        def close(self):
            pass

    class FakeSession:
        def post(self, url, json=None, timeout=None, stream=False):
            calls.append((url, stream))
            if url.endswith('/v3/tnrs/match_names'):
                return FakeResponse(
                    json_lib.dumps(tnrs_json),
                    tnrs_status,
                    json_data=tnrs_json,
                )
            if url.endswith('/v3/tree_of_life/induced_subtree'):
                return FakeResponse(
                    json_lib.dumps(induced_subtree_json),
                    induced_subtree_status,
                    json_data=induced_subtree_json,
                )
            raise AssertionError('Unexpected OpenTree URL: {}'.format(url))

    monkeypatch.setattr(root_mod.requests, 'Session', FakeSession)
    return calls


class TestMidpointRooting:
    def test_basic(self):
        tree = Tree('(A:1,B:3,(C:2,D:4):2);', parser=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)

    def test_already_rooted(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_symmetric_clock_exact_distances(self):
        """On a symmetric clock tree, all root-to-tip distances are equal."""
        tree = Tree('((A:2,B:2):1,(C:2,D:2):1);', parser=1)
        rooted = midpoint_rooting(tree)
        dists = [safe_get_distance(rooted, rooted, l) for l in rooted.leaves()]
        assert all(abs(d - 3.0) < 1e-6 for d in dists)

    def test_pairwise_distance_preservation(self):
        """Midpoint rooting must preserve all pairwise distances."""
        nwk = '((A:3,B:1):2,(C:4,D:2):1);'
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = midpoint_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert abs(original.get_distance(l1.name, l2.name) -
                               rooted.get_distance(l1, l2)) < 1e-6

    def test_nonzero_root_distance_does_not_crash(self):
        tree = Tree('((A:1,B:3):1,(C:1,D:1):1):2;', parser=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_path_overflow_scale_invariance(self):
        subtree = '(T0:1,T1:1):1'
        for index in range(2, 23):
            subtree = '({},T{}:1):1'.format(subtree, index)
        tree = Tree('({},T23:1);'.format(subtree), parser=1)
        huge_tree = tree.copy(method='deepcopy')
        scale = 10 ** 307
        for node in huge_tree.traverse():
            if not node.is_root:
                node.dist *= scale

        rooted = midpoint_rooting(tree)
        huge_rooted = midpoint_rooting(huge_tree)

        expected_root_sides = {
            frozenset(child.leaf_names())
            for child in rooted.get_children()
        }
        observed_root_sides = {
            frozenset(child.leaf_names())
            for child in huge_rooted.get_children()
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
        tree = Tree('(A:1,(B:1,(C:1,D:1):1):1);', parser=1)
        rooted = outgroup_rooting(tree, 'A')
        assert is_rooted(rooted)
        # A should be one of the subroot children
        subroot_children = rooted.get_children()
        outgroup_leaves = set()
        for child in subroot_children:
            outgroup_leaves.update(child.leaf_names())
        assert 'A' in outgroup_leaves

    def test_multiple_outgroups(self):
        tree = Tree('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);', parser=1)
        rooted = outgroup_rooting(tree, 'D,E')
        assert is_rooted(rooted)
        leaf_names = set(rooted.leaf_names())
        assert leaf_names == {'A', 'B', 'C', 'D', 'E'}

    def test_multiple_outgroups_with_spaces(self):
        tree = Tree('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);', parser=1)
        rooted = outgroup_rooting(tree, 'D, E')
        children = rooted.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {'D', 'E'} in child_leaf_sets

    def test_outgroup_not_found(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        with pytest.raises(ValueError, match='Outgroup label'):
            outgroup_rooting(tree, 'Z')

    def test_outgroup_partially_missing_raises(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        with pytest.raises(ValueError, match='Outgroup label'):
            outgroup_rooting(tree, 'A,Z')

    def test_outgroup_whole_tree_raises(self):
        tree = Tree('(A:1,B:1);', parser=1)
        with pytest.raises(ValueError, match='whole tree'):
            outgroup_rooting(tree, 'A,B')

    def test_outgroup_single_leaf_tree_raises(self):
        tree = Tree('A;', parser=1)
        with pytest.raises(ValueError, match='whole tree'):
            outgroup_rooting(tree, 'A')

    def test_nonzero_root_distance_does_not_crash(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1):2;', parser=1)
        rooted = outgroup_rooting(tree, 'A')
        children = rooted.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {'A'} in child_leaf_sets

    def test_pairwise_distance_preservation(self):
        """Outgroup rooting must preserve all pairwise distances."""
        nwk = '((A:2,B:3):1,(C:4,D:5):1);'
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = outgroup_rooting(tree, 'A')
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert abs(original.get_distance(l1.name, l2.name) -
                               rooted.get_distance(l1, l2)) < 1e-6

    def test_multiple_outgroup_exact_bipartition(self):
        """Outgroup rooting with multiple tips creates correct bipartition."""
        tree = Tree('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);', parser=1)
        rooted = outgroup_rooting(tree, 'D,E')
        children = rooted.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {'D', 'E'} in child_leaf_sets
        assert {'A', 'B', 'C'} in child_leaf_sets


class TestMadRooting:
    def test_basic(self):
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)

    def test_preserves_leaves(self):
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mad_rooting(tree)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_asymmetric_branch_lengths(self):
        tree = Tree('((A:10,B:1):1,(C:1,D:1):1);', parser=1)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_preserves_quoted_leaf_names_with_punctuation(self):
        tree = Tree()
        ingroup = tree.add_child(dist=1.0)
        ingroup.add_child(name='A,B', dist=1.0)
        ingroup.add_child(name='C:D', dist=2.0)
        tree.add_child(name="E(1)'", dist=3.0)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A,B', 'C:D', "E(1)'"}

    def test_requires_at_least_three_leaves(self):
        tree = Tree('(A:1,B:1);', parser=1)
        with pytest.raises(ValueError, match='at least 3 leaves'):
            mad_rooting(tree)

    def test_requires_at_least_three_positive_branches(self):
        tree = Tree('(A:0,B:1,C:1,D:0);', parser=1)
        with pytest.raises(ValueError, match='at least 3 positive branch lengths'):
            mad_rooting(tree)

    def test_requires_every_branch_length(self):
        tree = Tree('(A,B:1,C:2,D:3);', parser=1)

        with pytest.raises(ValueError, match='every branch length'):
            mad_rooting(tree)

    @pytest.mark.parametrize(
        'bad_length',
        [-1.0, float('inf'), float('nan')],
        ids=('negative', 'infinite', 'nan'),
    )
    def test_rejects_invalid_branch_lengths(self, bad_length):
        tree = Tree('(A:1,B:2,C:3,D:4);', parser=1)
        next(leaf for leaf in tree.leaves() if leaf.name == 'A').dist = bad_length

        with pytest.raises(ValueError, match='finite, non-negative'):
            mad_rooting(tree)

    def test_zero_distance_tips_are_retained_as_one_effective_tip(self):
        tree = Tree('(E:3,D:2,C:1,B:0,A:0);', parser=1)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_tied_root_edges_use_a_deterministic_split(self):
        expected_split = canonical_split(
            frozenset({'A'}),
            frozenset({'B', 'C', 'D'}),
        )
        observed_splits = list()
        for newick in (
            '(D:1,C:1,B:1,A:1);',
            '(B:1,D:1,A:1,C:1);',
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
            '((A:1,B:2):1,(C:3,(D:1,E:2):1):1);',
            parser=1,
        )
        observations = list()
        for scale in (10 ** -160, 1.0, 10 ** 160):
            scaled_tree = base_tree.copy(method='deepcopy')
            for node in scaled_tree.traverse():
                if not node.is_root:
                    node.dist = float(node.dist) * scale
            rooted = mad_rooting(scaled_tree)
            child_dist_by_taxa = {
                frozenset(child.leaf_names()): child.dist
                for child in rooted.get_children()
            }
            root_total = sum(child_dist_by_taxa.values())
            observations.append((
                frozenset(child_dist_by_taxa),
                child_dist_by_taxa[frozenset({'A', 'B'})] / root_total,
            ))

        assert all(
            split == frozenset((
                frozenset({'A', 'B'}),
                frozenset({'C', 'D', 'E'}),
            ))
            for split, _ in observations
        )
        assert [position for _, position in observations] == pytest.approx(
            [observations[1][1]] * 3,
            rel=10 ** -12,
            abs=10 ** -12,
        )

    def test_extreme_within_tree_length_range_does_not_overflow(self):
        rooted = mad_rooting(
            Tree(
                '(A:1e-160,B:1e-160,(C:1,D:1):1);',
                parser=1,
            )
        )

        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D'}
        assert all(
            math.isfinite(float(node.dist))
            for node in rooted.traverse()
            if not node.is_root
        )

    def test_heterogeneous_lengths_do_not_corrupt_path_delta_cancellation(self):
        tree = Tree(
            '(X1:1086.5437592639944,'
            '(X3:0.00085007102927072819,'
            '(X2:1.0785237796461825e-06,'
            'X0:0.00057934556794108716):75531.807545509451)'
            ':0.00763325784066344);',
            parser=1,
        )

        rooted = mad_rooting(tree)

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist
            for child in rooted.get_children()
        }
        assert set(child_dist_by_taxa) == {
            frozenset({'X0', 'X2'}),
            frozenset({'X1', 'X3'}),
        }
        assert child_dist_by_taxa[frozenset({'X0', 'X2'})] == pytest.approx(
            38033.662192492,
        )

    def test_deep_short_clade_does_not_lose_root_distance_precision(self):
        tree = Tree(
            '(X4:0.0063973377995121773,'
            'X3:9.9523327390797709e-06,'
            '(X0:87377874218.288742,'
            'X1:3.9277797271534321e-10):3.7230847581009794e-12,'
            'X2:1540148.2696119624);',
            parser=1,
        )

        rooted = mad_rooting(tree)

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist
            for child in rooted.get_children()
        }
        assert frozenset({'X0'}) in child_dist_by_taxa
        assert child_dist_by_taxa[frozenset({'X0'})] == pytest.approx(
            43689129622.58888,
        )

    def test_finite_root_halves_are_joined_without_intermediate_overflow(self):
        maximum_float = float.fromhex('0x1.fffffffffffffp+1023')
        tree = Tree(
            '((A:1,B:2):{0},(C:3,D:4):{0});'.format(maximum_float),
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
            '(A:5e-324,B:1,C:2,D:1e308);',
            parser=1,
        )

        with pytest.raises(ValueError, match='branch-length dynamic range'):
            mad_rooting(tree)

    def test_scoring_does_not_walk_tree_paths_for_every_pair(self):
        tree = Tree(
            '((A:1,B:2):1,(C:3,(D:1,E:2):1):1);',
            parser=1,
        )

        assert 'get_distance' not in mad_rooting.__code__.co_names
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)

    def test_balanced_320_tip_tree_completes_within_generous_budget(self):
        def balanced_newick(labels, depth=0):
            if len(labels) == 1:
                index = int(labels[0][1:])
                length = 0.7 + ((index * 17 + depth * 3) % 19) / 10
                return '{}:{}'.format(labels[0], length)
            midpoint = len(labels) // 2
            length = 0.5 + ((len(labels) * 11 + depth) % 13) / 10
            return '({},{}):{}'.format(
                balanced_newick(labels[:midpoint], depth + 1),
                balanced_newick(labels[midpoint:], depth + 1),
                length,
            )

        labels = ['T{}'.format(index) for index in range(320)]
        tree = Tree(balanced_newick(labels) + ';', parser=1)

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
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        original = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mad_rooting(tree)
        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist
            for child in rooted.get_children()
        }
        assert set(child_dist_by_taxa) == {
            frozenset({'A', 'B'}),
            frozenset({'C', 'D', 'E'}),
        }
        assert child_dist_by_taxa[frozenset({'A', 'B'})] == pytest.approx(
            1.584611134,
            abs=10 ** -8,
        )
        assert child_dist_by_taxa[frozenset({'C', 'D', 'E'})] == pytest.approx(
            0.415388866,
            abs=10 ** -8,
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
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mv_rooting(tree)
        assert is_rooted(rooted)

    def test_preserves_leaves(self):
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mv_rooting(tree)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_clock_like_tree(self):
        """On a clock-like tree, MV root should achieve near-zero variance."""
        import numpy as np
        tree = Tree('((A:2,B:2):1,(C:2,D:2):1);', parser=1)
        rooted = mv_rooting(tree)
        dists = [safe_get_distance(rooted, rooted, leaf) for leaf in rooted.leaves()]
        assert np.var(dists) < 1e-10

    def test_asymmetric_branch_lengths(self):
        tree = Tree('((A:10,B:1):1,(C:1,D:1):1);', parser=1)
        rooted = mv_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_wiki_exact_root_position(self):
        """MV rooting on wiki tree: root at x=5/12 from (C,(D,E)) side.

        Input: ((A:1,B:2):1,(C:3,(D:1,E:2):1):1);
        Expected root-to-tip distances:
          C: 3 + 5/12 = 41/12, D: 1+1+5/12 = 29/12, E: 2+1+5/12 = 41/12,
          A: 1 + 19/12 = 31/12, B: 2 + 19/12 = 43/12
        """
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mv_rooting(tree)
        dists = {l.name: safe_get_distance(rooted, rooted, l) for l in rooted.leaves()}
        assert abs(dists['C'] - 41/12) < 1e-4
        assert abs(dists['D'] - 29/12) < 1e-4
        assert abs(dists['E'] - 41/12) < 1e-4
        assert abs(dists['A'] - 31/12) < 1e-4
        assert abs(dists['B'] - 43/12) < 1e-4
        # Root children: (C,(D,E)) side=5/12, (A,B) side=19/12
        children = rooted.get_children()
        for c in children:
            if 'C' in c.leaf_names():
                assert abs(c.dist - 5/12) < 1e-4
            else:
                assert abs(c.dist - 19/12) < 1e-4

    def test_pairwise_distance_preservation(self):
        """MV rooting must preserve all pairwise distances."""
        nwk = '((A:1,B:2):1,(C:3,(D:1,E:2):1):1);'
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = mv_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert abs(original.get_distance(l1.name, l2.name) -
                               rooted.get_distance(l1, l2)) < 1e-4

    def test_three_taxa_exact(self):
        """MV rooting on 3-taxon tree: verify pairwise distances preserved."""
        tree = Tree('(A:1,B:3,C:5);', parser=1)
        rooted = mv_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    d = rooted.get_distance(l1, l2)
                    if {l1.name, l2.name} == {'A', 'B'}:
                        assert abs(d - 4) < 1e-6
                    elif {l1.name, l2.name} == {'A', 'C'}:
                        assert abs(d - 6) < 1e-6
                    elif {l1.name, l2.name} == {'B', 'C'}:
                        assert abs(d - 8) < 1e-6

    def test_nonzero_root_distance_does_not_crash(self):
        tree = Tree('((A:1,B:2):1,C:3):2;', parser=1)
        rooted = mv_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C'}

    def test_squared_length_overflow_scale_invariance(self):
        tree = Tree(
            '((A:1,B:2):1,(C:3,(D:1,E:2):1):1);',
            parser=1,
        )
        huge_tree = tree.copy(method='deepcopy')
        scale = 10 ** 154
        for node in huge_tree.traverse():
            if not node.is_root:
                node.dist *= scale

        rooted = mv_rooting(tree)
        huge_rooted = mv_rooting(huge_tree)

        assert is_rooted(huge_rooted)
        expected_root_sides = {
            frozenset(child.leaf_names())
            for child in rooted.get_children()
        }
        observed_root_sides = {
            frozenset(child.leaf_names())
            for child in huge_rooted.get_children()
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
        'rooter',
        [
            midpoint_rooting,
            mv_rooting,
            lambda tree: outgroup_rooting(tree, 'A'),
        ],
    )
    def test_split_key_rooters_reject_duplicate_leaf_names(self, rooter):
        tree = Tree('(A:1,A:2,(B:3,C:4):5);', parser=1)

        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            rooter(tree)

    @pytest.mark.parametrize(
        ('method_name', 'rooter'),
        [
            ('Midpoint', midpoint_rooting),
            ('MV', mv_rooting),
        ],
    )
    @pytest.mark.parametrize(
        'bad_length',
        [-1.0, float('inf'), float('nan')],
        ids=('negative', 'infinite', 'nan'),
    )
    def test_numeric_rooters_reject_invalid_nonroot_lengths(
            self, method_name, rooter, bad_length):
        tree = Tree('(A:1,B:2,(C:3,D:4):5);', parser=1)
        next(leaf for leaf in tree.leaves() if leaf.name == 'A').dist = bad_length

        with pytest.raises(
            ValueError,
            match='{} rooting requires finite, non-negative branch lengths'.format(
                method_name
            ),
        ):
            rooter(tree)

    @pytest.mark.parametrize('rooter', [midpoint_rooting, mv_rooting])
    def test_numeric_rooters_preserve_root_stem(self, rooter):
        tree = Tree('((A:1,B:2):3,(C:4,D:5):6);', parser=1)
        tree.dist = 7.5 * (10 ** 250)

        rooted = rooter(tree)

        assert rooted.dist == tree.dist

    @pytest.mark.parametrize(
        ('method_name', 'rooter'),
        [
            ('midpoint', midpoint_rooting),
            ('outgroup', lambda tree: outgroup_rooting(tree, 'A')),
            ('mad', mad_rooting),
            ('mv', mv_rooting),
            ('taxonomy-edge', lambda tree: root_mod._root_by_outgroup_set(tree, {'E', 'F'})),
        ],
    )
    def test_preserves_branch_root_and_tip_annotations(self, method_name, rooter):
        tree, expected_by_split = _annotated_reroot_tree()
        original_root_children = [set(child.leaf_names()) for child in tree.get_children()]

        rooted = rooter(tree)

        assert rooted is not tree
        _assert_reroot_annotations(rooted, expected_by_split)
        _assert_reroot_annotations(tree, expected_by_split)
        assert [set(child.leaf_names()) for child in tree.get_children()] == original_root_children

    def test_singleton_root_collapse_preserves_child_metadata(self):
        tree = Tree('(((A:1,B:1):1,C:1)INNER:2):3;', parser=1)
        tree.get_children()[0].props['root_tag'] = 'kept'
        source = Tree('((A:1,B:1):1,C:1);', parser=1)

        rooted = transfer_root(tree, source)

        assert rooted.name == 'INNER'
        assert rooted.props['root_tag'] == 'kept'
        assert rooted.dist == pytest.approx(5.0)

    def test_singleton_root_collapse_rejects_nonfinite_combined_length(self):
        tree = Tree(
            '(((A:1,B:1):1,C:1):1e308):1e308;',
            parser=1,
        )
        source = Tree('((A:1,B:1):1,C:1);', parser=1)

        with pytest.raises(ValueError, match='non-finite root branch length'):
            transfer_root(tree, source)

    @pytest.mark.parametrize(
        'rooter',
        [
            midpoint_rooting,
            lambda tree: outgroup_rooting(tree, 'A'),
            mv_rooting,
        ],
    )
    def test_rerooting_preserves_missing_branch_lengths(self, rooter):
        tree = Tree('(A,(B,(C,D)));', parser=1)
        rooted = rooter(tree)
        assert all(node.dist is None for node in rooted.traverse() if not node.is_root)

    def test_unchanged_root_preserves_one_sided_missing_root_length(self):
        tree = Tree('((A:1,B:1),(C:1,D:1):0);', parser=1)

        rooted = midpoint_rooting(tree)

        child_by_taxa = {
            frozenset(child.leaf_names()): child
            for child in rooted.get_children()
        }
        assert child_by_taxa[frozenset({'A', 'B'})].dist is None
        assert child_by_taxa[frozenset({'C', 'D'})].dist == pytest.approx(0.0)


class TestTransferRoot:
    def test_transfer(self):
        # tree_from is rooted with (A,B) | (C,D)
        tree_from = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        # tree_to is unrooted
        tree_to = Tree('(A:1,B:1,(C:1,D:1):1);', parser=1)
        result = transfer_root(tree_to, tree_from)
        assert is_rooted(result)
        assert set(result.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_pairwise_distance_preservation(self):
        """Transfer root must preserve pairwise distances from tree_to."""
        tree_from = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        tree_to = Tree('(A:2,B:3,(C:4,D:5):1);', parser=1)
        original = Tree('(A:2,B:3,(C:4,D:5):1);', parser=1)
        result = transfer_root(tree_to, tree_from)
        for l1 in result.leaves():
            for l2 in result.leaves():
                if l1.name != l2.name:
                    orig_d = original.get_distance(l1.name, l2.name)
                    new_d = result.get_distance(l1, l2)
                    assert abs(orig_d - new_d) < 1e-6, \
                        f'{l1.name}-{l2.name}: {orig_d} vs {new_d}'

    def test_preserves_internal_support_on_unrooted_splits(self):
        tree_from = Tree(
            '((A:1,B:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);',
            parser=1,
        )
        tree_to = Tree(
            '((E:1,F:1)40:0.5,'
            '(D:1,(C:1,(A:1,B:1)20:1)30:1)40:0.5);',
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
            canonical_split(frozenset({'A', 'B'}), frozenset({'C', 'D', 'E', 'F'}))
        ] == {20.0}
        assert support_by_split[
            canonical_split(frozenset({'A', 'B', 'C'}), frozenset({'D', 'E', 'F'}))
        ] == {30.0}
        assert support_by_split[
            canonical_split(frozenset({'E', 'F'}), frozenset({'A', 'B', 'C', 'D'}))
        ] == {40.0}

    def test_zero_subroot_length_in_source_does_not_crash(self):
        tree_from = Tree('((A:0,B:0):0,(C:0,D:0):0);', parser=1)
        tree_to = Tree('((A:3,B:3):3,(C:1,D:1):1);', parser=1)
        result = transfer_root(tree_to, tree_from)
        subroot_dists = sorted([child.dist for child in result.get_children()])
        assert subroot_dists == [1.0, 3.0]

    def test_nonzero_root_distance_in_target_does_not_crash(self):
        tree_from = Tree('((A:1,B:1):1,(C:1,D:1):1):0;', parser=1)
        tree_to = Tree('(A:1,B:1,(C:1,D:1):1):2;', parser=1)
        result = transfer_root(tree_to, tree_from)
        assert is_rooted(result)
        assert set(result.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_source_subroot_none_distance_does_not_redistribute(self):
        tree_from = Tree('((A:1,B:1),(C:1,D:1):1);', parser=1)
        tree_to = Tree('((A:1,B:1):2,(C:1,D:1):2);', parser=1)
        result = transfer_root(tree_to, tree_from)
        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist
            for child in result.get_children()
        }
        assert child_dist_by_taxa == {
            frozenset({'A', 'B'}): 2.0,
            frozenset({'C', 'D'}): 2.0,
        }

    @pytest.mark.parametrize('invalid_length', [float('nan'), float('inf'), -1.0])
    def test_invalid_source_root_length_does_not_redistribute(
            self, invalid_length):
        tree_from = Tree('((A:1,B:1):1,(C:1,D:1):3);', parser=1)
        tree_from.get_children()[0].dist = invalid_length
        tree_to = Tree('((A:1,B:1):2,(C:1,D:1):4);', parser=1)

        result = transfer_root(tree_to, tree_from)

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist
            for child in result.get_children()
        }
        assert child_dist_by_taxa == {
            frozenset({'A', 'B'}): 2.0,
            frozenset({'C', 'D'}): 4.0,
        }

    def test_singleton_root_in_target_does_not_create_unnamed_leaf(self):
        tree_from = Tree('((A:1,B:1):1,C:1);', parser=1)
        tree_to = Tree('(((A:1,B:1):1,C:1):1);', parser=1)
        result = transfer_root(tree_to, tree_from)
        leaf_names = list(result.leaf_names())
        assert set(leaf_names) == {'A', 'B', 'C'}
        assert None not in leaf_names

    def test_exact_transfer_rejects_different_leaf_sets(self):
        source = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        target = Tree('((A:1,B:1):1,C:1,D:1,X:1);', parser=1)
        with pytest.raises(ValueError, match='identical leaf labels'):
            transfer_root(target, source)

    def test_intersection_transfer_honors_root_length_redistribution(self):
        target = Tree('((A:1,B:1):2,C:1,D:1,X:1);', parser=1)
        source = Tree('((A:1,B:1):1,(C:1,D:1):3);', parser=1)
        redistributed = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode='intersection',
            redistribute_root_length=True,
        )
        unchanged = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode='intersection',
            redistribute_root_length=False,
        )
        assert sorted(child.dist for child in redistributed.get_children()) == pytest.approx([0.5, 1.5])
        assert sorted(child.dist for child in unchanged.get_children()) == pytest.approx([1.0, 1.0])

    def test_matching_intersection_root_still_honors_length_redistribution(self):
        target = Tree(
            '((A:1,B:1):2,(C:1,D:1,X:1):2);',
            parser=1,
        )
        source = Tree('((A:1,B:1):1,(C:1,D:1):3);', parser=1)

        redistributed = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode='intersection',
            redistribute_root_length=True,
        )
        unchanged = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode='intersection',
            redistribute_root_length=False,
        )

        redistributed_by_taxa = {
            frozenset(child.leaf_names()): child.dist
            for child in redistributed.get_children()
        }
        assert redistributed_by_taxa[frozenset({'A', 'B'})] == pytest.approx(1.0)
        assert redistributed_by_taxa[frozenset({'C', 'D', 'X'})] == pytest.approx(3.0)
        assert sorted(child.dist for child in unchanged.get_children()) == [2.0, 2.0]

    def test_intersection_transfer_skips_missing_source_root_length(self):
        target = Tree(
            '((A:1,B:1):2,(C:1,D:1,X:1):4);',
            parser=1,
        )
        source = Tree('((A:1,B:1),(C:1,D:1):3);', parser=1)

        rooted = root_mod.transfer_root_with_taxon_mode(
            target,
            source,
            taxon_mode='intersection',
            redistribute_root_length=True,
        )

        child_dist_by_taxa = {
            frozenset(child.leaf_names()): child.dist
            for child in rooted.get_children()
        }
        assert child_dist_by_taxa == {
            frozenset({'A', 'B'}): 2.0,
            frozenset({'C', 'D', 'X'}): 4.0,
        }

    def test_matching_exact_root_does_not_fill_missing_target_root_length(self):
        target = Tree('((A:1,B:1),(C:1,D:1):4);', parser=1)
        source = Tree('((A:1,B:1):1,(C:1,D:1):3);', parser=1)

        rooted = transfer_root(target, source)

        child_by_taxa = {
            frozenset(child.leaf_names()): child
            for child in rooted.get_children()
        }
        assert child_by_taxa[frozenset({'A', 'B'})].dist is None
        assert child_by_taxa[frozenset({'C', 'D'})].dist == pytest.approx(4.0)


class TestTaxonomyHttpBounds:
    def test_streaming_reader_stops_at_decoded_body_limit(self, monkeypatch):
        monkeypatch.setattr(root_mod, 'TAXONOMY_HTTP_MAX_BYTES', 8)
        monkeypatch.setattr(root_mod, 'TAXONOMY_HTTP_CHUNK_BYTES', 4)

        class StreamingResponse:
            headers = {}
            encoding = 'utf-8'

            def __init__(self):
                self.num_chunks_yielded = 0
                self.closed = False

            @property
            def text(self):
                raise AssertionError('Streaming responses must not access .text')

            def iter_content(self, chunk_size):
                assert chunk_size == 4
                for chunk in (b'1234', b'5678', b'9abc', b'defg'):
                    self.num_chunks_yielded += 1
                    yield chunk

            def close(self):
                self.closed = True

        response = StreamingResponse()
        with pytest.raises(ValueError, match='exceeds the size limit'):
            root_mod._bounded_response_text(response, 'test source')

        assert response.num_chunks_yielded == 3
        assert response.closed is True


class TestTaxonomyRooting:
    def test_malformed_opentree_json_falls_back_to_next_source(self, monkeypatch):
        install_fake_opentree(
            monkeypatch,
            tnrs_json=[],
            induced_subtree_json={'newick': ''},
        )
        reference = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        monkeypatch.setattr(
            root_mod,
            '_build_ncbi_reference_tree',
            lambda **kwargs: (reference, {'A', 'B', 'C', 'D'}, set()),
        )
        tree = Tree('((A:1,B:1):1,C:1,D:1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='opentree,ncbi', rank='no')
        assert {'A', 'B'} in [set(child.leaf_names()) for child in rooted.get_children()]

    def test_default_source_chain_constant(self):
        assert DEFAULT_TAXONOMY_SOURCE_CHAIN == 'ncbi,opentree,timetree'

    def test_ncbi_source_passes_args_to_ncbi_helpers(self, monkeypatch, tmp_path):
        observed = dict()

        def fake_resolve_ncbi_lineages(tree, taxid_tsv=None, rank='no', args=None, verbose=False):
            observed['lineages_download_dir'] = getattr(args, 'download_dir', None)
            return ({
                'A': [1, 10],
                'B': [1, 10],
                'C': [1, 20],
                'D': [1, 20],
            }, {})

        def fake_taxid2tree(lineages, taxid_counts, args=None):
            observed['tree_download_dir'] = getattr(args, 'download_dir', None)
            return Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)

        monkeypatch.setattr(root_mod, '_resolve_ncbi_lineages', fake_resolve_ncbi_lineages)
        monkeypatch.setattr(root_mod, 'taxid2tree', fake_taxid2tree)

        tree = Tree('(A:1,B:1,(C:1,D:1):1);', parser=1)
        args = make_args(download_dir=str(tmp_path / 'cache'))
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', rank='no', args=args)
        assert observed['lineages_download_dir'] == str(tmp_path / 'cache')
        assert observed['tree_download_dir'] == str(tmp_path / 'cache')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'A', 'B'} in child_leaf_sets
        assert {'C', 'D'} in child_leaf_sets

    def test_infers_ncbi_lineages_from_leaf_labels(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
            },
        )
        original = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        tree = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert abs(original.get_distance(l1.name, l2.name) - rooted.get_distance(l1, l2)) < 1e-6

    def test_ncbi_taxonomic_uses_taxonomy_query_fallbacks(self, monkeypatch):
        observed = dict(queries=[])

        class FakeNCBI:
            def __init__(self, *args, **kwargs):
                self.db = None

            def get_name_translator(self, names):
                observed['queries'].extend(list(names))
                mapping = dict()
                for name in names:
                    if name == 'Dictyostelium discoideum':
                        mapping[name] = [101]
                    elif name == 'Amoeba':
                        mapping[name] = [102]
                    elif name == 'Homo sapiens':
                        mapping[name] = [201]
                    elif name == 'Pan troglodytes':
                        mapping[name] = [202]
                return mapping

            def get_lineage(self, taxid):
                lineage_by_taxid = {
                    1: [1],
                    10: [1, 10],
                    20: [1, 20],
                    101: [1, 10, 101],
                    102: [1, 10, 102],
                    201: [1, 20, 201],
                    202: [1, 20, 202],
                }
                return lineage_by_taxid[int(taxid)]

            def get_taxid_translator(self, taxids):
                names = {
                    1: 'root',
                    10: 'Amoebozoa',
                    20: 'Primates',
                    101: 'Dictyostelium discoideum',
                    102: 'Amoeba',
                    201: 'Homo sapiens',
                    202: 'Pan troglodytes',
                }
                return {int(taxid): names[int(taxid)] for taxid in taxids if int(taxid) in names}

            def get_rank(self, taxids):
                ranks = {
                    1: 'no rank',
                    10: 'phylum',
                    20: 'order',
                    101: 'species',
                    102: 'genus',
                    201: 'species',
                    202: 'species',
                }
                return {int(taxid): ranks.get(int(taxid), 'no rank') for taxid in taxids}

        monkeypatch.setattr(constrain_mod.ete4, 'NCBITaxa', FakeNCBI)

        tree = Tree(
            '(Dictyostelium_cf_discoideum:1,Amoeba_sp_JDSRuffled:1,(Homo_sapiens:1,Pan_troglodytes:1):1);',
            parser=1,
        )
        args = make_args(species_parser='taxonomic')
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no', args=args)
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Dictyostelium_cf_discoideum', 'Amoeba_sp_JDSRuffled'} in child_leaf_sets
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert 'Dictyostelium discoideum' in observed['queries']
        assert 'Amoeba' in observed['queries']
        assert 'Dictyostelium cf discoideum' not in observed['queries']
        assert 'Amoeba sp JDSRuffled' not in observed['queries']

    def test_taxid_tsv_requires_exact_leaf_label_match(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                101: [1, 10, 101],
                200: [1, 20, 200],
                201: [1, 20, 201],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['A', 'B', 'C', 'X'], 'taxid': [100, 101, 200, 201]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(A:1,B:1,(C:1,D:1):1);', parser=1)
        with pytest.raises(ValueError, match='match the leaf labels'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')

    def test_taxid_tsv_numeric_leaf_label_mismatch_raises_valueerror(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                200: [1, 20, 200],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': [1, 3], 'taxid': [100, 200]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(1:1,2:1);', parser=1)
        with pytest.raises(ValueError, match='match the leaf labels'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')

    def test_taxid_tsv_numeric_leaf_labels_are_accepted(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                101: [1, 10, 101],
                200: [1, 20, 200],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': [1, 2, 3], 'taxid': [100, 101, 200]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(1:1,2:1,3:1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'1', '2'} in child_leaf_sets
        assert {'3'} in child_leaf_sets

    def test_taxid_tsv_leading_zero_leaf_labels_are_accepted(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                101: [1, 10, 101],
                200: [1, 20, 200],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['001', '002', '003'], 'taxid': [100, 101, 200]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(001:1,002:1,003:1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'001', '002'} in child_leaf_sets
        assert {'003'} in child_leaf_sets

    def test_taxid_tsv_na_literal_leaf_labels_are_accepted(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                101: [1, 10, 101],
                200: [1, 20, 200],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['NA', 'B', 'C'], 'taxid': [100, 101, 200]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(NA:1,B:1,C:1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'NA', 'B'} in child_leaf_sets
        assert {'C'} in child_leaf_sets

    def test_ambiguous_taxonomy_root_raises(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                100: [1, 100],
                200: [1, 200],
                300: [1, 300],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['A', 'B', 'C'], 'taxid': [100, 200, 300]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(A:1,B:1,C:1);', parser=1)
        with pytest.raises(ValueError, match='ambiguous'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')

    def test_incompatible_taxonomy_bipartition_raises(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
            },
        )
        tree = Tree('((Homo_sapiens:1,Arabidopsis_thaliana:1):1,(Pan_troglodytes:1,Oryza_sativa:1):1);', parser=1)
        with pytest.raises(ValueError, match='ncbi-derived root bipartition'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no')

    def test_ncbi_placeholder_taxa_are_excluded_when_not_interfering(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
                'Unknown': 32644,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                2787823: [1, 2787823],
                12908: [1, 2787823, 12908],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
                32644: [1, 2787823, 12908, 32644],
            },
            taxid_to_name={
                1: 'root',
                10: 'Primates',
                20: 'Plants',
                2787823: 'unclassified entries',
                12908: 'unclassified sequences',
                9606: 'Homo sapiens',
                9598: 'Pan troglodytes',
                3702: 'Arabidopsis thaliana',
                4530: 'Oryza sativa',
                32644: 'unidentified',
            },
            rank_by_taxid={32644: 'species'},
        )
        tree = Tree(
            '((Homo_sapiens:1,(Pan_troglodytes:1,(Unknown_amoeba:1,Unknown_amoeboid:1):1):1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets
        assert {'Homo_sapiens', 'Pan_troglodytes', 'Unknown_amoeba', 'Unknown_amoeboid'} in child_leaf_sets

    def test_ncbi_placeholder_taxa_on_root_path_raise(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
                'Unknown': 32644,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                2787823: [1, 2787823],
                12908: [1, 2787823, 12908],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
                32644: [1, 2787823, 12908, 32644],
            },
            taxid_to_name={
                1: 'root',
                10: 'Primates',
                20: 'Plants',
                2787823: 'unclassified entries',
                12908: 'unclassified sequences',
                9606: 'Homo sapiens',
                9598: 'Pan troglodytes',
                3702: 'Arabidopsis thaliana',
                4530: 'Oryza sativa',
                32644: 'unidentified',
            },
            rank_by_taxid={32644: 'species'},
        )
        tree = Tree(
            '(((Homo_sapiens:1,Pan_troglodytes:1):1,(Unknown_amoeba:1,Unknown_amoeboid:1):1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        with pytest.raises(ValueError, match='interfere'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no')

    def test_timetree_maps_species_names_back_to_original_labels(self, monkeypatch):
        calls = install_fake_timetree(
            monkeypatch,
            upload_html='<div id="prunetree-msg-box"></div>',
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,(Arabidopsis_thaliana:1,Oryza_sativa:1):20);',
        )
        tree = Tree(
            '(Homo_sapiens_gene1:1,Pan_troglodytes_gene1:1,(Arabidopsis_thaliana_gene1:1,Oryza_sativa_gene1:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='timetree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens_gene1', 'Pan_troglodytes_gene1'} in child_leaf_sets
        assert {'Arabidopsis_thaliana_gene1', 'Oryza_sativa_gene1'} in child_leaf_sets
        assert len(calls) == 2
        assert all(stream is True for _, stream in calls)

    def test_timetree_allows_monophyletic_duplicate_species_labels(self, monkeypatch):
        install_fake_timetree(
            monkeypatch,
            upload_html='<div id="prunetree-msg-box"></div>',
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,((Arabidopsis_thaliana:1,Oryza_sativa:1):20,Saccharomyces_cerevisiae:1):30);',
        )
        tree = Tree(
            '((Arabidopsis_thaliana_gene1:1,Oryza_sativa_gene1:1):1,(Saccharomyces_cerevisiae_gene1:1,((Homo_sapiens_gene1:1,Homo_sapiens_gene2:1):1,Pan_troglodytes_gene1:1):1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='timetree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens_gene1', 'Homo_sapiens_gene2', 'Pan_troglodytes_gene1'} in child_leaf_sets
        assert {'Arabidopsis_thaliana_gene1', 'Oryza_sativa_gene1', 'Saccharomyces_cerevisiae_gene1'} in child_leaf_sets

    def test_timetree_unresolved_names_are_excluded_when_not_interfering(self, monkeypatch):
        install_fake_timetree(
            monkeypatch,
            upload_html=(
                '<button id="prunetree-msg-btn" class="error-btn-enabled">Unresolved Names (1)</button>'
                '<div id="unresolved-names">Unknown amoeba (not found in NCBI taxonomy)<br/></div>'
            ),
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,(Arabidopsis_thaliana:1,Oryza_sativa:1):20);',
        )
        tree = Tree(
            '((Homo_sapiens:1,(Pan_troglodytes:1,Unknown_amoeba:1):1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='timetree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes', 'Unknown_amoeba'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_source_chain_falls_back_to_timetree(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_timetree(
            monkeypatch,
            upload_html='<div id="prunetree-msg-box"></div>',
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,(Arabidopsis_thaliana:1,Oryza_sativa:1):20);',
        )
        tree = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi,timetree', taxid_tsv=None, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_opentree_maps_species_names_back_to_original_labels(self, monkeypatch):
        calls = install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((((Homo_sapiens,Pan_troglodytes)Homininae)Primates,(((Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Embryophyta)Viridiplantae)Eukaryota);',
            },
        )
        tree = Tree(
            '(Homo_sapiens_gene1:1,Pan_troglodytes_gene1:1,(Arabidopsis_thaliana_gene1:1,Oryza_sativa_gene1:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='opentree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens_gene1', 'Pan_troglodytes_gene1'} in child_leaf_sets
        assert {'Arabidopsis_thaliana_gene1', 'Oryza_sativa_gene1'} in child_leaf_sets
        assert len(calls) == 2
        assert all(stream is True for _, stream in calls)

    def test_opentree_allows_monophyletic_duplicate_species_labels(self, monkeypatch):
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 5, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,((Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae,Saccharomyces_cerevisiae)Opisthokonta)Eukaryota;',
            },
        )
        tree = Tree(
            '((Arabidopsis_thaliana_gene1:1,Oryza_sativa_gene1:1):1,(Saccharomyces_cerevisiae_gene1:1,((Homo_sapiens_gene1:1,Homo_sapiens_gene2:1):1,Pan_troglodytes_gene1:1):1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='opentree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens_gene1', 'Homo_sapiens_gene2', 'Pan_troglodytes_gene1'} in child_leaf_sets
        assert {'Arabidopsis_thaliana_gene1', 'Oryza_sativa_gene1', 'Saccharomyces_cerevisiae_gene1'} in child_leaf_sets

    def test_opentree_unresolved_name_is_excluded_when_not_interfering(self, monkeypatch):
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': []},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        tree = Tree(
            '((Homo_sapiens:1,(Pan_troglodytes:1,Batrachochytrium_dendrobatidis:1):1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='opentree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes', 'Batrachochytrium_dendrobatidis'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_opentree_unresolved_name_on_root_path_raises(self, monkeypatch):
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': []},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        tree = Tree(
            '(((Homo_sapiens:1,Pan_troglodytes:1):1,Batrachochytrium_dendrobatidis:1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        with pytest.raises(ValueError, match='interfere'):
            taxonomy_rooting(tree, taxonomy_source='opentree', rank='no')

    def test_source_chain_falls_back_to_opentree(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        tree = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi,opentree', taxid_tsv=None, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_default_source_chain_falls_back_to_opentree(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        tree = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        rooted = taxonomy_rooting(tree, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_unknown_taxonomy_source_raises(self):
        tree = Tree('(A:1,B:1,(C:1,D:1):1);', parser=1)
        with pytest.raises(ValueError, match='Unknown taxonomy source'):
            taxonomy_rooting(tree, taxonomy_source='ncbi,unknown', rank='no')


class TestRootMain:
    def test_writes_preserved_nhx_properties_after_rerooting(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk(
            '(((A:2[&&NHX:tip_tag=tip_A],B:3)AB:5[&&NHX:edge_tag=ab],C:7)'
            'ABC:11[&&NHX:root_edge_tag=same],'
            '(D:13,(E:17,F:19)EF:23[&&NHX:edge_tag=ef])'
            'DEF:29[&&NHX:root_edge_tag=same])ROOT:0[&&NHX:root_tag=root_value];'
        )
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            format='1',
            outformat='1',
            method='outgroup',
            outgroup='A',
        )

        root_main(args)

        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        assert tree.props['root_tag'] == 'root_value'
        leaf_a = next(leaf for leaf in tree.leaves() if leaf.name == 'A')
        assert leaf_a.props['tip_tag'] == 'tip_A'
        all_taxa = frozenset(tree.leaf_names())
        internal_properties = {
            canonical_split(
                frozenset(node.leaf_names()),
                all_taxa - frozenset(node.leaf_names()),
            ): node.props.get('edge_tag')
            for node in tree.traverse()
            if not node.is_leaf and not node.is_root
        }
        assert internal_properties[
            canonical_split(frozenset({'A', 'B'}), frozenset({'C', 'D', 'E', 'F'}))
        ] == 'ab'
        assert internal_properties[
            canonical_split(frozenset({'E', 'F'}), frozenset({'A', 'B', 'C', 'D'}))
        ] == 'ef'

    def test_midpoint(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(A:1,B:5,(C:2,D:4):2);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='midpoint',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)

    def test_outgroup(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(A:1,(B:1,(C:1,D:1):1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='A',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_outgroup_requires_label(self, tmp_nwk):
        path = tmp_nwk('(A:1,(B:1,(C:1,D:1):1):1);')
        args = make_args(
            infile=path, outfile='-',
            method='outgroup', outgroup='',
        )
        with pytest.raises(ValueError, match='outgroup'):
            root_main(args)

    def test_outgroup_requires_label_none(self, tmp_nwk):
        path = tmp_nwk('(A:1,(B:1,(C:1,D:1):1):1);')
        args = make_args(
            infile=path, outfile='-',
            method='outgroup', outgroup=None,
        )
        with pytest.raises(ValueError, match='outgroup'):
            root_main(args)

    def test_unknown_method_raises(self, tmp_nwk):
        path = tmp_nwk('(A:1,(B:1,(C:1,D:1):1):1);')
        args = make_args(
            infile=path, outfile='-',
            method='unknown',
        )
        with pytest.raises(ValueError, match='Unknown rooting method'):
            root_main(args)

    def test_transfer(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile=tmp_outfile,
            method='transfer', infile2=path2, format2='auto',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)

    def test_transfer_single_leaf_tree(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('A;', 'tree1.nwk')
        path2 = tmp_nwk('A;', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile=tmp_outfile,
            method='transfer', infile2=path2, format='1', format2='1', outformat='1',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A'}

    def test_transfer_mismatched_raises(self, tmp_nwk):
        path1 = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,(C:1,E:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile='-',
            method='transfer', infile2=path2, format2='auto',
        )
        with pytest.raises(Exception, match='Leaf labels'):
            root_main(args)

    def test_transfer_requires_rooted_infile2(self, tmp_nwk):
        path1 = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('(A:1,B:1,C:1,D:1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile='-',
            method='transfer', infile2=path2, format2='auto',
        )
        with pytest.raises(ValueError, match='must be rooted'):
            root_main(args)

    def test_transfer_requires_bifurcating_root_infile2(self, tmp_nwk):
        path1 = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('(((A:1,B:1):1,C:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile='-',
            method='transfer', infile2=path2, format2='auto',
        )
        with pytest.raises(ValueError, match='exactly two children'):
            root_main(args)

    def test_transfer_incompatible_root_bipartition_raises(self, tmp_nwk):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,C:1):1,(B:1,D:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile='-',
            method='transfer', infile2=path2, format2='auto',
        )
        with pytest.raises(ValueError, match='No root bipartition'):
            root_main(args)

    def test_transfer_requires_infile2(self, tmp_nwk):
        path1 = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'tree1.nwk')
        args = make_args(
            infile=path1, outfile='-',
            method='transfer', infile2='', format2='auto',
        )
        with pytest.raises(ValueError, match='infile2'):
            root_main(args)

    def test_transfer_duplicate_leaf_names_raise(self, tmp_nwk):
        path1 = tmp_nwk('((A:1,A:2):1,B:1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,A:2):1,B:1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile='-',
            method='transfer', infile2=path2, format2='auto',
        )
        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            root_main(args)

    def test_transfer_empty_leaf_labels_raise(self, tmp_nwk):
        path1 = tmp_nwk('(A:1,(:1,B:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,:1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile='-',
            method='transfer', infile2=path2, format2='auto',
        )
        with pytest.raises(ValueError, match='Empty leaf labels'):
            root_main(args)

    def test_mad(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='mad',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_mad_with_too_few_leaves_raises(self, tmp_nwk):
        path = tmp_nwk('(A:1,B:1);')
        args = make_args(
            infile=path, outfile='-',
            method='mad',
        )
        with pytest.raises(ValueError, match='at least 3 leaves'):
            root_main(args)

    def test_mv(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='mv',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_taxonomy(self, monkeypatch, tmp_nwk, tmp_outfile):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
            },
        )
        path = tmp_nwk('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='taxonomy', taxonomy_source='ncbi', taxid_tsv=None, rank='no',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        child_leaf_sets = [set(child.leaf_names()) for child in tree.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_taxonomy_source_chain(self, monkeypatch, tmp_nwk, tmp_outfile):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_timetree(
            monkeypatch,
            upload_html='<div id="prunetree-msg-box"></div>',
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,(Arabidopsis_thaliana:1,Oryza_sativa:1):20);',
        )
        path = tmp_nwk('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='taxonomy', taxonomy_source='ncbi,timetree', taxid_tsv=None, rank='no',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        child_leaf_sets = [set(child.leaf_names()) for child in tree.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_taxonomy_opentree(self, monkeypatch, tmp_nwk, tmp_outfile):
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        path = tmp_nwk('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='taxonomy', taxonomy_source='opentree', taxid_tsv=None, rank='no',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        child_leaf_sets = [set(child.leaf_names()) for child in tree.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_taxonomy_defaults(self, monkeypatch, tmp_nwk, tmp_outfile):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        path = tmp_nwk('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='taxonomy',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        child_leaf_sets = [set(child.leaf_names()) for child in tree.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_wiki_outgroup_single(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method outgroup --outgroup a

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        Output: (a:0.5,(b:1,(c:1,((d:1,e:1):1,f:1):2):1):0.5):0;
        """
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='a',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # 'a' should be in one of the two subroot clades, by itself
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {'a'} in child_leaf_sets or any('a' in s and len(s) == 1 for s in child_leaf_sets)

    def test_wiki_outgroup_multiple(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method outgroup --outgroup a,b

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        Output: ((a:1,b:1):0.5,(c:1,((d:1,e:1):1,f:1):2):0.5):0;
        """
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='a,b',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # {a,b} should be one of the two subroot clades
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {'a', 'b'} in child_leaf_sets

    def test_wiki_root_transfer(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method transfer

        input1.nwk: (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        input2.nwk: ((((a:1,b:1):1,c:1),f:1):1,(d:3,e:3):1):0;
        Output: ((d:1,e:1):0.5,(f:1,(c:1,(b:1,a:1):1):2):0.5)Root:0;
        """
        path1 = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;', 'tree1.nwk')
        path2 = tmp_nwk('((((a:1,b:1):1,c:1),f:1):1,(d:3,e:3):1):0;', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile=tmp_outfile,
            method='transfer', infile2=path2, format2='auto',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # {d,e} should be in one subroot clade (like the source tree)
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {'d', 'e'} in child_leaf_sets

    def test_wiki_midpoint_rooting(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method midpoint

        Input:  ((((a:5,b:1):1,c:3):1,f:1):1,(d:1,e:1):1):0;
        Output: ((a:5,b:1):0.5,(c:3,(f:1,(d:1,e:1):2):1):0.5):0;
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk('((((a:5,b:1):1,c:3):1,f:1):1,(d:1,e:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='midpoint',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # Verify exact root-to-tip distances from wiki output
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['a'] - 5.5) < 1e-6
        assert abs(dists['b'] - 1.5) < 1e-6
        assert abs(dists['c'] - 3.5) < 1e-6
        assert abs(dists['f'] - 2.5) < 1e-6
        assert abs(dists['d'] - 4.5) < 1e-6
        assert abs(dists['e'] - 4.5) < 1e-6
        # Root children should both have dist 0.5
        children = tree.get_children()
        for c in children:
            assert abs(c.dist - 0.5) < 1e-6

    def test_wiki_outgroup_single_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki outgroup single: verify exact root-to-tip distances.

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        Output: (a:0.5,(b:1,(c:1,((d:1,e:1):1,f:1):2):1):0.5):0;
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='a',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['a'] - 0.5) < 1e-6
        assert abs(dists['b'] - 1.5) < 1e-6
        assert abs(dists['c'] - 2.5) < 1e-6
        assert abs(dists['f'] - 4.5) < 1e-6
        assert abs(dists['d'] - 5.5) < 1e-6
        assert abs(dists['e'] - 5.5) < 1e-6

    def test_wiki_outgroup_multiple_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki outgroup multiple: verify exact root-to-tip distances.

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        Output: ((a:1,b:1):0.5,(c:1,((d:1,e:1):1,f:1):2):0.5):0;
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='a,b',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['a'] - 1.5) < 1e-6
        assert abs(dists['b'] - 1.5) < 1e-6
        assert abs(dists['c'] - 1.5) < 1e-6
        assert abs(dists['f'] - 3.5) < 1e-6
        assert abs(dists['d'] - 4.5) < 1e-6
        assert abs(dists['e'] - 4.5) < 1e-6

    def test_wiki_transfer_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki transfer: verify exact root-to-tip distances.

        input1: (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        input2: ((((a:1,b:1):1,c:1):1,f:1):1,(d:3,e:3):1):0;
        Output: ((d:1,e:1):0.5,(f:1,(c:1,(b:1,a:1):1):2):0.5)Root:0;
        """
        # Note: Added explicit :1 to internal nodes that were missing dist in original wiki example
        path1 = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;', 'tree1.nwk')
        path2 = tmp_nwk('((((a:1,b:1):1,c:1):1,f:1):1,(d:3,e:3):1):0;', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile=tmp_outfile,
            method='transfer', infile2=path2, format2='auto',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['d'] - 1.5) < 1e-6
        assert abs(dists['e'] - 1.5) < 1e-6
        assert abs(dists['f'] - 1.5) < 1e-6
        assert abs(dists['c'] - 3.5) < 1e-6
        assert abs(dists['a'] - 4.5) < 1e-6
        assert abs(dists['b'] - 4.5) < 1e-6

    def test_wiki_mv_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki MV rooting: verify exact root-to-tip distances.

        Input: ((A:1,B:2):1,(C:3,(D:1,E:2):1):1);
        Output: ((C:3,(D:1,E:2):1):0.416667,(A:1,B:2):1.58333);
        Root at x=5/12 from (C,(D,E)) side.
        """
        path = tmp_nwk('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='mv',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['C'] - 41/12) < 1e-4
        assert abs(dists['D'] - 29/12) < 1e-4
        assert abs(dists['E'] - 41/12) < 1e-4
        assert abs(dists['A'] - 31/12) < 1e-4
        assert abs(dists['B'] - 43/12) < 1e-4
