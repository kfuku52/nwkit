import re
import random
import numpy
import pandas
import pytest
from argparse import Namespace
from ete4 import Tree

import nwkit.constrain as constrain
import nwkit.intersection as intersection
import nwkit.mark as mark
import nwkit.mcmctree as mcmctree
import nwkit.root as root
import nwkit.skim as skim
from nwkit.util import get_subtree_leaf_name_sets


def make_random_binary_tree(num_leaves, rng, leaf_prefix='L'):
    nodes = [
        f'{leaf_prefix}{i}:{rng.uniform(0.1, 2.0):.6f}'
        for i in range(num_leaves)
    ]
    while len(nodes) > 1:
        rng.shuffle(nodes)
        left = nodes.pop()
        right = nodes.pop()
        nodes.append(f'({left},{right}):{rng.uniform(0.1, 2.0):.6f}')
    nwk = re.sub(r':[0-9.]+$', '', nodes[0]) + ';'
    return Tree(nwk, parser=1)


def old_get_remove_names(arr1, arr2, match):
    if match == 'complete':
        return [a1 for a1 in arr1 if not any(a1 == a2 for a2 in arr2)]
    if match == 'prefix':
        out = []
        for a1 in arr1:
            is_hit = False
            for a2 in arr2:
                if a1.startswith(a2) or a2.startswith(a1):
                    is_hit = True
                    break
            if not is_hit:
                out.append(a1)
        return out
    if match == 'backward':
        out = []
        for a1 in arr1:
            is_hit = False
            for a2 in arr2:
                if a1.endswith(a2) or a2.endswith(a1):
                    is_hit = True
                    break
            if not is_hit:
                out.append(a1)
        return out
    raise ValueError(match)


def old_annotate_tree_attr(tree, args):
    for node in tree.traverse():
        node.add_props(
            is_target_leaf=False,
            is_descendant_all_target=False,
            is_target_only_mrca=False,
            is_target_only_mrca_clade=False,
            is_all_mrca=False,
            is_all_mrca_clade=False,
        )
    for node in tree.leaves():
        if re.fullmatch(args.pattern, node.name):
            node.props['is_target_leaf'] = True
    for node in tree.traverse(strategy='postorder'):
        if all([l.props.get('is_target_leaf') for l in node.leaves()]):
            node.props['is_descendant_all_target'] = True
    for node in tree.traverse():
        if node.is_root:
            continue
        if (
            (not node.up.props.get('is_descendant_all_target'))
            and node.props.get('is_descendant_all_target')
        ):
            node.props['is_target_only_mrca'] = True
            node.props['is_target_only_mrca_clade'] = True
            for clade_node in node.descendants():
                clade_node.props['is_target_only_mrca_clade'] = True
    target_leaves = [leaf for leaf in tree.leaves() if leaf.props.get('is_target_leaf')]
    if len(target_leaves) > 0:
        if len(target_leaves) == 1:
            all_mrca_node = target_leaves[0]
        else:
            all_mrca_node = tree.common_ancestor(target_leaves)
        all_mrca_node.props['is_all_mrca'] = True
        all_mrca_node.props['is_all_mrca_clade'] = True
        for clade_node in all_mrca_node.descendants():
            clade_node.props['is_all_mrca_clade'] = True
    return tree


def node_prop_map(tree, prop_keys):
    out = dict()
    for node in tree.traverse():
        sig = tuple(sorted(node.leaf_names()))
        out[sig] = tuple(bool(node.props.get(k)) for k in prop_keys)
    return out


def old_are_both_lineage_included(node, leaf_names):
    child_names_list = [list(child.leaf_names()) for child in node.get_children()]
    for child_names in child_names_list:
        if not any([child_name in leaf_names for child_name in child_names]):
            return False
    return True


def old_are_two_lineage_rank_differentiated(node, taxids, ta_leaf_names):
    child_taxids = list()
    for child in node.get_children():
        child_leaf_names = list(child.leaf_names())
        ch_taxids = [t for t, ln in zip(taxids, ta_leaf_names) if ln in child_leaf_names]
        child_taxids.append(ch_taxids)
    if len(set(child_taxids[0]) - set(child_taxids[1])) == 0:
        return False
    if len(set(child_taxids[1]) - set(child_taxids[0])) == 0:
        return False
    return True

def old_outgroup_rooting(tree, outgroup_str):
    for node in tree.traverse():
        if node.dist is None:
            node.dist = 0.0
    outgroup_list = outgroup_str.split(',')
    outgroup_nodes = [node for node in tree.traverse() if node.name in outgroup_list]
    if len(outgroup_nodes) == 0:
        raise SystemExit(1)
    elif len(outgroup_nodes) == 1:
        outgroup_node = outgroup_nodes[0]
    else:
        outgroup_node = tree.common_ancestor(outgroup_nodes)
    if outgroup_node is tree:
        non_outgroup_leaves = [node for node in tree.leaves() if node not in outgroup_nodes]
        tree.set_outgroup(non_outgroup_leaves[0])
        outgroup_node = tree.common_ancestor(outgroup_nodes)
    if outgroup_node is tree:
        raise SystemExit(1)
    tree.set_outgroup(outgroup_node)
    return tree


def old_populate_leaves(new_clade, taxid, lineages):
    for sp in lineages.keys():
        lineage = lineages[sp]
        if not taxid in lineage:
            continue
        new_leaf = Tree({'name': sp})
        new_leaf.add_props(ancestors=lineage)
        new_clade.add_child(new_leaf)
    return new_clade


def old_add_new_clade(clades, new_clade):
    flag_append = True
    new_leaf_set = set(new_clade.leaf_names())
    delete_index = list()
    for j in range(len(clades)):
        clade = clades[j]
        clade_leaf_set = set(clade.leaf_names())
        intersection_set = clade_leaf_set.intersection(new_leaf_set)
        if intersection_set == new_leaf_set:
            flag_append = False
            continue
        elif intersection_set == clade_leaf_set:
            for nnode in new_clade.leaves():
                for int_name in intersection_set:
                    if nnode.name == int_name:
                        nnode.delete()
                        break
            new_clade.add_child(clade)
            delete_index.append(j)
    for di in sorted(delete_index)[::-1]:
        del clades[di]
    if flag_append:
        new_clade_leaves = list(new_clade.leaf_names())
        if len(new_clade_leaves) != len(set(new_clade_leaves)):
            txt = 'Redundant leaves in the new clade: {}'.format(', '.join(new_clade_leaves))
            raise Exception(txt)
        clades.append(new_clade)
    return clades


def old_taxid2tree(lineages, taxid_counts):
    ncbi = constrain.ete4.NCBITaxa()
    is_multiple = taxid_counts[:, 1] > 1
    multi_counts = taxid_counts[is_multiple, :]
    clades = list()
    for i in numpy.arange(multi_counts.shape[0]):
        taxid = multi_counts[i, 0]
        ancestors = ncbi.get_lineage(taxid)
        new_clade = Tree()
        new_clade.add_props(ancestors=ancestors)
        new_clade = old_populate_leaves(new_clade, taxid, lineages)
        clades = old_add_new_clade(clades, new_clade)
    assert len(clades) == 1, 'Failed to merge clades into a single tree.'
    return clades[0]


def old_read_trait(args, tree):
    if args.trait is None:
        return pandas.DataFrame({'leaf_name': list(tree.leaf_names())})
    trait_df = pandas.read_csv(args.trait, sep='\t')
    leaf_names_list = list(tree.leaf_names())
    trait_df = trait_df[trait_df['leaf_name'].isin(leaf_names_list)]
    for leaf_name in leaf_names_list:
        if leaf_name not in trait_df['leaf_name'].values:
            trait_df = pandas.concat(
                [trait_df, pandas.DataFrame({'leaf_name': [leaf_name]})],
                ignore_index=True,
            )
    return trait_df


def old_sample_from_groups(trait_df, args):
    if args.prioritize_non_missing and args.group_by is not None:
        if args.filter_by is None:
            sampled_df = trait_df.groupby('group', group_keys=False).apply(
                lambda g: g.sample(frac=1).sort_values(by=[args.group_by]).head(args.retain_per_clade),
                include_groups=False,
            )
        else:
            if args.filter_mode == 'ascending':
                sampled_df = trait_df.groupby('group', group_keys=False).apply(
                    lambda g: g.sample(frac=1).sort_values(
                        by=[args.group_by, args.filter_by],
                        ascending=True,
                    ).head(args.retain_per_clade),
                    include_groups=False,
                )
            elif args.filter_mode == 'descending':
                sampled_df = trait_df.groupby('group', group_keys=False).apply(
                    lambda g: g.sample(frac=1).sort_values(
                        by=[args.group_by, args.filter_by],
                        ascending=False,
                    ).head(args.retain_per_clade),
                    include_groups=False,
                )
            else:
                raise ValueError(args.filter_mode)
    else:
        if args.filter_by is None:
            sampled_df = trait_df.groupby('group', group_keys=False).apply(
                lambda g: g.sample(frac=1).head(args.retain_per_clade),
                include_groups=False,
            )
        else:
            if args.filter_mode == 'ascending':
                sampled_df = trait_df.groupby('group', group_keys=False).apply(
                    lambda g: g.sample(frac=1).sort_values(
                        by=args.filter_by,
                        ascending=True,
                    ).head(args.retain_per_clade),
                    include_groups=False,
                )
            elif args.filter_mode == 'descending':
                sampled_df = trait_df.groupby('group', group_keys=False).apply(
                    lambda g: g.sample(frac=1).sort_values(
                        by=args.filter_by,
                        ascending=False,
                    ).head(args.retain_per_clade),
                    include_groups=False,
                )
            else:
                raise ValueError(args.filter_mode)
    sampled_df = sampled_df.merge(trait_df[['leaf_name', 'group']], on=['leaf_name'])
    sampled_df = sampled_df[trait_df.columns]
    return sampled_df


class TestIntersectionRefactorValidation:
    def test_get_remove_names_matches_naive_randomized(self):
        rng = random.Random(11)
        alphabet = 'abcd'
        for _ in range(120):
            arr1 = []
            arr2 = []
            for _ in range(rng.randint(0, 12)):
                length = rng.randint(0, 4)
                arr1.append(''.join(rng.choice(alphabet) for _ in range(length)))
            for _ in range(rng.randint(0, 10)):
                length = rng.randint(0, 4)
                arr2.append(''.join(rng.choice(alphabet) for _ in range(length)))
            for match in ['complete', 'prefix', 'backward']:
                expected = old_get_remove_names(arr1, arr2, match=match)
                actual = intersection.get_remove_names(arr1, arr2, match=match)
                assert actual == expected


class TestMarkRefactorValidation:
    def test_annotate_tree_attr_matches_old_logic_randomized(self):
        rng = random.Random(23)
        prop_keys = [
            'is_target_leaf',
            'is_descendant_all_target',
            'is_target_only_mrca',
            'is_target_only_mrca_clade',
            'is_all_mrca',
            'is_all_mrca_clade',
        ]
        for i in range(50):
            tree = make_random_binary_tree(num_leaves=rng.randint(5, 12), rng=rng, leaf_prefix=f'T{i}_')
            leaf_names = list(tree.leaf_names())
            mode = rng.randint(0, 3)
            if mode == 0:
                pattern = '.*'
            elif mode == 1:
                pattern = 'NO_MATCH_PATTERN'
            elif mode == 2:
                chosen = rng.sample(leaf_names, k=rng.randint(1, len(leaf_names)))
                pattern = '|'.join(re.escape(x) for x in chosen)
            else:
                chosen = rng.choice(leaf_names)
                pattern = re.escape(chosen.split('_')[0]) + '.*'
            args = Namespace(pattern=pattern)
            tree_old = old_annotate_tree_attr(tree.copy(), args)
            tree_new = mark.annotate_tree_attr(tree.copy(), args)
            expected = node_prop_map(tree_old, prop_keys)
            actual = node_prop_map(tree_new, prop_keys)
            assert actual == expected


class TestMcmctreeRefactorValidation:
    def test_helper_equivalence_with_and_without_cache(self):
        rng = random.Random(31)
        for i in range(40):
            tree = make_random_binary_tree(num_leaves=rng.randint(6, 12), rng=rng, leaf_prefix=f'M{i}_')
            subtree_sets = get_subtree_leaf_name_sets(tree)
            internal_nodes = [n for n in tree.traverse() if not n.is_leaf]
            node = rng.choice(internal_nodes)
            leaves = list(tree.leaf_names())
            selected = rng.sample(leaves, k=rng.randint(0, len(leaves)))
            old_inc = old_are_both_lineage_included(node, selected)
            new_inc_plain = mcmctree.are_both_lineage_included(node, selected, subtree_leaf_name_sets=None)
            new_inc_cached = mcmctree.are_both_lineage_included(node, selected, subtree_leaf_name_sets=subtree_sets)
            assert new_inc_plain == old_inc
            assert new_inc_cached == old_inc
            ta_leaf_names = rng.sample(leaves, k=rng.randint(2, len(leaves)))
            taxids = [rng.randint(1, 20) for _ in ta_leaf_names]
            old_diff = old_are_two_lineage_rank_differentiated(node, taxids, ta_leaf_names)
            new_diff_plain = mcmctree.are_two_lineage_rank_differentiated(
                node, taxids, ta_leaf_names, subtree_leaf_name_sets=None,
            )
            new_diff_cached = mcmctree.are_two_lineage_rank_differentiated(
                node, taxids, ta_leaf_names, subtree_leaf_name_sets=subtree_sets,
            )
            assert new_diff_plain == old_diff
            assert new_diff_cached == old_diff


class TestConstrainRefactorValidation:
    @staticmethod
    def _star_tree(leaf_names):
        leaves = sorted(set(leaf_names))
        assert len(leaves) >= 2
        nwk = '(' + ','.join(f'{x}:1' for x in leaves) + ');'
        return Tree(nwk, parser=1)

    @staticmethod
    def _clade_signatures(clades):
        return sorted([tuple(sorted(c.leaf_names())) for c in clades])

    def test_add_new_clade_matches_old_logic_randomized(self):
        rng = random.Random(37)
        universe = [f'X{i}' for i in range(9)]
        for _ in range(80):
            existing_sets = []
            for _ in range(rng.randint(1, 4)):
                k = rng.randint(2, 6)
                existing_sets.append(set(rng.sample(universe, k=k)))
            new_set = set(rng.sample(universe, k=rng.randint(2, 6)))
            old_clades = [self._star_tree(s) for s in existing_sets]
            new_clades = [self._star_tree(s) for s in existing_sets]
            old_new = self._star_tree(new_set)
            new_new = self._star_tree(new_set)
            old_exc = None
            new_exc = None
            try:
                old_out = old_add_new_clade(old_clades, old_new)
            except Exception as e:
                old_exc = type(e)
                old_out = None
            try:
                new_out = constrain.add_new_clade(new_clades, new_new)
            except Exception as e:
                new_exc = type(e)
                new_out = None
            assert old_exc == new_exc
            if old_exc is None:
                assert self._clade_signatures(new_out) == self._clade_signatures(old_out)

    def test_taxid2tree_matches_old_logic_with_fake_ncbi(self, monkeypatch):
        rng = random.Random(41)

        class FakeNCBI:
            def get_lineage(self, taxid):
                taxid = int(taxid)
                return [1, taxid // 10, taxid]

        monkeypatch.setattr(constrain.ete4, 'NCBITaxa', FakeNCBI)
        for _ in range(8):
            lineages = {}
            num_species = rng.randint(12, 24)
            for i in range(num_species):
                family = 10 + (i // 6)
                genus = 100 + (i // 2)
                species = 1000 + i
                lineages[f'SP{i}'] = [1, family, genus, species]
            taxid_counts = constrain.get_taxid_counts(lineages)
            old_tree = old_taxid2tree(lineages, taxid_counts)
            new_tree = constrain.taxid2tree(lineages, taxid_counts)
            assert set(old_tree.leaf_names()) == set(new_tree.leaf_names())
            old_clades = {tuple(sorted(n.leaf_names())) for n in old_tree.traverse()}
            new_clades = {tuple(sorted(n.leaf_names())) for n in new_tree.traverse()}
            assert new_clades == old_clades


class TestSkimRefactorValidation:
    def test_read_trait_matches_old_logic(self, tmp_path):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        trait_path = tmp_path / 'trait.tsv'
        pandas.DataFrame(
            {
                'leaf_name': ['A', 'C', 'E'],
                'trait': ['x', 'y', 'z'],
                'score': [10, 20, 30],
            }
        ).to_csv(trait_path, sep='\t', index=False)
        args = Namespace(trait=str(trait_path))
        old_df = old_read_trait(args, tree).sort_values('leaf_name').reset_index(drop=True)
        new_df = skim.read_trait(args, tree).sort_values('leaf_name').reset_index(drop=True)
        pandas.testing.assert_frame_equal(new_df, old_df, check_dtype=False)

    @pytest.mark.parametrize(
        'prioritize_non_missing,filter_by,filter_mode',
        [
            (True, None, 'ascending'),
            (True, 'score', 'ascending'),
            (True, 'score', 'descending'),
            (False, 'score', 'ascending'),
            (False, 'score', 'descending'),
        ],
    )
    def test_sample_from_groups_matches_old_on_deterministic_inputs(
        self, prioritize_non_missing, filter_by, filter_mode,
    ):
        trait_df = pandas.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C', 'D', 'E', 'F'],
                'group': [1, 1, 2, 2, 3, 3],
                'trait': [0, 1, 0, 1, 0, 1],
                'score': [10, 20, 30, 40, 50, 60],
            }
        )
        args = Namespace(
            prioritize_non_missing=prioritize_non_missing,
            group_by='trait',
            retain_per_clade=1,
            filter_by=filter_by,
            filter_mode=filter_mode,
        )
        numpy.random.seed(7)
        old_df = old_sample_from_groups(trait_df.copy(), args).sort_values('leaf_name').reset_index(drop=True)
        numpy.random.seed(7)
        new_df = skim.sample_from_groups(trait_df.copy(), args).sort_values('leaf_name').reset_index(drop=True)
        pandas.testing.assert_frame_equal(new_df, old_df, check_dtype=False)


class TestRootRefactorValidation:
    def test_collect_leaf_distance_stats_matches_bruteforce(self):
        rng = random.Random(53)
        for _ in range(30):
            tree = make_random_binary_tree(num_leaves=rng.randint(5, 10), rng=rng, leaf_prefix='R')
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

    def test_outgroup_rooting_matches_old_logic_randomized(self):
        rng = random.Random(67)
        for _ in range(40):
            tree = make_random_binary_tree(num_leaves=rng.randint(6, 12), rng=rng, leaf_prefix='O')
            leaf_names = list(tree.leaf_names())
            out_size = rng.randint(1, min(3, len(leaf_names) - 1))
            outgroup_names = rng.sample(leaf_names, k=out_size)
            outgroup_str = ','.join(outgroup_names)
            old_tree = old_outgroup_rooting(tree.copy(), outgroup_str)
            new_tree = root.outgroup_rooting(tree.copy(), outgroup_str)
            old_splits = {frozenset(child.leaf_names()) for child in old_tree.get_children()}
            new_splits = {frozenset(child.leaf_names()) for child in new_tree.get_children()}
            assert old_splits == new_splits
            for l1 in leaf_names:
                for l2 in leaf_names:
                    if l1 == l2:
                        continue
                    old_d = old_tree.get_distance(l1, l2)
                    new_d = new_tree.get_distance(l1, l2)
                    assert old_d == pytest.approx(new_d, abs=1e-9)
