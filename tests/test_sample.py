import pandas as pd
import pytest
from ete4 import Tree

from nwkit.sample import (
    _leaf_path_edges,
    _path_gain,
    apply_filters,
    parse_filter_spec,
    sample_main,
    select_max_pd,
    sort_candidates,
)
from tests.helpers import make_args


def make_sample_args(**kwargs):
    defaults = {
        'trait': None,
        'n': 1,
        'method': 'max-pd',
        'filter': [],
        'rank': [],
        'allow_fewer': False,
        'report': None,
    }
    defaults.update(kwargs)
    return make_args(**defaults)


class TestFilterAndRankParsing:
    def test_parse_filter_spec_accepts_threshold_expression(self):
        assert parse_filter_spec('busco_complete_pct:ge:80') == (
            'busco_complete_pct',
            'ge',
            '80',
        )

    def test_parse_filter_spec_rejects_invalid_operator(self):
        with pytest.raises(ValueError, match='Invalid --filter operator'):
            parse_filter_spec('busco_complete_pct:min:80')

    def test_apply_filters_supports_numeric_bounds(self):
        dataframe = pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C'],
                'busco_complete_pct': [90, 75, 85],
                'num_seq': [50000, 10000, 300000],
            }
        )
        filtered = apply_filters(
            dataframe,
            ['busco_complete_pct:ge:80', 'num_seq:le:200000'],
        )
        assert filtered['leaf_name'].tolist() == ['A']

    def test_sort_candidates_uses_multiple_rank_columns(self):
        dataframe = pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C'],
                'busco_complete_pct': [90, 95, 95],
                'num_seq': [50000, 100000, 75000],
            }
        )
        ranked = sort_candidates(
            dataframe,
            ['busco_complete_pct:desc', 'num_seq:asc'],
        )
        assert ranked['leaf_name'].tolist() == ['C', 'B', 'A']


class TestSampleMain:
    def test_select_max_pd_matches_reference_greedy_order(self):
        tree = Tree('(((A:1.5,B:0.5):2,C:1):1,(D:3,(E:1,F:1):2):1);', parser=1)
        candidate_order = ['C', 'E', 'A', 'F', 'D', 'B']
        leaf_by_name = {leaf.name: leaf for leaf in tree.leaves()}
        path_edges_by_leaf = {
            leaf_name: _leaf_path_edges(leaf_by_name[leaf_name])
            for leaf_name in candidate_order
        }

        covered_edges = set()
        expected = []
        pd_total = 0.0
        for _ in range(4):
            best_leaf = None
            best_gain = None
            for leaf_name in candidate_order:
                if any(row['leaf_name'] == leaf_name for row in expected):
                    continue
                gain = _path_gain(path_edges_by_leaf[leaf_name], covered_edges)
                if best_leaf is None or gain > best_gain:
                    best_leaf = leaf_name
                    best_gain = gain
            for edge, _length in path_edges_by_leaf[best_leaf]:
                covered_edges.add(edge)
            pd_total += best_gain
            expected.append({'leaf_name': best_leaf, 'pd_gain': best_gain, 'pd_total': pd_total})

        assert select_max_pd(candidate_order, path_edges_by_leaf, 4) == expected

    def test_sample_without_trait_uses_tree_only_max_pd(self, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(((A:1,B:1):1,C:1):1,(D:1,E:1):1);', encoding='utf-8')
        out_tree = tmp_path / 'sampled.nwk'
        out_table = tmp_path / 'sampled.tsv'

        args = make_sample_args(
            infile=str(tree_path),
            outfile=str(out_tree),
            n=2,
            report=str(out_table),
        )
        sample_main(args)

        selected = pd.read_csv(out_table, sep='\t')
        assert selected['leaf_name'].tolist() == ['A', 'D']
        out_text = out_tree.read_text(encoding='utf-8')
        assert 'A' in out_text
        assert 'D' in out_text
        for excluded in ['B', 'C', 'E']:
            assert excluded not in out_text

    def test_sample_filters_quality_size_and_uses_rank_as_pd_tiebreaker(self, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        trait_path = tmp_path / 'trait.tsv'
        out_tree = tmp_path / 'sampled.nwk'
        out_table = tmp_path / 'sampled.tsv'
        tree_path.write_text('((A:1,B:1):1,(C:1,D:1):1,E:1);', encoding='utf-8')
        pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C', 'D', 'E'],
                'busco_complete_pct': [90, 79, 95, 95, 85],
                'num_seq': [50000, 10000, 150000, 300000, 100000],
            }
        ).to_csv(trait_path, sep='\t', index=False)

        args = make_sample_args(
            infile=str(tree_path),
            outfile=str(out_tree),
            trait=str(trait_path),
            n=2,
            filter=['busco_complete_pct:ge:80', 'num_seq:le:200000'],
            rank=['num_seq:asc', 'busco_complete_pct:desc'],
            report=str(out_table),
        )
        sample_main(args)

        selected = pd.read_csv(out_table, sep='\t')
        assert selected['leaf_name'].tolist() == ['A', 'C']
        assert selected['num_seq'].tolist() == [50000, 150000]
        assert 'B' not in selected['leaf_name'].tolist()
        assert 'D' not in selected['leaf_name'].tolist()

    def test_ranked_method_selects_ranked_candidates_after_filters(self, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        trait_path = tmp_path / 'trait.tsv'
        out_tree = tmp_path / 'sampled.nwk'
        out_table = tmp_path / 'sampled.tsv'
        tree_path.write_text('((A:1,B:1):1,(C:1,D:1):1,E:1);', encoding='utf-8')
        pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C', 'D', 'E'],
                'busco_complete_pct': [90, 79, 95, 95, 85],
                'num_seq': [50000, 10000, 150000, 300000, 100000],
            }
        ).to_csv(trait_path, sep='\t', index=False)

        args = make_sample_args(
            infile=str(tree_path),
            outfile=str(out_tree),
            trait=str(trait_path),
            n=2,
            method='ranked',
            filter=['busco_complete_pct:ge:80', 'num_seq:le:200000'],
            rank=['num_seq:asc'],
            report=str(out_table),
        )
        sample_main(args)

        selected = pd.read_csv(out_table, sep='\t')
        assert selected['leaf_name'].tolist() == ['A', 'E']

    def test_sample_requires_enough_filtered_candidates_by_default(self, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        trait_path = tmp_path / 'trait.tsv'
        tree_path.write_text('(A:1,B:1,C:1);', encoding='utf-8')
        pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C'],
                'busco_complete_pct': [90, 70, 60],
            }
        ).to_csv(trait_path, sep='\t', index=False)

        args = make_sample_args(
            infile=str(tree_path),
            outfile=str(tmp_path / 'sampled.nwk'),
            trait=str(trait_path),
            n=2,
            filter=['busco_complete_pct:ge:80'],
        )
        with pytest.raises(ValueError, match='fewer than --n'):
            sample_main(args)
