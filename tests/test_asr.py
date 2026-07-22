import math

import pandas as pd
import pytest
from ete4 import Tree

import nwkit.asr as asr
from nwkit.asr import _er_transition_matrix, asr_main
from tests.helpers import make_args


def _write_trait(tmp_path, rows, name='traits.tsv'):
    path = tmp_path / name
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
    return str(path)


class TestAsrMain:
    def test_rejects_colliding_primary_and_model_outputs(self, tmp_path):
        output = tmp_path / 'same.tsv'
        args = make_args(
            outfile=str(output),
            model_out=str(output),
            tree_out=None,
            stochastic_map_out=None,
        )
        with pytest.raises(ValueError, match='Output paths must be distinct'):
            asr_main(args)

    def test_probabilities_report_internal_nodes_and_missing_tips(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'habitat': 'x'},
                {'leaf_name': 'B', 'habitat': 'x'},
                {'leaf_name': 'C', 'habitat': 'y'},
                {'leaf_name': 'D', 'habitat': ''},
            ],
        )
        outfile = tmp_path / 'asr.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column='habitat',
            states='x,y',
            missing_values=None,
            model='ER',
            rate=0.1,
            rate_bounds=None,
            root_prior='equal',
            target='intnode,missing_tip',
            output='probabilities',
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert {'branch_id', 'parent', 'node_class', 'name', 'map_state', 'map_probability', 'p_x', 'p_y'}.issubset(table.columns)
        assert len(table.index) == 4
        assert set(table['node_class']) == {'root', 'intnode', 'leaf'}
        for _, row in table.iterrows():
            assert abs((row['p_x'] + row['p_y']) - 1.0) < 1e-9
            assert row['map_state'] in ['x', 'y']
            assert abs(row['map_probability'] - max(row['p_x'], row['p_y'])) < 1e-9
        d_row = table.loc[table['name'] == 'D'].iloc[0]
        assert bool(d_row['is_imputed']) is True
        assert d_row['p_y'] > d_row['p_x']

    def test_map_output_can_report_all_nodes(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'red'},
                {'leaf_name': 'C', 'state': 'blue'},
            ],
        )
        outfile = tmp_path / 'asr_map.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.2,
            rate_bounds=None,
            root_prior='equal',
            target='all',
            output='map',
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert list(table.columns) == [
            'branch_id',
            'parent',
            'node_class',
            'name',
            'observed_state',
            'is_imputed',
            'state',
            'probability',
        ]
        assert len(table.index) == 5
        a_row = table.loc[table['name'] == 'A'].iloc[0]
        assert a_row['state'] == 'red'
        assert abs(a_row['probability'] - 1.0) < 1e-9

    def test_default_output_includes_observed_tips(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'red'},
                {'leaf_name': 'C', 'state': 'blue'},
            ],
        )
        outfile = tmp_path / 'asr_default.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.2,
            rate_bounds=None,
            root_prior='equal',
            output='probabilities',
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert len(table.index) == 5
        tip_rows = table.loc[table['node_class'] == 'leaf'].set_index('name')
        assert set(tip_rows.index) == {'A', 'B', 'C'}
        assert tip_rows.loc['A', 'observed_state'] == 'red'
        assert bool(tip_rows.loc['A', 'is_imputed']) is False
        assert tip_rows.loc['A', 'p_red'] == pytest.approx(1.0)
        assert tip_rows.loc['C', 'observed_state'] == 'blue'
        assert tip_rows.loc['C', 'p_blue'] == pytest.approx(1.0)

    def test_omitted_trait_leaf_is_imputed(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'red'},
            ],
            name='partial_traits.tsv',
        )
        outfile = tmp_path / 'imputed.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.2,
            rate_bounds=None,
            root_prior='equal',
            target='missing_tip',
            output='probabilities',
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert table['name'].tolist() == ['C']
        assert bool(table.loc[0, 'is_imputed']) is True
        assert abs((table.loc[0, 'p_red'] + table.loc[0, 'p_blue']) - 1.0) < 1e-9

    def test_rejects_observed_state_not_listed_in_states(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(A:1,B:1);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'green'},
            ],
        )
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / 'out.tsv'),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.2,
            rate_bounds=None,
            root_prior='equal',
            target='all',
            output='probabilities',
        )
        with pytest.raises(ValueError, match='--states'):
            asr_main(args)

    def test_zero_rate_rejects_conflicting_tip_states(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(A:1,B:1);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'blue'},
            ],
        )
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / 'out.tsv'),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.0,
            rate_bounds=None,
            root_prior='equal',
            target='all',
            output='probabilities',
        )
        with pytest.raises(ValueError, match='zero likelihood'):
            asr_main(args)

    def test_er_transition_matrix_matches_two_state_formula(self):
        matrix = _er_transition_matrix(branch_length=1.5, rate=0.2, num_states=2)
        decay = math.exp(-2.0 * 0.2 * 1.5)
        assert matrix[0, 0] == pytest.approx(0.5 + 0.5 * decay)
        assert matrix[0, 1] == pytest.approx(0.5 - 0.5 * decay)
        assert matrix[1, 0] == pytest.approx(0.5 - 0.5 * decay)
        assert matrix[1, 1] == pytest.approx(0.5 + 0.5 * decay)

    def test_inside_likelihood_reuses_transition_matrix_for_equal_branch_lengths(self, monkeypatch):
        tree = Tree('((A:1,B:1):1,C:1);', parser=1)
        states = ['red', 'blue']
        likelihood_by_leaf = {
            'A': pd.Series([1.0, 0.0]).to_numpy(),
            'B': pd.Series([1.0, 0.0]).to_numpy(),
            'C': pd.Series([0.0, 1.0]).to_numpy(),
        }
        rate_matrix = asr._build_rate_matrix('ER', states, [0.2])
        call_count = {'value': 0}
        original_transition_matrix = asr._transition_matrix

        def counted_transition_matrix(rate_matrix_arg, branch_length_arg):
            call_count['value'] += 1
            return original_transition_matrix(rate_matrix_arg, branch_length_arg)

        monkeypatch.setattr(asr, '_transition_matrix', counted_transition_matrix)
        asr._compute_inside_likelihoods(tree, likelihood_by_leaf, rate_matrix)
        assert call_count['value'] == 1

    def test_model_out_reports_fixed_er_metadata(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'red'},
                {'leaf_name': 'C', 'state': 'blue'},
            ],
        )
        model_out = tmp_path / 'model.tsv'
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / 'asr.tsv'),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.2,
            rate_bounds=None,
            root_prior='equal',
            target='all',
            output='map',
            model_out=str(model_out),
        )
        asr_main(args)
        table = pd.read_csv(model_out, sep='\t')
        row = table.iloc[0]
        assert row['model'] == 'ER'
        assert bool(row['rate_estimated']) is False
        assert row['states'] == 'red,blue'
        assert row['rate'] == pytest.approx(0.2)
        assert row['rate_red_to_blue'] == pytest.approx(0.2)
        assert row['rate_blue_to_red'] == pytest.approx(0.2)

    def test_sym_and_ard_models_estimate_rates(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'blue'},
                {'leaf_name': 'C', 'state': 'green'},
                {'leaf_name': 'D', 'state': 'red'},
            ],
        )
        for model in ['SYM', 'ARD']:
            outfile = tmp_path / '{}.tsv'.format(model.lower())
            model_out = tmp_path / '{}_model.tsv'.format(model.lower())
            args = make_args(
                infile=infile,
                outfile=str(outfile),
                trait=trait,
                state_column='state',
                states='red,blue,green',
                missing_values=None,
                model=model,
                rate=None,
                rate_bounds='1e-4,10',
                root_prior='equal',
                target='intnode',
                output='probabilities',
                model_out=str(model_out),
            )
            asr_main(args)
            model_table = pd.read_csv(model_out, sep='\t')
            assert model_table.loc[0, 'model'] == model
            assert bool(model_table.loc[0, 'rate_estimated']) is True
            rate_columns = [
                column for column in model_table.columns
                if column.startswith('rate_') and column != 'rate_bounds' and column != 'rate_estimated'
            ]
            assert all(float(model_table.loc[0, column]) >= 0.0 for column in rate_columns)
            out_table = pd.read_csv(outfile, sep='\t')
            probability_columns = [column for column in out_table.columns if column.startswith('p_')]
            for _, row in out_table.iterrows():
                assert abs(sum(float(row[column]) for column in probability_columns) - 1.0) < 1e-6

    def test_ambiguous_states_use_multi_hot_tip_likelihood(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red|blue'},
                {'leaf_name': 'B', 'state': 'red'},
                {'leaf_name': 'C', 'state': 'blue'},
            ],
        )
        outfile = tmp_path / 'ambiguous.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.2,
            rate_bounds=None,
            root_prior='empirical',
            target='tip',
            output='probabilities',
            ambiguous_separator='|',
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep='\t')
        a_row = table.loc[table['name'] == 'A'].iloc[0]
        assert a_row['observed_state'] == 'red|blue'
        assert a_row['p_red'] > 0.0
        assert a_row['p_blue'] > 0.0

    def test_tree_out_writes_nhx_annotations(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'red'},
                {'leaf_name': 'C', 'state': 'blue'},
            ],
        )
        tree_out = tmp_path / 'annotated.nwk'
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / 'asr.tsv'),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.2,
            rate_bounds=None,
            root_prior='equal',
            target='intnode',
            output='map',
            tree_out=str(tree_out),
            tree_outformat='auto',
            tree_annotation='all',
        )
        asr_main(args)
        annotated = tree_out.read_text()
        assert '&&NHX' in annotated
        assert 'asr_state=' in annotated
        assert 'asr_probability=' in annotated
        assert 'asr_p_red=' in annotated

    def test_stochastic_map_out_is_reproducible_with_seed(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'red'},
                {'leaf_name': 'C', 'state': 'blue'},
            ],
        )
        first_out = tmp_path / 'map1.tsv'
        second_out = tmp_path / 'map2.tsv'
        common = dict(
            infile=infile,
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.5,
            rate_bounds=None,
            root_prior='equal',
            target='intnode',
            output='map',
            n_sim=20,
            seed=7,
        )
        asr_main(make_args(outfile=str(tmp_path / 'asr1.tsv'), stochastic_map_out=str(first_out), **common))
        asr_main(make_args(outfile=str(tmp_path / 'asr2.tsv'), stochastic_map_out=str(second_out), **common))
        first = pd.read_csv(first_out, sep='\t')
        second = pd.read_csv(second_out, sep='\t')
        pd.testing.assert_frame_equal(first, second)
        assert {'from_state', 'to_state', 'mean_count', 'posterior_frequency', 'num_simulations'}.issubset(first.columns)
        assert set(first['num_simulations']) == {20}

    def test_stochastic_map_reuses_uniformization_context_for_equal_branch_lengths(self, tmp_nwk, tmp_path, monkeypatch):
        infile = tmp_nwk('((A:1,B:1):1,C:1);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'red'},
                {'leaf_name': 'C', 'state': 'blue'},
            ],
        )
        call_lengths = list()
        original_builder = asr._build_uniformization_context

        def counted_builder(rate_matrix, branch_length):
            call_lengths.append(float(branch_length))
            return original_builder(rate_matrix, branch_length)

        monkeypatch.setattr(asr, '_build_uniformization_context', counted_builder)
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / 'asr.tsv'),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.5,
            rate_bounds=None,
            root_prior='equal',
            target='intnode',
            output='map',
            stochastic_map_out=str(tmp_path / 'smap.tsv'),
            n_sim=5,
            seed=7,
            threads=1,
        )
        asr_main(args)
        assert call_lengths == [1.0]

    def test_threaded_stochastic_map_matches_single_thread_with_seed(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'red'},
                {'leaf_name': 'C', 'state': 'blue'},
            ],
        )
        single_out = tmp_path / 'single.tsv'
        threaded_out = tmp_path / 'threaded.tsv'
        common = dict(
            infile=infile,
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.5,
            rate_bounds=None,
            root_prior='equal',
            target='intnode',
            output='map',
            n_sim=30,
            seed=11,
        )
        asr_main(make_args(outfile=str(tmp_path / 'asr_single.tsv'), stochastic_map_out=str(single_out), threads=1, **common))
        asr_main(make_args(outfile=str(tmp_path / 'asr_threaded.tsv'), stochastic_map_out=str(threaded_out), threads=4, **common))
        single = pd.read_csv(single_out, sep='\t')
        threaded = pd.read_csv(threaded_out, sep='\t')
        pd.testing.assert_frame_equal(single, threaded)

    def test_stochastic_map_rejects_non_positive_threads(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(A:1,B:1);', 'tree.nwk')
        trait = _write_trait(
            tmp_path,
            [
                {'leaf_name': 'A', 'state': 'red'},
                {'leaf_name': 'B', 'state': 'blue'},
            ],
        )
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / 'asr.tsv'),
            trait=trait,
            state_column='state',
            states='red,blue',
            missing_values=None,
            model='ER',
            rate=0.5,
            rate_bounds=None,
            root_prior='equal',
            target='all',
            output='map',
            stochastic_map_out=str(tmp_path / 'smap.tsv'),
            n_sim=5,
            seed=1,
            threads=0,
        )
        with pytest.raises(ValueError, match='--threads'):
            asr_main(args)
