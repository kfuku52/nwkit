import math

import pandas as pd
import pytest

from nwkit.asr import _er_transition_matrix, asr_main
from tests.helpers import make_args


def _write_trait(tmp_path, rows, name='traits.tsv'):
    path = tmp_path / name
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
    return str(path)


class TestAsrMain:
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
        assert {'branch_id', 'parent', 'node_type', 'name', 'map_state', 'map_probability', 'p_x', 'p_y'}.issubset(table.columns)
        assert len(table.index) == 4
        assert set(table['node_type']) == {'root', 'intnode', 'tip'}
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
            'node_type',
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
