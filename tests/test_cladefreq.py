import pandas as pd
import pytest

from nwkit.cladefreq import cladefreq_main
from tests.helpers import make_args


def _write_tree_collection(tmp_path, trees, name='trees.nwk'):
    path = tmp_path / name
    path.write_text('\n'.join(trees) + '\n')
    return str(path)


class TestCladefreqMain:
    @pytest.mark.parametrize('threads', [1, 2])
    def test_rejects_unrooted_input_in_every_execution_mode(self, tmp_path, threads):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '(A:1,B:1,(C:1,D:1):1);',
            ],
            name='contains-unrooted.nwk',
        )

        with pytest.raises(ValueError, match=r'Input tree 2 is not rooted'):
            cladefreq_main(
                make_args(
                    infile=infile,
                    outfile=str(tmp_path / 'cladefreq.tsv'),
                    reference=None,
                    reference_format='auto',
                    weight_tsv=None,
                    support_scale='percent',
                    threads=threads,
                )
            )

    def test_threaded_cladefreq_matches_single_thread(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
                '((A:1,D:1):1,(B:1,C:1):1);',
            ],
            name='threaded.nwk',
        )
        single_out = tmp_path / 'single.tsv'
        threaded_out = tmp_path / 'threaded.tsv'
        common = dict(
            infile=infile,
            reference=None,
            reference_format='auto',
            weight_tsv=None,
            support_scale='percent',
        )

        cladefreq_main(make_args(outfile=str(single_out), threads=1, **common))
        cladefreq_main(make_args(outfile=str(threaded_out), threads=2, **common))

        pd.testing.assert_frame_equal(
            pd.read_csv(single_out, sep='\t'),
            pd.read_csv(threaded_out, sep='\t'),
        )

    def test_reports_internal_clade_frequencies(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
        )
        outfile = tmp_path / 'cladefreq.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            reference=None,
            reference_format='auto',
            weight_tsv=None,
            support_scale='percent',
        )
        cladefreq_main(args)
        table = pd.read_csv(outfile, sep='\t')
        ab_row = table.loc[table['descendant_taxa'] == 'A,B'].iloc[0]
        ac_row = table.loc[table['descendant_taxa'] == 'A,C'].iloc[0]
        assert abs(ab_row['frequency'] - (200.0 / 3.0)) < 1e-4
        assert abs(ac_row['frequency'] - (100.0 / 3.0)) < 1e-4

    def test_reference_tree_flags_present_clades(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
            name='with_reference.nwk',
        )
        reference = tmp_path / 'reference.nwk'
        reference.write_text('((A:1,C:1)42:1,(B:1,D:1)84:1);')
        outfile = tmp_path / 'cladefreq.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            reference=str(reference),
            reference_format='auto',
            weight_tsv=None,
            support_scale='percent',
        )
        cladefreq_main(args)
        table = pd.read_csv(outfile, sep='\t')
        ac_row = table.loc[table['descendant_taxa'] == 'A,C'].iloc[0]
        ab_row = table.loc[table['descendant_taxa'] == 'A,B'].iloc[0]
        assert bool(ac_row['in_reference']) is True
        assert abs(ac_row['reference_support'] - 42.0) < 1e-6
        assert bool(ab_row['in_reference']) is False

    def test_reference_tree_rejects_duplicate_leaf_labels(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            ['((A:1,B:1):1,(C:1,D:1):1);'],
        )
        reference = tmp_path / 'duplicate-reference.nwk'
        reference.write_text(
            '(((A:1,A:1):1,B:1):1,(C:1,D:1):1);'
        )

        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            cladefreq_main(
                make_args(
                    infile=infile,
                    outfile=str(tmp_path / 'cladefreq.tsv'),
                    reference=str(reference),
                    reference_format='auto',
                    weight_tsv=None,
                    support_scale='percent',
                )
            )

    def test_reference_tree_rejects_unrooted_representation(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            ['((A:1,B:1):1,(C:1,D:1):1);'],
        )
        reference = tmp_path / 'unrooted-reference.nwk'
        reference.write_text(
            '(A:1,B:1,(C:1,D:1):1);'
        )

        with pytest.raises(ValueError, match='reference.*must be rooted'):
            cladefreq_main(
                make_args(
                    infile=infile,
                    outfile=str(tmp_path / 'cladefreq.tsv'),
                    reference=str(reference),
                    reference_format='auto',
                    weight_tsv=None,
                    support_scale='percent',
                )
            )

    def test_weighted_frequencies_use_tree_weights(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
            name='weighted.nwk',
        )
        weight_tsv = tmp_path / 'weights.tsv'
        pd.DataFrame({'weight': [3.0, 1.0]}).to_csv(weight_tsv, sep='\t', index=False)
        outfile = tmp_path / 'cladefreq.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            reference=None,
            reference_format='auto',
            weight_tsv=str(weight_tsv),
            support_scale='proportion',
        )
        cladefreq_main(args)
        table = pd.read_csv(outfile, sep='\t')
        ab_row = table.loc[table['descendant_taxa'] == 'A,B'].iloc[0]
        assert abs(ab_row['frequency'] - 0.75) < 1e-6
        assert abs(ab_row['weight_sum'] - 3.0) < 1e-6

    def test_rejects_an_unrepresentable_total_weight(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1);',
            ],
            name='overflowing-weights.nwk',
        )
        weight_tsv = tmp_path / 'weights.tsv'
        pd.DataFrame({'weight': [1e308, 1e308]}).to_csv(
            weight_tsv,
            sep='\t',
            index=False,
        )

        with pytest.raises(ValueError, match='sum of tree weights is too large'):
            cladefreq_main(
                make_args(
                    infile=infile,
                    outfile=str(tmp_path / 'cladefreq.tsv'),
                    reference=None,
                    reference_format='auto',
                    weight_tsv=str(weight_tsv),
                    support_scale='proportion',
                )
            )

    def test_rejects_nan_weights(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
            name='weighted_nan.nwk',
        )
        weight_tsv = tmp_path / 'weights.tsv'
        pd.DataFrame({'weight': [float('nan'), 1.0]}).to_csv(weight_tsv, sep='\t', index=False)
        outfile = tmp_path / 'cladefreq.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            reference=None,
            reference_format='auto',
            weight_tsv=str(weight_tsv),
            support_scale='percent',
        )
        with pytest.raises(ValueError, match='weight'):
            cladefreq_main(args)
