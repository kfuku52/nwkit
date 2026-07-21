import re
from argparse import Namespace

import pandas as pd
import pytest
from ete4 import Tree

from nwkit.draw import _get_species_overlap_node_types, draw_main


def make_draw_args(**kwargs):
    defaults = {
        'infile': '-',
        'format': 'auto',
        'quoted_node_names': True,
        'outfile': None,
        'image_format': 'auto',
        'species_parser': 'legacy',
        'species_regex': r'^([^_]+_[^_]+)(?:_|$)',
        'species_map_tsv': None,
        'species_overlap_node_plot': 'auto',
        'ladderize': False,
        'trait': None,
        'group_by': None,
        'trait_palette': 'tab10',
        'support_labels': True,
        'support_min': None,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


def extract_svg_text_positions(svg_text):
    pattern = re.compile(r'<text[^>]*y="([^"]+)"[^>]*>([^<]+)</text>')
    return {content: float(y) for y, content in pattern.findall(svg_text)}


class TestDrawMain:
    def test_species_overlap_node_types_use_descendant_species_once(self):
        tree = Tree('(((Homo_sapiens_G1:1,Homo_sapiens_G2:1):1,Mus_musculus_G1:1):1,Danio_rerio_G1:1);', parser=1)
        args = make_draw_args()

        node_type_by_node, parsed = _get_species_overlap_node_types(
            tree=tree,
            args=args,
            require_all_tip_labels=True,
        )

        assert parsed is True
        assert node_type_by_node[tree.common_ancestor(['Homo_sapiens_G1', 'Homo_sapiens_G2'])] == 'duplication'
        assert node_type_by_node[tree.common_ancestor(['Homo_sapiens_G1', 'Mus_musculus_G1'])] == 'speciation'
        assert node_type_by_node[tree] == 'speciation'

    def test_draw_writes_svg_with_species_overlap_legend(self, tmp_nwk, tmp_path):
        nwk = '((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,(Homo_sapiens_G2:1,Danio_rerio_G1:1):1);'
        infile = tmp_nwk(nwk)
        outfile = tmp_path / 'tree.svg'
        args = make_draw_args(infile=infile, outfile=str(outfile))

        draw_main(args)

        assert outfile.exists()
        text = outfile.read_text(encoding='utf-8')
        assert 'Speciation node' in text
        assert 'Duplication node' in text
        assert 'Homo_sapiens_G1' in text
        assert 'Helvetica' in text
        assert 'width="259.2pt"' in text

    def test_draw_writes_support_values_on_branches(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1)0.95:1,(C:1,D:1)0.88:1);')
        outfile = tmp_path / 'tree.svg'
        args = make_draw_args(infile=infile, outfile=str(outfile), format='0', species_overlap_node_plot='no')

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        positions = extract_svg_text_positions(text)
        assert '0.95' in positions
        assert '0.88' in positions

    def test_draw_skips_species_overlap_legend_when_labels_do_not_match_regex(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        outfile = tmp_path / 'tree.svg'
        args = make_draw_args(infile=infile, outfile=str(outfile))

        draw_main(args)

        assert outfile.exists()
        text = outfile.read_text(encoding='utf-8')
        assert 'Speciation node' not in text
        assert 'Duplication node' not in text
        assert 'A</text>' in text

    def test_draw_skips_species_overlap_legend_for_unrooted_tree_in_auto_mode(self, tmp_nwk, tmp_path, capsys):
        infile = tmp_nwk('(Homo_sapiens_G1:1,Mus_musculus_G1:1,Danio_rerio_G1:1);')
        outfile = tmp_path / 'unrooted_auto.svg'
        args = make_draw_args(infile=infile, outfile=str(outfile), species_overlap_node_plot='auto')

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        captured = capsys.readouterr()
        assert 'Speciation node' not in text
        assert 'Duplication node' not in text
        assert 'unrooted' in captured.err

    def test_draw_raises_for_unrooted_tree_in_yes_mode(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(Homo_sapiens_G1:1,Mus_musculus_G1:1,Danio_rerio_G1:1);')
        outfile = tmp_path / 'unrooted_yes.svg'
        args = make_draw_args(infile=infile, outfile=str(outfile), species_overlap_node_plot='yes')

        with pytest.raises(ValueError, match='requires a rooted tree'):
            draw_main(args)

    def test_draw_handles_topology_without_branch_lengths(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A,B),(C,D));')
        outfile = tmp_path / 'topology.svg'
        args = make_draw_args(
            infile=infile,
            outfile=str(outfile),
            species_overlap_node_plot='no',
        )

        draw_main(args)

        assert outfile.exists()
        assert outfile.stat().st_size > 0

    def test_draw_uses_tight_tip_label_spacing(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        outfile = tmp_path / 'spacing.svg'
        args = make_draw_args(infile=infile, outfile=str(outfile), species_overlap_node_plot='no')

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        positions = extract_svg_text_positions(text)
        assert abs(positions['B'] - positions['A'] - 8.5) < 0.2
        assert abs(positions['C'] - positions['B'] - 8.5) < 0.2
        assert abs(positions['D'] - positions['C'] - 8.5) < 0.2

    def test_draw_colors_tip_labels_from_trait_table(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame({'leaf_name': ['A', 'B', 'C'], 'group': ['x', 'x', 'y']}).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'trait.svg'
        args = make_draw_args(
            infile=infile,
            outfile=str(outfile),
            species_overlap_node_plot='no',
            trait=str(trait_path),
            group_by='group',
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert 'x' in text
        assert 'y' in text
        assert '#1f77b4' in text

    def test_draw_can_filter_support_labels(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1)0.95:1,(C:1,D:1)0.88:1);')
        outfile = tmp_path / 'support_filter.svg'
        args = make_draw_args(
            infile=infile,
            outfile=str(outfile),
            format='0',
            species_overlap_node_plot='no',
            support_min=0.9,
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert '0.95' in text
        assert '0.88' not in text

    def test_draw_rejects_unknown_trait_leaf_names(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:1);')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame({'leaf_name': ['A', 'Z'], 'group': ['x', 'y']}).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'trait.svg'
        args = make_draw_args(
            infile=infile,
            outfile=str(outfile),
            species_overlap_node_plot='no',
            trait=str(trait_path),
            group_by='group',
            unmatched='error',
        )

        with pytest.raises(ValueError, match='--trait and tree tips differ'):
            draw_main(args)

    def test_draw_accepts_numeric_leaf_names_in_trait_table(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((1:1,2:1):1,3:1);')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame({'leaf_name': [1, 2], 'group': ['x', 'y']}).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'trait_numeric.svg'
        args = make_draw_args(
            infile=infile,
            outfile=str(outfile),
            species_overlap_node_plot='no',
            trait=str(trait_path),
            group_by='group',
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert '1</text>' in text
        assert '2</text>' in text
        assert 'x' in text
        assert 'y' in text

    def test_draw_accepts_leading_zero_leaf_names_in_trait_table(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((001:1,002:1):1,003:1);')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame({'leaf_name': ['001', '002'], 'group': ['x', 'y']}).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'trait_leading_zero.svg'
        args = make_draw_args(
            infile=infile,
            outfile=str(outfile),
            species_overlap_node_plot='no',
            trait=str(trait_path),
            group_by='group',
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert '001</text>' in text
        assert '002</text>' in text
        assert 'x' in text
        assert 'y' in text

    def test_draw_accepts_na_literal_leaf_names_in_trait_table(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((NA:1,B:1):1,C:1);')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame({'leaf_name': ['NA', 'B'], 'group': ['x', 'y']}).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'trait_na_literal.svg'
        args = make_draw_args(
            infile=infile,
            outfile=str(outfile),
            species_overlap_node_plot='no',
            trait=str(trait_path),
            group_by='group',
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert 'NA</text>' in text
        assert 'B</text>' in text
        assert 'x' in text
        assert 'y' in text
