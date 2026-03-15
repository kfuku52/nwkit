import re
from argparse import Namespace

import pytest

from nwkit.draw import draw_main


def make_draw_args(**kwargs):
    defaults = {
        'infile': '-',
        'format': 'auto',
        'quoted_node_names': True,
        'outfile': None,
        'image_format': 'auto',
        'species_regex': r'^([^_]+_[^_]+)(?:_|$)',
        'species_overlap_node_plot': 'auto',
        'ladderize': False,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


def extract_svg_text_positions(svg_text):
    pattern = re.compile(r'<text[^>]*y="([^"]+)"[^>]*>([^<]+)</text>')
    return {content: float(y) for y, content in pattern.findall(svg_text)}


class TestDrawMain:
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
