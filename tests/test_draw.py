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
        'figure_width': 3.6,
        'figure_height': None,
        'label_panel_width': None,
        'tip_image_manifest': None,
        'tip_image_root': None,
        'tip_image_size': 18.0,
        'tip_image_gap': 4.0,
        'font_size': 8.0,
        'font_family': 'Helvetica',
        'branch_color': '#000000',
        'terminal_branch_color': None,
        'tip_label_position': 'aligned',
        'root_marker': 'none',
        'root_marker_color': '#0072B2',
        'root_marker_size': None,
        'tip_badge_property': None,
        'tip_badge_missing_label': None,
        'node_pie_properties': None,
        'node_pie_target': 'root,intnode',
        'node_pie_leaf_filter': [],
        'node_label_property': None,
        'node_label_target': 'intnode',
        'node_label_filter': [],
        'node_label_decimals': 2,
        'node_label_prefix': '',
        'property_color': [],
        'legend': True,
        'transparent': False,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


def write_test_png(path, color=(20, 120, 200, 255)):
    Image = pytest.importorskip('PIL.Image')
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.new('RGBA', (12, 8), color).save(path)


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

    def test_draw_uses_explicit_panel_dimensions_and_root_marker(self, tmp_nwk, tmp_path, monkeypatch):
        infile = tmp_nwk('((A:1,B:1):1,C:2);')
        outfile = tmp_path / 'panel.svg'
        from matplotlib.axes import Axes

        marker_areas = []
        original_scatter = Axes.scatter

        def record_scatter(axis, *args, **kwargs):
            marker_areas.append(kwargs.get('s'))
            return original_scatter(axis, *args, **kwargs)

        monkeypatch.setattr(Axes, 'scatter', record_scatter)
        args = make_draw_args(
            infile=infile,
            outfile=str(outfile),
            species_overlap_node_plot='no',
            figure_width=1.2,
            figure_height=2.0,
            label_panel_width=0.3,
            root_marker='diamond',
            root_marker_color='#0072B2',
            root_marker_size=4.0,
            transparent=True,
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert 'width="86.4pt"' in text
        assert 'height="144pt"' in text
        assert '#0072b2' in text.lower()
        assert marker_areas == [16.0]

    def test_draw_rejects_nonpositive_root_marker_size(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:2);')
        args = make_draw_args(
            infile=infile,
            outfile=str(tmp_path / 'panel.svg'),
            species_overlap_node_plot='no',
            root_marker='diamond',
            root_marker_size=0.0,
        )

        with pytest.raises(ValueError, match='root-marker-size'):
            draw_main(args)

    def test_draw_displays_nhx_tip_badges(self, tmp_path):
        tree = Tree('((A:1,B:1):1,C:2);', parser=1)
        tree['A'].add_props(state='X')
        tree['B'].add_props(state='Y')
        tree['C'].add_props(state='X')
        infile = tmp_path / 'annotated.nhx'
        infile.write_text(
            tree.write(props=['state'], parser=1, format_root_node=True),
            encoding='utf-8',
        )
        outfile = tmp_path / 'badges.svg'
        args = make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            tip_badge_property='state',
            property_color=['X=#BFD7EA', 'Y=#F7D59C'],
            legend=False,
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert '>X</text>' in text
        assert '>Y</text>' in text
        assert '#bfd7ea' in text.lower()
        assert '#f7d59c' in text.lower()

    def test_draw_can_label_missing_tip_properties(self, tmp_path):
        tree = Tree('((A:1,B:1):1,C:2);', parser=1)
        tree['A'].add_props(state='X')
        tree['B'].add_props(state='')
        infile = tmp_path / 'partial.nhx'
        infile.write_text(
            tree.write(props=['state'], parser=1, format_root_node=True),
            encoding='utf-8',
        )
        outfile = tmp_path / 'missing_badges.svg'
        args = make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            tip_badge_property='state',
            tip_badge_missing_label='missing',
            property_color=['state:X=#BFD7EA', 'state:missing=#F2F2F2'],
            legend=False,
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert text.count('>missing</text>') == 2
        assert '#f2f2f2' in text.lower()

    def test_draw_displays_filtered_nhx_probability_nodes(self, tmp_path):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        ab = tree.common_ancestor(['A', 'B'])
        cd = tree.common_ancestor(['C', 'D'])
        for node, probability_y in ((tree, 0.55), (ab, 0.80), (cd, 0.60)):
            node.add_props(
                asr_p_X=1.0 - probability_y,
                asr_p_Y=probability_y,
                asr_probability=max(probability_y, 1.0 - probability_y),
            )
        infile = tmp_path / 'asr.nhx'
        properties = ['asr_p_X', 'asr_p_Y', 'asr_probability']
        infile.write_text(
            tree.write(props=properties, parser=1, format_root_node=True),
            encoding='utf-8',
        )
        outfile = tmp_path / 'asr.svg'
        args = make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            node_pie_properties='asr_p_X,asr_p_Y',
            node_pie_target='root,intnode',
            node_label_property='asr_p_Y',
            node_label_target='intnode',
            node_label_filter=['asr_probability:ge:0.68'],
            node_label_prefix='P(Y)=',
            property_color=['X=#56B4E9', 'Y=#E69F00'],
            legend=False,
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8')
        assert '>P(Y)=0.80</text>' in text
        assert '>P(Y)=0.60</text>' not in text
        assert '#56b4e9' in text.lower()
        assert '#e69f00' in text.lower()

    def test_draw_leaf_pie_filter_keeps_internal_and_matching_leaf_pies(self, tmp_path):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        for node in tree.traverse():
            node.add_props(asr_p_X=0.4, asr_p_Y=0.6)
        for leaf in tree.leaves():
            leaf.add_props(asr_observed_state='X')
        tree['D'].add_props(asr_observed_state='')
        infile = tmp_path / 'leaf_filter.nhx'
        infile.write_text(
            tree.write(
                props=['asr_p_X', 'asr_p_Y', 'asr_observed_state'],
                parser=1,
                format_root_node=True,
            ),
            encoding='utf-8',
        )
        outfile = tmp_path / 'leaf_filter.svg'
        args = make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            node_pie_properties='asr_p_X,asr_p_Y',
            node_pie_target='all',
            node_pie_leaf_filter=['asr_observed_state:eq:'],
            property_color=['X=#56B4E9', 'Y=#E69F00'],
            legend=False,
        )

        draw_main(args)

        text = outfile.read_text(encoding='utf-8').lower()
        assert text.count('#56b4e9') == 4
        assert text.count('#e69f00') == 4

    def test_draw_embeds_manifest_tip_images_and_expands_auto_height(self, tmp_path, capsys):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('((A:1,B:1):1,C:2);')
        image_dir = tmp_path / 'assets' / 'images'
        write_test_png(image_dir / 'a.png', color=(220, 60, 60, 255))
        write_test_png(image_dir / 'b.png', color=(60, 180, 80, 255))
        write_test_png(image_dir / 'c.png', color=(60, 90, 220, 255))
        manifest = tmp_path / 'assets' / 'manifest.tsv'
        manifest.write_text(
            'leaf_name\tlocal_path\n'
            'A\timages/a.png\n'
            'B\timages/b.png\n'
            'C\timages/c.png\n'
        )
        plain_outfile = tmp_path / 'plain.svg'
        image_outfile = tmp_path / 'with-images.svg'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(plain_outfile),
            species_overlap_node_plot='no',
        ))
        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(image_outfile),
            species_overlap_node_plot='no',
            tip_image_manifest=str(manifest),
            tip_image_size=24.0,
        ))

        plain_text = plain_outfile.read_text(encoding='utf-8')
        image_text = image_outfile.read_text(encoding='utf-8')
        plain_height = float(re.search(r'height="([0-9.]+)pt"', plain_text).group(1))
        image_height = float(re.search(r'height="([0-9.]+)pt"', image_text).group(1))
        assert image_text.count('<image') == 3
        assert image_height > plain_height
        assert 'A</text>' in image_text
        assert 'Loaded tip images for 3 tree tip(s)' in capsys.readouterr().err

    def test_draw_resolves_tip_images_from_explicit_root(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('(A:1,B:1);')
        image_root = tmp_path / 'image-root'
        write_test_png(image_root / 'a.png')
        write_test_png(image_root / 'b.png')
        manifest_dir = tmp_path / 'metadata'
        manifest_dir.mkdir()
        manifest = manifest_dir / 'manifest.tsv'
        manifest.write_text(
            'leaf_name\tlocal_path\n'
            'A\ta.png\n'
            'B\tb.png\n'
        )
        outfile = tmp_path / 'rooted-assets.svg'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            tip_image_manifest=str(manifest),
            tip_image_root=str(image_root),
        ))

        assert outfile.read_text(encoding='utf-8').count('<image') == 2

    def test_draw_rejects_duplicate_tip_image_rows(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('(A:1,B:1);')
        image = tmp_path / 'image.png'
        write_test_png(image)
        manifest = tmp_path / 'manifest.tsv'
        manifest.write_text(
            'leaf_name\tlocal_path\n'
            'A\timage.png\n'
            'A\timage.png\n'
        )

        with pytest.raises(ValueError, match="Duplicated 'leaf_name'"):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'duplicate.svg'),
                species_overlap_node_plot='no',
                tip_image_manifest=str(manifest),
            ))

    def test_draw_applies_unmatched_error_to_tip_image_manifest(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('(A:1,B:1);')
        image = tmp_path / 'image.png'
        write_test_png(image)
        manifest = tmp_path / 'manifest.tsv'
        manifest.write_text('leaf_name\tlocal_path\nA\timage.png\n')

        with pytest.raises(ValueError, match='tip-image-manifest and tree tips differ'):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'missing.svg'),
                species_overlap_node_plot='no',
                tip_image_manifest=str(manifest),
                unmatched='error',
            ))

    def test_draw_rejects_fixed_height_that_overlaps_tip_images(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('(A:1,B:1);')
        image = tmp_path / 'image.png'
        write_test_png(image)
        manifest = tmp_path / 'manifest.tsv'
        manifest.write_text(
            'leaf_name\tlocal_path\n'
            'A\timage.png\n'
            'B\timage.png\n'
        )

        with pytest.raises(ValueError, match='non-overlapping tip images'):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'short.svg'),
                species_overlap_node_plot='no',
                tip_image_manifest=str(manifest),
                tip_image_size=36.0,
                figure_height=0.5,
            ))

    def test_draw_rasterizes_svg_tip_images(self, tmp_path):
        pytest.importorskip('cairosvg')
        infile = tmp_path / 'tree.nwk'
        infile.write_text('(A:1,B:1);')
        silhouette = tmp_path / 'silhouette.svg'
        silhouette.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="40" height="20">'
            '<path d="M2 18 L20 2 L38 18 Z" fill="#111111"/>'
            '</svg>'
        )
        manifest = tmp_path / 'manifest.tsv'
        manifest.write_text(
            'leaf_name\tlocal_path\n'
            'A\tsilhouette.svg\n'
            'B\tsilhouette.svg\n'
        )
        outfile = tmp_path / 'silhouettes.svg'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            tip_image_manifest=str(manifest),
        ))

        assert outfile.read_text(encoding='utf-8').count('<image') == 2
