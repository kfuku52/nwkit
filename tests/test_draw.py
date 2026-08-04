import json
import math
import re
from argparse import Namespace

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from ete4 import Tree
from matplotlib.figure import Figure

from nwkit.draw import (
    _add_scale_bar,
    _get_species_overlap_node_types,
    _load_tip_image,
    _prepare_tip_label_text,
    _read_tip_image_manifest,
    _wrap_tip_label,
    draw_main,
)
from nwkit.draw_layouts import make_tree_layout
from nwkit.draw_quality import DrawingArtist, evaluate_drawing


def make_draw_args(**kwargs):
    defaults = {
        'infile': '-',
        'format': 'auto',
        'quoted_node_names': True,
        'outfile': None,
        'image_format': 'auto',
        'layout': 'rectangular',
        'subtree_packing': 'standard',
        'spiral_turns': None,
        'fan_span': 180.0,
        'unrooted_method': 'equal-angle',
        'daylight_iterations': 5,
        'max_visible_tips': None,
        'collapse_label': None,
        'collapse_property_aggregation': 'none',
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
        'branch_width': 0.8,
        'terminal_branch_color': None,
        'branch_color_property': None,
        'branch_width_property': None,
        'branch_width_range': '0.4,2.5',
        'scale_bar': 'none',
        'depth_guide': 'none',
        'branch_length_unit': '',
        'tip_labels': True,
        'tip_label_position': 'aligned',
        'tip_label_wrap': 'none',
        'tip_spacing': 'uniform',
        'tip_label_font_style': 'plain',
        'tip_track': [],
        'tip_track_type': 'auto',
        'tip_track_size': 5.0,
        'tip_track_palette': 'viridis',
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
        'legend_columns': 'auto',
        'legend_position': 'auto',
        'collision_policy': 'resolve',
        'layout_report': None,
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
    def test_scale_bar_label_is_above_bar(self):
        figure, axes = plt.subplots(figsize=(3.0, 2.0))
        axes.set_xlim(0.0, 3.0)
        artist = _add_scale_bar(
            figure=figure,
            axes=axes,
            size=1.0,
            label='1 substitutions/site',
            color='#000000',
            font_family='DejaVu Sans',
            font_size=8.0,
            anchor_x=0.1,
            anchor_y=0.02,
        )

        try:
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            label_bounds = artist.txt_label.get_window_extent(renderer)
            bar_bounds = artist.size_bar.get_window_extent(renderer)
        finally:
            plt.close(figure)

        assert label_bounds.y0 >= bar_bounds.y1

    def test_collision_solver_moves_lower_priority_annotation(self):
        figure, axes = plt.subplots(figsize=(2.0, 2.0))
        fixed = axes.annotate('fixed', xy=(0.5, 0.5), xytext=(0.0, 0.0),
                              textcoords='offset points')
        movable = axes.annotate('movable', xy=(0.5, 0.5), xytext=(0.0, 0.0),
                                textcoords='offset points')
        artists = [
            DrawingArtist(fixed, kind='node_label', priority=100),
            DrawingArtist(
                movable,
                kind='node_label',
                priority=10,
                movable=True,
            ),
        ]

        try:
            report = evaluate_drawing(figure, artists, [], policy='resolve')
        finally:
            plt.close(figure)

        assert report['initial_collision_count'] == 1
        assert report['final_collision_count'] == 0
        assert report['moved_artist_count'] == 1

    @pytest.mark.parametrize(
        'layout',
        [
            'slanted',
            'cladogram',
            'circular',
            'fan',
            'radial',
            'unrooted',
            'spiral',
            'fractal',
        ],
    )
    def test_draw_writes_svg_with_modern_layout(self, tmp_nwk, tmp_path, layout):
        infile = tmp_nwk('(((A:1,B:1):1,C:1):1,(D:1,E:8):1);')
        outfile = tmp_path / '{}.svg'.format(layout)

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            image_format='svg',
            layout=layout,
            tip_label_position='branch-end',
            species_overlap_node_plot='no',
            figure_width=5.0,
            figure_height=4.0,
        ))

        svg = outfile.read_text()
        assert '<svg' in svg
        for leaf_name in ('A', 'B', 'C', 'D', 'E'):
            assert '>{}<'.format(leaf_name) in svg

    def test_tidy_layout_preserves_branch_axis_and_compacts_tree(self):
        tree = Tree('(((A:1,B:1):1,C:1):1,(D:1,E:8):1);', parser=1)
        rectangular = make_tree_layout(tree, layout='rectangular')
        tidy = make_tree_layout(
            tree,
            layout='rectangular',
            subtree_packing='tidy',
        )

        assert tidy.name == 'rectangular'
        assert tidy.metadata['subtree_packing'] == 'tidy'
        assert tidy.xcoord == rectangular.xcoord
        rectangular_span = max(rectangular.ycoord.values()) - min(rectangular.ycoord.values())
        tidy_span = max(tidy.ycoord.values()) - min(tidy.ycoord.values())
        assert tidy_span < rectangular_span

    def test_slanted_layout_preserves_phylogram_depth_with_straight_edges(self):
        tree = Tree('((A:1,B:2):3,C:4);', parser=1)
        rectangular = make_tree_layout(tree, layout='rectangular')
        slanted = make_tree_layout(tree, layout='slanted')

        assert slanted.xcoord == rectangular.xcoord
        assert all(len(path) == 2 for path in slanted.edge_paths.values())

    def test_cladogram_aligns_tips_and_ignores_branch_lengths(self):
        tree = Tree('(((A:1,B:9):2,C:7):3,D:4);', parser=1)
        cladogram = make_tree_layout(tree, layout='cladogram')

        tip_positions = {cladogram.xcoord[leaf] for leaf in tree.leaves()}
        assert len(tip_positions) == 1
        assert cladogram.xcoord[tree] == 0.0
        assert all(len(path) == 2 for path in cladogram.edge_paths.values())

    @pytest.mark.parametrize('layout', ['circular', 'fan'])
    def test_circular_layouts_use_arc_then_radial_paths(self, layout):
        tree = Tree('(((A:1,B:1):1,C:1):1,(D:1,E:1):1);', parser=1)
        drawing = make_tree_layout(tree, layout=layout)

        assert drawing.spatial is True
        assert drawing.equal_aspect is True
        assert any(len(path) > 3 for path in drawing.edge_paths.values())

    def test_fan_span_controls_tip_angular_span(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        fan = make_tree_layout(tree, layout='fan', fan_span=60.0)

        angles = sorted(fan.label_angles.values())
        assert angles[-1] - angles[0] == pytest.approx(60.0)

    def test_fan_defaults_to_right_facing_semicircle(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        fan = make_tree_layout(tree, layout='fan')

        angles = sorted(fan.label_angles.values())
        assert angles[0] == pytest.approx(-90.0)
        assert angles[-1] == pytest.approx(90.0)

    def test_fan_span_360_matches_circular_geometry(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        circular = make_tree_layout(tree, layout='circular')
        fan = make_tree_layout(tree, layout='fan', fan_span=360.0)

        for node in tree.traverse():
            assert fan.xcoord[node] == pytest.approx(circular.xcoord[node])
            assert fan.ycoord[node] == pytest.approx(circular.ycoord[node])

    def test_radial_layout_uses_straight_edges(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        radial = make_tree_layout(tree, layout='radial')

        assert radial.spatial is True
        assert all(len(path) == 2 for path in radial.edge_paths.values())

    def test_unrooted_layout_suppresses_degree_two_root_on_joined_edge(self):
        tree = Tree('((A:1,B:1):2,(C:1,D:1):4);', parser=1)
        unrooted = make_tree_layout(tree, layout='unrooted')
        first, second = tree.get_children()

        joined_length = math.hypot(
            unrooted.xcoord[second] - unrooted.xcoord[first],
            unrooted.ycoord[second] - unrooted.ycoord[first],
        )
        first_to_root = math.hypot(
            unrooted.xcoord[tree] - unrooted.xcoord[first],
            unrooted.ycoord[tree] - unrooted.ycoord[first],
        )
        assert joined_length == pytest.approx(6.0)
        assert first_to_root == pytest.approx(2.0)
        assert unrooted.root_path == []

    def test_equal_daylight_unrooted_layout_preserves_every_branch_length(self):
        tree = Tree(
            '(((A:1,B:2):1,C:3):2,((D:1,E:2):1,(F:1,G:1):2):4);',
            parser=1,
        )
        equal_angle = make_tree_layout(tree, layout='unrooted')
        daylight = make_tree_layout(
            tree,
            layout='unrooted',
            unrooted_method='equal-daylight',
            daylight_iterations=8,
        )

        assert daylight.name == 'unrooted-daylight'
        assert any(
            daylight.xcoord[node] != pytest.approx(equal_angle.xcoord[node])
            or daylight.ycoord[node] != pytest.approx(equal_angle.ycoord[node])
            for node in tree.traverse()
        )
        for node in tree.traverse():
            if node.is_root:
                continue
            distance = math.hypot(
                daylight.xcoord[node] - daylight.xcoord[node.up],
                daylight.ycoord[node] - daylight.ycoord[node.up],
            )
            assert distance == pytest.approx(float(node.dist))

    @pytest.mark.parametrize('fan_span', [0.0, 360.1])
    def test_fan_layout_rejects_invalid_span(self, fan_span):
        tree = Tree('(A:1,B:1);', parser=1)

        with pytest.raises(ValueError, match='--fan-span'):
            make_tree_layout(tree, layout='fan', fan_span=fan_span)

    @pytest.mark.parametrize('layout_name', ['circular', 'radial', 'fractal'])
    def test_label_aware_spatial_layout_allocates_more_angle_to_tall_labels(
        self,
        layout_name,
    ):
        tree = Tree('(A:1,B:1,C:1);', parser=1)
        sizes = {
            tree['A']: (0.4, 0.1),
            tree['B']: (0.4, 1.0),
            tree['C']: (0.4, 0.1),
        }
        drawing = make_tree_layout(
            tree,
            layout=layout_name,
            aspect_ratio=1.4,
            label_size_by_leaf=sizes,
            tip_spacing='label-aware',
        )

        angles = [
            math.radians(drawing.label_angles[tree[name]]) % (2.0 * math.pi)
            for name in ('A', 'B', 'C')
        ]
        gaps = [
            (angles[(index + 1) % 3] - angles[index]) % (2.0 * math.pi)
            for index in range(3)
        ]
        assert drawing.name == layout_name
        assert drawing.spatial is True
        assert gaps[0] > gaps[2]
        assert gaps[1] > gaps[2]

    @pytest.mark.parametrize('layout_name', ['rectangular', 'slanted', 'cladogram'])
    def test_label_aware_cartesian_layout_uses_variable_tip_rows(
        self,
        layout_name,
    ):
        tree = Tree('(A:1,B:1,C:1);', parser=1)
        sizes = {
            tree['A']: (0.4, 0.1),
            tree['B']: (0.4, 1.0),
            tree['C']: (0.4, 0.1),
        }
        uniform = make_tree_layout(tree, layout=layout_name)
        aware = make_tree_layout(
            tree,
            layout=layout_name,
            label_size_by_leaf=sizes,
            tip_spacing='label-aware',
        )

        assert aware.ycoord[tree['B']] - aware.ycoord[tree['A']] > (
            uniform.ycoord[tree['B']] - uniform.ycoord[tree['A']]
        )
        assert aware.ycoord[tree['C']] - aware.ycoord[tree['B']] > (
            uniform.ycoord[tree['C']] - uniform.ycoord[tree['B']]
        )

    def test_label_aware_tidy_layout_uses_variable_leaf_boxes(self):
        tree = Tree('(A:1,B:1,C:1);', parser=1)
        uniform = make_tree_layout(
            tree,
            layout='rectangular',
            subtree_packing='tidy',
        )
        aware = make_tree_layout(
            tree,
            layout='rectangular',
            subtree_packing='tidy',
            label_size_by_leaf={
                tree['A']: (0.4, 0.1),
                tree['B']: (0.4, 1.0),
                tree['C']: (0.4, 0.1),
            },
            tip_spacing='label-aware',
        )

        uniform_span = max(uniform.ycoord.values()) - min(uniform.ycoord.values())
        aware_span = max(aware.ycoord.values()) - min(aware.ycoord.values())
        assert aware_span > uniform_span

    def test_tidy_multiline_label_reserves_height_only_at_terminal_extent(self):
        def drawing_for(newick):
            tree = Tree(newick, parser=1)
            drawing = make_tree_layout(
                tree,
                layout='rectangular',
                subtree_packing='tidy',
                label_size_by_leaf={
                    tree['A']: (2.0, 2.0),
                    tree['B']: (0.2, 0.1),
                    tree['C']: (0.2, 0.1),
                },
                terminal_extent_by_leaf={
                    tree['A']: 2.0,
                    tree['B']: 0.2,
                    tree['C']: 0.2,
                },
                tip_spacing='label-aware',
            )
            return tree, drawing

        separated_tree, separated = drawing_for('(A:10,B:1,C:1);')
        overlapping_tree, overlapping = drawing_for('(A:1,B:1,C:1);')
        separated_gap = (
            separated.ycoord[separated_tree['B']]
            - separated.ycoord[separated_tree['A']]
        )
        overlapping_gap = (
            overlapping.ycoord[overlapping_tree['B']]
            - overlapping.ycoord[overlapping_tree['A']]
        )

        assert separated_gap < overlapping_gap
        assert separated_gap == pytest.approx(1.0)

    @pytest.mark.parametrize('tip_spacing', ['uniform', 'label-aware'])
    def test_tidy_parents_are_centered_on_direct_children(self, tip_spacing):
        tree = Tree('(((A:1,B:2):1,C:1):1,(D:1,(E:2,F:1):1):2);', parser=1)
        drawing = make_tree_layout(
            tree,
            layout='rectangular',
            subtree_packing='tidy',
            label_size_by_leaf={
                leaf: (
                    0.4,
                    1.2 if leaf.name == 'B' else (
                        0.6 if leaf.name == 'E' else 0.1
                    ),
                )
                for leaf in tree.leaves()
            },
            tip_spacing=tip_spacing,
        )

        for node in tree.traverse():
            children = node.get_children()
            if not children:
                continue
            expected = sum(drawing.ycoord[child] for child in children) / len(children)
            assert drawing.ycoord[node] == pytest.approx(expected)
            assert min(drawing.ycoord[child] for child in children) <= drawing.ycoord[node]
            assert drawing.ycoord[node] <= max(
                drawing.ycoord[child] for child in children
            )

    def test_label_aware_circular_layout_preserves_root_distance_as_radius(self):
        tree = Tree('(((A:1,B:2):3,C:2):1,(D:5,E:1):2);', parser=1)
        rectangular = make_tree_layout(tree, layout='rectangular')
        circular = make_tree_layout(
            tree,
            layout='circular',
            label_size_by_leaf={
                leaf: (0.5, 0.1 + index * 0.1)
                for index, leaf in enumerate(tree.leaves())
            },
            tip_spacing='label-aware',
        )

        assert circular.name == 'circular'
        for node in tree.traverse():
            radius = math.hypot(circular.xcoord[node], circular.ycoord[node])
            assert radius == pytest.approx(rectangular.xcoord[node])

    @pytest.mark.parametrize(
        'layout_name',
        [
            'rectangular',
            'slanted',
            'cladogram',
            'circular',
            'fan',
            'radial',
            'unrooted',
            'spiral',
            'fractal',
        ],
    )
    def test_every_layout_accepts_label_aware_tip_spacing(self, layout_name):
        tree = Tree('((A:1,B:2):1,(C:1,D:3):2);', parser=1)
        drawing = make_tree_layout(
            tree,
            layout=layout_name,
            label_size_by_leaf={
                leaf: (0.4, 0.1 + index * 0.25)
                for index, leaf in enumerate(tree.leaves())
            },
            tip_spacing='label-aware',
        )

        assert drawing.leaf_order
        assert set(drawing.xcoord) == set(tree.traverse())
        assert set(drawing.ycoord) == set(tree.traverse())

    @pytest.mark.parametrize('layout_name', ['rectangular', 'circular', 'fan', 'spiral'])
    def test_tidy_subtree_packing_composes_with_supported_layouts(
        self,
        layout_name,
    ):
        tree = Tree('(((A:1,B:2):1,C:1):1,(D:1,(E:2,F:1):1):2);', parser=1)
        drawing = make_tree_layout(
            tree,
            layout=layout_name,
            subtree_packing='tidy',
            label_size_by_leaf={
                leaf: (0.4, 1.0 if leaf.name == 'B' else 0.1)
                for leaf in tree.leaves()
            },
            terminal_extent_by_leaf={leaf: 0.4 for leaf in tree.leaves()},
            tip_spacing='label-aware',
        )

        assert drawing.name == layout_name
        assert drawing.metadata['subtree_packing'] == 'tidy'
        assert drawing.leaf_order == list(tree.leaves())
        assert all(math.isfinite(value) for value in drawing.xcoord.values())
        assert all(math.isfinite(value) for value in drawing.ycoord.values())

    @pytest.mark.parametrize(
        'layout_name',
        ['slanted', 'cladogram', 'radial', 'unrooted', 'fractal'],
    )
    def test_tidy_subtree_packing_rejects_unverified_geometries(
        self,
        layout_name,
    ):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)

        with pytest.raises(ValueError, match='supported only with rectangular'):
            make_tree_layout(
                tree,
                layout=layout_name,
                subtree_packing='tidy',
            )

    def test_tidy_circular_packing_reserves_the_wraparound_seam(self):
        tree = Tree('(((A:1,B:1):1,C:1):1,(D:1,(E:1,F:1):1):1);', parser=1)
        drawing = make_tree_layout(
            tree,
            layout='circular',
            subtree_packing='tidy',
            label_size_by_leaf={
                leaf: (0.4, 1.5 if leaf.name in {'A', 'F'} else 0.1)
                for leaf in tree.leaves()
            },
            terminal_extent_by_leaf={leaf: 0.4 for leaf in tree.leaves()},
            tip_spacing='label-aware',
        )
        angles = [
            math.radians(drawing.label_angles[leaf]) % (2.0 * math.pi)
            for leaf in drawing.leaf_order
        ]
        cyclic_gaps = [
            (angles[(index + 1) % len(angles)] - angles[index]) % (2.0 * math.pi)
            for index in range(len(angles))
        ]

        assert min(cyclic_gaps) > 0.0
        assert cyclic_gaps[-1] > min(cyclic_gaps[1:-1])

    @pytest.mark.parametrize('layout_name', ['rectangular', 'circular', 'fan', 'spiral'])
    def test_tidy_subtree_packing_is_deterministic(self, layout_name):
        tree = Tree('(((A:1,B:2):1,C:1):1,(D:1,(E:2,F:1):1):2);', parser=1)
        kwargs = {
            'layout': layout_name,
            'subtree_packing': 'tidy',
            'label_size_by_leaf': {
                leaf: (0.4, 1.0 if leaf.name == 'B' else 0.1)
                for leaf in tree.leaves()
            },
            'terminal_extent_by_leaf': {leaf: 0.4 for leaf in tree.leaves()},
            'tip_spacing': 'label-aware',
        }

        first = make_tree_layout(tree, **kwargs)
        second = make_tree_layout(tree, **kwargs)

        assert first.xcoord == second.xcoord
        assert first.ycoord == second.ycoord
        assert first.edge_paths == second.edge_paths

    def test_label_aware_spiral_changes_tip_allocation(self):
        tree = Tree('(A:1,B:1,C:1,D:1);', parser=1)
        sizes = {
            tree['A']: (0.4, 0.1),
            tree['B']: (0.4, 1.0),
            tree['C']: (0.4, 0.1),
            tree['D']: (0.4, 0.1),
        }
        uniform = make_tree_layout(tree, layout='spiral')
        aware = make_tree_layout(
            tree,
            layout='spiral',
            label_size_by_leaf=sizes,
            tip_spacing='label-aware',
        )

        assert any(
            aware.xcoord[leaf] != pytest.approx(uniform.xcoord[leaf])
            or aware.ycoord[leaf] != pytest.approx(uniform.ycoord[leaf])
            for leaf in tree.leaves()
        )

    def test_label_aware_unrooted_changes_angles_and_preserves_edge_lengths(self):
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):2);', parser=1)
        uniform = make_tree_layout(tree, layout='unrooted')
        aware = make_tree_layout(
            tree,
            layout='unrooted',
            label_size_by_leaf={
                leaf: (0.4, 1.2 if leaf.name == 'D' else 0.1)
                for leaf in tree.leaves()
            },
            tip_spacing='label-aware',
        )

        assert any(
            aware.xcoord[leaf] != pytest.approx(uniform.xcoord[leaf])
            or aware.ycoord[leaf] != pytest.approx(uniform.ycoord[leaf])
            for leaf in tree.leaves()
        )
        uniform_center = min(
            tree.traverse(),
            key=lambda node: math.hypot(
                uniform.xcoord[node],
                uniform.ycoord[node],
            ),
        )
        aware_center = min(
            tree.traverse(),
            key=lambda node: math.hypot(
                aware.xcoord[node],
                aware.ycoord[node],
            ),
        )
        assert aware_center is uniform_center
        for node, path in aware.edge_paths.items():
            length = sum(
                math.hypot(end[0] - start[0], end[1] - start[1])
                for start, end in zip(path, path[1:])
            )
            assert length == pytest.approx(float(node.dist))

    @pytest.mark.parametrize('removed_layout', ['packed', 'packed-phylogram', 'tidy'])
    def test_removed_packed_layout_names_have_no_compatibility_alias(
        self,
        removed_layout,
    ):
        tree = Tree('(A:1,B:1);', parser=1)

        with pytest.raises(ValueError, match="Unsupported '--layout'"):
            make_tree_layout(tree, layout=removed_layout)

    def test_invalid_tip_spacing_is_rejected(self):
        tree = Tree('(A:1,B:1);', parser=1)

        with pytest.raises(ValueError, match='--tip-spacing'):
            make_tree_layout(tree, tip_spacing='packed')

    def test_invalid_subtree_packing_is_rejected(self):
        tree = Tree('(A:1,B:1);', parser=1)

        with pytest.raises(ValueError, match='--subtree-packing'):
            make_tree_layout(tree, subtree_packing='dense')

    def test_tip_label_wrap_prefers_delimiters_and_hard_wraps(self):
        assert _wrap_tip_label('Arabidopsis_thaliana', 12) == 'Arabidopsis_\nthaliana'
        assert _wrap_tip_label('ABCDEFGHIJ', 4) == 'ABCD\nEFGH\nIJ'

    def test_auto_tip_label_wraps_long_names_but_leaves_short_names(self):
        tree = Tree('(Short:1,VeryLongSpeciesNameWithoutDelimiter:1);', parser=1)
        text_by_leaf, size_by_leaf = _prepare_tip_label_text(
            leaf_order=list(tree.leaves()),
            wrap='auto',
            font_size=8.0,
            font_family='DejaVu Sans',
            layout_name='fractal',
            panel_width_in=7.2,
            panel_height_in=5.0,
        )

        assert text_by_leaf[tree['Short']] == 'Short'
        assert '\n' in text_by_leaf[tree['VeryLongSpeciesNameWithoutDelimiter']]
        assert tree['VeryLongSpeciesNameWithoutDelimiter'].name == (
            'VeryLongSpeciesNameWithoutDelimiter'
        )
        assert (
            size_by_leaf[tree['VeryLongSpeciesNameWithoutDelimiter']][0]
            < 1.0
        )

    def test_taxonomy_wrap_keeps_binomial_together(self):
        tree = Tree('(Arabidopsis_thaliana_accession_Col_0:1);', parser=1)
        text_by_leaf, _ = _prepare_tip_label_text(
            leaf_order=list(tree.leaves()),
            wrap='taxonomy',
            font_size=8.0,
            font_family='DejaVu Sans',
            layout_name='fractal',
            panel_width_in=1.0,
            panel_height_in=1.0,
        )

        displayed = text_by_leaf[tree['Arabidopsis_thaliana_accession_Col_0']]
        assert displayed.startswith('Arabidopsis_thaliana\n')

    def test_draw_rejects_invalid_tip_label_wrap(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(A:1,B:1);')

        with pytest.raises(ValueError, match='--tip-label-wrap'):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'invalid-wrap.svg'),
                species_overlap_node_plot='no',
                tip_label_wrap='zero',
            ))

    def test_spiral_layout_warps_connectors_into_sampled_curves(self):
        tree = Tree('(((A:1,B:1):1,C:1):1,(D:1,E:8):1);', parser=1)
        spiral = make_tree_layout(tree, layout='spiral', spiral_turns=2.5)

        assert spiral.spatial is True
        assert spiral.equal_aspect is True
        assert any(len(path) > 3 for path in spiral.edge_paths.values())
        assert all(math.isfinite(value) for value in spiral.xcoord.values())
        assert all(math.isfinite(value) for value in spiral.ycoord.values())

    def test_fractal_layout_fits_requested_rectangular_aspect(self):
        tree = Tree('(((A:1,B:1):1,C:1):1,(D:1,E:8):1);', parser=1)
        fractal = make_tree_layout(tree, layout='fractal', aspect_ratio=1.8)

        node_x_span = max(fractal.xcoord.values()) - min(fractal.xcoord.values())
        node_y_span = max(fractal.ycoord.values()) - min(fractal.ycoord.values())
        assert node_x_span / node_y_span == pytest.approx(1.8)

    def test_spiral_layout_rejects_nonpositive_turn_count(self):
        tree = Tree('(A:1,B:1);', parser=1)

        with pytest.raises(ValueError, match='--spiral-turns'):
            make_tree_layout(tree, layout='spiral', spiral_turns=0)

    def test_draw_can_hide_tip_labels_for_dense_overview(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((Alpha:1,Beta:1):1,(Gamma:1,Delta:1):1);')
        outfile = tmp_path / 'unlabelled-fractal.svg'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            image_format='svg',
            layout='fractal',
            tip_labels=False,
            tip_label_position='branch-end',
            species_overlap_node_plot='no',
        ))

        svg = outfile.read_text()
        assert 'Alpha' not in svg
        assert 'Beta' not in svg

    def test_draw_rejects_unreasonably_large_tip_image_size(
        self,
        tmp_nwk,
        tmp_path,
    ):
        infile = tmp_nwk('(A:1,B:1);')

        with pytest.raises(
            ValueError,
            match='--tip-image-size must be no greater than 512 points',
        ):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'oversized.svg'),
                species_overlap_node_plot='no',
                tip_image_size=1e308,
            ))

    def test_draw_closes_figure_when_save_fails(
        self,
        monkeypatch,
        tmp_nwk,
        tmp_path,
    ):
        infile = tmp_nwk('(A:1,B:1);')
        open_figures = set(plt.get_fignums())

        def fail_save(*args, **kwargs):
            raise OSError('simulated renderer failure')

        monkeypatch.setattr(Figure, 'savefig', fail_save)

        with pytest.raises(OSError, match='simulated renderer failure'):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'failed.svg'),
                species_overlap_node_plot='no',
            ))

        assert set(plt.get_fignums()) == open_figures

    def test_draw_restores_process_wide_matplotlib_settings(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(A:1,B:1);')
        outfile = tmp_path / 'isolated.svg'
        keys = ('font.size', 'font.family', 'font.sans-serif', 'svg.fonttype', 'svg.hashsalt')
        before = {
            key: list(matplotlib.rcParams[key])
            if isinstance(matplotlib.rcParams[key], list)
            else matplotlib.rcParams[key]
            for key in keys
        }

        draw_main(make_draw_args(
            infile=infile,
            outfile=str(outfile),
            species_overlap_node_plot='no',
            font_size=17.0,
            font_family='DejaVu Sans',
        ))

        after = {
            key: list(matplotlib.rcParams[key])
            if isinstance(matplotlib.rcParams[key], list)
            else matplotlib.rcParams[key]
            for key in keys
        }
        assert after == before

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

    def test_draw_writes_scale_bar_and_layout_quality_report(
        self,
        tmp_nwk,
        tmp_path,
    ):
        infile = tmp_nwk('((A:1,B:2):1,(C:3,D:4):2);')
        outfile = tmp_path / 'scaled.svg'
        report_path = tmp_path / 'layout.json'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            image_format='svg',
            species_overlap_node_plot='no',
            tip_spacing='label-aware',
            scale_bar='auto',
            branch_length_unit='substitutions/site',
            layout_report=str(report_path),
        ))

        svg = outfile.read_text(encoding='utf-8')
        report = json.loads(report_path.read_text(encoding='utf-8'))
        assert 'substitutions/site' in svg
        assert report['branch_lengths_encoded'] is True
        assert report['tip_spacing'] == 'label-aware'
        assert report['subtree_packing'] == 'standard'
        assert report['scale_bar'] > 0.0
        assert report['scale_bar_position'] == 'bottom-reserved'
        assert report['scale_bar_label_position'] == 'above'
        assert not any(
            'scale_bar' in collision_kind
            for collision_kind in report['final_collisions_by_kind']
        )
        assert report['visible_tip_count'] == 4
        assert 'final_collisions_by_kind' in report

    def test_draw_reports_composed_tidy_packing(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(((A:1,B:2):1,C:1):1,(D:1,E:2):1);')
        report_path = tmp_path / 'circular-tidy.json'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(tmp_path / 'circular-tidy.svg'),
            image_format='svg',
            layout='circular',
            subtree_packing='tidy',
            tip_label_position='branch-end',
            tip_spacing='label-aware',
            species_overlap_node_plot='no',
            layout_report=str(report_path),
        ))

        report = json.loads(report_path.read_text(encoding='utf-8'))
        assert report['layout_requested'] == 'circular'
        assert report['layout'] == 'circular'
        assert report['subtree_packing'] == 'tidy'

    def test_draw_rejects_scale_bar_for_topology_only_layout(
        self,
        tmp_nwk,
        tmp_path,
    ):
        infile = tmp_nwk('((A:1,B:2):1,C:3);')

        with pytest.raises(ValueError, match='branch-length-preserving'):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'fractal.svg'),
                species_overlap_node_plot='no',
                layout='fractal',
                tip_label_position='branch-end',
                scale_bar='auto',
            ))

    @pytest.mark.parametrize(
        ('layout_name', 'guide_type', 'encoding'),
        [
            ('slanted', 'axis-grid', 'depth-projection'),
            ('radial', 'concentric-rings', 'depth-projection'),
            ('spiral', 'spiral-depth-key', 'warped-depth'),
        ],
    )
    def test_draw_writes_layout_specific_depth_guide(
        self,
        tmp_nwk,
        tmp_path,
        layout_name,
        guide_type,
        encoding,
    ):
        infile = tmp_nwk('(((A:1,B:2):1,C:3):1,D:5);')
        outfile = tmp_path / '{}.svg'.format(layout_name)
        report_path = tmp_path / '{}.json'.format(layout_name)

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            image_format='svg',
            layout=layout_name,
            tip_label_position='branch-end',
            species_overlap_node_plot='no',
            depth_guide='1',
            branch_length_unit='substitutions/site',
            layout_report=str(report_path),
            figure_width=5.0,
            figure_height=4.0,
        ))

        svg = outfile.read_text(encoding='utf-8')
        report = json.loads(report_path.read_text(encoding='utf-8'))
        assert 'root-to-node distance' in svg.lower()
        assert report['depth_guide_interval'] == pytest.approx(1.0)
        assert report['depth_guide_type'] == guide_type
        assert report['branch_length_encoding'] == encoding
        assert report['branch_lengths_encoded'] is True
        if layout_name == 'radial':
            assert report['depth_guide_in_panel_labels'] is True
        assert not any(
            'depth_guide' in collision_kind
            for collision_kind in report['final_collisions_by_kind']
        )

    def test_draw_rejects_depth_guide_for_incompatible_layout(
        self,
        tmp_nwk,
        tmp_path,
    ):
        infile = tmp_nwk('(A:1,B:2);')

        with pytest.raises(ValueError, match='slanted, radial, and spiral'):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'tree.svg'),
                species_overlap_node_plot='no',
                depth_guide='auto',
            ))

    def test_draw_rejects_depth_guide_without_positive_branch_lengths(
        self,
        tmp_nwk,
        tmp_path,
    ):
        infile = tmp_nwk('(A,B);')

        with pytest.raises(ValueError, match='positive input branch lengths'):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'tree.svg'),
                layout='slanted',
                tip_label_position='branch-end',
                species_overlap_node_plot='no',
                depth_guide='auto',
            ))

    def test_draw_auto_collapses_only_the_rendering_copy(
        self,
        tmp_nwk,
        tmp_path,
    ):
        infile = tmp_nwk(
            '(((A:1,B:1):1,(C:1,D:1):1):1,'
            '((E:1,F:1):1,(G:1,H:1):1):1);'
        )
        outfile = tmp_path / 'collapsed.svg'
        report_path = tmp_path / 'collapsed.json'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            max_visible_tips=4,
            layout_report=str(report_path),
        ))

        svg = outfile.read_text(encoding='utf-8')
        report = json.loads(report_path.read_text(encoding='utf-8'))
        assert 'A…B (n=2)' in svg
        assert report['input_tip_count'] == 8
        assert report['visible_tip_count'] == 4
        assert len(report['collapsed_clades']) == 4
        with open(infile, encoding='utf-8') as handle:
            assert handle.read().count('A') == 1

    def test_draw_maps_tip_tracks_and_branch_styles_from_nhx(self, tmp_path):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        for index, node in enumerate(tree.traverse()):
            if not node.is_root:
                node.add_props(
                    regime='foreground' if index % 2 else 'background',
                    signal=index + 1,
                )
        for index, leaf in enumerate(tree.leaves()):
            leaf.add_props(
                state='X' if index % 2 else 'Y',
                score=index / 3.0,
            )
        infile = tmp_path / 'layers.nhx'
        infile.write_text(
            tree.write(
                props=['regime', 'signal', 'state', 'score'],
                parser=1,
                format_root_node=True,
            ),
            encoding='utf-8',
        )
        outfile = tmp_path / 'layers.svg'
        report_path = tmp_path / 'layers.json'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            tip_track=['state', 'score'],
            branch_color_property='regime',
            branch_width_property='signal',
            property_color=[
                'regime:foreground=#D55E00',
                'regime:background=#0072B2',
                'state:X=#009E73',
                'state:Y=#CC79A7',
            ],
            layout_report=str(report_path),
        ))

        svg = outfile.read_text(encoding='utf-8').lower()
        report = json.loads(report_path.read_text(encoding='utf-8'))
        assert '#d55e00' in svg
        assert '#0072b2' in svg
        assert '#009e73' in svg
        assert '#cc79a7' in svg
        assert report['artist_counts']['tip_track'] == 8
        assert report['tip_track_properties'] == ['state', 'score']

    def test_draw_taxonomy_typography_italicizes_exact_binomial(
        self,
        tmp_nwk,
        tmp_path,
    ):
        infile = tmp_nwk('(Arabidopsis_thaliana:1,Sample_1:1);')
        outfile = tmp_path / 'taxonomy.svg'

        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(outfile),
            species_overlap_node_plot='no',
            tip_label_font_style='taxonomy',
        ))

        svg = outfile.read_text(encoding='utf-8')
        assert "font: italic 8px 'Helvetica'" in svg

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

    def test_draw_uses_first_duplicate_tip_image_row(self, tmp_path):
        first_image = tmp_path / 'first.png'
        second_image = tmp_path / 'second.png'
        write_test_png(first_image)
        write_test_png(second_image)
        manifest = tmp_path / 'manifest.tsv'
        manifest.write_text(
            'leaf_name\tlocal_path\n'
            'A\tfirst.png\n'
            'A\tsecond.png\n'
        )

        path_by_leaf = _read_tip_image_manifest(
            str(manifest),
            tree=Tree('(A:1,B:1);'),
            unmatched='ignore',
        )

        assert path_by_leaf['A'] == str(first_image.resolve())

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
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" '
            '"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="40" height="20">'
            '<path d="M2 18 L20 2 L38 18 Z" fill="#111111"/>'
            '</svg>'
        )
        source_before = silhouette.read_bytes()
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
        assert silhouette.read_bytes() == source_before

    def test_draw_rejects_utf16_svg_before_entity_parsing(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('(A:1,B:1);')
        silhouette = tmp_path / 'utf16-entity.svg'
        silhouette.write_bytes((
            '<?xml version="1.0" encoding="UTF-16"?>'
            '<!DOCTYPE svg [<!ENTITY local "expanded">]>'
            '<svg xmlns="http://www.w3.org/2000/svg" width="40" height="20">'
            '<text>&local;</text>'
            '</svg>'
        ).encode('utf-16'))
        manifest = tmp_path / 'manifest.tsv'
        manifest.write_text(
            'leaf_name\tlocal_path\n'
            'A\tutf16-entity.svg\n'
            'B\tutf16-entity.svg\n'
        )

        with pytest.raises(RuntimeError, match='UTF-16|encoding'):
            draw_main(make_draw_args(
                infile=str(infile),
                outfile=str(tmp_path / 'utf16.svg'),
                species_overlap_node_plot='no',
                tip_image_manifest=str(manifest),
            ))

    def test_tip_image_loader_downsamples_before_retaining_array(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'large.png'
        Image.new('RGB', (400, 200), (20, 40, 60)).save(source)

        loaded = _load_tip_image(str(source), max_edge_px=40)

        assert loaded.shape == (20, 40, 4)
