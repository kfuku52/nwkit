import json
import math
from pathlib import Path
import random

import matplotlib.pyplot as plt
import pytest
from ete4 import Tree
from matplotlib.transforms import Bbox

import nwkit.draw as draw_module
from nwkit.draw import (
    _branch_style_maps,
    _property_color,
    _tip_track_colors,
    draw_main,
)
from nwkit.draw_layouts import (
    _graph_has_crossing,
    _undirected_graph_edges,
    _unrooted_graph,
    make_tree_layout,
)
from nwkit.draw_prep import collapse_tree_for_drawing
from nwkit.draw_quality import (
    DrawingArtist,
    _rectangle_union_area,
    evaluate_drawing,
)
from tests.test_draw import make_draw_args, write_test_png


def _random_binary_tree(tip_count, seed):
    generator = random.Random(seed)
    pool = []
    for index in range(tip_count):
        leaf = Tree()
        leaf.name = 'L{}'.format(index)
        pool.append(leaf)
    while len(pool) > 1:
        first = pool.pop(generator.randrange(len(pool)))
        second = pool.pop(generator.randrange(len(pool)))
        parent = Tree()
        parent.add_child(first, dist=generator.uniform(0.1, 2.0))
        parent.add_child(second, dist=generator.uniform(0.1, 2.0))
        pool.append(parent)
    return pool[0]


def _deep_caterpillar(depth):
    root = Tree()
    node = root
    for index in range(depth):
        node.add_child(name='L{}'.format(index), dist=1.0)
        node = node.add_child(dist=1.0)
    node.name = 'last'
    return root


def test_equal_daylight_randomized_layouts_preserve_lengths_without_crossings():
    for seed in range(12):
        tree = _random_binary_tree(40, 9200 + seed)
        layout = make_tree_layout(
            tree,
            layout='unrooted',
            unrooted_method='equal-daylight',
        )
        adjacency, _ = _unrooted_graph(tree, use_topology_depth=False)
        edges = _undirected_graph_edges(adjacency)

        assert not _graph_has_crossing(edges, layout.xcoord, layout.ycoord)
        for first, neighbors in adjacency.items():
            for second, expected_length in neighbors:
                observed_length = math.hypot(
                    layout.xcoord[first] - layout.xcoord[second],
                    layout.ycoord[first] - layout.ycoord[second],
                )
                assert observed_length == pytest.approx(expected_length, abs=1e-10)
        assert layout.metadata['daylight_iterations_actual'] <= 5


@pytest.mark.parametrize(
    'layout_name',
    ['rectangular', 'cladogram', 'fractal', 'circular'],
)
def test_deep_caterpillar_layouts_do_not_depend_on_python_recursion(layout_name):
    tree = _deep_caterpillar(1200)

    layout = make_tree_layout(tree, layout=layout_name)

    assert len(layout.xcoord) == 2401
    assert len(layout.leaf_order) == 1201


def test_deep_caterpillar_tidy_packing_does_not_depend_on_python_recursion():
    tree = _deep_caterpillar(1200)

    layout = make_tree_layout(
        tree,
        layout='rectangular',
        subtree_packing='tidy',
    )

    assert len(layout.xcoord) == 2401
    assert len(layout.leaf_order) == 1201


def test_deep_caterpillar_collapse_uses_iterative_copy():
    tree = _deep_caterpillar(1200)

    drawing_tree, collapsed = collapse_tree_for_drawing(tree, 100)

    assert len(list(drawing_tree.leaves())) <= 100
    assert collapsed


def test_draw_rejects_output_report_alias_before_mutation(tmp_nwk, tmp_path):
    infile = tmp_nwk('(A:1,B:1);')
    destination = tmp_path / 'same.svg'
    destination.write_text('sentinel', encoding='utf-8')

    with pytest.raises(ValueError, match='different paths'):
        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(destination),
            layout_report=str(destination),
            species_overlap_node_plot='no',
        ))

    assert destination.read_text(encoding='utf-8') == 'sentinel'


def test_draw_rejects_report_aliasing_input_before_mutation(tmp_nwk, tmp_path):
    infile = tmp_nwk('(A:1,B:1);')
    original = Path(infile).read_text(encoding='utf-8')

    with pytest.raises(ValueError, match='must not overwrite'):
        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(tmp_path / 'tree.svg'),
            layout_report=str(infile),
            species_overlap_node_plot='no',
        ))

    assert Path(infile).read_text(encoding='utf-8') == original


def test_draw_rejects_output_aliasing_tip_image_asset(tmp_nwk, tmp_path):
    infile = tmp_nwk('(A:1,B:1);')
    image_path = tmp_path / 'tip.png'
    write_test_png(image_path)
    original = image_path.read_bytes()
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text(
        'leaf_name\tlocal_path\nA\ttip.png\n',
        encoding='utf-8',
    )

    with pytest.raises(ValueError, match='must not overwrite'):
        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(image_path),
            tip_image_manifest=str(manifest),
            species_overlap_node_plot='no',
        ))

    assert image_path.read_bytes() == original


def test_noop_max_visible_tips_can_be_combined_with_trait(tmp_nwk, tmp_path):
    infile = tmp_nwk('(A:1,B:1);')
    trait = tmp_path / 'trait.tsv'
    trait.write_text('leaf_name\tgroup\nA\tx\nB\ty\n', encoding='utf-8')
    outfile = tmp_path / 'tree.svg'

    draw_main(make_draw_args(
        infile=str(infile),
        outfile=str(outfile),
        max_visible_tips=100,
        trait=str(trait),
        group_by='group',
        species_overlap_node_plot='no',
    ))

    assert outfile.is_file()


def test_collapsed_property_summary_distinguishes_missing_mixed_and_mean():
    tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
    values = {
        'A': {'state': 'X', 'score': 10},
        'B': {'score': 20},
        'C': {'state': 'Y', 'score': 30},
        'D': {},
    }
    for leaf in tree.leaves():
        leaf.add_props(**values[leaf.name])

    default_tree, default_report = collapse_tree_for_drawing(tree, 2)
    mean_tree, mean_report = collapse_tree_for_drawing(
        tree,
        2,
        property_aggregation='mean',
    )

    assert [leaf.props['state'] for leaf in default_tree.leaves()] == [
        'partial', 'partial'
    ]
    default_leaves = list(default_tree.leaves())
    mean_leaves = list(mean_tree.leaves())
    assert default_leaves[0].props['score'] == 'mixed'
    assert mean_leaves[0].props['score'] == pytest.approx(15.0)
    assert mean_leaves[1].props['score'] == 'partial'
    assert default_report[0]['property_summary']['score']['status'] == 'mixed'
    assert mean_report[0]['property_summary']['score']['status'] == 'mean'


def test_implicit_branch_colors_are_stable_distinct_and_match_legend_mapping():
    tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
    values = ['z', 'x', 'y']
    nodes = [node for node in tree.traverse() if not node.is_root]
    for index, node in enumerate(nodes):
        node.add_props(regime=values[index % len(values)])

    colors, _ = _branch_style_maps(
        tree,
        base_color='#000000',
        base_width=0.8,
        color_property='regime',
        width_property=None,
        width_range='0.4,2.5',
        property_colors={},
        palette='tab10',
    )
    expected = {
        value: _property_color(
            'regime',
            value,
            {},
            fallback_index=index,
            palette='tab10',
        )
        for index, value in enumerate(sorted(values))
    }

    assert len(set(expected.values())) == 3
    assert all(colors[node] == expected[node.props['regime']] for node in nodes)


def test_tip_track_legends_cover_constant_missing_and_many_categories():
    constant_tree = Tree('(A:1,B:1,C:1);', parser=1)
    constant_tree['A'].add_props(score=2)
    constant_tree['B'].add_props(score=2)
    _, constant_legend = _tip_track_colors(
        constant_tree,
        ['score'],
        'auto',
        {},
        'tab10',
        'viridis',
    )
    assert [label for label, _ in constant_legend] == [
        'score: 2', 'score: missing'
    ]

    category_tree = Tree(
        '(' + ','.join('L{}:1'.format(index) for index in range(13)) + ');',
        parser=1,
    )
    for index, leaf in enumerate(category_tree.leaves()):
        leaf.add_props(group='g{}'.format(index))
    _, category_legend = _tip_track_colors(
        category_tree,
        ['group'],
        'categorical',
        {},
        'tab20',
        'viridis',
    )
    assert len(category_legend) == 13


@pytest.mark.parametrize('layout_name', ['slanted', 'radial'])
def test_scale_bar_rejects_depth_projection_layouts(
    tmp_nwk,
    tmp_path,
    layout_name,
):
    infile = tmp_nwk('((A:1,B:2):1,C:3);')

    with pytest.raises(ValueError, match='directly measurable'):
        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(tmp_path / '{}.svg'.format(layout_name)),
            layout=layout_name,
            tip_label_position='branch-end',
            scale_bar='1',
            species_overlap_node_plot='no',
        ))


def test_scale_bar_rejects_value_larger_than_tree_depth(tmp_nwk, tmp_path):
    infile = tmp_nwk('(A:1,B:2);')

    with pytest.raises(ValueError, match='must not exceed'):
        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(tmp_path / 'tree.svg'),
            scale_bar='100',
            species_overlap_node_plot='no',
        ))


def test_layout_report_discloses_overflow_and_rendering_metadata(
    tmp_nwk,
    tmp_path,
    capsys,
):
    infile = tmp_nwk('({}:1,B:1);'.format('LONG_' + ('X' * 500)))
    report_path = tmp_path / 'layout.json'

    draw_main(make_draw_args(
        infile=str(infile),
        outfile=str(tmp_path / 'tree.svg'),
        layout_report=str(report_path),
        figure_width=1.5,
        figure_height=1.5,
        tip_label_position='branch-end',
        species_overlap_node_plot='no',
    ))

    report = json.loads(report_path.read_text(encoding='utf-8'))
    assert report['fits_within_figure'] is False
    assert report['maximum_overflow_points'] > 0.0
    assert report['branch_collision_check_complete'] is True
    assert report['nwkit_version']
    assert report['font_size_points'] == 8.0
    assert report['output_format'] == 'svg'
    assert 'exceed the fixed figure boundary' in capsys.readouterr().err


def test_new_matplotlib_figures_are_closed_when_rendering_raises(
    tmp_nwk,
    tmp_path,
    monkeypatch,
):
    infile = tmp_nwk('(A:1,B:1);')
    existing = set(plt.get_fignums())

    def fail_quality_check(*args, **kwargs):
        raise RuntimeError('quality failure')

    monkeypatch.setattr(draw_module, 'evaluate_drawing', fail_quality_check)
    with pytest.raises(RuntimeError, match='quality failure'):
        draw_main(make_draw_args(
            infile=str(infile),
            outfile=str(tmp_path / 'tree.svg'),
            species_overlap_node_plot='no',
        ))

    assert set(plt.get_fignums()) == existing


def test_rectangle_union_area_uses_exact_sweep_line_union():
    rectangles = [
        Bbox.from_extents(0, 0, 2, 2),
        Bbox.from_extents(1, 1, 3, 3),
        Bbox.from_extents(10, 10, 11, 11),
    ]

    assert _rectangle_union_area(rectangles) == pytest.approx(8.0)


def test_collision_audit_remains_complete_above_500_artists():
    figure, axes = plt.subplots(figsize=(4.0, 12.0))
    axes.set_xlim(0, 1)
    axes.set_ylim(0, 600)
    artists = [
        DrawingArtist(
            axes.text(0.8, index, 'x', fontsize=4),
            kind='tip_label',
            owner=index,
        )
        for index in range(600)
    ]
    line = axes.plot([0.1, 0.1], [0, 600])[0]
    try:
        report = evaluate_drawing(
            figure,
            artists,
            [('tree', line)],
            policy='ignore',
        )
    finally:
        plt.close(figure)

    assert report['branch_collision_check_complete'] is True


@pytest.mark.parametrize('layout_name', ['radial', 'circular'])
def test_spatial_root_stub_occupies_an_empty_angular_gap(layout_name):
    tree = Tree('(A:1,B:1,C:1,D:1);', parser=1)
    layout = make_tree_layout(tree, layout=layout_name)
    root_start, root_end = layout.root_path
    stub_angle = math.atan2(
        root_start[1] - root_end[1],
        root_start[0] - root_end[0],
    )
    child_angles = [
        math.atan2(
            layout.ycoord[child] - layout.ycoord[tree],
            layout.xcoord[child] - layout.xcoord[tree],
        )
        for child in tree.get_children()
    ]

    assert all(
        abs(math.sin(stub_angle - child_angle)) > 1e-6
        for child_angle in child_angles
    )
