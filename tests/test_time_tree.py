import numpy as np
import pytest
from ete4 import Tree

from nwkit.time_tree import (
    _credible_interval,
    _hpd_interval,
    annotate_mcmctree_calibrations,
    parse_mcmctree_calibration,
    paml_node_mapping,
    posterior_sample_tree,
    read_mcmctree_posterior,
    summarize_mcmctree_posterior,
)


def test_mcmctree_interval_indices_match_paml_order_statistics():
    values = np.arange(10, dtype=float)

    assert _hpd_interval(values, 0.8) == (0.0, 8.0)
    assert _credible_interval(values, 0.8, 'equal-tail') == (1.0, 9.0)


@pytest.mark.parametrize(
    ('text', 'kind', 'lower', 'upper'),
    [
        ('@12.5', 'point', 12.5, 12.5),
        ('B(10, 20, 0.025, 0.025)', 'bounded', 10.0, 20.0),
        ('L(10, 0.1, 1, 0.025)', 'lower', 10.0, None),
        ('U(20, 0.025)', 'upper', None, 20.0),
        ('>.06<.08', 'bounded', 0.06, 0.08),
        ('>.06', 'lower', 0.06, None),
        ('<.08', 'upper', None, 0.08),
    ],
)
def test_parse_mcmctree_calibration_variants(text, kind, lower, upper):
    record = parse_mcmctree_calibration(text)

    assert record['type'] == kind
    assert record.get('lower') == lower
    assert record.get('upper') == upper


@pytest.mark.parametrize(
    ('text', 'message'),
    [
        ('B(20,10,0.025,0.025)', 'must not exceed'),
        ('B(10,20,0,0.025)', 'tail probabilities'),
        ('L(10,0,-1,0.025)', 'scales must be positive'),
        ('@1e999', 'must be finite'),
    ],
)
def test_parse_mcmctree_calibration_rejects_invalid_values(text, message):
    with pytest.raises(ValueError, match=message):
        parse_mcmctree_calibration(text)


def test_paml_node_mapping_is_tips_then_internal_preorder():
    tree = Tree('(((A,B),C),(D,E));', parser=1)
    mapping = paml_node_mapping(tree)

    assert [mapping[index].name for index in range(1, 6)] == ['A', 'B', 'C', 'D', 'E']
    assert mapping[6].is_root
    assert set(mapping[7].leaf_names()) == {'A', 'B', 'C'}
    assert set(mapping[8].leaf_names()) == {'A', 'B'}
    assert set(mapping[9].leaf_names()) == {'D', 'E'}


def test_read_and_summarize_posterior_reconstructs_dated_tree(tmp_path):
    tree = Tree("((A,B)'B(2,6,0.025,0.025)',(C,D));", parser=1)
    posterior_path = tmp_path / 'mcmc.txt'
    posterior_path.write_text(
        'Gen\tt_n5\tt_n6\tt_n7\tlnL\n'
        '1\t10\t4\t3\t-1\n'
        '2\t12\t5\t4\t-1\n'
        '3\t11\t6\t5\t-1\n'
        '4\t13\t7\t6\t-1\n'
    )

    annotate_mcmctree_calibrations(tree)
    posterior = read_mcmctree_posterior(posterior_path, tree, burnin=1, thin=2)
    summarize_mcmctree_posterior(tree, posterior, point='median', ci='equal-tail', level=0.8)

    assert posterior.sample_count == 2
    assert posterior.node_ids == (5, 6, 7)
    assert tree.props['age'] == pytest.approx(12.5)
    assert tree.common_ancestor(['A', 'B']).props['age'] == pytest.approx(6.0)
    assert tree.common_ancestor(['A', 'B']).props['calibration_type'] == 'bounded'
    assert tree['A'].dist == pytest.approx(6.0)
    assert tree.common_ancestor(['A', 'B']).dist == pytest.approx(6.5)
    assert tree.props['age_ci_kind'] == 'equal-tail'
    assert np.isfinite(float(tree.props['age_ci_low']))

    sample = posterior_sample_tree(tree, posterior, 0)
    assert sample['A'].dist == pytest.approx(5.0)
    assert sample.common_ancestor(['A', 'B']).dist == pytest.approx(7.0)


def test_posterior_rejects_chronologically_invalid_sample(tmp_path):
    tree = Tree('((A,B),(C,D));', parser=1)
    path = tmp_path / 'mcmc.txt'
    path.write_text('Gen t_n5 t_n6 t_n7\n1 4 5 2\n')

    with pytest.raises(ValueError, match='older than its parent'):
        read_mcmctree_posterior(path, tree)


def test_posterior_rejects_nonbinary_topology(tmp_path):
    tree = Tree('(A,B,C);', parser=1)
    path = tmp_path / 'mcmc.txt'
    path.write_text('Gen t_n4\n1 4\n')

    with pytest.raises(ValueError, match='rooted and binary'):
        read_mcmctree_posterior(path, tree)
