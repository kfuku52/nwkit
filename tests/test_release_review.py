"""Regression checks for the September 2026 integration review."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nwkit.asr import _symmetric_transition_matrices
from nwkit.cli import main
from nwkit.shift_native_model import ShiftLayout, ShiftTree
from nwkit.util import read_tree
from tests.test_shift_native_model import TREE, dense_reference


@pytest.mark.parametrize("components", [1, 2])
def test_long_symmetric_branches_preserve_each_stationary_component(components):
    generator = (np.ones((4, 4)) - 4 * np.eye(4)) / 3
    matrix = np.kron(np.eye(components), generator)
    expected = np.kron(np.eye(components), np.full((4, 4), 0.25))
    transitions = _symmetric_transition_matrices(matrix, [0.0, 1e20])
    np.testing.assert_array_equal(transitions[0.0], np.eye(len(matrix)))
    np.testing.assert_allclose(transitions[1e20], expected, atol=1e-13)
    np.testing.assert_allclose(transitions[1e20].sum(axis=1), 1.0, atol=1e-13)


def test_disconnected_slow_component_keeps_its_finite_time_dynamics():
    generator = (np.ones((4, 4)) - 4 * np.eye(4)) / 3
    matrix = np.kron(np.diag([1.0, 1e-20]), generator)
    transition = _symmetric_transition_matrices(matrix, [1e20])[1e20]
    expected = np.zeros((8, 8))
    expected[:4, :4] = 0.25
    expected[4:, 4:] = 0.25 + (np.eye(4) - 0.25) * np.exp(-4 / 3)
    np.testing.assert_allclose(transition, expected, atol=1e-13)
    np.testing.assert_allclose(transition.sum(axis=1), 1.0, atol=1e-13)


def test_shift_iterator_preserves_requested_regimes():
    tree = ShiftTree.build(read_tree(TREE, "auto", True, quiet=True))
    layout = ShiftLayout.build(tree, (branch for branch in [1, 7]))
    assert layout.shifts == (1, 7)
    expected, _ = dense_reference(tree, layout, 0.7, 1, np.zeros(8), "OUfixedRoot")
    np.testing.assert_allclose(layout.design(tree, 0.7), expected, atol=1e-14)
    with pytest.raises(ValueError, match="distinct non-root"):
        ShiftLayout.build(tree, iter([1, 1]))


def test_simulation_reads_utf8_parameters_with_non_utf8_locale(tmp_path, monkeypatch):
    tree = tmp_path / "tree.nwk"
    tree.write_text(TREE, encoding="utf-8")
    parameters = tmp_path / "parameters.json"
    names = ["体長", "体重"]
    parameters.write_text(
        json.dumps(
            {
                "trait_names": names,
                "alpha": 0.7,
                "process_tip_covariance": [[1.0, 0.2], [0.2, 1.0]],
                "regime_optima": [[0.0, 0.0]],
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    original = Path.read_text

    def locale_read(path, encoding=None, errors=None, **kwargs):
        return original(path, encoding=encoding or "cp1252", errors=errors, **kwargs)

    monkeypatch.setattr(Path, "read_text", locale_read)
    output, truth = tmp_path / "traits.tsv", tmp_path / "truth.json"
    main(
        [
            "shift-simulate",
            "--infile",
            str(tree),
            "--input-rooted",
            "yes",
            "--parameters",
            str(parameters),
            "--outfile",
            str(output),
            "--truth-out",
            str(truth),
        ]
    )
    assert set(names) <= set(pd.read_csv(output, sep="\t").columns)
    assert json.loads(truth.read_text(encoding="utf-8"))["trait_names"] == names
