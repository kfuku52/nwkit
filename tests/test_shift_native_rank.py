"""Structural screening shares work only for identical covariance/missingness."""

import numpy as np
import pytest

from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_search import observable_layout
from tests.test_shift_native_search import sample_data


@pytest.mark.parametrize("different_mask", [False, True])
@pytest.mark.parametrize("alphas", [1.0, [0.0, float("inf")]])
def test_structural_screen_keeps_trait_specific_masks_and_alpha(
    monkeypatch, different_mask, alphas
):
    base = sample_data()
    values = base.values.copy()
    if different_mask:
        values[0, 1] = np.nan
    data = ShiftData.build(base.tree, values, base.trait_names)
    # The first terminal shift is unobserved in the second trait when masked.
    first_tip = data.tree.branch_ids[data.tree.compiled.leaf_indices[0]]
    layout = ShiftLayout.build(data.tree, [first_tip])
    expected = True
    for j, alpha in enumerate(np.broadcast_to(alphas, (2,))):
        design = layout.design(data.tree, alpha)[np.isfinite(data.values[:, j])]
        norms = np.linalg.norm(design, axis=0)
        expected &= bool(
            np.all(norms > 0)
            and np.linalg.matrix_rank(design / norms) == design.shape[1]
        )
    calls = []
    original = ShiftLayout.design

    def counted(self, tree, alpha_height, **kwargs):
        calls.append(alpha_height)
        return original(self, tree, alpha_height, **kwargs)

    monkeypatch.setattr(ShiftLayout, "design", counted)
    assert observable_layout(data, layout, alphas) == expected
    assert len(calls) == (1 if not different_mask and np.ndim(alphas) == 0 else 2)
