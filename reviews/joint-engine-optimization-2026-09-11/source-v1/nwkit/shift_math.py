"""Scale-aware checks and branch geometry for OU shift inference."""

import math
from typing import Any


def remaining_heights(tree):
    """Measure time forward to tips without subtracting large root depths."""
    heights: dict[Any, float] = {}
    for node in tree.traverse("postorder"):
        heights[node] = (
            0.0
            if node.is_leaf
            else max(child.dist + heights[child] for child in node.children)
        )
    return heights


def close_in_units(actual, expected, *, operands=(), rtol=1e-7):
    """Allow relative error and rounding of source operands, not a unitful floor."""
    values = (actual, expected, *operands)
    if not all(math.isfinite(x) for x in values):
        return False
    tolerance = 32 * max(math.ulp(float(x)) for x in values)
    return math.isclose(actual, expected, rel_tol=rtol, abs_tol=tolerance)


def ancestral_effect_scales(tree, ids, intercept, effects, field):
    """Source and intermediate magnitudes along each ancestry only."""
    scales: dict[Any, float] = {}
    partial: dict[Any, float] = {}
    for node in tree.traverse("preorder"):
        delta = effects.get(ids[node], {}).get(field)
        partial[node] = intercept if node.is_root else partial[node.up] + (delta or 0.0)
        scales[node] = max(
            abs(intercept) if node.is_root else scales[node.up],
            abs(delta) if delta is not None else 0.0,
            abs(partial[node]),
        )
    return scales


def observation_variances(standard_errors):
    """Square finite nonnegative SEs, rejecting loss of representable variance."""
    values = []
    for error in standard_errors:
        error = float(error)
        if not math.isfinite(error) or error < 0:
            raise ValueError("Standard errors must be finite and non-negative.")
        variance = error * error
        if not math.isfinite(variance) or (error > 0 and variance == 0):
            raise ValueError(
                "A squared standard error is outside floating-point range."
            )
        values.append(variance)
    return values
