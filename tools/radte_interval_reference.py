"""Independent dense GLS reference for branch-only log-duration contrasts.

No RADTE objective, gradient, precision builder or interval function is used.
The reference is exact for a linear Gaussian log-duration model with known K.
"""

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.stats import norm, t


def gaussian_contrast_interval(
    observations, design, covariance, contrast, *, level=0.95, sd=None
):
    y, x, k, contrast = map(np.asarray, (observations, design, covariance, contrast))
    if not 0 < level < 1 or (sd is not None and (not np.isfinite(sd) or sd <= 0)):
        raise ValueError("Interior level and positive supplied SD required.")
    whitening = cholesky(k, lower=True)
    wx = solve_triangular(whitening, x, lower=True)
    wy = solve_triangular(whitening, y, lower=True)
    beta, _, rank, _ = np.linalg.lstsq(wx, wy, rcond=None)
    df = len(y) - rank
    if rank != x.shape[1] or (sd is None and df <= 0):
        raise ValueError("Identified coefficients and positive residual df required.")
    residual = wy - wx @ beta
    variance = residual @ residual / df if sd is None else sd**2
    se = np.sqrt(variance * (contrast @ np.linalg.solve(wx.T @ wx, contrast)))
    center = float(contrast @ beta)
    critical = t.ppf((1 + level) / 2, df) if sd is None else norm.ppf((1 + level) / 2)
    return {
        "estimate": center,
        "lower": float(center - critical * se),
        "upper": float(center + critical * se),
        "df": int(df),
        "variance": float(variance),
    }
