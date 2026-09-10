"""Research-only extension of the known-error variance grid with a tail bound.

For C(v) = v A + D (A positive semidefinite, D positive definite),
log det C(v) is nondecreasing. Maximizing over any mean cannot give likelihood
above -0.5 * (d log(2 pi) + log det C(v)). A grid tail can be excluded only
when this bound is below every candidate's incumbent for every sample.
"""

import numpy as np
from scipy.linalg import solve_triangular

from nwkit.shift_calibration import CalibratedSearch


class TailCertifiedSearch(CalibratedSearch):
    """Same alpha grid; extend the geometric variance grid until its tail is ruled out.

    This is a research backend, not the production CLI. Certification concerns
    the upper tail, not spacing between variance grid points or continuous alpha.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self.known_error or np.any(self.variances <= 0):
            raise ValueError("Tail research requires strictly positive known errors")
        self.extension_count = 0
        self.tail_log_likelihood_upper_bound = np.inf

    def _extend_variance_grid(self):
        variance = float(self.variance_grid[-1] * 10**0.25)
        if not np.isfinite(variance):
            raise ValueError("Variance tail cannot be represented; rescale inputs")
        for ai, alpha in enumerate(self.grid):
            K, weights = self.geometry(alpha)
            X = np.einsum(
                "mij,mjk->mik", self.loads * weights[:, None, :], self.transform
            )
            projected = np.einsum("di,mij->mdj", self.q, X)
            covariance = variance * K + np.diag(self.variances)
            L = np.linalg.cholesky(self.q @ covariance @ self.q.T)
            W = (
                solve_triangular(
                    L, projected.transpose(1, 0, 2).reshape(self.d, -1), lower=True
                )
                .reshape(self.d, -1, 2)
                .transpose(1, 0, 2)
            )
            self.cache.append(
                (ai, variance, L, W, self._mean_solver(W), 2 * np.log(np.diag(L)).sum())
            )
        self.variance_grid = np.r_[self.variance_grid, variance]
        self.extension_count += 1

    def profile(self, Z):
        Z = np.asarray(Z, dtype=float)
        if Z.ndim == 1:
            Z = Z[:, None]
        best = np.full((len(self.models), Z.shape[1]), -np.inf)
        at = np.zeros(best.shape, dtype=int)
        start = 0
        for extension in range(65):
            for index in range(start, len(self.cache)):
                score, _, _ = self._at(Z, self.cache[index])
                improved = score > best + 1e-10
                at[improved] = index
                best[improved] = score[improved]
            bounds = [
                -0.5 * (self.d * np.log(2 * np.pi) + item[5])
                for item in self.cache
                if item[1] == self.variance_grid[-1]
            ]
            upper = float(max(bounds))
            self.tail_log_likelihood_upper_bound = upper
            if upper < float(best.min()) - 1e-8:
                return best, at
            if extension == 64:
                break
            start = len(self.cache)
            self._extend_variance_grid()
        raise ValueError("Upper variance tail not certified after 64 extensions")
