"""Complete-search bootstrap calibration of small-tree OU contrast likelihoods."""

import numpy as np
from scipy.linalg import helmert, solve_triangular

from nwkit.shift_candidates import enumerate_candidates
from nwkit.util import assign_branch_ids

ALPHA_HEIGHT_GRID = np.r_[0.0, np.inf, np.geomspace(0.001, 1000.0, 25)]


def validate_calibration_options(replicates, level):
    if (
        isinstance(replicates, (bool, np.bool_))
        or not isinstance(replicates, (int, np.integer))
        or replicates < 19
        or not 0 < level < 1
        or 1 / (replicates + 1) > level
    ):
        raise ValueError(
            "Calibration needs integer >=19 replicates and Monte Carlo resolution no larger than the test level."
        )


class CalibratedSearch:
    def __init__(
        self,
        tree,
        max_shifts=2,
        convergence=True,
        variances=None,
        limit=5000,
        alpha_grid=None,
    ):
        self.tree = tree
        leaves = list(tree.leaves())
        self.n = len(leaves)
        self.d = self.n - 1
        if self.n < 4 or max_shifts >= self.d:
            raise ValueError(
                "Calibrated selection requires 4..16 tips and positive residual degrees of freedom."
            )
        self.q = helmert(self.n, full=False)
        ids = assign_branch_ids(tree)
        bid = {i: node for node, i in ids.items()}
        self.height = max(tree.get_distance(tree, t) for t in leaves)
        self.distance = np.array(
            [[tree.get_distance(a, b) / self.height for b in leaves] for a in leaves]
        )
        self.models, self.counts = enumerate_candidates(
            tree, max_shifts, limit, convergence=convergence
        )
        self.dim = np.array([len(m["groups"]) - 1 for m in self.models])
        self.k = np.array([len(m["shift_branch_ids"]) for m in self.models])
        self._build_mean_design(tree, ids, bid, leaves)
        self.variances = (
            np.zeros(self.n)
            if variances is None
            else np.asarray(variances, dtype=float)
        )
        if (
            self.variances.shape != (self.n,)
            or not np.all(np.isfinite(self.variances))
            or np.any(self.variances < 0)
        ):
            raise ValueError(
                "Known observation variances must be finite, nonnegative and tip-aligned."
            )
        self.known_error = bool(np.any(self.variances > 0))
        self.grid = (
            ALPHA_HEIGHT_GRID.copy()
            if alpha_grid is None
            else np.asarray(alpha_grid, dtype=float)
        )
        if (
            self.grid.ndim != 1
            or not len(self.grid)
            or np.any(np.isnan(self.grid))
            or np.any(self.grid < 0)
        ):
            raise ValueError("Invalid alpha-height grid.")
        self.variance_grid = (
            np.r_[
                0,
                np.median(self.variances[self.variances > 0])
                * np.geomspace(1e-6, 1e6, 49),
            ]
            if self.known_error
            else np.array([np.nan])
        )
        if self.known_error and (
            not np.all(np.isfinite(self.variance_grid))
            or np.any(self.variance_grid[1:] <= 0)
        ):
            raise ValueError(
                "Variance grid is not representable; rescale the trait and its standard errors together."
            )
        self.cache = []
        for ai, alpha in enumerate(self.grid):
            K, weights = self.geometry(alpha)
            X = np.einsum(
                "mij,mjk->mik", self.loads * weights[:, None, :], self.transform
            )
            projected = np.einsum("di,mij->mdj", self.q, X)
            for variance in self.variance_grid:
                if (
                    self.known_error
                    and variance == 0
                    and np.count_nonzero(self.variances == 0) > 1
                ):
                    # Rank(Q D Q') < n-1 with two or more exact observations.
                    # Do not let roundoff turn a singular boundary into a fit.
                    continue
                covariance = (
                    variance * K + np.diag(self.variances) if self.known_error else K
                )
                try:
                    L = np.linalg.cholesky(self.q @ covariance @ self.q.T)
                except np.linalg.LinAlgError:
                    if self.known_error and variance == 0:
                        continue
                    raise
                W = (
                    solve_triangular(
                        L, projected.transpose(1, 0, 2).reshape(self.d, -1), lower=True
                    )
                    .reshape(self.d, -1, 2)
                    .transpose(1, 0, 2)
                )
                solver = self._mean_solver(W)
                self.cache.append(
                    (ai, variance, L, W, solver, 2 * np.log(np.diag(L)).sum())
                )
        families = [np.flatnonzero(self.k == 0), np.flatnonzero(self.k <= 1)]
        if convergence:
            families.append(np.flatnonzero(self.dim <= 1))
        families.append(np.arange(len(self.models)))
        self.families: list[np.ndarray] = []
        for family in families:
            if not any(np.array_equal(family, previous) for previous in self.families):
                self.families.append(family)

    def _build_mean_design(self, tree, ids, bid, leaves):
        self.loads = np.zeros((len(self.models), self.n, 2))
        self.age = np.zeros((len(self.models), 2))
        self.transform = np.zeros((len(self.models), 2, 2))
        for m, model in enumerate(self.models):
            aliases = {b: g for g, group in enumerate(model["groups"]) for b in group}
            base = aliases[0]
            groupids = {
                g: i
                for i, g in enumerate(
                    g for g in range(len(model["groups"])) if g != base
                )
            }
            for j, branch in enumerate(model["shift_branch_ids"]):
                node = bid[branch]
                descendants = set(node.leaves())
                self.loads[m, :, j] = [tip in descendants for tip in leaves]
                self.age[m, j] = (
                    self.height - tree.get_distance(tree, node.up)
                ) / self.height
                ancestor = node.up
                while ids[ancestor] not in aliases:
                    ancestor = ancestor.up
                for sign, group in [(1, aliases[branch]), (-1, aliases[ids[ancestor]])]:
                    if group != base:
                        self.transform[m, j, groupids[group]] += sign

    def geometry(self, alpha):
        if alpha == 0:
            return 1 - self.distance / 2, self.age
        if np.isinf(alpha):
            return np.eye(self.n), (self.age > 0).astype(float)
        K = (
            np.exp(-alpha * self.distance)
            * (-np.expm1(-alpha * (2 - self.distance)))
            / (-np.expm1(-2 * alpha))
        )
        return K, -np.expm1(-alpha * self.age) / -np.expm1(-alpha)

    def _mean_solver(self, W):
        # QR avoids squaring the condition number through normal equations.
        solver = np.zeros((len(self.models), 2, self.d))
        for dimension in (1, 2):
            selected = self.dim == dimension
            if not np.any(selected):
                continue
            orthogonal, triangular = np.linalg.qr(W[selected, :, :dimension])
            solver[selected, :dimension, :] = np.linalg.solve(
                triangular, orthogonal.transpose(0, 2, 1)
            )
        return solver

    def _at(self, Z, item):
        _, _, L, W, solver, logdet = item
        w = solve_triangular(L, Z, lower=True)
        beta = np.einsum("mid,db->mib", solver, w)
        # Explicit residuals avoid subtraction of nearly equal quadratic forms.
        residual = w[None, :, :] - np.einsum("mdi,mib->mdb", W, beta)
        rss = np.einsum("mdb,mdb->mb", residual, residual)
        threshold = 100 * np.finfo(float).eps ** 2 * np.sum(w * w, axis=0)
        if not self.known_error and np.any(rss <= threshold[None, :]):
            raise ValueError(
                "Zero residual variance: calibrated Gaussian selection is undefined."
            )
        if self.known_error:
            score = -0.5 * (self.d * np.log(2 * np.pi) + logdet + rss)
        else:
            score = -0.5 * (
                self.d * (np.log(2 * np.pi) + 1 + np.log(rss / self.d)) + logdet
            )
        if not np.all(np.isfinite(score)):
            raise ValueError(
                "Nonfinite contrast likelihood; rescale the trait and its standard errors together."
            )
        return score, beta, rss

    def profile(self, Z):
        Z = np.asarray(Z, dtype=float)
        if Z.ndim == 1:
            Z = Z[:, None]
        best = np.full((len(self.models), Z.shape[1]), -np.inf)
        at = np.zeros(best.shape, dtype=int)
        for index, item in enumerate(self.cache):
            score, _, _ = self._at(Z, item)
            improved = score > best + 1e-10
            at[improved] = index
            best[improved] = score[improved]
        if self.known_error:
            selected_variances = np.array([item[1] for item in self.cache])[at]
            if np.any(selected_variances == self.variance_grid[-1]):
                raise ValueError(
                    "Process variance reached the upper grid boundary; review inputs before inference."
                )
        return best, at

    def alpha_profile(self, z, winner):
        profile = np.full(len(self.grid), -np.inf)
        for item in self.cache:
            score, _, _ = self._at(np.asarray(z)[:, None], item)
            profile[item[0]] = max(profile[item[0]], score[winner, 0])
        return profile

    def fit(self, values, seed=1, replicates=199, level=0.05):
        validate_calibration_options(replicates, level)
        values = np.asarray(values, dtype=float)
        if values.shape != (self.n,) or not np.all(np.isfinite(values)):
            raise ValueError("Finite tip-aligned observations are required.")
        z = self.q @ (values - values.mean())
        best, at = self.profile(z)
        best, at = best[:, 0], at[:, 0]
        rng = np.random.default_rng(seed)
        tests = []
        for stage, family in enumerate(self.families):
            winner = int(family[np.argmax(best[family])])
            item = self.cache[at[winner]]
            _, variance, L, W, _, _ = item
            _, coefficients, rss = self._at(z[:, None], item)
            beta = coefficients[winner, :, 0]
            scale = (
                1.0
                if self.known_error
                else rss[winner, 0] / (self.d - self.dim[winner])
            )
            if stage == len(self.families) - 1:
                break
            mean = L @ W[winner] @ beta
            observed = float(2 * (best.max() - best[winner]))
            p = self._calibrated_probability(
                mean, L, scale, family, observed, rng, replicates
            )
            tests.append(
                {
                    "stage": stage,
                    "candidate_count": len(family),
                    "statistic": observed,
                    "p_value": p,
                }
            )
            if p > level:
                break
        alpha_height = self.grid[item[0]]
        process_variance = float(
            variance if self.known_error else rss[winner, 0] / self.d
        )
        profile = self.alpha_profile(z, winner)
        support = profile >= profile.max() - 1.920729410347062
        limit_supported = bool(
            np.any(support & ((self.grid == 0) | np.isinf(self.grid)))
        )
        status = (
            "brownian_limit"
            if alpha_height == 0
            else ("independent_limit" if np.isinf(alpha_height) else "finite")
        )
        K, weights = self.geometry(alpha_height)
        effects = self.transform[winner] @ beta
        mean_effects = weights[winner] * effects
        mean = self.loads[winner] @ mean_effects
        covariance = process_variance * K + np.diag(self.variances)
        inverse_one = (
            (self.variances == 0).astype(float)
            if process_variance == 0 and np.any(self.variances == 0)
            else np.linalg.solve(covariance, np.ones(self.n))
        )
        intercept = float(inverse_one @ (values - mean) / inverse_one.sum())
        identifiable = bool(
            status == "finite" and not limit_supported and process_variance > 0
        )
        optimum_effects = (
            effects / -np.expm1(-alpha_height) if identifiable else np.full(2, np.nan)
        )
        return {
            "winner": winner,
            "model": self.models[winner],
            "tests": tests,
            "alpha": None
            if np.isinf(alpha_height)
            else float(alpha_height / self.height),
            "alpha_height": None if np.isinf(alpha_height) else float(alpha_height),
            "alpha_status": status,
            "alpha_limit_supported": limit_supported,
            "alpha_profile": [
                {
                    "alpha_height": None if np.isinf(a) else float(a),
                    "limit": "infinity" if np.isinf(a) else None,
                    "log_likelihood": float(ll),
                    "supported": bool(s),
                }
                for a, ll, s in zip(self.grid, profile, support, strict=True)
            ],
            "process_tip_variance": process_variance,
            "sigma2": None
            if np.isinf(alpha_height)
            else float(
                process_variance
                / self.height
                * (
                    1
                    if alpha_height == 0
                    else 2 * alpha_height / -np.expm1(-2 * alpha_height)
                )
            ),
            "intercept": intercept,
            "predicted": (intercept + mean).tolist(),
            "mean_effects": mean_effects.tolist(),
            "optimum_effects": [
                float(e) if identifiable else None for e in optimum_effects
            ],
            "optimum_identifiable": identifiable,
            "contrast_log_likelihood": float(best[winner]),
            "candidate_count": len(self.models),
            "calibration_replicates": replicates,
            "calibration_level": level,
            "seed": seed,
        }

    def _calibrated_probability(
        self, mean, L, scale, family, observed, rng, replicates
    ):
        # Keep the original row-major RNG sequence, but never build a
        # candidate-by-contrast-by-all-replicates residual tensor.
        noise = rng.normal(size=(self.d, replicates))
        exceedances = 0
        for start in range(0, replicates, 64):
            sim = mean[:, None] + np.sqrt(scale) * L @ noise[:, start : start + 64]
            sb, _ = self.profile(sim)
            statistics = 2 * (sb.max(axis=0) - sb[family].max(axis=0))
            exceedances += np.count_nonzero(statistics >= observed - 1e-10)
        return float((1 + exceedances) / (replicates + 1))
