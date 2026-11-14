"""Independent dense OU contrast profiles for small scientific audits.

Branch innovations and mean recursion define this oracle; it does not import
the production covariance, mean-design, likelihood or optimization routines.
Finite profiles refine every local maximum of a logarithmic scouting grid.
This is a numerical reference, not a certificate of a continuous global optimum.
"""

import numpy as np
from scipy.linalg import null_space, solve_triangular
from scipy.optimize import minimize_scalar

from nwkit.util import assign_branch_ids


class DenseOUReference:
    def __init__(self, tree, model, variances=None):
        self.tree = tree
        self.nodes = list(tree.traverse("preorder"))
        self.tips = list(tree.leaves())
        self.n = len(self.tips)
        self.d = self.n - 1
        self.height = tree.get_distance(tree, self.tips[0])
        self.q = null_space(np.ones((1, self.n))).T
        self.variances = (
            np.zeros(self.n)
            if variances is None
            else np.asarray(variances, dtype=float)
        )
        self.known_error = bool(np.any(self.variances))
        self.ids = assign_branch_ids(tree)
        aliases = {b: g for g, group in enumerate(model["groups"]) for b in group}
        baseline = aliases[0]
        nonbaseline = [g for g in range(len(model["groups"])) if g != baseline]
        self.dimension = len(nonbaseline)
        self.groups = {}
        self.targets = {}
        for node in self.nodes:
            group = aliases.get(
                self.ids[node], baseline if node.is_root else self.groups[node.up]
            )
            self.groups[node] = group
            self.targets[node] = np.array([float(group == g) for g in nonbaseline])

    def branch_moments(self, alpha_height):
        mean = {self.tree: np.zeros(self.dimension)}
        innovations = {self.tree: np.zeros(len(self.nodes))}
        a = float(alpha_height)
        for index, node in enumerate(self.nodes[1:], 1):
            elapsed = node.dist / self.height
            if a == 0:
                decay, response, variance = 1.0, elapsed, elapsed
            elif np.isinf(a):
                decay, response, variance = 0.0, 1.0, 1.0
            else:
                decay = np.exp(-a * elapsed)
                response = -np.expm1(-a * elapsed) / -np.expm1(-a)
                variance = -np.expm1(-2 * a * elapsed) / -np.expm1(-2 * a)
            mean[node] = decay * mean[node.up] + response * self.targets[node]
            innovations[node] = decay * innovations[node.up]
            innovations[node][index] += np.sqrt(variance)
        design = np.array([mean[tip] for tip in self.tips])
        loads = np.array([innovations[tip] for tip in self.tips])
        return design, loads @ loads.T

    def at(self, values, alpha_height, process_variance=None):
        """Profile mean coefficients, and analytic process scale without errors."""
        design, K = self.branch_moments(alpha_height)
        if self.known_error:
            if process_variance is None or process_variance < 0:
                raise ValueError(
                    "Known-error reference needs a nonnegative process variance"
                )
            if process_variance == 0 and np.count_nonzero(self.variances == 0) > 1:
                raise ValueError("Singular zero-process covariance")
            covariance = process_variance * K + np.diag(self.variances)
        else:
            covariance = K
        L = np.linalg.cholesky(self.q @ covariance @ self.q.T)
        response = solve_triangular(L, self.q @ np.asarray(values), lower=True)
        X = solve_triangular(L, self.q @ design, lower=True)
        coefficients, _, rank, _ = np.linalg.lstsq(X, response, rcond=None)
        if rank != self.dimension:
            raise ValueError("Rank-deficient candidate mean design")
        residual = response - X @ coefficients
        rss = float(residual @ residual)
        logdet = 2 * np.log(np.diag(L)).sum()
        if self.known_error:
            ll = -0.5 * (self.d * np.log(2 * np.pi) + logdet + rss)
            variance = process_variance
        else:
            if rss <= 0:
                raise ValueError("Zero residual process variance")
            variance = rss / self.d
            ll = -0.5 * (self.d * (np.log(2 * np.pi * variance) + 1) + logdet)
        return {
            "log_likelihood": float(ll),
            "process_tip_variance": float(variance),
            "coefficients": coefficients,
            "mean_contrasts": self.q @ design @ coefficients,
        }

    def at_profiled_variance(self, values, alpha_height):
        if not self.known_error:
            return self.at(values, alpha_height)
        scale = float(np.median(self.variances[self.variances > 0]))
        candidates = []
        if np.count_nonzero(self.variances == 0) <= 1:
            candidates.append(self.at(values, alpha_height, 0))
        grid = np.linspace(
            np.log(scale) - 12 * np.log(10), np.log(scale) + 12 * np.log(10), 97
        )
        fits = [self.at(values, alpha_height, np.exp(x)) for x in grid]
        candidates.extend(fits)
        for index in range(1, len(grid) - 1):
            if fits[index]["log_likelihood"] >= max(
                fits[index - 1]["log_likelihood"], fits[index + 1]["log_likelihood"]
            ):
                fit = minimize_scalar(
                    lambda x: (
                        -self.at(values, alpha_height, np.exp(x))["log_likelihood"]
                    ),
                    bounds=(grid[index - 1], grid[index + 1]),
                    method="bounded",
                    options={"xatol": 1e-10},
                )
                candidates.append(self.at(values, alpha_height, np.exp(fit.x)))
        if fits[-1]["log_likelihood"] >= fits[-2]["log_likelihood"]:
            raise ValueError("Continuous reference variance exceeds its scouting range")
        return max(candidates, key=lambda fit: fit["log_likelihood"])

    def profile(self, values, grid_points=161):
        alphas = [0.0, np.inf, *np.geomspace(1e-8, 1e8, grid_points)]
        fits = [self.at_profiled_variance(values, a) for a in alphas]
        for index in range(3, len(alphas) - 1):
            if fits[index]["log_likelihood"] >= max(
                fits[index - 1]["log_likelihood"], fits[index + 1]["log_likelihood"]
            ):
                result = minimize_scalar(
                    lambda x: (
                        -self.at_profiled_variance(values, np.exp(x))["log_likelihood"]
                    ),
                    bounds=(np.log(alphas[index - 1]), np.log(alphas[index + 1])),
                    method="bounded",
                    options={"xatol": 1e-10},
                )
                alphas.append(float(np.exp(result.x)))
                fits.append(self.at_profiled_variance(values, alphas[-1]))
        maximum = max(fit["log_likelihood"] for fit in fits)
        winner = next(
            i for i, fit in enumerate(fits) if fit["log_likelihood"] >= maximum - 1e-10
        )
        return {
            **fits[winner],
            "alpha_height": alphas[winner],
            "scouting_grid_points": grid_points,
            "continuous_global_optimum_certified": False,
        }
