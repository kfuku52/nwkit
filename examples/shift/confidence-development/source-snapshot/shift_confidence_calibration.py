"""Research-only Berger--Boos envelope for the no-error, no-shift family.

The confidence set is obtained by inverting finite Monte Carlo likelihood-ratio
tests on the declared alpha grid. No continuous-alpha or later-stage guarantee
is implied. This module does not change the production selection rule.
"""

import hashlib

import numpy as np
from scipy.linalg import solve_triangular

from nwkit.shift_calibration import validate_calibration_options


def restricted_envelope(
    shift_probabilities, confidence_probabilities, beta, level=0.05
):
    shift = np.asarray(shift_probabilities, dtype=float)
    confidence = np.asarray(confidence_probabilities, dtype=float)
    if (
        shift.ndim != 1
        or not len(shift)
        or confidence.shape != shift.shape
        or not np.all(np.isfinite(shift))
        or not np.all(np.isfinite(confidence))
        or np.any((shift <= 0) | (shift > 1))
        or np.any((confidence <= 0) | (confidence > 1))
        or not 0 <= beta < level < 1
    ):
        raise ValueError("Invalid probabilities or confidence-set error allocation")
    retained = confidence > beta
    maximum = float(shift[retained].max()) if np.any(retained) else 0.0
    p_value = min(1.0, maximum + beta)
    return dict(
        beta=beta,
        confidence_level=1 - beta,
        retained_indices=np.flatnonzero(retained).tolist(),
        retained_count=int(retained.sum()),
        empty_confidence_set=not bool(retained.any()),
        restricted_maximum=maximum,
        p_value=p_value,
        reject=p_value <= level,
        scope="No-error no-shift family and declared finite alpha grid only",
    )


class NullConfidenceBank:
    """Full-search and confidence-test draws, reusable within one paired block."""

    def __init__(self, search):
        if search.known_error or len(search.cache) != len(search.grid):
            raise ValueError("Confidence prototype currently requires no known errors")
        if len(search.families[0]) != 1 or search.dim[search.families[0][0]] != 0:
            raise ValueError("The first family must be the no-shift model")
        self.search = search
        self.null_index = int(search.families[0][0])

    def null_profile(self, Z):
        """Null-only profiled Gaussian density, without the mean-search helpers."""
        Z = np.asarray(Z, dtype=float)
        if Z.ndim == 1:
            Z = Z[:, None]
        if Z.ndim != 2 or Z.shape[0] != self.search.d or not np.all(np.isfinite(Z)):
            raise ValueError("Finite contrast-aligned values are required")
        result = []
        for item in self.search.cache:
            w = solve_triangular(item[2], Z, lower=True)
            rss = np.sum(w * w, axis=0)
            if np.any(rss <= 0):
                raise ValueError("Null likelihood has zero residual process variance")
            result.append(
                -0.5
                * (
                    self.search.d
                    * (np.log(2 * np.pi) + 1 + np.log(rss / self.search.d))
                    + item[5]
                )
            )
        return np.array(result)

    def simulate(self, seed, replicates=999):
        validate_calibration_options(replicates, 0.05)
        noise = np.random.default_rng(seed).normal(size=(self.search.d, replicates))
        shift_rows, confidence_rows = [], []
        for index, item in enumerate(self.search.cache):
            shifts, confidences = [], []
            for start in range(0, replicates, 64):
                Z = item[2] @ noise[:, start : start + 64]
                best, _ = self.search.profile(Z)
                null = self.null_profile(Z)
                shifts.extend(2 * (best.max(axis=0) - best[self.null_index]))
                confidences.extend(2 * (null.max(axis=0) - null[index]))
            shift_rows.append(shifts)
            confidence_rows.append(confidences)
        self.shift_statistics = np.asarray(shift_rows)
        self.confidence_statistics = np.asarray(confidence_rows)
        self.replicates = replicates
        self.seed = seed
        return dict(
            seed=seed,
            replicates=replicates,
            generator_points=len(self.search.cache),
            fitted_bootstrap_datasets=len(self.search.cache) * replicates,
            noise_sha256=hashlib.sha256(noise.tobytes()).hexdigest(),
            shift_statistics_sha256=hashlib.sha256(
                self.shift_statistics.tobytes()
            ).hexdigest(),
            confidence_statistics_sha256=hashlib.sha256(
                self.confidence_statistics.tobytes()
            ).hexdigest(),
        )

    def evaluate(self, values, betas=(0.005, 0.01)):
        values = np.asarray(values, dtype=float)
        if values.shape != (self.search.n,) or not np.all(np.isfinite(values)):
            raise ValueError("Finite tip-aligned values are required")
        z = self.search.q @ (values - values.mean())
        best, at = self.search.profile(z)
        observed = float(2 * (best[:, 0].max() - best[self.null_index, 0]))
        null = self.null_profile(z)[:, 0]
        confidence_statistics = 2 * (null.max() - null)
        shift_counts = np.count_nonzero(
            self.shift_statistics >= observed - 1e-10, axis=1
        )
        confidence_counts = np.count_nonzero(
            self.confidence_statistics >= confidence_statistics[:, None] - 1e-10, axis=1
        )
        shift_p = (1 + shift_counts) / (self.replicates + 1)
        confidence_p = (1 + confidence_counts) / (self.replicates + 1)
        best_index = int(np.argmax(best[:, 0]))
        alpha_index = int(at[best_index, 0])
        return dict(
            shift_statistic=observed,
            null_log_likelihood=best[self.null_index, 0],
            full_log_likelihood=best[best_index, 0],
            full_winner=best_index,
            full_model=self.search.models[best_index],
            full_alpha_index=alpha_index,
            null_grid_log_likelihood=null.tolist(),
            confidence_statistics=confidence_statistics.tolist(),
            shift_exceedances=shift_counts.tolist(),
            confidence_exceedances=confidence_counts.tolist(),
            shift_probabilities=shift_p.tolist(),
            confidence_probabilities=confidence_p.tolist(),
            full_envelope_p=float(shift_p.max()),
            full_envelope_reject=bool(shift_p.max() <= 0.05),
            plugin_p=float(shift_p[at[self.null_index, 0]]),
            candidates=[restricted_envelope(shift_p, confidence_p, b) for b in betas],
        )
