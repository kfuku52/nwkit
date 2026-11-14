"""Cached covariance profiles for ranking proposed native shift layouts."""

import math

import numpy as np

from nwkit.gaussian_whitening import TreeWhitening
from nwkit.shift_native_model import _normalized_ou_duration, covariance_geometry
from nwkit.shift_native_screen import _effect_design, descendant_design


class NativeQuickProfile:
    def __init__(self, data, null_fit, branches, fixed_alpha=None):
        self.data = data
        self.branches = tuple(branches)
        branches = self.branches
        self.columns = {branch: i for i, branch in enumerate(branches)}
        self.node_indices = {b: i for i, b in enumerate(data.tree.branch_ids)}
        self.age = self._ages()
        design, _ = descendant_design(data.tree, branches)
        dimension = len(data.trait_names)
        alphas = [np.asarray([fit.alpha_height for fit in null_fit["fits"]])]
        if fixed_alpha is None:
            alphas.extend(np.full(dimension, value) for value in (0.3, 3.0))
        else:
            alphas = [
                np.broadcast_to(np.asarray(fixed_alpha, dtype=float), (dimension,))
            ]
        self.profiles = [
            tuple(
                self._trait_profile(j, alpha[j], fit, design)
                for j, fit in enumerate(null_fit["fits"])
            )
            for alpha in alphas
        ]
        self.evaluations = 0

    def _ages(self):
        tree = self.data.tree
        return {
            branch: tree.remaining_times[tree.compiled.parents[i]]
            for branch, i in self.node_indices.items()
            if i
        }

    def _trait_profile(self, trait, alpha, fit, design):
        data = self.data
        if not data.tree.exact_ultrametric:
            design = _effect_design(data.tree, design, alpha, self.branches)
        mask = np.isfinite(data.values[:, trait])
        observed = tuple(
            i
            for i, keep in zip(data.tree.compiled.leaf_indices, mask, strict=True)
            if keep
        )
        slopes, innovations, root_variance = covariance_geometry(
            data.tree, alpha, fit.process_variance, fit.root_model
        )
        factor = TreeWhitening.build(
            data.tree.compiled,
            observed,
            slopes,
            innovations,
            data.variances[mask, trait] + fit.measurement_variance,
            root_variance=root_variance,
        )
        white = factor.apply(
            np.column_stack(
                (np.ones(sum(mask)), data.values[mask, trait], design[mask])
            )
        )
        intercept = white[:, 0] / np.linalg.norm(white[:, 0])
        projected = (
            white[:, 1:] - intercept[:, None] * (intercept @ white[:, 1:])[None, :]
        )
        # A single orthogonal reduction preserves the joint response/design
        # geometry for every later layout. Keep the response column so that its
        # residual outside the candidate span is retained, without normal
        # equations or a subtraction of nearly equal squared norms.
        if projected.shape[0] > projected.shape[1]:
            projected = np.linalg.qr(projected, mode="r")
        constant = -0.5 * (
            sum(mask) * math.log(2 * math.pi) + factor.log_determinant
        ) - sum(mask) * math.log(data.scales[trait])
        return alpha, projected[:, 0], projected[:, 1:], constant, int(sum(mask))

    def _transform(self, layout):
        groups = layout.node_groups(self.data.tree)
        transform = np.zeros((len(layout.shifts), len(layout.groups) - 1))
        for row, branch in enumerate(layout.shifts):
            node = self.node_indices[branch]
            for sign, group in [
                (1, groups[node]),
                (-1, groups[self.data.tree.compiled.parents[node]]),
            ]:
                if group:
                    transform[row, group - 1] += sign
        return transform

    def score(self, layout):
        self.evaluations += 1
        columns = [self.columns[b] for b in layout.shifts]
        ages = np.array([self.age[b] for b in layout.shifts])
        transform = self._transform(layout)
        totals = []
        for profiles in self.profiles:
            total = 0.0
            for alpha, response, matrix, constant, observed_count in profiles:
                selected = matrix[:, columns]
                if self.data.tree.exact_ultrametric and not math.isinf(alpha):
                    selected = selected * _normalized_ou_duration(alpha, ages)
                predictors = selected @ transform
                residual = response
                if predictors.shape[1]:
                    norms = np.linalg.norm(predictors, axis=0)
                    if np.any(norms == 0):
                        total = -math.inf
                        break
                    normalized = predictors / norms
                    singular = np.linalg.svd(normalized, compute_uv=False)
                    # Preserve the original n-tip rank threshold after reducing
                    # the row count; compression must not relax identifiability.
                    tolerance = (
                        singular[0]
                        * max(observed_count, predictors.shape[1])
                        * np.finfo(normalized.dtype).eps
                    )
                    if np.count_nonzero(singular > tolerance) < predictors.shape[1]:
                        total = -math.inf
                        break
                    q, _ = np.linalg.qr(normalized)
                    residual = response - q @ (q.T @ response)
                total += constant - 0.5 * float(residual @ residual)
            totals.append(total)
        return max(totals)
