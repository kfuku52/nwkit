"""Native diagonal-trait OU models with shared, fixed shift/regime layouts.

Each trait has its own covariance parameters and mean coefficients. Layouts
are shared across traits; this is not a full evolutionary trait covariance.
Scaled regime offsets have well-defined Brownian and independent-tip limits.
"""

import math
from dataclasses import dataclass

import numpy as np

from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_whitening import TreeWhitening, tree_gls
from nwkit.util import assign_branch_ids


@dataclass(frozen=True)
class ShiftTree:
    compiled: CompiledTree
    branch_ids: tuple[int, ...]
    height: float
    times: np.ndarray
    levels: tuple[np.ndarray, ...]
    indices_by_branch: tuple[int, ...]
    tip_intervals: tuple[tuple[int, int], ...]
    remaining_times: np.ndarray
    exact_ultrametric: bool

    @classmethod
    def build(cls, tree):
        compiled = CompiledTree.from_tree(tree)
        if len(compiled.leaf_indices) < 3:
            raise ValueError("Native shift inference needs at least three tips.")
        if any(not str(node.name) for node in tree.leaves()):
            raise ValueError("Native shift inference needs nonempty tip names.")
        if any(len(children) not in (0, 2) for children in compiled.children):
            raise ValueError("Native shift inference needs a binary rooted tree.")
        times = np.zeros(len(compiled.nodes))
        depths = times.copy()
        edge_depths = np.zeros(len(compiled.nodes), dtype=int)
        levels: dict[int, list[int]] = {}
        for index, node in enumerate(compiled.nodes[1:], 1):
            times[index] = float(node.dist)
            if not math.isfinite(times[index]) or times[index] <= 0:
                raise ValueError("Native shift branches must be finite and positive.")
            depths[index] = depths[compiled.parents[index]] + times[index]
            edge_depths[index] = edge_depths[compiled.parents[index]] + 1
            levels.setdefault(int(edge_depths[index]), []).append(index)
        heights = depths[list(compiled.leaf_indices)]
        height = float(max(heights))
        if not np.isfinite(depths).all() or not height > 0:
            raise ValueError("Native shift tree height must be finite and positive.")
        if np.ptp(heights) > height * 1e-8:
            raise ValueError(
                "Native shift inference needs an ultrametric tree; no repair is applied."
            )
        ids = assign_branch_ids(tree)
        branch_ids = tuple(ids[node] for node in compiled.nodes)
        intervals = {node: (i, i + 1) for i, node in enumerate(compiled.leaf_indices)}
        remaining = np.zeros(len(compiled.nodes))
        for i in compiled.postorder:
            children = compiled.children[i]
            if children:
                intervals[i] = (intervals[children[0]][0], intervals[children[-1]][1])
            if i:
                parent = compiled.parents[i]
                remaining[parent] = max(
                    remaining[parent], times[i] / height + remaining[i]
                )
        return cls(
            compiled,
            branch_ids,
            height,
            times / height,
            tuple(np.asarray(levels[level], dtype=int) for level in sorted(levels)),
            tuple(int(i) for i in np.argsort(branch_ids)),
            tuple(intervals[i] for i in range(len(compiled.nodes))),
            remaining,
            bool(np.ptp(heights) == 0),
        )

    @property
    def leaf_names(self):
        return tuple(
            str(self.compiled.nodes[i].name) for i in self.compiled.leaf_indices
        )


@dataclass(frozen=True)
class ShiftLayout:
    shifts: tuple[int, ...]
    groups: tuple[tuple[int, ...], ...]

    @classmethod
    def build(cls, tree: ShiftTree, shifts=(), groups=None):
        if any(
            isinstance(b, bool) or not isinstance(b, (int, np.integer)) for b in shifts
        ):
            raise ValueError("Shift branches must be integer branch IDs.")
        shifts = tuple(sorted(int(b) for b in shifts))
        if len(set(shifts)) != len(shifts) or any(
            b == 0 or b not in tree.branch_ids for b in shifts
        ):
            raise ValueError("Shift branches must be distinct non-root branch IDs.")
        if groups is None:
            groups = ((b,) for b in (0, *shifts))
        groups = tuple(tuple(sorted(group)) for group in groups)
        if any(not group for group in groups):
            raise ValueError("Regime groups cannot be empty.")
        members = [b for group in groups for b in group]
        if any(
            isinstance(b, bool) or not isinstance(b, (int, np.integer)) for b in members
        ):
            raise ValueError("Regime groups must contain integer branch IDs.")
        if sorted(members) != [0, *shifts]:
            raise ValueError(
                "Regime groups must partition the background and all shifts."
            )
        result = cls(shifts, tuple(sorted(groups)))
        labels = result.node_groups(tree)
        for branch in shifts:
            index = tree.indices_by_branch[branch]
            if labels[index] == labels[tree.compiled.parents[index]]:
                raise ValueError("A shift cannot retain its immediate parent's regime.")
        return result

    def node_groups(self, tree):
        aliases = {
            branch: index for index, group in enumerate(self.groups) for branch in group
        }
        labels = np.zeros(len(tree.branch_ids), dtype=int)
        # Preorder makes each subtree contiguous. Apply ancestral regimes first,
        # then overwrite their descendant intervals with each nested shift.
        for index, group in sorted(
            (tree.indices_by_branch[branch], group)
            for branch, group in aliases.items()
        ):
            last_tip = tree.tip_intervals[index][1] - 1
            stop = tree.compiled.leaf_indices[last_tip] + 1
            labels[index:stop] = group
        return labels

    def design(self, tree, alpha_height, *, node_groups=None):
        if math.isnan(alpha_height) or alpha_height < 0:
            raise ValueError("Alpha times tree height must be nonnegative.")
        labels = self.node_groups(tree) if node_groups is None else node_groups
        if not tree.exact_ultrametric:
            # Trees within the input rounding tolerance still use their actual
            # branch times. The constant-age interval shortcut is inapplicable.
            slopes, weights = mean_geometry(tree, alpha_height)
            matrix = np.zeros((len(tree.branch_ids), len(self.groups)))
            matrix[:, 0] = 1
            parents = np.asarray(tree.compiled.parents)
            for indices in tree.levels:
                matrix[indices, 1:] = (
                    slopes[indices, None] * matrix[parents[indices], 1:]
                )
                shifted = indices[labels[indices] != 0]
                matrix[shifted, labels[shifted]] += weights[shifted]
            return matrix[list(tree.compiled.leaf_indices)]
        matrix = np.zeros((len(tree.compiled.leaf_indices), len(self.groups)))
        matrix[:, 0] = 1.0
        # On an ultrametric tree each shift contributes the same age-weighted
        # jump to its contiguous descendant-tip interval. Nested shifts subtract
        # the inherited regime, exactly reproducing branchwise OU propagation.
        for branch in self.shifts:
            index = tree.indices_by_branch[branch]
            parent = tree.compiled.parents[index]
            age = tree.remaining_times[parent]
            weight = (
                age
                if alpha_height == 0
                else (
                    1.0
                    if math.isinf(alpha_height)
                    else -math.expm1(-alpha_height * age) / -math.expm1(-alpha_height)
                )
            )
            first, last = tree.tip_intervals[index]
            if labels[index]:
                matrix[first:last, labels[index]] += weight
            if labels[parent]:
                matrix[first:last, labels[parent]] -= weight
        return matrix


def mean_geometry(tree, alpha_height):
    alpha_height = float(alpha_height)
    if math.isnan(alpha_height) or alpha_height < 0:
        raise ValueError("Alpha times tree height must be nonnegative.")
    if alpha_height == 0:
        return np.ones(len(tree.times)), tree.times.copy()
    if math.isinf(alpha_height):
        return np.zeros(len(tree.times)), np.ones(len(tree.times))
    return np.exp(-alpha_height * tree.times), -np.expm1(
        -alpha_height * tree.times
    ) / -np.expm1(-alpha_height)


def covariance_geometry(tree, alpha_height, process_variance, root_model):
    if root_model not in {"OUfixedRoot", "OUrandomRoot"}:
        raise ValueError("Unknown native OU root model.")
    process_variance = float(process_variance)
    if not math.isfinite(process_variance) or process_variance < 0:
        raise ValueError("Process tip variance must be finite and nonnegative.")
    slopes, _ = mean_geometry(tree, alpha_height)
    if alpha_height == 0:
        if root_model == "OUrandomRoot":
            raise ValueError("Stationary-root OU is undefined at alpha zero.")
        return slopes, process_variance * tree.times, 0.0
    denominator = (
        1.0
        if math.isinf(alpha_height) or root_model == "OUrandomRoot"
        else -np.expm1(-2 * alpha_height)
    )
    fractions = (
        np.ones(len(tree.times))
        if math.isinf(alpha_height)
        else -np.expm1(-2 * alpha_height * tree.times)
    )
    innovations = process_variance * (fractions / denominator)
    root_variance = process_variance if root_model == "OUrandomRoot" else 0.0
    return slopes, innovations, root_variance


@dataclass(frozen=True)
class ShiftData:
    tree: ShiftTree
    trait_names: tuple[str, ...]
    values: np.ndarray
    variances: np.ndarray
    centers: np.ndarray
    scales: np.ndarray

    @classmethod
    def build(cls, tree, values, trait_names, variances=None):
        tree = tree if isinstance(tree, ShiftTree) else ShiftTree.build(tree)
        names = tuple(trait_names)
        if (
            not names
            or len(set(names)) != len(names)
            or any(not isinstance(x, str) or not x for x in names)
        ):
            raise ValueError("Trait names must be distinct nonempty strings.")
        values = np.asarray(values, dtype=float)
        if values.shape != (len(tree.leaf_names), len(names)) or np.isinf(values).any():
            raise ValueError(
                "Trait values must match tips by traits; missing coordinates use NaN."
            )
        errors = (
            np.zeros_like(values)
            if variances is None
            else np.asarray(variances, dtype=float)
        )
        if (
            errors.shape != values.shape
            or not np.isfinite(errors).all()
            or np.any(errors < 0)
        ):
            raise ValueError(
                "Known observation variances must be finite, nonnegative and tip/trait aligned."
            )
        if np.any(np.sum(np.isfinite(values), axis=0) < 3):
            raise ValueError(
                "Each native shift trait needs at least three observed tips."
            )
        centers = np.nanmean(values, axis=0)
        centered = values - centers
        scales = np.nanmax(np.abs(centered), axis=0)
        scales = np.maximum(scales, np.sqrt(np.max(errors, axis=0)))
        if np.any(scales == 0) or not np.isfinite(scales).all():
            raise ValueError(
                "A constant error-free trait has no regular Gaussian variance fit."
            )
        normalized_errors = errors / scales / scales
        if not np.isfinite(normalized_errors).all() or np.any(
            (errors > 0) & (normalized_errors == 0)
        ):
            raise ValueError("Trait/error scale is not representable.")
        return cls(tree, names, centered / scales, normalized_errors, centers, scales)


@dataclass(frozen=True)
class NativeTraitFit:
    alpha_height: float
    process_variance: float
    measurement_variance: float
    coefficients: np.ndarray
    coefficient_covariance: np.ndarray
    predicted: np.ndarray
    log_likelihood: float
    quadratic: float
    num_observations: int
    num_mean_parameters: int
    root_model: str


def evaluate_trait(
    data,
    layout,
    trait,
    alpha_height,
    process_variance,
    measurement_variance=0.0,
    *,
    root_model="OUfixedRoot",
    _node_groups=None,
):
    """Fit the mean at supplied *normalized* covariance parameters.

    Output stays normalized for optimization. ``restore_trait_fit`` restores
    original trait units, likelihood density units, and physical tree time.
    """
    if not math.isfinite(measurement_variance) or measurement_variance < 0:
        raise ValueError(
            "Additional measurement variance must be finite and nonnegative."
        )
    mask = np.isfinite(data.values[:, trait])
    observed = tuple(
        i for i, keep in zip(data.tree.compiled.leaf_indices, mask, strict=True) if keep
    )
    design = layout.design(data.tree, alpha_height, node_groups=_node_groups)
    slopes, innovations, root_variance = covariance_geometry(
        data.tree, alpha_height, process_variance, root_model
    )
    factor = TreeWhitening.build(
        data.tree.compiled,
        observed,
        slopes,
        innovations,
        data.variances[mask, trait] + measurement_variance,
        root_variance=root_variance,
    )
    coefficients, likelihood, quadratic, covariance = tree_gls(
        factor, data.values[mask, trait], design[mask]
    )
    return NativeTraitFit(
        float(alpha_height),
        float(process_variance),
        float(measurement_variance),
        coefficients,
        covariance,
        design @ coefficients,
        likelihood,
        quadratic,
        len(observed),
        design.shape[1],
        root_model,
    )


def restore_trait_fit(data, trait, fit):
    scale, center = data.scales[trait], data.centers[trait]
    alpha = fit.alpha_height
    status = (
        "brownian_limit"
        if alpha == 0
        else ("independent_limit" if math.isinf(alpha) else "finite")
    )
    coefficients = fit.coefficients * scale
    coefficients[0] += center
    sigma2 = None
    if not math.isinf(alpha):
        ratio = 1.0 if alpha == 0 else 2 * alpha / -np.expm1(-2 * alpha)
        if fit.root_model == "OUrandomRoot":
            ratio = 2 * alpha
        sigma2 = fit.process_variance * scale**2 * ratio / data.tree.height
    return {
        "trait": data.trait_names[trait],
        "alpha": None if math.isinf(alpha) else alpha / data.tree.height,
        "alpha_height": None if math.isinf(alpha) else alpha,
        "alpha_status": status,
        "sigma2": sigma2,
        "process_tip_variance": fit.process_variance * scale**2,
        "measurement_variance": fit.measurement_variance * scale**2,
        "scaled_regime_coefficients": coefficients.tolist(),
        "coefficient_covariance": (fit.coefficient_covariance * scale**2).tolist(),
        "predicted": (center + scale * fit.predicted).tolist(),
        "log_likelihood": fit.log_likelihood - fit.num_observations * math.log(scale),
        "num_observations": fit.num_observations,
        "num_mean_parameters": fit.num_mean_parameters,
    }
