"""Generative native shift models, including joint evolutionary covariance.

Simulation is unconditional on observed trait values. The supplied regime means,
root treatment, process, observation errors and missing mask define the model.
"""

from dataclasses import dataclass

import numpy as np

from nwkit.shift_joint_model import joint_covariance_geometry
from nwkit.shift_native_model import (
    ShiftData,
    ShiftLayout,
    ShiftTree,
    _restore_trait_variance,
    mean_geometry,
    restore_trait_covariance,
)
from nwkit.vector_gaussian import VectorProcess, VectorTransition
from nwkit.vector_simulation import simulate_vector_process


@dataclass(frozen=True)
class ShiftSimulation:
    tree: ShiftTree
    layout: ShiftLayout
    trait_names: tuple
    alpha_height: np.ndarray
    covariance_coordinate: np.ndarray
    coefficients: np.ndarray
    root_model: str
    sampling_variances: np.ndarray
    measurement_covariance: np.ndarray
    missing: np.ndarray

    def validate(self):
        self._validated_geometry()

    def _validated_geometry(self):
        n, p = len(self.tree.leaf_names), len(self.trait_names)
        if (
            np.shape(self.alpha_height) != (p,)
            or np.asarray(self.missing).dtype != bool
        ):
            raise ValueError(
                "Simulation alpha/missing mask must match trait dimensions and use boolean missingness."
            )
        ShiftLayout.build(self.tree, self.layout.shifts, self.layout.groups)
        if (
            not p
            or len(set(self.trait_names)) != p
            or any(not isinstance(x, str) or not x for x in self.trait_names)
        ):
            raise ValueError("Simulation traits must have distinct nonempty names.")
        for value, shape, label in (
            (self.coefficients, (len(self.layout.groups), p), "Regime coefficients"),
            (self.sampling_variances, (n, p), "Sampling variances"),
            (self.measurement_covariance, (p, p), "Measurement covariance"),
        ):
            if np.shape(value) != shape or not np.isfinite(value).all():
                raise ValueError(f"{label} must be finite with shape {shape}.")
        if np.any(self.sampling_variances < 0) or np.shape(self.missing) != (n, p):
            raise ValueError(
                "Simulation variances/missing mask must be tip/trait aligned."
            )
        if not np.allclose(
            self.measurement_covariance,
            self.measurement_covariance.T,
            atol=1e-13,
            rtol=1e-10,
        ):
            raise ValueError("Measurement covariance must be symmetric.")
        eigenvalues = np.linalg.eigvalsh(self.measurement_covariance)
        if eigenvalues.min() < -1e-12 * max(1.0, np.max(np.abs(eigenvalues))):
            raise ValueError("Measurement covariance must be positive semidefinite.")
        return joint_covariance_geometry(
            self.tree, self.alpha_height, self.covariance_coordinate, self.root_model
        )


def simulation_from_fit(data, result):
    """Construct original-unit generating parameters from an in-memory fit."""
    if "joint_fit" in result:
        joint = result["joint_fit"]
        alpha = joint.alpha_height.copy()
        covariance = restore_trait_covariance(joint.covariance_coordinate, data.scales)
        coefficients = joint.coefficients * data.scales[None]
        noise = joint.measurement_variance
    else:
        alpha = np.array([fit.alpha_height for fit in result["fits"]])
        variance = np.array([fit.process_variance for fit in result["fits"]])
        covariance = np.diag(_restore_trait_variance(variance, data.scales))
        coefficients = np.column_stack(
            [fit.coefficients * data.scales[j] for j, fit in enumerate(result["fits"])]
        )
        noise = np.array([fit.measurement_variance for fit in result["fits"]])
    coefficients[0] += data.centers
    return ShiftSimulation(
        data.tree,
        result["layout"],
        data.trait_names,
        alpha,
        covariance,
        coefficients,
        result["root_model"],
        _restore_trait_variance(data.variances, data.scales),
        np.diag(_restore_trait_variance(noise, data.scales)),
        ~np.isfinite(data.values),
    )


def simulate_shift(spec, replicates=1, *, seed=1):
    """Return observed (replicate,tip,trait) and latent (replicate,node,trait)."""
    slopes, innovations, root, _, _ = spec._validated_geometry()
    if (
        isinstance(replicates, bool)
        or not isinstance(replicates, (int, np.integer))
        or replicates < 1
    ):
        raise ValueError("Simulation replicate count must be a positive integer.")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("Simulation seed must be a nonnegative integer.")
    process_seed, error_seed = np.random.SeedSequence(int(seed)).spawn(2)
    p = len(spec.trait_names)
    weights = np.column_stack(
        [mean_geometry(spec.tree, float(alpha))[1] for alpha in spec.alpha_height]
    )
    groups = spec.layout.node_groups(spec.tree)
    baseline = spec.coefficients[0]
    offsets = spec.coefficients.copy()
    offsets[0] = 0
    transitions = {
        node: VectorTransition(
            np.diag(slopes[i]),
            weights[i] * offsets[groups[i]],
            innovations[i],
        )
        for i, node in enumerate(spec.tree.compiled.nodes)
        if i
    }
    # Propagate deviations from the background optimum so small changes can
    # accumulate over several branches before returning to original units.
    process = VectorProcess(spec.tree.compiled.tree, p, transitions, np.zeros(p), root)
    _, latent = simulate_vector_process(
        process,
        replicates,
        seed=int(process_seed.generate_state(1)[0]),
        compiled_tree=spec.tree.compiled,
    )
    values = latent[:, spec.tree.compiled.leaf_indices, :].copy()
    if np.any(spec.measurement_covariance) or np.any(spec.sampling_variances):
        rng = np.random.default_rng(error_seed)
        zero = np.zeros(p)
        for i, sampling in enumerate(spec.sampling_variances):
            covariance = spec.measurement_covariance + np.diag(sampling)
            values[:, i] += rng.multivariate_normal(
                zero, covariance, size=replicates, check_valid="raise"
            )
    latent += baseline
    values += baseline
    values[:, spec.missing] = np.nan
    return values, latent


def simulate_joint_data(data, result, rng):
    spec = simulation_from_fit(data, result)
    values, _ = simulate_shift(spec, seed=int(rng.integers(0, 2**63 - 1)))
    return ShiftData.build(
        data.tree, values[0], data.trait_names, spec.sampling_variances
    )
