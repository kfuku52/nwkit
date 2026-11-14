"""Generative native shift models, including joint evolutionary covariance.

Simulation is unconditional on observed trait values. The supplied regime means,
root treatment, process, observation errors and missing mask define the model.
"""

from dataclasses import dataclass

import numpy as np

from nwkit.shift_joint_model import joint_covariance_geometry
from nwkit.shift_native_model import ShiftData, ShiftLayout, ShiftTree, mean_geometry
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
        joint_covariance_geometry(
            self.tree, self.alpha_height, self.covariance_coordinate, self.root_model
        )


def simulation_from_fit(data, result):
    """Construct original-unit generating parameters from an in-memory fit."""
    if "joint_fit" in result:
        joint = result["joint_fit"]
        alpha = joint.alpha_height.copy()
        covariance = (
            joint.covariance_coordinate * data.scales[:, None] * data.scales[None]
        )
        coefficients = joint.coefficients * data.scales[None]
        noise = joint.measurement_variance * data.scales**2
    else:
        alpha = np.array([fit.alpha_height for fit in result["fits"]])
        covariance = np.diag(
            [
                fit.process_variance * data.scales[j] ** 2
                for j, fit in enumerate(result["fits"])
            ]
        )
        coefficients = np.column_stack(
            [fit.coefficients * data.scales[j] for j, fit in enumerate(result["fits"])]
        )
        noise = np.array(
            [
                fit.measurement_variance * data.scales[j] ** 2
                for j, fit in enumerate(result["fits"])
            ]
        )
    coefficients[0] += data.centers
    return ShiftSimulation(
        data.tree,
        result["layout"],
        data.trait_names,
        alpha,
        covariance,
        coefficients,
        result["root_model"],
        data.variances * data.scales[None] ** 2,
        np.diag(noise),
        ~np.isfinite(data.values),
    )


def simulate_shift(spec, replicates=1, *, seed=1):
    """Return observed (replicate,tip,trait) and latent (replicate,node,trait)."""
    spec.validate()
    if (
        isinstance(replicates, bool)
        or not isinstance(replicates, (int, np.integer))
        or replicates < 1
    ):
        raise ValueError("Simulation replicate count must be a positive integer.")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("Simulation seed must be a nonnegative integer.")
    process_seed, error_seed = np.random.SeedSequence(int(seed)).spawn(2)
    slopes, innovations, root, _, _ = joint_covariance_geometry(
        spec.tree, spec.alpha_height, spec.covariance_coordinate, spec.root_model
    )
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
            (1 - slopes[i]) * baseline + weights[i] * offsets[groups[i]],
            innovations[i],
        )
        for i, node in enumerate(spec.tree.compiled.nodes)
        if i
    }
    process = VectorProcess(spec.tree.compiled.tree, p, transitions, baseline, root)
    _, latent = simulate_vector_process(
        process, replicates, seed=int(process_seed.generate_state(1)[0])
    )
    values = latent[:, spec.tree.compiled.leaf_indices, :].copy()
    rng = np.random.default_rng(error_seed)
    for i in range(len(spec.tree.leaf_names)):
        covariance = spec.measurement_covariance + np.diag(spec.sampling_variances[i])
        values[:, i] += rng.multivariate_normal(
            np.zeros(p), covariance, size=replicates, check_valid="raise"
        )
    values[:, spec.missing] = np.nan
    return values, latent


def simulate_joint_data(data, result, rng):
    spec = simulation_from_fit(data, result)
    values, _ = simulate_shift(spec, seed=int(rng.integers(0, 2**63 - 1)))
    return ShiftData.build(
        data.tree, values[0], data.trait_names, spec.sampling_variances
    )
