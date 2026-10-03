"""Bounded conditional MSC fits and explicitly local identifiability diagnostics."""

import math
from dataclasses import dataclass
from itertools import product
from typing import Any

import numpy as np
from scipy.optimize import minimize

from nwkit.mul_msc_model import retime_candidate, score_gene


def species_topology_signature(gene, parser):
    signatures: dict[Any, tuple] = {}
    for node in gene.traverse("postorder"):
        signatures[node] = (
            ("tip", parser.parse(node.name).species_label)
            if node.is_leaf
            else ("node", *sorted(signatures[c] for c in node.children))
        )
    return signatures[gene]


@dataclass(frozen=True)
class FitSettings:
    mode: str
    age_bounds: tuple[float, float] | None = None
    ne_bounds: tuple[float, float] | None = None
    fixed_age: float | None = None
    fixed_scale: float | None = None
    grid_points: int = 5
    starts: int = 3
    maxiter: int = 200
    max_evaluations: int = 5000
    max_states: int = 100000
    max_assignments: int = 10000

    def __post_init__(self):
        if self.mode not in {"age", "ne", "joint"}:
            raise ValueError("MSC fit mode must be age, ne or joint.")
        for bounds, enabled, name in (
            (self.age_bounds, self.mode in {"age", "joint"}, "age"),
            (self.ne_bounds, self.mode in {"ne", "joint"}, "Ne"),
        ):
            if enabled != (bounds is not None):
                raise ValueError(f"MSC {name} bounds must match the fit mode.")
            if bounds is not None and (
                len(bounds) != 2
                or not 0 < bounds[0] < bounds[1]
                or any(not math.isfinite(v) for v in bounds)
            ):
                raise ValueError(f"MSC {name} bounds require finite 0 < lower < upper.")
        for value, needed, name in (
            (self.fixed_age, self.mode == "ne", "age"),
            (self.fixed_scale, self.mode == "age", "time scale"),
        ):
            if needed and (value is None or not math.isfinite(value) or value <= 0):
                raise ValueError(f"MSC fixed {name} must be finite and positive.")
        if self.ne_bounds is not None and not math.isfinite(2 * self.ne_bounds[1]):
            raise ValueError("MSC Ne scale must be finite.")
        if (
            self.grid_points < 3
            or min(
                self.starts,
                self.maxiter,
                self.max_evaluations,
                self.max_states,
                self.max_assignments,
            )
            < 1
        ):
            raise ValueError(
                "MSC fitting requires >=3 grid points and positive limits."
            )


class MscObjective:
    def __init__(self, candidate, genes, parser, settings, *, weights=None):
        self.candidate, self.parser, self.settings = candidate, parser, settings
        supplied = (
            np.ones(len(genes)) if weights is None else np.asarray(weights, dtype=float)
        )
        if (
            supplied.shape != (len(genes),)
            or np.any(~np.isfinite(supplied))
            or np.any(supplied < 0)
            or not np.any(supplied > 0)
        ):
            raise ValueError(
                "MSC family weights must be finite, nonnegative and nonempty."
            )
        groups = {}
        self.indices = []
        self.genes: list[Any] = []
        counts: list[float] = []
        for gene, weight in zip(genes, supplied, strict=True):
            signature = species_topology_signature(gene, parser)
            if signature not in groups:
                groups[signature] = len(self.genes)
                self.genes.append(gene)
                counts.append(0.0)
            index = groups[signature]
            counts[index] += float(weight)
            self.indices.append(index)
        self.weights = np.asarray(counts)
        self.input_weights = supplied
        self.parameters, self.bounds = [], []
        if settings.age_bounds is not None:
            if candidate.age_bounds is None:
                raise ValueError("MSC fitting needs temporally bounded candidates.")
            self.parameters.append("hybridization_age")
            self.bounds.append(candidate.age_bounds)
        if settings.ne_bounds is not None:
            self.parameters.append("effective_population_size")
            self.bounds.append(tuple(math.log(v) for v in settings.ne_bounds))
        self.cache = {}

    def decode(self, coordinates):
        age, scale = self.settings.fixed_age, self.settings.fixed_scale
        values = {}
        for name, position, (lower, upper) in zip(
            self.parameters, coordinates, self.bounds, strict=True
        ):
            value = lower + float(position) * (upper - lower)
            if name == "effective_population_size":
                value = math.exp(value)
                scale = 2 * value
            else:
                age = value
            values[name] = value
        if (
            age is None
            or scale is None
            or not math.isfinite(age)
            or not math.isfinite(scale)
        ):
            raise ValueError("MSC decoded time/Ne scale must be finite.")
        return age, scale, values

    def evaluate(self, coordinates):
        key = tuple(float(value) for value in coordinates)
        if any(not math.isfinite(v) or not 0 <= v <= 1 for v in key):
            raise ValueError("MSC optimizer coordinates must be finite and in [0,1].")
        if key not in self.cache:
            if len(self.cache) >= self.settings.max_evaluations:
                raise ValueError("MSC fitting exceeds --msc-max-evaluations.")
            age, scale, _ = self.decode(key)
            candidate = retime_candidate(self.candidate, age)
            rows = [
                score_gene(
                    gene,
                    candidate,
                    self.parser,
                    time_scale=scale,
                    max_states=self.settings.max_states,
                    max_assignments=self.settings.max_assignments,
                )
                for gene in self.genes
            ]
            logs = np.asarray([row[0] for row in rows])
            total = math.fsum(
                float(w) * float(p) for w, p in zip(self.weights, logs, strict=True)
            )
            if not math.isfinite(total):
                raise ArithmeticError("MSC fitted log likelihood must be finite.")
            self.cache[key] = (total, logs, rows)
        return self.cache[key]

    def loss(self, coordinates):
        return -self.evaluate(coordinates)[0]


def optimize_coordinates(objective, *, fixed=None):
    fixed = {} if fixed is None else fixed
    free = [i for i in range(len(objective.parameters)) if i not in fixed]
    knots = np.linspace(0, 1, objective.settings.grid_points)

    def expand(values):
        result = dict(zip(free, values, strict=True)) | fixed
        return np.asarray([result[i] for i in range(len(objective.parameters))])

    grid = sorted(
        (
            (objective.loss(expand(values)), tuple(values))
            for values in product(knots, repeat=len(free))
        ),
        key=lambda pair: (pair[0], pair[1]),
    )
    if not free:
        return expand(()), [], grid[0][0]
    attempts, successful = [], []
    for _, start in grid[: objective.settings.starts]:
        result = minimize(
            lambda values: objective.loss(expand(values)),
            start,
            method="L-BFGS-B",
            bounds=[(0, 1)] * len(free),
            options={
                "maxiter": objective.settings.maxiter,
                "ftol": 1e-11,
                "gtol": 1e-6,
                "maxls": 30,
            },
        )
        attempts.append(
            {
                "start": list(start),
                "success": bool(result.success),
                "message": str(result.message),
                "iterations": int(result.nit),
                "loss": float(result.fun),
            }
        )
        if result.success and math.isfinite(result.fun):
            successful.append((float(result.fun), tuple(result.x)))
    if not successful:
        raise ArithmeticError(
            "No MSC optimizer start converged; no grid/fixed fallback."
        )
    best = min(successful)
    if best[0] > grid[0][0] + 1e-7:
        raise ArithmeticError("MSC converged fit is worse than its coarse grid.")
    return expand(best[1]), attempts, grid[0][0]


def local_diagnostics(objective, optimum):
    scores = [
        objective.evaluate(point)[0]
        for point in product((0.0, 0.5, 1.0), repeat=len(optimum))
    ]
    if max(scores) - min(scores) <= 1e-8:
        return {
            "status": "flat",
            "rank": 0,
            "singular_values": [],
            "step_ranks": [],
            "boundary": False,
            "meaning": "Likelihood is flat across the evaluated search grid; fitted parameters are withheld.",
        }
    boundary = any(min(value, 1 - value) <= 1e-5 for value in optimum)
    if boundary:
        return {
            "status": "boundary",
            "rank": None,
            "singular_values": [],
            "step_ranks": [],
            "boundary": True,
            "meaning": "Optimum touches supplied/temporal search bounds; no interior point estimate or interval.",
        }
    ranks = []
    singular = np.asarray([], dtype=float)
    for step in (1e-4, 2e-4):
        step = min(step, min(min(v, 1 - v) for v in optimum) / 2)
        columns = []
        for i in range(len(optimum)):
            left, right = optimum.copy(), optimum.copy()
            left[i] -= step
            right[i] += step
            columns.append(
                (objective.evaluate(right)[1] - objective.evaluate(left)[1])
                / (2 * step)
            )
        jacobian = np.column_stack(columns) * np.sqrt(objective.weights)[:, None]
        singular = np.linalg.svd(jacobian, compute_uv=False)
        cutoff = max(1e-9, (float(singular[0]) if len(singular) else 0) * 1e-7)
        ranks.append(int(np.count_nonzero(singular > cutoff)))
    status = (
        "locally-distinguishable"
        if min(ranks) == len(optimum)
        else "locally-unidentified"
    )
    return {
        "status": status,
        "rank": min(ranks),
        "singular_values": singular.tolist(),
        "step_ranks": ranks,
        "boundary": False,
        "meaning": "Two-step finite-difference rank of observed-pattern log probabilities, in normalized search coordinates; local sample diagnostic, not proof of global identifiability or calibrated uncertainty.",
    }


def fit_candidate(candidate, genes, parser, settings, *, weights=None):
    objective = MscObjective(candidate, genes, parser, settings, weights=weights)
    optimum, attempts, grid_loss = optimize_coordinates(objective)
    total, _, grouped = objective.evaluate(optimum)
    diagnostics = local_diagnostics(objective, optimum)
    profiles = []
    for i, name in enumerate(objective.parameters):
        positions = sorted(
            set(np.linspace(0, 1, settings.grid_points)) | {float(optimum[i])}
        )
        for position in positions:
            solution, profile_attempts, _ = optimize_coordinates(
                objective, fixed={i: float(position)}
            )
            score = objective.evaluate(solution)[0]
            if score > total + 1e-6:
                raise ArithmeticError(
                    "MSC nuisance profile exceeds fitted maximum; increase grid/starts."
                )
            _, _, values = objective.decode(solution)
            profiles.append(
                {
                    "mul.tree": candidate.id,
                    "parameter": name,
                    "value": values[name],
                    "log_likelihood": score,
                    "delta_log_likelihood": max(0.0, total - score),
                    "nuisance_parameters": {
                        k: v for k, v in values.items() if k != name
                    },
                    "optimizer_attempts": profile_attempts,
                }
            )
    age, scale, values = objective.decode(optimum)
    estimates = (
        values
        if diagnostics["status"] == "locally-distinguishable"
        else {name: None for name in values}
    )
    return {
        "id": candidate.id,
        "log_likelihood": total,
        "candidate": retime_candidate(candidate, age),
        "time_scale": scale,
        "estimates": estimates,
        "diagnostics": diagnostics,
        "numerical_solution": values,
        "normalized_coordinates": optimum.tolist(),
        "optimizer_attempts": attempts,
        "coarse_grid_log_likelihood": -grid_loss,
        "evaluations": len(objective.cache),
        "profiles": profiles,
        "gene_rows": [
            {
                "mul.tree": candidate.id,
                "gene.tree": i + 1,
                "log_likelihood": grouped[index][0],
                "num_assignments": grouped[index][1],
                "coalescent_states": grouped[index][2],
                "family_weight": float(objective.input_weights[i]),
            }
            for i, index in enumerate(objective.indices)
        ],
    }
