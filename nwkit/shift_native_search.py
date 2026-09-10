"""Shared native search accounting and exhaustive fixed-layout reference."""

import itertools
import math
from dataclasses import dataclass

import numpy as np

from nwkit.shift_candidates import set_partitions
from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_model import ShiftLayout


def layout_complexity(layout):
    """Prespecified nested families: number of locations plus free regime offsets."""
    return len(layout.shifts) + len(layout.groups) - 1


def observable_layout(data, layout, alpha_height=1.0):
    # This structural screen depends on the observation mask, not trait values.
    alphas = np.broadcast_to(
        np.asarray(alpha_height, dtype=float), (len(data.trait_names),)
    )
    checked = set()
    for trait in range(len(data.trait_names)):
        observed = np.isfinite(data.values[:, trait])
        key = (float(alphas[trait]), observed.tobytes())
        if key in checked:
            continue
        design = layout.design(data.tree, alphas[trait])
        rows = design[observed]
        norms = np.linalg.norm(rows, axis=0)
        if (
            len(rows) <= rows.shape[1]
            or np.any(norms == 0)
            or np.linalg.matrix_rank(rows / norms) < rows.shape[1]
        ):
            return False
        checked.add(key)
    return True


def _bell_number(size):
    row = [1]
    for _ in range(size):
        following = [row[-1]]
        following.extend(0 for _ in row)
        for i, value in enumerate(row):
            following[i + 1] = following[i] + value
        row = following
    return row[0]


def enumerate_native_layouts(data, max_shifts=2, *, convergence=False, limit=5000):
    branches = tuple(b for b in data.tree.branch_ids if b)
    if (
        isinstance(max_shifts, bool)
        or not isinstance(max_shifts, int)
        or not 0 <= max_shifts < len(data.tree.leaf_names) - 1
    ):
        raise ValueError(
            "Maximum shifts must be a nonnegative integer smaller than tips minus one."
        )
    if limit < 1:
        raise ValueError("Exhaustive candidate limit must be positive.")
    upper = sum(
        math.comb(len(branches), k) * (_bell_number(k + 1) if convergence else 1)
        for k in range(max_shifts + 1)
    )
    if upper > max(10000, 20 * limit):
        raise ValueError(
            "Exhaustive traversal budget exceeded before enumeration; use native heuristic search or explicitly increase the candidate limit."
        )
    candidates, excluded, considered = [], 0, 0
    for size in range(max_shifts + 1):
        for shifts in itertools.combinations(branches, size):
            groups = set_partitions((0, *shifts)) if convergence else [None]
            for group in groups:
                considered += 1
                try:
                    layout = ShiftLayout.build(data.tree, shifts, group)
                except ValueError:
                    excluded += 1
                    continue
                if not observable_layout(data, layout):
                    excluded += 1
                    continue
                candidates.append(layout)
                if len(candidates) > limit:
                    raise ValueError(
                        "Exhaustive candidate limit exceeded before fitting."
                    )
    candidates.sort(
        key=lambda layout: (layout_complexity(layout), layout.shifts, layout.groups)
    )
    return candidates, {
        "considered": considered,
        "structurally_excluded": excluded,
        "eligible": len(candidates),
        "complete_discrete_enumeration": True,
    }


@dataclass
class NativeSearchResult:
    records: list[dict]
    best_by_complexity: dict[int, dict]
    metadata: dict

    def families(self):
        families, best = [], None
        for complexity, result in sorted(self.best_by_complexity.items()):
            if best is None or result["log_likelihood"] > best["log_likelihood"] + 1e-9:
                best = result
            families.append((complexity, best))
        return families


class NativeLayoutEvaluator:
    """Retain scalar candidate records and at most one full fit per complexity."""

    def __init__(self, data, fit_arguments):
        self.data = data
        self.fit_arguments = fit_arguments
        self.records: list[dict] = []
        self.scores: dict[ShiftLayout, float] = {}
        self.best: dict[int, dict] = {}

    def evaluate(self, layout):
        if layout in self.scores:
            return self.scores[layout]
        fixed_alpha = self.fit_arguments.get("alpha_height")
        if fixed_alpha is not None and not observable_layout(
            self.data, layout, fixed_alpha
        ):
            self.scores[layout] = -math.inf
            self.records.append(
                {
                    "shift_branch_ids": list(layout.shifts),
                    "groups": [list(g) for g in layout.groups],
                    "complexity": layout_complexity(layout),
                    "log_likelihood": None,
                    "evaluations": 0,
                    "status": "structurally_excluded_at_fixed_alpha",
                }
            )
            return -math.inf
        result = fit_native_layout(self.data, layout, **self.fit_arguments)
        if any(
            trait["optimizer"]["nuisance_variance_at_numerical_bound"]
            for trait in result["traits"]
        ):
            raise ValueError(
                "Native search reached a numerical variance bound; no unconstrained variance optimum is established."
            )
        unresolved = [
            mode
            for trait in result["traits"]
            for mode in trait["optimizer"]["alpha_candidates"]
            if not mode["success"]
            and "rank deficient" not in str(mode.get("message", ""))
        ]
        if unresolved:
            raise ValueError(
                "A native search covariance mode failed; complete the fit before selection: "
                + str(unresolved[0])
            )
        score = result["log_likelihood"]
        self.scores[layout] = score
        complexity = layout_complexity(layout)
        self.records.append(
            {
                "shift_branch_ids": list(layout.shifts),
                "groups": [list(g) for g in layout.groups],
                "complexity": complexity,
                "log_likelihood": score,
                "evaluations": sum(
                    r["optimizer"]["evaluations"] for r in result["traits"]
                ),
                "complete_covariance_modes": not unresolved,
            }
        )
        if (
            complexity not in self.best
            or score > self.best[complexity]["log_likelihood"] + 1e-9
        ):
            self.best[complexity] = result
        return score

    def finish(self, metadata):
        return NativeSearchResult(self.records, self.best, metadata)


def exhaustive_native_search(
    data, *, max_shifts=2, convergence=False, limit=5000, fit_arguments=None
):
    layouts, metadata = enumerate_native_layouts(
        data, max_shifts, convergence=convergence, limit=limit
    )
    evaluator = NativeLayoutEvaluator(data, fit_arguments or {})
    for layout in layouts:
        evaluator.evaluate(layout)
    return evaluator.finish(
        {
            "strategy": "exhaustive",
            **metadata,
            "continuous_global_optimum_certified": False,
        }
    )
