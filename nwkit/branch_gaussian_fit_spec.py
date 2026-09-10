"""Explicit parameter sharing and finite bounds for fixed branch assignments."""

import math
from dataclasses import dataclass, replace

import numpy as np

from nwkit.branch_gaussian import OUBranch
from nwkit.branch_gaussian_input import BranchGaussianAssignment, _read_rows

_FIELDS = {"sigma2": "variance_rate", "alpha": "alpha", "theta": "optimum"}


@dataclass(frozen=True)
class BranchFitParameter:
    group: str
    parameter: str
    regimes: tuple[str, ...]
    initial: float
    lower: float
    upper: float
    time_scale: float = 1.0

    def coordinate(self, value):
        if self.parameter == "sigma2":
            return (math.log(value) - math.log(self.lower)) / (
                math.log(self.upper) - math.log(self.lower)
            )
        if self.parameter == "alpha" and self.lower > 0:
            return (math.log(value) - math.log(self.lower)) / (
                math.log(self.upper) - math.log(self.lower)
            )
        if self.parameter == "alpha":
            if value == 0:
                return 0.0
            return float(
                np.logaddexp(0.0, math.log(value) + math.log(self.time_scale))
            ) / float(
                np.logaddexp(0.0, math.log(self.upper) + math.log(self.time_scale))
            )
        return (value - self.lower) / (self.upper - self.lower)

    def value(self, coordinate):
        if coordinate <= 0:
            return self.lower
        if coordinate >= 1:
            return self.upper
        if self.parameter == "sigma2" or (self.parameter == "alpha" and self.lower > 0):
            return math.exp(
                (1.0 - coordinate) * math.log(self.lower)
                + coordinate * math.log(self.upper)
            )
        if self.parameter == "alpha":
            transformed = coordinate * float(
                np.logaddexp(0.0, math.log(self.upper) + math.log(self.time_scale))
            )
            log_value = (
                transformed
                + math.log(-math.expm1(-transformed))
                - math.log(self.time_scale)
            )
            return math.exp(log_value)
        return (1.0 - coordinate) * self.lower + coordinate * self.upper


def _parameter_row(row, definitions, context):
    regime, name, group = (row[key] for key in ("regime", "parameter", "group"))
    if regime not in definitions or name not in _FIELDS or not group:
        raise ValueError(
            f"Invalid fit regime, parameter or empty group in {context}; "
            "only sigma2, alpha and theta can be estimated."
        )
    diffusion = definitions[regime].diffusion
    if diffusion is None or (name != "sigma2" and not isinstance(diffusion, OUBranch)):
        raise ValueError(
            f"Parameter {name} is not used by regime {regime} in {context}."
        )
    try:
        lower, upper = float(row["lower"]), float(row["upper"])
    except ValueError as exc:
        raise ValueError(f"Non-numeric parameter bounds in {context}.") from exc
    if (
        not all(math.isfinite(x) for x in (lower, upper, upper - lower))
        or lower >= upper
    ):
        raise ValueError(
            f"Fit bounds must be finite, representably increasing in {context}."
        )
    if (name == "sigma2" and lower <= 0) or (name == "alpha" and lower < 0):
        raise ValueError(
            "Estimated sigma2 needs a positive lower bound; alpha needs a nonnegative lower bound."
        )
    initial = getattr(diffusion, _FIELDS[name])
    if not lower <= initial <= upper:
        raise ValueError(
            f"Initial {name} for regime {regime} is outside its fit bounds."
        )
    return BranchFitParameter(group, name, (regime,), initial, lower, upper)


def load_branch_fit_spec(path, assignment, tree):
    """One row per estimated regime/parameter; identical group names tie values."""
    regimes = assignment.regime_by_branch_id
    if regimes is None:
        raise ValueError("--branch-fit requires --branch-regimes and --regime-models.")
    definitions = {
        regime: assignment.models_by_branch_id[identifier]
        for identifier, regime in regimes.items()
    }
    columns = {"regime", "parameter", "group", "lower", "upper"}
    rows = _read_rows(path, required=columns, allowed=columns)
    time_scale = max(
        (float(node.dist) for node in tree.traverse() if not node.is_root), default=0.0
    )
    groups: dict[str, BranchFitParameter] = {}
    seen = set()
    for line, row in rows:
        parameter = _parameter_row(row, definitions, f"{path}, line {line}")
        if parameter.parameter == "alpha":
            if not math.isfinite(time_scale) or time_scale <= 0:
                raise ValueError(
                    "Alpha estimation requires a positive finite branch length."
                )
            parameter = replace(parameter, time_scale=time_scale)
            if (
                parameter.lower == 0
                and np.logaddexp(0.0, math.log(parameter.upper) + math.log(time_scale))
                == 0
            ):
                raise ValueError(
                    "Alpha bounds are not numerically resolvable on this tree."
                )
        if (
            parameter.lower > 0
            and parameter.parameter != "theta"
            and math.log(parameter.upper) == math.log(parameter.lower)
        ):
            raise ValueError(
                "Fit bounds are indistinguishable on the logarithmic parameter scale."
            )
        key = (row["regime"], parameter.parameter)
        if key in seen:
            raise ValueError(f"Duplicate estimated regime/parameter: {key}.")
        seen.add(key)
        old = groups.get(parameter.group)
        if old is not None:
            if replace(old, regimes=parameter.regimes) != parameter:
                raise ValueError(
                    f"Shared group {parameter.group} must have the same parameter, "
                    "initial value and bounds in every regime."
                )
            parameter = replace(
                parameter, regimes=tuple(sorted((*old.regimes, *parameter.regimes)))
            )
        groups[parameter.group] = parameter
    if not 1 <= len(groups) <= 20:
        raise ValueError(
            "--branch-fit must specify between 1 and 20 free parameter groups."
        )
    return tuple(groups[key] for key in sorted(groups))


def fitted_assignment(assignment, parameters, coordinates):
    """Replace only specified diffusion parameters, preserving all end jumps."""
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (len(parameters),) or np.any(~np.isfinite(coordinates)):
        raise ValueError("Invalid branch-fit parameter coordinates.")
    if np.any(coordinates < 0) or np.any(coordinates > 1):
        raise ValueError("Branch-fit coordinates are outside their bounds.")
    changes: dict[str, dict[str, float]] = {}
    for parameter, coordinate in zip(parameters, coordinates, strict=True):
        for regime in parameter.regimes:
            changes.setdefault(regime, {})[_FIELDS[parameter.parameter]] = (
                parameter.value(float(coordinate))
            )
    regimes = assignment.regime_by_branch_id
    if regimes is None:
        raise ValueError("Parameter estimation requires named branch regimes.")
    models = {}
    for identifier, model in assignment.models_by_branch_id.items():
        updates = changes.get(regimes[identifier])
        models[identifier] = (
            replace(model, diffusion=replace(model.diffusion, **updates))
            if updates
            else model
        )
    return BranchGaussianAssignment(models, regimes)
