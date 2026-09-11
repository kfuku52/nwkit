"""Fixed-parameter scalar BM, OU and Gaussian jumps assigned by branch ID."""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral
from typing import Any

import numpy as np

from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
)
from nwkit.util import assign_branch_ids


def _finite(value: float, name: str, *, nonnegative: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be numeric, not boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite.") from exc
    if not math.isfinite(result) or (nonnegative and result < 0.0):
        raise ValueError(
            f"{name} must be finite" + (" and non-negative." if nonnegative else ".")
        )
    return result


@dataclass(frozen=True, slots=True)
class BrownianBranch:
    """Brownian diffusion with variance rate in trait² per tree-time unit."""

    variance_rate: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "variance_rate",
            _finite(self.variance_rate, "variance_rate", nonnegative=True),
        )


@dataclass(frozen=True, slots=True)
class OUBranch:
    """OU diffusion; alpha=0 is exactly Brownian with the same variance rate."""

    alpha: float
    variance_rate: float = 1.0
    optimum: float = 0.0

    def __post_init__(self) -> None:
        for name in ("alpha", "variance_rate", "optimum"):
            object.__setattr__(
                self,
                name,
                _finite(getattr(self, name), name, nonnegative=name != "optimum"),
            )


@dataclass(frozen=True, slots=True)
class GaussianJump:
    """One prescribed independent additive Gaussian event, at the branch end.

    Variance is per event, not per time. This is not a Poisson jump mixture.
    """

    mean: float = 0.0
    variance: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "mean", _finite(self.mean, "jump mean"))
        object.__setattr__(
            self, "variance", _finite(self.variance, "jump variance", nonnegative=True)
        )


def _product_ratio(*factors: float, denominator: float = 1.0) -> float:
    """Multiply finite nonnegative factors without intermediate range loss."""
    if any(value == 0.0 for value in factors):
        return 0.0
    mantissa, exponent = 1.0, 0
    for value in factors:
        part, power = math.frexp(value)
        mantissa *= part
        exponent += power
    divisor, power = math.frexp(denominator)
    try:
        return math.ldexp(mantissa / divisor, exponent - power)
    except OverflowError as exc:
        raise ValueError("Gaussian transition overflows floating-point range.") from exc


def _variance_product(*factors: float, denominator: float = 1.0) -> float:
    variance = _product_ratio(*factors, denominator=denominator)
    if variance == 0.0 and all(value > 0.0 for value in factors):
        raise ValueError("Positive diffusion variance underflows floating-point range.")
    return variance


def _diffusion_transition(diffusion, length):
    if diffusion is None:
        return GaussianTransition(1.0, 0.0, 0.0)
    rate = diffusion.variance_rate
    if isinstance(diffusion, BrownianBranch) or diffusion.alpha == 0.0:
        return GaussianTransition(1.0, 0.0, _variance_product(rate, length))
    exponent = diffusion.alpha * length
    attenuation = -math.expm1(-exponent)
    innovation_fraction = -math.expm1(-2.0 * exponent)
    # Avoid forming stationary variance or rounding a tiny rate*time product
    # before multiplying by another parameter that restores its scale.
    if exponent <= 1.0:
        factor = innovation_fraction / (2.0 * exponent) if exponent else 1.0
        variance = _variance_product(rate, length, factor)
    else:
        variance = _variance_product(
            rate, 0.5 * innovation_fraction, denominator=diffusion.alpha
        )
    if exponent < 1e-8:
        factor = attenuation / exponent if exponent else 1.0
        intercept = math.copysign(
            _product_ratio(diffusion.alpha, length, abs(diffusion.optimum), factor),
            diffusion.optimum,
        )
    else:
        intercept = attenuation * diffusion.optimum
    return GaussianTransition(math.exp(-exponent), intercept, variance)


@dataclass(frozen=True, slots=True)
class BranchGaussianModel:
    """Diffusion along a branch followed by an optional independent end jump."""

    diffusion: BrownianBranch | OUBranch | None = None
    jump: GaussianJump | None = None

    def __post_init__(self) -> None:
        if self.diffusion is not None and not isinstance(
            self.diffusion, (BrownianBranch, OUBranch)
        ):
            raise ValueError("diffusion must be BrownianBranch, OUBranch, or None.")
        if self.jump is not None and not isinstance(self.jump, GaussianJump):
            raise ValueError("jump must be GaussianJump or None.")
        if self.diffusion is None and self.jump is None:
            raise ValueError("A branch model requires a diffusion or a jump.")

    def transition(self, length: float) -> GaussianTransition:
        """Compile an affine transition; zero length still permits an end jump."""
        length = _finite(length, "Branch length", nonnegative=True)
        transition = _diffusion_transition(self.diffusion, length)
        if self.jump is None:
            return transition
        return GaussianTransition(
            transition.slope,
            transition.intercept + self.jump.mean,
            transition.variance + self.jump.variance,
        )


def build_branch_gaussian_process(
    tree: Any,
    models_by_branch_id: Mapping[int, BranchGaussianModel],
    *,
    root: GaussianRootPrior,
) -> GaussianTreeProcess:
    """Assign every non-root incoming branch exactly once, with an explicit root.

    IDs follow ``nwkit.util.assign_branch_ids`` (level order). Root ID 0 must
    not occur in the map. Recompute assignments after changing tree topology
    or child ordering. Returned processes use the existing likelihood,
    conditioning, covariance and simulation APIs without conversion.
    """
    if not isinstance(root, GaussianRootPrior):
        raise ValueError("root must be an explicit GaussianRootPrior.")
    if not tree.is_root:
        raise ValueError(
            "tree must be a root node; detach a subtree before assigning IDs."
        )
    if not isinstance(models_by_branch_id, Mapping):
        raise ValueError("Branch models must be a mapping from integer IDs to models.")
    identifiers = assign_branch_ids(tree)
    nodes = {
        identifier: node for node, identifier in identifiers.items() if not node.is_root
    }
    supplied = dict(models_by_branch_id)
    if any(isinstance(key, bool) or not isinstance(key, Integral) for key in supplied):
        raise ValueError("Branch IDs must be integers.")
    if set(supplied) != set(nodes):
        raise ValueError(
            "Branch models must cover every non-root branch ID exactly; "
            f"missing={sorted(set(nodes) - set(supplied))}, "
            f"extra={sorted(set(supplied) - set(nodes))}."
        )
    transitions = {}
    for identifier, node in nodes.items():
        model = supplied[identifier]
        if not isinstance(model, BranchGaussianModel):
            raise ValueError(f"Branch {identifier} requires a BranchGaussianModel.")
        transitions[node] = model.transition(node.dist)
    return GaussianTreeProcess(tree, transitions, root, model="branch-heterogeneous")
