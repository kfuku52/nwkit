"""Rebuild fitted scalar ASR models as shared Gaussian tree processes."""

from nwkit.evolution import build_evolutionary_process
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
    brownian_transition,
    ou_transition,
)
from nwkit.regime_drift_asr import build_regime_drift_process
from nwkit.regime_gaussian_asr import build_regime_ou_process


def fitted_scalar_process(tree, model, fit, *, root_prior=None, regime_assignment=None):
    """Return the exact affine-Gaussian process represented by an ASR fit."""

    model = str(model).upper()
    if model == "BM":
        return build_evolutionary_process(
            tree,
            model="brownian",
            root_mode="flat",
            variance_scale=fit.sigma2,
            allow_zero=True,
        )
    if model in {"LAMBDA", "KAPPA", "DELTA", "EB", "ACDC"}:
        return build_evolutionary_process(
            tree,
            model=model.lower(),
            parameter=fit.evolution_parameter,
            root_mode="flat",
            variance_scale=fit.sigma2,
            allow_zero=True,
        )
    if model == "BM-DRIFT":
        transitions = {
            node: GaussianTransition(
                1.0,
                fit.drift * float(node.dist),
                fit.sigma2 * float(node.dist),
            )
            for node in tree.traverse()
            if not node.is_root
        }
        return GaussianTreeProcess(
            tree,
            transitions,
            GaussianRootPrior("flat", 0.0, None),
            "brownian-drift",
        )
    if model == "BMS":
        if regime_assignment is None:
            raise ValueError("A fitted BMS process requires its regime assignment.")
        transitions = {
            node: brownian_transition(
                fit.sigma2_by_regime[regime_assignment.by_node[node]] * float(node.dist)
            )
            for node in tree.traverse()
            if not node.is_root
        }
        return GaussianTreeProcess(
            tree,
            transitions,
            GaussianRootPrior("flat", 0.0, None),
            "bms",
        )
    if model == "BMS-DRIFT":
        if regime_assignment is None:
            raise ValueError(
                "A fitted BMS-DRIFT process requires its regime assignment."
            )
        return build_regime_drift_process(
            tree,
            regime_assignment,
            sigma2_by_regime=fit.sigma2_by_regime,
            drift_by_regime=fit.drift_by_regime,
        )
    if model == "OU":
        mode = root_prior or getattr(fit, "root_prior", "stationary")
        transitions = {
            node: ou_transition(float(node.dist), fit.alpha, fit.sigma2, fit.theta)[0]
            for node in tree.traverse()
            if not node.is_root
        }
        if mode == "stationary":
            root = GaussianRootPrior(
                "stationary", fit.theta, fit.sigma2 / (2.0 * fit.alpha)
            )
        elif mode == "fixed":
            root = GaussianRootPrior("fixed", fit.root_mean, 0.0)
        elif mode == "gaussian":
            root = GaussianRootPrior("gaussian", fit.root_mean, fit.root_variance)
        else:
            raise ValueError(f"Unsupported fitted OU root prior: {mode}.")
        return GaussianTreeProcess(tree, transitions, root, "ou", fit.alpha)
    if model in {"OUM", "OUMA", "OUMV", "OUMVA"}:
        if regime_assignment is None:
            raise ValueError(
                "A fitted OU regime process requires its regime assignment."
            )
        return build_regime_ou_process(
            tree,
            regime_assignment,
            alpha_by_regime=fit.alpha_by_regime,
            sigma2_by_regime=fit.sigma2_by_regime,
            theta_by_regime=fit.theta_by_regime,
        )
    raise ValueError(
        f"Gaussian process diagnostics are not implemented for --model {model}."
    )
