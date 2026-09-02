"""CLI-facing simulation diagnostics for fitted scalar continuous ASR models."""

import math

import numpy as np
import pandas as pd

from nwkit.asr_diagnostics import gaussian_posterior_predictive, parametric_bootstrap
from nwkit.compiled_tree import CompiledTree
from nwkit.continuous_asr_process import fitted_scalar_process
from nwkit.gaussian_inference import (
    condition_gaussian_tree,
    sample_gaussian_posterior,
    simulate_gaussian_process,
)
from nwkit.util import assign_branch_ids, get_node_class


def _positive_count(value, option, default):
    value = default if value is None else value
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{option} must be a positive integer.")
    try:
        number = float(value)
        result = int(number)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{option} must be a positive integer.") from exc
    if not math.isfinite(number) or result != number or result <= 0:
        raise ValueError(f"{option} must be a positive integer.")
    return result


def posterior_samples_table(tree, samples, trait):
    """Convert joint all-node draws to a stable long-form table."""

    compiled = CompiledTree.from_tree(tree)
    if tuple(samples.nodes) != compiled.nodes:
        raise ValueError("Posterior samples do not follow the input tree order.")
    ids = assign_branch_ids(tree)
    num_samples, num_nodes = samples.values.shape
    node_indices = np.tile(np.arange(num_nodes), num_samples)
    sample_indices = np.repeat(np.arange(num_samples), num_nodes)
    nodes = compiled.nodes
    return pd.DataFrame(
        {
            "sample": sample_indices,
            "branch_id": [ids[nodes[index]] for index in node_indices],
            "parent": [
                -1 if nodes[index].is_root else ids[nodes[index].up]
                for index in node_indices
            ],
            "node_class": [get_node_class(nodes[index]) for index in node_indices],
            "name": [
                "" if nodes[index].name in (None, "") else str(nodes[index].name)
                for index in node_indices
            ],
            "trait": trait,
            "value": samples.values.reshape(-1),
        }
    )


def _simulated_tip_data(process, observed, errors, seed):
    seed_sequence = np.random.SeedSequence(seed)
    process_seed, error_seed = seed_sequence.spawn(2)
    conditioned = condition_gaussian_tree(process, observed, standard_errors=errors)
    root_values = (
        conditioned.marginals[process.tree].mean
        if process.root.mode == "flat"
        else None
    )
    simulated = simulate_gaussian_process(
        process,
        num_samples=1,
        seed=int(process_seed.generate_state(1)[0]),
        root_values=root_values,
        compiled_tree=conditioned.compiled_tree,
    )
    node_index = {node: index for index, node in enumerate(simulated.nodes)}
    error_rng = np.random.default_rng(error_seed)
    result = dict.fromkeys(observed)
    for leaf in process.tree.leaves():
        name = str(leaf.name)
        if observed.get(name) is None:
            continue
        value = float(simulated.values[0, node_index[leaf]])
        error = 0.0 if errors is None else float(errors[name])
        if error > 0.0:
            value += float(error_rng.normal(0.0, error))
        result[name] = value
    return result


def _refit(tree, simulated, errors, args, settings, original_fit, regime_assignment):
    model = settings.model
    if model == "BM":
        from nwkit.continuous_asr import compute_bm_marginals

        return compute_bm_marginals(
            tree,
            simulated,
            sigma2=None if original_fit.sigma2_estimated else original_fit.sigma2,
            standard_errors=errors,
            _tree_validated=True,
        )[1]
    if model == "BMS":
        from nwkit.regime_continuous_asr import compute_bms_marginals

        return compute_bms_marginals(
            tree,
            simulated,
            regime_assignment,
            sigma2_by_regime=None
            if original_fit.sigma2_estimated
            else original_fit.sigma2_by_regime,
            standard_errors=errors,
            _tree_validated=True,
        )[1]
    if model == "BMS-DRIFT":
        from nwkit.asr_regimes import read_regime_parameters
        from nwkit.regime_drift_asr import compute_regime_drift_marginals

        fixed = read_regime_parameters(
            getattr(args, "regime_parameters", None),
            regime_assignment.regimes,
            ("sigma2", "drift"),
        )
        return compute_regime_drift_marginals(
            tree,
            simulated,
            regime_assignment,
            sigma2=getattr(args, "sigma2", None),
            drift=getattr(args, "drift", None),
            regime_parameters=None if fixed is None else fixed[0],
            regime_parameters_source=None if fixed is None else fixed[1],
            standard_errors=errors,
            _tree_validated=True,
        )[1]
    if model in {"LAMBDA", "KAPPA", "DELTA", "EB", "ACDC"}:
        from nwkit.transformed_continuous_asr import compute_transformed_bm_marginals

        return compute_transformed_bm_marginals(
            tree,
            simulated,
            model=model.lower(),
            sigma2=None if original_fit.sigma2_estimated else original_fit.sigma2,
            evolution_parameter=(
                None
                if original_fit.evolution_parameter_estimated
                else original_fit.evolution_parameter
            ),
            evolution_parameter_bounds=original_fit.evolution_parameter_bounds,
            standard_errors=errors,
            _tree_validated=True,
        )[1]
    if model == "BM-DRIFT":
        from nwkit.nonstationary_continuous_asr import compute_bm_drift_marginals

        return compute_bm_drift_marginals(
            tree,
            simulated,
            sigma2=None if original_fit.sigma2_estimated else original_fit.sigma2,
            drift=None if original_fit.drift_estimated else original_fit.drift,
            standard_errors=errors,
            _tree_validated=True,
        )[1]
    if model == "OU":
        if settings.root_prior == "stationary":
            from nwkit.ou_asr import compute_ou_marginals

            return compute_ou_marginals(
                tree,
                simulated,
                alpha=None if original_fit.alpha_estimated else original_fit.alpha,
                sigma2=None if original_fit.sigma2_estimated else original_fit.sigma2,
                theta=None if original_fit.theta_estimated else original_fit.theta,
                alpha_bounds=original_fit.alpha_bounds,
                standard_errors=errors,
                _tree_validated=True,
            )[1]
        from nwkit.general_ou_asr import compute_general_ou_marginals

        return compute_general_ou_marginals(
            tree,
            simulated,
            root_prior=settings.root_prior,
            root_mean=original_fit.root_mean,
            root_variance=(
                original_fit.root_variance
                if settings.root_prior == "gaussian"
                else None
            ),
            alpha=None if original_fit.alpha_estimated else original_fit.alpha,
            sigma2=None if original_fit.sigma2_estimated else original_fit.sigma2,
            theta=None if original_fit.theta_estimated else original_fit.theta,
            alpha_bounds=original_fit.alpha_bounds,
            standard_errors=errors,
        )[1]
    if model in {"OUM", "OUMA", "OUMV", "OUMVA"}:
        from nwkit.asr_regimes import read_regime_parameters
        from nwkit.regime_gaussian_asr import (
            compute_regime_ou_marginals,
            regime_parameter_columns,
        )

        fixed = read_regime_parameters(
            getattr(args, "regime_parameters", None),
            regime_assignment.regimes,
            regime_parameter_columns(model),
        )
        return compute_regime_ou_marginals(
            tree,
            simulated,
            regime_assignment,
            model=model,
            alpha=getattr(args, "alpha", None),
            sigma2=getattr(args, "sigma2", None),
            theta=getattr(args, "theta", None),
            regime_parameters=None if fixed is None else fixed[0],
            regime_parameters_source=None if fixed is None else fixed[1],
            alpha_bounds=original_fit.alpha_bounds,
            standard_errors=errors,
            _tree_validated=True,
        )[1]
    raise ValueError(f"Parametric bootstrap is unavailable for --model {model}.")


def _fit_row(fit, model):
    row = {"model": model, "fit_model_status": getattr(fit, "fit_status", "")}
    for name in (
        "sigma2",
        "alpha",
        "theta",
        "drift",
        "evolution_parameter",
        "log_likelihood",
        "restricted_log_likelihood",
    ):
        value = getattr(fit, name, None)
        if value is not None:
            row[name] = value
    for name in (
        "sigma2_by_regime",
        "alpha_by_regime",
        "theta_by_regime",
        "drift_by_regime",
    ):
        values = getattr(fit, name, None)
        if values is not None:
            for regime, value in values.items():
                safe = str(regime).encode("utf-8").hex()
                row[f"{name.removesuffix('_by_regime')}_regime_{safe}"] = value
    return row


def write_continuous_diagnostics(
    tree,
    observed,
    errors,
    args,
    settings,
    fit,
    regime_assignment,
):
    """Write each requested scalar Gaussian diagnostic."""

    sample_out = getattr(args, "posterior_samples_out", None)
    predictive_out = getattr(args, "posterior_predictive_out", None)
    bootstrap_out = getattr(args, "bootstrap_out", None)
    if all(
        value in (None, "") for value in (sample_out, predictive_out, bootstrap_out)
    ):
        return
    process = fitted_scalar_process(
        tree,
        settings.model,
        fit,
        root_prior=settings.root_prior,
        regime_assignment=regime_assignment,
    )
    seed_sequence = np.random.SeedSequence(getattr(args, "seed", None))
    sample_seed, predictive_seed, bootstrap_seed = seed_sequence.spawn(3)
    if sample_out not in (None, ""):
        draws = _positive_count(
            getattr(args, "posterior_samples", None), "--posterior-samples", 1000
        )
        samples = sample_gaussian_posterior(
            process,
            observed,
            num_samples=draws,
            standard_errors=errors,
            seed=int(sample_seed.generate_state(1)[0]),
        )
        posterior_samples_table(tree, samples, args.state_column).to_csv(
            sample_out, sep="\t", index=False
        )
    if predictive_out not in (None, ""):
        simulations = _positive_count(
            getattr(args, "posterior_predictive_simulations", None),
            "--posterior-predictive-simulations",
            1000,
        )
        gaussian_posterior_predictive(
            process,
            observed,
            standard_errors=errors,
            num_simulations=simulations,
            seed=int(predictive_seed.generate_state(1)[0]),
        ).to_csv(predictive_out, sep="\t", index=False)
    if bootstrap_out not in (None, ""):
        simulations = _positive_count(
            getattr(args, "bootstrap_simulations", None),
            "--bootstrap-simulations",
            100,
        )
        rows, failures = parametric_bootstrap(
            lambda replicate_seed: _simulated_tip_data(
                process, observed, errors, replicate_seed
            ),
            lambda values: _refit(
                tree,
                values,
                errors,
                args,
                settings,
                fit,
                regime_assignment,
            ),
            lambda fitted: _fit_row(fitted, settings.model),
            num_simulations=simulations,
            seed=int(bootstrap_seed.generate_state(1)[0]),
        )
        if len(failures):
            rows = pd.concat((rows, failures), ignore_index=True, sort=False)
        rows.sort_values("replicate").to_csv(bootstrap_out, sep="\t", index=False)
