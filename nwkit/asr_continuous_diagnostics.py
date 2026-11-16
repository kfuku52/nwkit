"""CLI-facing simulation diagnostics for fitted scalar continuous ASR models."""

import math

import numpy as np
import pandas as pd

from nwkit.asr_diagnostics import gaussian_posterior_predictive, parametric_bootstrap
from nwkit.compiled_tree import CompiledTree
from nwkit.continuous_asr_process import fitted_scalar_process
from nwkit.gaussian_inference import (
    condition_gaussian_tree,
    prepare_gaussian_single_draw,
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
    nodes = compiled.nodes
    return pd.DataFrame(
        {
            "sample": np.repeat(np.arange(num_samples), num_nodes),
            "branch_id": np.tile([ids[node] for node in nodes], num_samples),
            "parent": np.tile(
                [-1 if node.is_root else ids[node.up] for node in nodes],
                num_samples,
            ),
            "node_class": np.tile(
                np.asarray([get_node_class(node) for node in nodes], dtype=object),
                num_samples,
            ),
            "name": np.tile(
                np.asarray(
                    [
                        "" if node.name in (None, "") else str(node.name)
                        for node in nodes
                    ],
                    dtype=object,
                ),
                num_samples,
            ),
            "trait": trait,
            "value": samples.values.reshape(-1),
        }
    )


def _tip_bootstrap_simulator(process, observed, errors):
    compiled = CompiledTree.from_tree(process.tree)
    root_values = None
    if process.root.mode == "flat":
        conditioned = condition_gaussian_tree(
            process, observed, standard_errors=errors, compiled_tree=compiled
        )
        root_values = conditioned.marginals[process.tree].mean
    records = [
        (name, index, 0.0 if errors is None else float(errors[name]))
        for name, index in compiled.leaf_index_by_name.items()
        if observed.get(name) is not None
    ]
    draw = prepare_gaussian_single_draw(process, compiled, root_values)
    names = [item[0] for item in records]
    indices = np.asarray([item[1] for item in records], dtype=int)
    deviations = np.asarray([item[2] for item in records], dtype=float)
    noisy = deviations > 0.0
    positive_deviations = deviations[noisy]

    def simulate(seed):
        process_seed, error_seed = np.random.SeedSequence(seed).spawn(2)
        rng = np.random.default_rng(int(process_seed.generate_state(1)[0]))
        values = draw(rng)[indices]
        error_rng = np.random.default_rng(error_seed)
        values[noisy] += error_rng.normal(0.0, positive_deviations)
        result = dict.fromkeys(observed)
        result.update(zip(names, values.tolist(), strict=True))
        return result

    return simulate


def _simulated_tip_data(process, observed, errors, seed, *, include_latent=False):
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
    return (result, simulated) if include_latent else result


def _refit(tree, simulated, errors, args, settings, original_fit, regime_assignment):
    if settings.model == "BRANCH-GAUSSIAN":
        from nwkit.branch_gaussian_asr import compute_branch_marginals

        return compute_branch_marginals(
            tree,
            simulated,
            errors,
            assignment=original_fit.branch_assignment,
            fit_spec=original_fit.fit_spec,
            root=original_fit.process.root,
            compute_posterior=False,
        )[1]
    model = settings.model
    if model == "BM":
        from nwkit.continuous_asr import compute_bm_marginals

        return compute_bm_marginals(
            tree,
            simulated,
            sigma2=None if original_fit.sigma2_estimated else original_fit.sigma2,
            standard_errors=errors,
            _tree_validated=True,
            _compute_posterior=False,
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
            _compute_posterior=False,
        )[1]
    if model == "BMS-DRIFT":
        from nwkit.regime_drift_asr import compute_regime_drift_marginals

        fixed = (
            None
            if original_fit.regime_parameters_source is None
            else {
                regime: {
                    "sigma2": original_fit.sigma2_by_regime[regime],
                    "drift": original_fit.drift_by_regime[regime],
                }
                for regime in regime_assignment.regimes
            }
        )
        return compute_regime_drift_marginals(
            tree,
            simulated,
            regime_assignment,
            sigma2=getattr(args, "sigma2", None),
            drift=getattr(args, "drift", None),
            regime_parameters=fixed,
            regime_parameters_source=original_fit.regime_parameters_source,
            standard_errors=errors,
            _tree_validated=True,
            _compute_posterior=False,
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
            _compute_posterior=False,
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
            _compute_posterior=False,
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
                _compute_posterior=False,
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
            _compute_posterior=False,
        )[1]
    if model in {"OUM", "OUMA", "OUMV", "OUMVA"}:
        from nwkit.regime_gaussian_asr import (
            compute_regime_ou_marginals,
            regime_parameter_columns,
        )

        fitted_maps = {
            "alpha": original_fit.alpha_by_regime,
            "sigma2": original_fit.sigma2_by_regime,
            "theta": original_fit.theta_by_regime,
        }
        fixed = (
            None
            if original_fit.regime_parameters_source is None
            else {
                regime: {
                    column: fitted_maps[column][regime]
                    for column in regime_parameter_columns(model)
                }
                for regime in regime_assignment.regimes
            }
        )
        return compute_regime_ou_marginals(
            tree,
            simulated,
            regime_assignment,
            model=model,
            alpha=getattr(args, "alpha", None),
            sigma2=getattr(args, "sigma2", None),
            theta=getattr(args, "theta", None),
            regime_parameters=fixed,
            regime_parameters_source=original_fit.regime_parameters_source,
            alpha_bounds=original_fit.alpha_bounds,
            standard_errors=errors,
            _tree_validated=True,
            _compute_posterior=False,
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

    if settings.model in {"MV-BM", "MV-OU", "MV-OU-DIAG", "MV-OU-FULL"}:
        from nwkit.asr_multivariate_diagnostics import write_multivariate_diagnostics

        write_multivariate_diagnostics(tree, observed, errors, args, settings, fit)
        return

    sample_out = getattr(args, "posterior_samples_out", None)
    predictive_out = getattr(args, "posterior_predictive_out", None)
    bootstrap_out = getattr(args, "bootstrap_out", None)
    cross_validation_out = getattr(args, "cross_validation_out", None)
    interval_out = getattr(args, "bootstrap_intervals_out", None)

    def refit_process(training):
        refitted = _refit(
            tree, training, errors, args, settings, fit, regime_assignment
        )
        return fitted_scalar_process(
            tree,
            settings.model,
            refitted,
            root_prior=settings.root_prior,
            regime_assignment=regime_assignment,
        )

    if cross_validation_out not in (None, ""):
        from nwkit.asr_cross_validation import gaussian_cross_validation

        gaussian_cross_validation(
            tree,
            observed,
            refit_process,
            errors=errors,
            mode=getattr(args, "cross_validation_unit", None) or "tip",
            level=getattr(args, "ci_level", None) or 0.95,
        ).to_csv(cross_validation_out, sep="\t", index=False)
    if all(
        value in (None, "")
        for value in (sample_out, predictive_out, bootstrap_out, interval_out)
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
    if interval_out not in (None, ""):
        from nwkit.asr_bootstrap_intervals import bootstrap_prediction_intervals

        bootstrap_prediction_intervals(
            process,
            observed,
            lambda seed: _simulated_tip_data(
                process, observed, errors, seed, include_latent=True
            ),
            refit_process,
            errors=errors,
            num_simulations=_positive_count(
                getattr(args, "bootstrap_interval_simulations", None),
                "--bootstrap-interval-simulations",
                100,
            ),
            seed=int(seed_sequence.spawn(1)[0].generate_state(1)[0]),
            level=getattr(args, "ci_level", None) or 0.95,
        ).to_csv(interval_out, sep="\t", index=False)
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
            _tip_bootstrap_simulator(process, observed, errors),
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
