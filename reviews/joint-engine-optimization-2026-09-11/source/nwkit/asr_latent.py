"""CLI contracts and outputs for fixed-parameter latent-history ASR."""

import json
import sys

import pandas as pd

from nwkit.latent_gaussian import positive_integer
from nwkit.util import assign_branch_ids, is_missing_table_value

LATENT_MODELS = frozenset({"JUMP-BM", "MM-BM", "MM-OU"})
LATENT_OPTIONS = (
    "jump_rate",
    "jump_sd",
    "latent_regime_config",
    "regime_column",
    "history_samples",
    "latent_history_out",
)


def regime_input_columns(args, model):
    column = getattr(args, "regime_column", None)
    return (column,) if model in {"MM-BM", "MM-OU"} and column else ()


def validate_latent_options(args, model):
    def supplied(name):
        return getattr(args, name, None) not in (None, "")

    if model not in LATENT_MODELS:
        if any(supplied(name) for name in LATENT_OPTIONS):
            raise ValueError("Latent-history options require JUMP-BM, MM-BM, or MM-OU.")
        return
    unsupported = (
        "figure_out",
        "tree_ensemble",
        "model_comparison_out",
        "replicate_observations",
        "posterior_predictive_out",
        "bootstrap_out",
        "cross_validation_out",
        "bootstrap_intervals_out",
    )
    if any(supplied(name) for name in unsupported):
        raise ValueError(
            "Latent-history models currently support summary, tree, model, history, "
            "and posterior-sample outputs only; Gaussian diagnostics, figures, "
            "ensembles, replicates, and IC comparison are unsupported."
        )
    required: tuple[str, ...]
    invalid: tuple[str, ...]
    if model == "JUMP-BM":
        required = ("sigma2", "jump_rate", "jump_sd")
        invalid = ("latent_regime_config", "regime_column")
    else:
        required = ("latent_regime_config", "regime_column")
        invalid = ("sigma2", "jump_rate", "jump_sd")
    if any(not supplied(name) for name in required) or any(
        supplied(name) for name in invalid
    ):
        raise ValueError(
            f"{model} requires only these fixed-process inputs: {', '.join(required)}."
        )
    raw_count = getattr(args, "history_samples", None)
    count = positive_integer(
        1000 if raw_count is None else raw_count, "--history-samples"
    )
    if count < 2:
        raise ValueError("--history-samples must be at least two.")


def _fit(tree, frame, observed, errors, args, model):
    options = dict(
        standard_errors=errors,
        history_samples=1000 if args.history_samples is None else args.history_samples,
        seed=getattr(args, "seed", None),
    )
    if model == "JUMP-BM":
        from nwkit.jump_asr import fit_jump_bm

        return fit_jump_bm(
            tree,
            observed,
            sigma2=args.sigma2,
            jump_rate=args.jump_rate,
            jump_sd=args.jump_sd,
            **options,
        )
    from nwkit.markov_gaussian_asr import fit_markov_gaussian

    if args.regime_column not in frame.columns:
        raise ValueError(f"Missing regime column: {args.regime_column}.")
    with open(args.latent_regime_config, encoding="utf-8") as stream:
        config = json.load(stream)
    regimes = {
        str(name): None
        if is_missing_table_value(value, getattr(args, "missing_values", None))
        else str(value)
        for name, value in zip(
            frame["leaf_name"], frame[args.regime_column], strict=True
        )
        if str(name) in observed
    }
    return fit_markov_gaussian(tree, observed, regimes, config, model=model, **options)


def latent_model_table(fit, args, level):
    return pd.DataFrame(
        [
            {
                "trait_type": "continuous",
                "trait": args.state_column,
                "model": fit.model,
                "root_prior": "flat",
                "regime_root_prior": "equal" if fit.model.startswith("MM-") else "",
                "estimation_method": "fixed_parameters_importance_integration",
                "log_likelihood": fit.log_likelihood,
                "continuous_conditional_log_likelihood": fit.continuous_log_likelihood,
                "discrete_log_likelihood": fit.discrete_log_likelihood,
                "likelihood_kind": "monte_carlo_flat_root_integrated",
                "rankable": False,
                "fit_status": fit.fit_status,
                "history_samples": len(fit.weights),
                "importance_ess": fit.effective_sample_size,
                "likelihood_relative_mcse": fit.relative_mcse,
                "max_history_weight": float(fit.weights.max()),
                "num_observed": fit.num_observed,
                "ci_level": level,
                "interval_kind": "latent_history_mixture_conditional_on_parameters",
                "history_uncertainty_included": True,
                "parameter_uncertainty_included": False,
                "tree_uncertainty_included": False,
                "fixed_process_parameters": json.dumps(fit.parameters, sort_keys=True),
                "seed": getattr(args, "seed", None),
            }
        ]
    )


def latent_history_table(tree, fit):
    """Importance-weighted histories; segment durations remain in tree units."""
    ids = assign_branch_ids(tree)
    rows = []
    for index, (weight, history) in enumerate(
        zip(fit.weights, fit.latent, strict=True)
    ):
        for node in fit.nodes:
            row = {"history": index, "weight": weight, "branch_id": ids[node]}
            if fit.model == "JUMP-BM":
                if node.is_root:
                    continue
                row["jump_count"] = history[node]
            else:
                states = fit.parameters["states"]
                row["regime"] = states[history["node_states"][node]]
                row["segments"] = json.dumps(
                    [
                        [states[state], duration]
                        for state, duration in history["segments"].get(node, ())
                    ]
                )
            rows.append(row)
    return pd.DataFrame(rows)


def run_latent_asr(tree, frame, observed, errors, args, settings, targets):
    from nwkit.asr import _should_output_node, _write_table
    from nwkit.asr_continuous_diagnostics import posterior_samples_table
    from nwkit.continuous_asr_io import continuous_output_table, write_continuous_tree

    fit = _fit(tree, frame, observed, errors, args, settings.model)
    selected = [
        node for node in tree.traverse() if _should_output_node(node, observed, targets)
    ]
    table = continuous_output_table(
        tree,
        selected,
        observed,
        errors,
        fit.marginals,
        trait=args.state_column,
        ci_level=settings.ci_level,
    )
    _write_table(table, args.outfile)
    if getattr(args, "model_out", None):
        _write_table(latent_model_table(fit, args, settings.ci_level), args.model_out)
    if getattr(args, "latent_history_out", None):
        _write_table(latent_history_table(tree, fit), args.latent_history_out)
    if getattr(args, "posterior_samples_out", None):
        count = 1000 if args.posterior_samples is None else args.posterior_samples
        # Separate stream from history proposals.
        import numpy as np

        sequence = np.random.SeedSequence(getattr(args, "seed", None)).spawn(2)[1]
        samples = fit.sample(
            observed, errors, count, int(sequence.generate_state(1)[0])
        )
        _write_table(
            posterior_samples_table(tree, samples, args.state_column),
            args.posterior_samples_out,
        )
    write_continuous_tree(
        tree, observed, errors, fit.marginals, args, settings.ci_level
    )
    sys.stderr.write(
        f"Latent-history ASR: {fit.fit_status}; ESS={fit.effective_sample_size:.1f}/"
        f"{len(fit.weights)}, relative likelihood MCSE={fit.relative_mcse:.3g}. "
        "Repeat with independent seeds and more histories to assess convergence; "
        "intervals exclude parameter/tree uncertainty.\n"
    )
