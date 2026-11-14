"""Individual-keyed ASR input and outputs for jointly estimated Sigma and W."""

import sys
from copy import copy
from dataclasses import dataclass

import numpy as np
import pandas as pd

from nwkit.asr_input import AsrSettings, asr_trait_columns
from nwkit.continuous_asr import GaussianMarginal
from nwkit.continuous_asr_io import _summary
from nwkit.evolution import build_evolutionary_process
from nwkit.individual_covariance import conditional_vector, fit_individual_covariance
from nwkit.multivariate_asr import (
    MultivariateGaussianMarginal,
    _trait_id,
    multivariate_covariance_table,
    multivariate_output_table,
    write_multivariate_tree,
)
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.trait_input import numeric_trait_value
from nwkit.util import read_tip_table


@dataclass(frozen=True)
class IndividualData:
    names: tuple
    traits: tuple
    keys: tuple
    matrix: np.ndarray
    species: np.ndarray
    individuals: np.ndarray
    coordinates: np.ndarray
    values: np.ndarray


def validate_individual_options(args):
    within = getattr(args, "within_species_covariance", None)
    if within is None:
        if any(
            getattr(args, name, None) is not None
            for name in ("covariance_method", "individual_out")
        ):
            raise ValueError(
                "--covariance-method and --individual-out require --within-species-covariance."
            )
        return None
    if args.model != "MV-BM":
        raise ValueError("--within-species-covariance requires --model MV-BM.")
    if getattr(args, "outfile", None) in (None, ""):
        raise ValueError("--outfile must be a file path or '-' for standard output.")
    if getattr(args, "trait_type", "auto") == "discrete":
        raise ValueError("Individual covariance requires continuous traits.")
    settings = AsrSettings.from_args(args, "continuous")
    unsupported = (
        "standard_error_column",
        "measurement_covariance",
        "replicate_observations",
        "compare_models",
        "model_comparison_out",
        "tree_ensemble",
        "tree_ensemble_out",
        "figure_out",
        "posterior_samples_out",
        "posterior_samples",
        "posterior_predictive_out",
        "posterior_predictive_simulations",
        "bootstrap_out",
        "bootstrap_simulations",
        "bootstrap_intervals_out",
        "bootstrap_interval_simulations",
        "cross_validation_out",
        "cross_validation",
        "profile_ci_level",
        "latent_regime_config",
        "latent_history_out",
    )
    supplied = [
        "--" + name.replace("_", "-")
        for name in unsupported
        if getattr(args, name, None) not in (None, "")
    ]
    if supplied:
        raise ValueError(
            "Joint individual covariance does not support: " + ", ".join(supplied)
        )
    return settings


def read_individual_data(args, names, traits):
    columns = ("individual_id", "trait", "value")
    table, _, _ = read_tip_table(
        args.trait,
        tree_leaf_names=names,
        required_columns=columns,
        duplicate_leaf_names="allow",
        preserve_columns=columns,
        missing_values="",  # Apply trait missing markers only to values, never IDs.
        unmatched=getattr(args, "unmatched", "warn"),
    )
    table = table[table.leaf_name.isin(names)].copy()
    if table.individual_id.isna().any() or any(
        not str(x).strip() for x in table.individual_id
    ):
        raise ValueError(
            "individual_id must be nonempty; IDs are scoped within each species."
        )
    if table.duplicated(["leaf_name", "individual_id", "trait"]).any():
        raise ValueError(
            "Duplicate (leaf_name, individual_id, trait) rows are not independent individuals."
        )
    if not set(table.trait).issubset(traits):
        raise ValueError("Individual TSV contains traits absent from --state-column.")
    keys = tuple(sorted(set(zip(table.leaf_name, table.individual_id, strict=True))))
    key_index = {key: i for i, key in enumerate(keys)}
    matrix = np.full((len(keys), len(traits)), np.nan)
    trait_index = {trait: i for i, trait in enumerate(traits)}
    for row in table.itertuples(index=False):
        matrix[
            key_index[(row.leaf_name, row.individual_id)], trait_index[row.trait]
        ] = numeric_trait_value(
            row.value, row.trait, getattr(args, "missing_values", None)
        )
    individual, coordinate = np.where(np.isfinite(matrix))
    name_index = {name: i for i, name in enumerate(names)}
    species = np.array([name_index[keys[i][0]] for i in individual], dtype=int)
    return IndividualData(
        tuple(names),
        tuple(traits),
        keys,
        matrix,
        species,
        individual,
        coordinate,
        matrix[individual, coordinate],
    )


def _node_marginals(tree, data, fit, process):
    # A single dense tree covariance supports both tips and internal nodes.
    nodes = tuple(tree.traverse())
    if 8 * len(nodes) ** 2 > 256 * 1024**2:
        raise ValueError(
            "Joint ASR node covariance exceeds the 256 MiB reconstruction budget."
        )
    covariance = process.covariance(nodes)
    index = {node: i for i, node in enumerate(nodes)}
    leaves = {str(node.name): node for node in tree.leaves()}
    observation_nodes = np.array([index[leaves[data.names[s]]] for s in data.species])
    posterior = {}
    for node in nodes:
        row = index[node]
        mean, variance = conditional_vector(
            fit,
            data.coordinates,
            covariance[row, observation_nodes],
            covariance[row, row],
        )
        posterior[node] = MultivariateGaussianMarginal(mean, variance)
    return posterior


def _species_summaries(data, fit):
    observed, errors, counts = {}, {}, {}
    for name in data.names:
        rows = [i for i, key in enumerate(data.keys) if key[0] == name]
        count = np.sum(np.isfinite(data.matrix[rows]), axis=0)
        total = np.nansum(data.matrix[rows], axis=0)
        observed[name] = tuple(
            float(total[t] / count[t]) if count[t] else None
            for t in range(len(data.traits))
        )
        errors[name] = tuple(
            float(np.sqrt(fit.within[t, t]) / np.sqrt(count[t])) if count[t] else None
            for t in range(len(data.traits))
        )
        counts[name] = count
    return observed, errors, counts


def _individual_table(data, fit, covariance, ci_level):
    rows = []
    name_index = {name: i for i, name in enumerate(data.names)}
    for individual, (name, identifier) in enumerate(data.keys):
        s = name_index[name]
        mean, variance = conditional_vector(
            fit,
            data.coordinates,
            covariance[s, data.species],
            covariance[s, s],
            same_individual=data.individuals == individual,
        )
        for t, trait in enumerate(data.traits):
            value = data.matrix[individual, t]
            missing = not np.isfinite(value)
            # Observed individual coordinates are conditioned on exactly; W
            # describes biological dispersion, not additional known instrument SE.
            marginal = GaussianMarginal(
                float(mean[t]) if missing else float(value),
                float(variance[t, t]) if missing else 0.0,
            )
            rows.append(
                {
                    "leaf_name": name,
                    "individual_id": identifier,
                    "trait": trait,
                    "observed_value": "" if missing else value,
                    "is_imputed": missing,
                    "estimand": "individual_value",
                    **_summary(marginal, ci_level),
                }
            )
    return pd.DataFrame(rows)


def _model_table(fit, data, args, ci_level):
    row = {
        "model": "MV-BM",
        "observation_model": "independent_individuals_with_common_W",
        "trait": ",".join(data.traits),
        "root_prior": "flat",
        "estimation_method": fit.method,
        "within_species_covariance": args.within_species_covariance,
        "log_likelihood": fit.log_likelihood if fit.method == "ML" else "",
        "restricted_log_likelihood": fit.log_likelihood if fit.method == "REML" else "",
        "likelihood_kind": "profile_ml"
        if fit.method == "ML"
        else "flat_root_integrated",
        "num_observed_coordinates": len(data.values),
        "num_observed_individuals": len(np.unique(data.individuals)),
        "num_observed_species": len(np.unique(data.species)),
        "residual_df": len(data.values) - len(data.traits),
        "fit_status": fit.fit_status,
        "optimizer_success": fit.optimizer_success,
        "optimizer_starts": fit.optimizer_starts,
        "optimizer_converged_starts": fit.optimizer_converged_starts,
        "identifiability_singular_value_ratio": fit.identifiability_ratio,
        "sigma_eigenvalue_ratio_scaled": fit.sigma_eigenvalue_ratio,
        "within_eigenvalue_ratio_scaled": fit.within_eigenvalue_ratio,
        "sigma_interpretation": "diffusion_covariance_per_branch_length",
        "within_interpretation": "individual_covariance_around_latent_species_mean",
        "ci_level": ci_level,
        "interval_kind": "conditional_on_fitted_Sigma_W",
        "root_mean_uncertainty_included": True,
        "parameter_uncertainty_included": False,
        "tree_uncertainty_included": False,
    }
    for a, first in enumerate(data.traits):
        row[f"root_mean_{_trait_id(first)}"] = fit.root_mean[a]
        for b, second in enumerate(data.traits):
            row[f"sigma_{_trait_id(first)}_to_{_trait_id(second)}"] = fit.sigma[a, b]
            row[f"within_{_trait_id(first)}_to_{_trait_id(second)}"] = fit.within[a, b]
    return pd.DataFrame([row])


def run_individual_asr(tree, args, settings):
    from nwkit.asr import _parse_targets, _should_output_node, _write_table

    traits = asr_trait_columns(args.state_column, "MV-BM")
    data = read_individual_data(args, sorted(tree.leaf_names()), traits)
    if 8 * sum(1 for _ in tree.traverse()) ** 2 > 256 * 1024**2:
        raise ValueError(
            "Joint ASR node covariance exceeds the 256 MiB reconstruction budget."
        )
    process = build_evolutionary_process(tree, allow_zero=True)
    covariance = process.tip_covariance(data.names)
    fit = fit_individual_covariance(
        covariance,
        data.values,
        data.species,
        data.individuals,
        data.coordinates,
        dimension=len(traits),
        within=args.within_species_covariance,
        method=getattr(args, "covariance_method", None) or "REML",
    )
    posterior = _node_marginals(tree, data, fit, process)
    observed, errors, counts = _species_summaries(data, fit)
    targets = _parse_targets(args.target)
    selected = [
        node for node in tree.traverse() if _should_output_node(node, observed, targets)
    ]
    table = multivariate_output_table(
        tree, selected, observed, posterior, traits, settings.ci_level, errors=errors
    )
    table = table.rename(
        columns={
            "observed_value": "sample_mean",
            "observed_se": "sample_mean_se_estimate",
        }
    )
    table["estimand"] = [
        "latent_species_mean" if kind == "leaf" else "ancestral_mean"
        for kind in table.node_class
    ]
    table["num_individuals"] = [
        int(counts[name][traits.index(trait)])
        if name in counts and kind == "leaf"
        else 0
        for name, trait, kind in zip(
            table.name, table.trait, table.node_class, strict=True
        )
    ]
    outputs = {"outfile": table}
    if getattr(args, "model_out", None):
        outputs["model_out"] = _model_table(fit, data, args, settings.ci_level)
    if getattr(args, "individual_out", None):
        outputs["individual_out"] = _individual_table(
            data, fit, covariance, settings.ci_level
        )
    if getattr(args, "covariance_out", None):
        outputs["covariance_out"] = multivariate_covariance_table(
            tree, selected, posterior, traits
        )
    paths = {
        key: getattr(args, key)
        for key in (*outputs, "tree_out")
        if getattr(args, key, None) not in (None, "", "-")
    }
    validate_output_targets(paths.values())
    with output_transaction(paths.values()) as staged:
        for key, output in outputs.items():
            if key in paths:
                _write_table(output, staged[paths[key]])
        if "tree_out" in paths:
            tree_args = copy(args)
            tree_args.tree_out = staged[paths["tree_out"]]
            tree_args.tree_annotation = settings.tree_annotation
            write_multivariate_tree(
                tree,
                observed,
                posterior,
                traits,
                tree_args,
                settings.ci_level,
                errors=errors,
            )
        if args.outfile == "-":
            _write_table(table, "-")
    sys.stderr.write(
        f"Joint individual MV-BM: {fit.method}; fit_status={fit.fit_status}; intervals condition on Sigma and W (parameter uncertainty excluded).\n"
    )
