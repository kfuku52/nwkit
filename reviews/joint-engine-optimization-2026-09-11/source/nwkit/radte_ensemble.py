"""Propagate paired input-tree samples without breaking species-age dependence.

The reference fit remains the reported point estimate. Ensemble percentiles
describe conditional refits, not an automatically combined Bayesian posterior.
"""

from types import SimpleNamespace

import numpy as np

from nwkit.radte_inputs import read_inputs
from nwkit.util import read_tree_strings


def read_ensembles(args):
    genes = (
        read_tree_strings(args.gene_tree_ensemble) if args.gene_tree_ensemble else None
    )
    species = (
        read_tree_strings(args.species_tree_ensemble)
        if args.species_tree_ensemble
        else None
    )
    lengths = [len(values) for values in [genes, species] if values is not None]
    if not lengths or min(lengths) < 2:
        raise ValueError("Input ensembles require at least two Newick trees.")
    if len(set(lengths)) != 1:
        raise ValueError(
            "Gene and species ensembles must have equal sample counts; their order defines paired samples."
        )
    return genes, species, lengths[0]


def sample_arguments(args, genes, species, index):
    values = {
        key: value for key, value in vars(args).items() if not key.startswith("_")
    }
    values.update(
        uncertainty="none", gene_tree_ensemble=None, species_tree_ensemble=None
    )
    if genes is not None:
        values["generax_nhx" if args.generax_nhx else "gene_tree"] = genes[index]
    if species is not None:
        values["species_tree"] = species[index]
        # Each entire species chronogram is conditioned on as one joint draw.
        values["species_node_bounds_tsv"] = None
        values["species_node_intervals_tsv"] = None
    values["seed"] = args.seed + index + 1
    return SimpleNamespace(**values)


def validate_sample(reference, sample, args):
    if {n.name for n in sample.gene.leaves()} != {
        n.name for n in reference.gene.leaves()
    }:
        raise ValueError("Gene-tree samples must have exactly the reference gene tips.")
    expected = reference.species_table.set_index("species_event_id")
    actual = sample.species_table.set_index("species_event_id")
    if set(expected.index) != set(actual.index):
        raise ValueError(
            "Species-tree samples must retain the reference species topology and tips."
        )
    if args.species_tree_ensemble and args.species_node_bounds_tsv:
        actual = actual.loc[expected.index]
        if np.any(actual.age < expected.age_min - 1e-8) or np.any(
            actual.age > expected.age_max + 1e-8
        ):
            raise ValueError(
                "Species-age sample lies outside the explicitly supplied hard bounds."
            )


def ensemble_intervals(reference, fit, args):
    from nwkit.radte import run_dating

    genes, species, count = read_ensembles(args)
    samples, sample_ids, failures, methods, clades = [], [], [], [], []
    diagnostics = []
    within_audit: list[dict] = []
    for index in range(count):
        try:
            options = sample_arguments(args, genes, species, index)
            chronology = read_inputs(options)
            validate_sample(reference, chronology, args)
            result, problem, _, _ = run_dating(chronology, options)
            if not problem.feasible(result.parameters):
                raise ValueError("Input-ensemble fit violated chronology constraints.")
            by_group = dict(
                zip(chronology.groups, result.ages * chronology.scale, strict=True)
            )
            samples.append(
                [
                    by_group.get(key, np.nan) / reference.scale
                    for key in reference.groups
                ]
            )
            sample_ids.append(index + 1)
            observed = set(
                zip(
                    chronology.events.gene_clade_id,
                    chronology.events.event_type,
                    strict=True,
                )
            )
            clades.append(
                [
                    (gid, event) in observed
                    for gid, event in zip(
                        reference.events.gene_clade_id,
                        reference.events.event_type,
                        strict=True,
                    )
                ]
            )
            methods.append(
                "branch-marginal"
                if problem.likelihood is None
                else "sequence-marginal"
                if getattr(problem, "marginal", False)
                else "sequence-joint-map"
            )
            from nwkit.radte_components import conditional_intervals

            fit.conditional_intervals.extend(
                conditional_intervals(
                    chronology,
                    result,
                    problem,
                    options,
                    index + 1,
                    methods[-1],
                    audit=within_audit,
                )
            )
            diagnostics.append(result.diagnostics)
        except (ValueError, FloatingPointError) as exc:
            failures.append(dict(sample=index + 1, reason=str(exc)))
    fit.diagnostics.append(f"input_ensemble_successes={len(samples)}/{count}")
    metadata = dict(
        sample_count=count,
        successes=len(samples),
        failures=failures,
        sample_methods=methods,
        sample_ids=sample_ids,
        sample_diagnostics=diagnostics,
        within_fit_diagnostics=within_audit,
        species_chronogram_samples=species is not None,
        pairing="same-Newick-order",
        within_fit_uncertainty=getattr(args, "ensemble_within_uncertainty", "none")
        != "none",
        within_fit_method=getattr(args, "ensemble_within_uncertainty", "none"),
        components_are_combined_posterior=False,
    )
    if not samples:
        fit.interval_status = "unavailable-all-input-ensemble-fits-failed"
        return metadata
    fit.samples = np.asarray(samples)
    fit.sample_ids = sample_ids
    fit.sample_presence = np.isfinite(fit.samples).mean(axis=0)
    fit.sample_clade_presence = np.mean(clades, axis=0)
    if len(samples) < max(20, int(np.ceil(0.9 * count))):
        fit.interval_status = "unavailable-insufficient-successful-input-samples"
        return metadata
    lower = np.full(len(reference.groups), np.nan)
    upper = lower.copy()
    tail = (1 - args.interval_level) / 2
    for j in range(len(reference.groups)):
        values = fit.samples[np.isfinite(fit.samples[:, j]), j]
        if len(values) >= max(20, int(np.ceil(0.9 * len(samples)))):
            lower[j], upper[j] = np.quantile(values, [tail, 1 - tail])
    fit.interval_lower, fit.interval_upper = lower, upper
    fit.interval_status = "input-ensemble-percentiles-conditional-refits"
    return metadata
