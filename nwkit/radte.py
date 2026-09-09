"""Native Reconciliation-Assisted Divergence Time Estimation entry point."""

import hashlib
import json
import sys
import time
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from nwkit import __version__
from nwkit.clade_index import CladeIndex
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.radte_inputs import read_inputs
from nwkit.radte_model import fit_dates, laplace_intervals
from nwkit.radte_sequence import (
    QuadraticLikelihood,
    SequenceLikelihood,
    build_quadratic,
)
from nwkit.radte_sequence_fit import default_sequence_model, fit_sequence_model
from nwkit.radte_studentized import studentized_intervals
from nwkit.radte_uncertainty import bootstrap_intervals, profile_intervals
from nwkit.util import (
    _serialize_newick_node_name,
    copy_tree_iteratively,
    read_tree,
    validate_outputs_do_not_replace_inputs,
)

RADTE_SUFFIXES = {
    "tree": ".dated.nwk",
    "nodes": ".nodes.tsv",
    "species": ".species.tsv",
    "events": ".events.tsv",
    "groups": ".shared-ages.tsv",
    "samples": ".age-samples.tsv",
    "conditional_intervals": ".conditional-intervals.tsv",
    "uncertainty_components": ".uncertainty-components.tsv",
    "likelihood": ".likelihood.json",
    "paml_trace": ".mcmctree-trace.tsv",
    "manifest": ".manifest.json",
}
RADTE_INPUTS = (
    "gene_tree",
    "generax_nhx",
    "species_tree",
    "notung_parsable",
    "reconciliation",
    "species_node_bounds_tsv",
    "species_node_intervals_tsv",
    "species_map_tsv",
    "alignment",
    "gene_tree_ensemble",
    "species_tree_ensemble",
    "likelihood_summary",
)


def radte_paths(prefix):
    if not isinstance(prefix, str) or prefix.strip() in {"", "-"}:
        raise ValueError("--out-prefix must be a nonempty filesystem prefix.")
    return {name: str(prefix) + suffix for name, suffix in RADTE_SUFFIXES.items()}


def _hash_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_options(args):
    if args.seed < 0 or args.starts < 1 or args.maxiter < 1:
        raise ValueError(
            "Seed must be nonnegative; starts and maxiter must be positive."
        )
    if not 0 <= args.rate_correlation < 1:
        raise ValueError("--rate-correlation must satisfy 0 <= rho < 1.")
    if not 0 < args.interval_level < 1:
        raise ValueError("--interval-level must be between zero and one.")
    if args.rate_sd is not None and args.rate_sd < 0:
        raise ValueError("--rate-sd must be nonnegative.")
    if args.max_age is not None and args.max_age <= 0:
        raise ValueError("--max-age must be positive.")
    validate_backend_options(args)
    ensemble = bool(args.gene_tree_ensemble or args.species_tree_ensemble)
    if ensemble != (args.uncertainty == "input-ensemble"):
        raise ValueError(
            "Tree ensembles require --uncertainty input-ensemble, and that method requires an ensemble input."
        )
    if getattr(args, "ensemble_within_uncertainty", "none") != "none" and not ensemble:
        raise ValueError("--ensemble-within-uncertainty requires an input ensemble.")
    if args.gene_tree_ensemble and args.likelihood_summary:
        raise ValueError(
            "Gene-tree ensembles cannot reuse a single-tree likelihood summary; supply the alignment instead."
        )
    if args.alignment and args.likelihood_summary:
        raise ValueError("Choose --alignment or --likelihood-summary, not both.")
    if args.inference != "auto" and not (args.alignment or args.likelihood_summary):
        raise ValueError(
            "--inference is a sequence option; tree-only inference already integrates branch rates."
        )
    if args.inference == "marginal" and args.likelihood == "exact":
        raise ValueError(
            "--inference marginal requires a quadratic likelihood; choose auto or quadratic."
        )
    sequence_options = [
        "substitution_model",
        "kappa",
        "gamma_shape",
        "gamma_categories",
        "likelihood",
        "gtr_exchangeabilities",
    ]
    if not args.alignment and any(
        getattr(args, name) is not None for name in sequence_options
    ):
        raise ValueError(
            "Sequence model/likelihood options require --alignment; a summary already specifies its likelihood."
        )
    for name in RADTE_INPUTS:
        if getattr(args, name, None) == "-":
            raise ValueError(
                f"--{name.replace('_', '-')} requires a file path for reproducible RADTE inputs."
            )


def validate_backend_options(args):
    paml_options = [
        key
        for key, value in vars(args).items()
        if key.startswith("mcmctree_") and value is not None
    ]
    if args.backend == "native" and paml_options:
        raise ValueError("MCMCTree controls require --backend mcmctree.")
    if args.backend != "mcmctree":
        return
    if not args.alignment or args.likelihood_summary:
        raise ValueError(
            "MCMCTree requires --alignment and cannot use a native likelihood summary."
        )
    if (
        args.uncertainty != "none"
        or args.rate_sd is not None
        or args.rate_correlation != 0
        or args.inference != "auto"
    ):
        raise ValueError(
            "MCMCTree supplies its own clock, inference and posterior intervals; native options cannot be combined."
        )
    if any(
        getattr(args, key) is not None
        for key in [
            "kappa",
            "gtr_exchangeabilities",
            "gamma_shape",
            "gamma_categories",
            "likelihood",
        ]
    ):
        raise ValueError(
            "Native sequence-model controls cannot be used with the MCMCTree reference."
        )
    for key in paml_options:
        value = getattr(args, key)
        if key.endswith("_prior"):
            numbers = [float(item) for item in value.split()]
            if len(numbers) != 3 or not np.isfinite(numbers).all() or min(numbers) <= 0:
                raise ValueError(
                    f"--{key.replace('_', '-')} requires three finite positive numbers."
                )
            setattr(args, key, " ".join(format(item, ".15g") for item in numbers))
        elif key not in {"mcmctree_bin", "mcmctree_clock", "mcmctree_likelihood"}:
            if value < (0 if key == "mcmctree_burnin" else 1):
                raise ValueError(f"Invalid positive MCMCTree control: {key}")


def _edge_ids(c):
    index = CladeIndex(c.gene)
    return [index.clade_id_for_node(n) for n in c.edges]


def read_likelihood_summary(path, c):
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    if data.get("schema") != "nwkit-radte-quadratic-v1":
        raise ValueError("Unsupported or unavailable likelihood summary.")
    edge_ids = data["edge_ids"]
    expected = _edge_ids(c)
    if len(edge_ids) != len(set(edge_ids)) or set(edge_ids) != set(expected):
        raise ValueError(
            "Likelihood summary must match the rooted gene topology and all tip labels."
        )
    mapping = np.asarray(data["mapping"], dtype=int)
    center = np.asarray(data["center"], dtype=float)
    gradient = np.asarray(data["gradient"], dtype=float)
    hessian = np.asarray(data["hessian"], dtype=float)
    nll = float(data["nll"])
    size = len(expected) - 1
    if (
        center.shape != (size,)
        or gradient.shape != (size,)
        or hessian.shape != (size, size)
        or mapping.shape != (len(expected),)
        or set(mapping) != set(range(size))
        or not all(
            np.isfinite(x).all() for x in [center, gradient, hessian, np.array(nll)]
        )
        or not np.allclose(hessian, hessian.T, atol=1e-9)
        or np.linalg.eigvalsh(hessian)[0] <= 0
    ):
        raise ValueError(
            "Invalid likelihood summary dimensions, values, or positive-definite Hessian."
        )
    lookup = {key: i for i, key in enumerate(edge_ids)}
    mapping = mapping[[lookup[key] for key in expected]]
    roots = [i for i, n in enumerate(c.edges) if n.up is c.gene]
    counts = np.bincount(mapping)
    if mapping[roots[0]] != mapping[roots[1]] or sorted(counts) != [1] * (size - 1) + [
        2
    ]:
        raise ValueError("Likelihood summary must combine exactly the two root edges.")
    return QuadraticLikelihood(center, gradient, hessian, nll, mapping), data


def _likelihood_data(likelihood, c, metadata):
    if not isinstance(likelihood, QuadraticLikelihood):
        return dict(
            schema="nwkit-radte-likelihood-unavailable-v1",
            reason="exact-or-tree-only-mode",
        )
    return dict(
        schema="nwkit-radte-quadratic-v1",
        edge_ids=_edge_ids(c),
        mapping=likelihood.mapping.tolist(),
        center=likelihood.center.tolist(),
        gradient=likelihood.gradient.tolist(),
        hessian=likelihood.hessian.tolist(),
        nll=likelihood.nll,
        metadata=metadata,
    )


def _fit(c, args, likelihood):
    return fit_dates(
        c,
        rho=args.rate_correlation,
        likelihood=likelihood,
        rate_sd=args.rate_sd,
        starts=args.starts,
        maxiter=args.maxiter,
        seed=args.seed,
        inference=args.inference,
    )


def run_dating(c, args):
    likelihood = None
    metadata = {}
    diagnostics = []
    if args.likelihood_summary:
        likelihood, data = read_likelihood_summary(args.likelihood_summary, c)
        metadata = dict(data.get("metadata", {}))
    elif args.alignment:
        metadata = dict(
            model=args.substitution_model or default_sequence_model(args.alignment),
            kappa=args.kappa if args.kappa is not None else 2.0,
            gamma_shape=args.gamma_shape if args.gamma_shape is not None else 1.0,
            gamma_categories=args.gamma_categories
            if args.gamma_categories is not None
            else 4,
        )
        if args.kappa is not None and metadata["model"] != "hky":
            raise ValueError("--kappa applies only to the HKY sequence model.")
        exchange = None
        if args.gtr_exchangeabilities is not None:
            if metadata["model"] != "gtr":
                raise ValueError(
                    "--gtr-exchangeabilities requires --substitution-model gtr."
                )
            exchange = np.array(
                [float(value) for value in args.gtr_exchangeabilities.split(",")]
            )
        exact = SequenceLikelihood(
            c, args.alignment, exchangeabilities=exchange, **metadata
        )
        initial_lengths = np.array([n.dist for n in c.edges])
        metadata["prefit"] = fit_sequence_model(
            exact,
            initial_lengths,
            fit_kappa=metadata["model"] == "hky" and args.kappa is None,
            fit_gtr=metadata["model"] == "gtr" and exchange is None,
            fit_gamma=args.gamma_shape is None,
            maxiter=args.maxiter,
        )
        metadata.update(
            kappa=exact.kappa,
            gamma_shape=exact.gamma_shape,
            frequencies=exact.pi.tolist(),
            gamma_rates=exact.rates.tolist(),
            exchangeabilities=exact.exchangeabilities.tolist()
            if exact.model == "gtr"
            else None,
            alignment_sha256=_hash_file(args.alignment),
        )
        initial_lengths = getattr(exact, "initial_lengths", initial_lengths)
        likelihood = exact
        if args.likelihood != "exact":
            try:
                likelihood = build_quadratic(exact, initial_lengths, args.maxiter)
            except ValueError as exc:
                if args.likelihood == "quadratic" or args.inference == "marginal":
                    raise
                diagnostics.append("quadratic_unavailable_using_exact: " + str(exc))
    fit, problem = _fit(c, args, likelihood)
    if isinstance(likelihood, QuadraticLikelihood):
        lengths = fit.rates * c.scale * c.durations(fit.ages)
        if likelihood.exact is None:
            combined = np.bincount(likelihood.mapping, weights=lengths)
            if np.max(abs(np.log(combined) - likelihood.center)) > 0.5:
                raise ValueError(
                    "Summary-only solution left the local log-length trust region; rerun with --alignment."
                )
            diagnostics.append("summary_only_approximation_not_independently_verified")
        else:
            checks = (
                problem.approximation_checks(fit.parameters)
                if getattr(problem, "marginal", False)
                else [likelihood.check(lengths)]
            )
            valid = all(check[0] for check in checks)
            error = max(check[1] for check in checks)
            gradient_error = max(check[2] for check in checks)
            metadata["approximation_nll_error"] = error
            metadata["approximation_standardized_score_error"] = gradient_error
            if not valid:
                if args.likelihood == "quadratic" or args.inference == "marginal":
                    raise ValueError(
                        "Quadratic likelihood failed exact-value/gradient validation; use --likelihood exact or auto."
                    )
                diagnostics.append("quadratic_failed_validation_refitted_exact")
                likelihood = likelihood.exact
                fit, problem = _fit(c, args, likelihood)
    fit.diagnostics.extend(diagnostics)
    if args.uncertainty == "laplace":
        laplace_intervals(fit, problem, args.interval_level)
    elif args.uncertainty == "studentized":
        studentized_intervals(fit, problem, args.interval_level)
    elif args.uncertainty == "profile":
        profile_intervals(
            fit,
            problem,
            level=args.interval_level,
            starts=args.starts,
            maxiter=args.maxiter,
            seed=args.seed,
        )
    elif args.uncertainty == "bootstrap":
        bootstrap_intervals(
            fit,
            problem,
            replicates=args.bootstrap_replicates,
            level=args.interval_level,
            rho=args.rate_correlation,
            starts=args.starts,
            maxiter=args.maxiter,
            seed=args.seed,
            rate_sd=args.rate_sd,
        )
    elif args.uncertainty == "input-ensemble":
        from nwkit.radte_ensemble import ensemble_intervals

        fit.ensemble_metadata = ensemble_intervals(c, fit, args)
    return fit, problem, _likelihood_data(likelihood, c, metadata), metadata


def dated_newick(c, fit):
    from ete4.parser.newick import PARSERS

    output = copy_tree_iteratively(c.gene)
    by_id = dict(
        zip(c.events.gene_clade_id, fit.ages[c.group_by_node] * c.scale, strict=True)
    )
    index = CladeIndex(output)
    for node in output.traverse():
        node.dist = (
            0.0
            if node is output
            else by_id[index.clade_id_for_node(node.up)]
            - by_id[index.clade_id_for_node(node)]
        )
    # ETE's default six significant figures can destroy equality of shared
    # ages along paths. Preserve all double precision in time branch lengths.
    parser = {
        key: [dict(item) for item in fields] for key, fields in PARSERS[1].items()
    }
    for fields in parser.values():
        fields[1]["write"] = lambda value: format(float(value), ".17g")
    parser["leaf"][0]["write"] = lambda value: _serialize_newick_node_name(
        value, is_internal=False
    )
    parser["internal"][0]["write"] = lambda value: _serialize_newick_node_name(
        value, is_internal=True
    )
    return (
        "[&R]"
        + output.write(
            parser=parser, format_root_node=True, props=["support", "D", "S", "H"]
        )
        + "\n"
    )


def result_tables(c, fit, *, bound_policy="hard"):
    nodes = c.events.copy()
    nodes["estimated_age"] = fit.ages[c.group_by_node] * c.scale
    nodes["interval_lower"] = (
        np.nan
        if fit.interval_lower is None
        else fit.interval_lower[c.group_by_node] * c.scale
    )
    nodes["interval_upper"] = (
        np.nan
        if fit.interval_upper is None
        else fit.interval_upper[c.group_by_node] * c.scale
    )
    nodes["interval_status"] = fit.interval_status
    nodes["age_identifiability"] = (
        "absolute-scale-unidentified"
        if "absolute_age_scale_unidentified_within_hard_bounds" in fit.diagnostics
        else "nonunique-optimum"
        if "nonunique_age_optimum" in fit.diagnostics
        else "no-nonuniqueness-detected"
    )
    nodes["sample_event_presence"] = (
        np.nan if fit.sample_presence is None else fit.sample_presence[c.group_by_node]
    )
    nodes["sample_clade_presence"] = (
        np.nan if fit.sample_clade_presence is None else fit.sample_clade_presence
    )
    edge_by_node = {node: i for i, node in enumerate(c.edges)}
    nodes["estimated_rate"] = [
        np.nan if n is c.gene else fit.rates[edge_by_node[n]] for n in c.nodes
    ]
    nodes["bound_policy"] = bound_policy
    nodes["constraint_status"] = "retained"
    if bound_policy == "PAML-soft-prior":
        unselected = (nodes.event_type == "duplication") & (
            nodes.parent_gene_clade_id != ""
        )
        nodes.loc[unselected, "constraint_status"] = "not-selected"
    nodes["within_original_bounds"] = (
        nodes.estimated_age >= nodes.age_min - c.scale * 1e-8
    ) & (nodes.estimated_age <= nodes.age_max + c.scale * 1e-8)
    rows = []
    for group, members in nodes[nodes.event_type == "speciation"].groupby(
        "shared_age_id", sort=True
    ):
        ages = members.estimated_age.to_numpy()
        rows.append(
            dict(
                shared_age_id=group,
                species_name=members.species_name.iloc[0],
                member_count=len(members),
                estimated_age=ages[0],
                max_member_age_difference=float(np.ptp(ages)),
                interval_lower=members.interval_lower.iloc[0],
                interval_upper=members.interval_upper.iloc[0],
                interval_status=fit.interval_status,
            )
        )
    groups = pd.DataFrame(
        rows,
        columns=[
            "shared_age_id",
            "species_name",
            "member_count",
            "estimated_age",
            "max_member_age_difference",
            "interval_lower",
            "interval_upper",
            "interval_status",
        ],
    )
    species = c.species_table.copy()
    group_ids = {key: i for i, key in enumerate(c.groups)}
    species["estimated_age"] = [
        fit.ages[group_ids["S:" + sid]] * c.scale for sid in species.species_event_id
    ]
    species_indices = np.array(
        [group_ids["S:" + sid] for sid in species.species_event_id]
    )
    species["interval_lower"] = (
        np.nan
        if fit.interval_lower is None
        else fit.interval_lower[species_indices] * c.scale
    )
    species["interval_upper"] = (
        np.nan
        if fit.interval_upper is None
        else fit.interval_upper[species_indices] * c.scale
    )
    species["interval_status"] = fit.interval_status
    represented = set(nodes.loc[nodes.event_type == "speciation", "shared_age_id"])
    species["estimation_status"] = [
        "represented-gene-event" if "S:" + sid in represented else "calibration-only"
        for sid in species.species_event_id
    ]
    if fit.ensemble_metadata and fit.ensemble_metadata.get(
        "species_chronogram_samples"
    ):
        species.loc[
            species.estimation_status == "calibration-only", "estimation_status"
        ] = "input-chronogram-samples"
    species.loc[
        species.estimation_status == "calibration-only",
        ["interval_lower", "interval_upper"],
    ] = np.nan
    if bound_policy == "PAML-soft-prior":
        absent = species.estimation_status == "calibration-only"
        species.loc[absent, "estimated_age"] = np.nan
        species.loc[absent, "estimation_status"] = "not-sampled"
    if fit.samples is None:
        samples = pd.DataFrame(columns=["replicate", "shared_age_id", "estimated_age"])
    else:
        samples = pd.DataFrame(
            [
                dict(
                    replicate=i + 1 if fit.sample_ids is None else fit.sample_ids[i],
                    shared_age_id=key,
                    estimated_age=float(age * c.scale),
                )
                for i, sample in enumerate(fit.samples)
                for key, age in zip(c.groups, sample, strict=True)
                if np.isfinite(age)
            ]
        )
    return dict(
        nodes=nodes, species=species, events=c.events, groups=groups, samples=samples
    )


def radte_main(args):
    start = time.monotonic()
    validate_options(args)
    requested_options = {
        key: value
        for key, value in vars(args).items()
        if key != "handler" and not key.startswith("_")
    }
    paths = radte_paths(args.out_prefix)
    figure_out = getattr(args, "figure_out", None)
    if figure_out:
        from nwkit.result_plot import figure_format

        figure_format(figure_out)
        paths["figure"] = figure_out
    inputs = [
        ("--" + name.replace("_", "-"), getattr(args, name, None))
        for name in RADTE_INPUTS
    ]
    validate_outputs_do_not_replace_inputs(
        inputs, list(paths.items()), label="RADTE output"
    )
    validate_output_targets(paths.values(), follow_symlinks=False)
    hashes = {name: _hash_file(path) for name, path in inputs if path}
    c = read_inputs(args)
    paml_trace = pd.DataFrame(columns=["chain", "Gen"])
    problem = None
    if args.backend == "mcmctree":
        from nwkit.radte_paml import run_mcmctree

        fit, paml_trace, model_metadata = run_mcmctree(c, args)
        likelihood_data = dict(
            schema="nwkit-radte-likelihood-unavailable-v1", reason="external-backend"
        )
    else:
        fit, problem, likelihood_data, model_metadata = run_dating(c, args)
        if not problem.feasible(fit.parameters):
            raise ValueError("Dated output failed chronology validation.")
    tables = result_tables(
        c, fit, bound_policy="hard" if args.backend == "native" else "PAML-soft-prior"
    )
    from nwkit.radte_components import component_tables

    tables.update(component_tables(c, fit))
    tables["paml_trace"] = paml_trace
    if args.backend == "native" and not tables["nodes"].within_original_bounds.all():
        raise ValueError("Dated output failed calibration validation.")
    text = dated_newick(c, fit)
    model = "marginal-lognormal"
    if args.backend == "mcmctree":
        model = "mcmctree-posterior-mean"
    elif problem is not None and problem.likelihood is not None:
        model = (
            "sequence-marginal-quadratic"
            if getattr(problem, "marginal", False)
            else "sequence-empirical-bayes-map"
        )
    manifest = dict(
        schema="nwkit-radte-run-v1",
        status="complete",
        version=__version__,
        method=model,
        clock=("paml-" + str(args.mcmctree_clock or 2))
        if args.backend == "mcmctree"
        else ("independent" if args.rate_correlation == 0 else "lineage-ar1"),
        rate_correlation=args.rate_correlation,
        log_rate_sd=fit.log_rate_sd,
        log_rate_mean=fit.log_rate_mean,
        objective=fit.objective,
        shared_speciation_ages=True,
        experimental_native_estimator=args.backend == "native",
        calibration_policy="hard-all-events"
        if args.backend == "native"
        else "PAML-soft-root-and-speciation-priors",
        uncertainty=fit.interval_status,
        interval_level=args.interval_level,
        conditional_on="rooted reconciliation, species calibration domain, substitution model",
        diagnostics=fit.diagnostics,
        optimizer_attempts=fit.attempts,
        sequence_model=model_metadata,
        input_ensemble=fit.ensemble_metadata,
        seed=args.seed,
        input_sha256=hashes,
        elapsed_seconds=time.monotonic() - start,
        python=sys.version.split()[0],
        numpy=np.__version__,
        scipy=scipy.__version__,
        options=requested_options,
    )
    # One transaction publishes all artifacts, including empty optional outputs,
    # so reruns cannot leave stale intervals or a previous sequence summary.
    Path(args.out_prefix).absolute().parent.mkdir(parents=True, exist_ok=True)
    with output_transaction(list(paths.values()), follow_symlinks=False) as staged:
        staged.write_text(paths["tree"], lambda handle: handle.write(text))
        for key, frame in tables.items():
            staged.write_text(
                paths[key], partial(frame.to_csv, sep="\t", index=False, na_rep="NA")
            )
        staged.write_text(
            paths["likelihood"],
            lambda handle: json.dump(
                likelihood_data, handle, indent=2, allow_nan=False
            ),
        )
        if figure_out:
            from nwkit.result_plot import save_result_figure
            from nwkit.result_plot_cli import result_figure_options
            from nwkit.result_plot_data import dating_plot_data

            data = dating_plot_data(
                read_tree(text, "auto", True),
                c.species,
                tables["nodes"],
                tables["species"],
                manifest,
            )
            save_result_figure(
                data,
                staged[figure_out],
                image_format=figure_format(figure_out),
                **result_figure_options(args),
            )
        manifest["output_sha256"] = {
            key: _hash_file(staged[path])
            for key, path in paths.items()
            if key != "manifest"
        }
        staged.write_text(
            paths["manifest"],
            lambda handle: json.dump(manifest, handle, indent=2, allow_nan=False),
        )
    print(f"Dated tree: {Path(paths['tree']).resolve()}", file=sys.stderr)
    print(
        f"Shared speciation ages: enforced; interval status: {fit.interval_status}",
        file=sys.stderr,
    )
    print("Method: " + model, file=sys.stderr)
    for diagnostic in fit.diagnostics:
        print("Diagnostic: " + diagnostic, file=sys.stderr)
