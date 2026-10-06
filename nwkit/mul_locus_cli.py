"""Staged outputs for the experimental finite-grid locus Monte Carlo model."""

import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import pandas as pd

from nwkit import __version__
from nwkit.mul_locus import signature_size
from nwkit.mul_locus_mc import (
    build_bank,
    calibrate,
    category_bound,
    integration_alpha,
    make_tasks,
    parameter_record,
    validate_model,
)
from nwkit.mul_msc import dated_text, read_msc_inputs
from nwkit.mul_msc_fit import species_topology_signature
from nwkit.mul_reconcile_model import validate_binary
from nwkit.output_transaction import output_transaction
from nwkit.util import validate_outputs_do_not_replace_inputs


def json_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate locus JSON key: {key}")
        result[key] = value
    return result


def finite_json(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [finite_json(item) for item in value]
    return value


def bank_record(bank, args=None):
    from nwkit.mul_locus_integral import IntegratedBank

    if isinstance(bank, IntegratedBank):
        return integrated_record(bank, args)
    return {
        "candidate": bank.candidate,
        "h2": bank.h2,
        "grid": bank.grid,
        "parameters": parameter_record(bank),
        "population_tree": dated_text(bank.population, args),
        "population_leaf_species": {
            leaf.name: leaf.props.get("mul_species", leaf.name)
            for leaf in bank.population.leaves()
        },
        "samples": bank.samples,
        "attempts": bank.attempts,
        "counts": [
            {"signature": key, "hits": count} for key, count in bank.counts.items()
        ],
        "strata": None
        if bank.strata is None
        else [
            {
                **{key: value for key, value in stratum.items() if key != "counts"},
                "counts": [
                    {"signature": key, "hits": count}
                    for key, count in stratum["counts"].items()
                ],
            }
            for stratum in bank.strata
        ],
    }


def integrated_record(bank, args=None):
    from dataclasses import asdict

    return {
        "candidate": bank.candidate,
        "h2": bank.h2,
        "grid": bank.grid,
        "parameters": parameter_record(bank),
        "population_tree": dated_text(bank.population, args),
        "population_leaf_species": {
            leaf.name: leaf.props.get("mul_species", leaf.name)
            for leaf in bank.population.leaves()
        },
        "samples": bank.samples,
        "attempts": bank.attempts,
        "method": bank.method,
        "interval_method": bank.interval_method,
        "strata": [
            {
                "condition": s["condition"],
                "weight": s["weight"],
                "samples": s["samples"],
                "selected": asdict(s["selected"]),
                "patterns": [
                    {"signature": key, **asdict(value)}
                    for key, value in s["patterns"].items()
                ],
            }
            for s in bank.strata
        ],
    }


def mul_locus_main(args):
    msc_only = (
        "species_time_unit",
        "effective_population_size",
        "hybridization_age",
        "max_coalescent_states",
        "max_coalescent_assignments",
        "msc_fit",
        "hybridization_age_bounds",
        "population_size_bounds",
        "msc_grid_points",
        "msc_fit_starts",
        "msc_maxiter",
        "msc_max_evaluations",
        "msc_profile_out",
    )
    if any(getattr(args, name, None) is not None for name in msc_only):
        raise ValueError("locus-mc uses its explicit JSON, not MSC-only options.")
    if args.locus_model is None or args.tree_out is not None:
        raise ValueError(
            "locus-mc requires --locus-model; no inferred point --tree-out."
        )
    replicates = args.locus_bootstrap if args.locus_bootstrap is not None else 0
    null_calibration = getattr(args, "locus_null_calibration", None)
    if replicates < 0 or (
        (args.locus_calibration_out is not None or null_calibration is not None)
        and replicates == 0
    ):
        raise ValueError("Calibration output requires positive bootstrap replicates.")
    species, genes, parser, paths = read_msc_inputs(
        args, extra_outputs=("locus_calibration_out",)
    )
    validate_outputs_do_not_replace_inputs(
        [("--locus-model", args.locus_model)], list(paths.items())
    )
    with open(args.locus_model) as handle:
        model = json.load(handle, object_pairs_hook=json_pairs)
    parameters = validate_model(model, species)
    observations = []
    for gene in genes:
        validate_binary(gene, "Gene tree")
        observation = species_topology_signature(gene, parser)
        if signature_size(observation) > model["max_observed_tips"]:
            raise ValueError(
                "Input family violates the declared size selection; no filtering."
            )
        if any(
            parser.parse(leaf.name).species_label not in model["detection"]
            for leaf in gene.leaves()
        ):
            raise ValueError("Gene tree has an unmatched species.")
        observations.append(observation)
    tasks, excluded = make_tasks(
        species, args.h1, args.h2, parameters, args.max_candidates
    )
    alpha = integration_alpha(model, len(tasks))
    categories = category_bound(len(model["detection"]), model["max_observed_tips"])
    if args.cpus == 1:
        banks = [build_bank(task, model) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=min(args.cpus, len(tasks))) as pool:
            banks = list(pool.map(partial(build_bank, model=model), tasks))
    scorer_args = {}
    if model.get("integration") in ("detection-rb", "hybrid-rb"):
        from nwkit.mul_locus_integral import score_integrated_bank

        scorer_args["scorer"] = score_integrated_bank
    fit, calibration = calibrate(
        banks,
        observations,
        model,
        alpha,
        replicates,
        null_calibration=null_calibration or "plug-in",
        **scorer_args,
    )
    write_results(
        args,
        model,
        observations,
        banks,
        fit,
        calibration,
        excluded,
        paths,
        categories,
        alpha,
    )


def write_results(
    args,
    model,
    observations,
    banks,
    fit,
    calibration,
    excluded,
    paths,
    categories,
    alpha,
):
    ranked = sorted(
        fit["rows"], key=lambda r: (-r["log_likelihood"], r["mul.tree"], r["grid"])
    )
    best = ranked[0]
    ambiguous = [r for r in ranked if r["mc_upper"] >= best["mc_lower"]]
    table = []
    for row in ranked:
        table.append(
            {
                **{
                    key: row[key]
                    for key in (
                        "mul.tree",
                        "h2.node",
                        "grid",
                        "log_likelihood",
                        "mc_lower",
                        "mc_upper",
                    )
                },
                "parameters": json.dumps(row["parameters"], sort_keys=True),
                "status": "mc-zero-estimate"
                if not math.isfinite(row["log_likelihood"])
                else "mc-grid-evaluated",
            }
        )
    table.extend(
        {
            "mul.tree": r["mul.tree"],
            "h2.node": r["h2.node"],
            "grid": r["grid"],
            "status": "excluded",
            "reason": r["reason"],
        }
        for r in excluded
    )
    details = []
    for bank, row in zip(banks, fit["rows"], strict=True):
        for i, observation in enumerate(observations, 1):
            pattern = next(p for p in row["patterns"] if p["signature"] == observation)
            details.append(
                {
                    "mul.tree": bank.candidate,
                    "grid": bank.grid,
                    "gene.tree": i,
                    "hits": pattern.get("hits"),
                    "samples": bank.samples,
                    "integration": model.get("integration", "selected-histogram"),
                    "probability_lower": pattern["probability_lower"],
                    "probability_upper": pattern["probability_upper"],
                    "probability_estimate": pattern["probability_estimate"],
                }
            )
    metadata = {
        "schema": "nwkit-mul-locus-mc-v1",
        "method": "linear-DL-daughter-bounded-MLC-finite-grid-MC-v1",
        "nwkit_version": __version__,
        "model": model,
        "hypothesis_scope": {
            "h1": args.h1,
            "h2": args.h2,
            "species_tree": dated_text(
                next(b.population for b in banks if b.candidate == 0), args
            ),
        },
        "observations": observations,
        "num_gene_trees": len(observations),
        "mc_category_union_bound": categories,
        "per_probability_alpha": alpha,
        "best_numerical_grid_point": {k: best[k] for k in ("mul.tree", "grid")},
        "mc_overlapping_grid_points": [
            {k: r[k] for k in ("mul.tree", "grid")} for r in ambiguous
        ],
        "ranking_status": "mc-resolved" if len(ambiguous) == 1 else "mc-unresolved",
        "scores": fit,
        "excluded": excluded,
        "calibration": calibration,
        "banks": [bank_record(bank, args) for bank in banks],
        "limitations": [
            "Experimental small-family generative model; finite-grid maximization is not continuous MLE.",
            "One ancestral origin locus, explicit finite ancestral DL stem, fixed species ages and shared diploid Ne/rates.",
            "One disomic direct-parent allotetraploid; no autopolyploidy, multiple events, hemiplasy or subgenome exchange.",
            "Gene copies are distinct loci, not allelic samples; detection occurs after full gene genealogy sampling.",
            "MC intervals bound bank sampling error simultaneously over the finite observation universe, not biological uncertainty.",
            "No pseudocounts or history truncation; resource or missing-support failures abort all outputs.",
            "Null calibration repeats the full supplied candidate/grid search; its method-specific restrictions are recorded in calibration.limitations.",
            "Grid-supremum calibration covers only the supplied finite null points, not off-grid nuisance parameters or model misspecification.",
            "Fixed rooted gene trees; gene-tree bootstrap stability is separate from event calibration.",
        ],
    }
    tables = {
        "outfile": pd.DataFrame(finite_json(table)),
        "report": pd.DataFrame(details),
        "check_out": pd.DataFrame(details),
    }
    if calibration is not None:
        tables["locus_calibration_out"] = pd.DataFrame(finite_json(calibration["rows"]))
    with output_transaction(paths.values()) as staged:
        for name, path in paths.items():
            if name in tables:
                staged.write_text(
                    path,
                    lambda handle, name=name: tables[name].to_csv(
                        handle, sep="\t", index=False
                    ),
                )
            else:
                staged.write_text(
                    path,
                    lambda handle: json.dump(
                        finite_json(metadata), handle, indent=2, allow_nan=False
                    ),
                )
    if args.outfile == "-":
        tables["outfile"].to_csv(sys.stdout, sep="\t", index=False)
    sys.stderr.write(
        f"Experimental locus MC: {len(banks)} grid banks, {metadata['ranking_status']}; finite-grid estimates, not calibrated biological support.\n"
    )
    if calibration is not None:
        label = (
            "Finite-grid supremum MC"
            if calibration["method"] == "finite-grid-supremum-Monte-Carlo-test"
            else "Plug-in null bootstrap"
        )
        scope = (
            "finite supplied null grid only"
            if label.startswith("Finite")
            else "not uniform composite-null coverage"
        )
        sys.stderr.write(
            f"{label} ({calibration['replicates']} replicates per null point): "
            f"P={calibration['p_value']:.8g}, MC-overlap interval "
            f"[{calibration['mc_p_lower']:.8g}, {calibration['mc_p_upper']:.8g}]; "
            f"{scope}.\n"
        )
