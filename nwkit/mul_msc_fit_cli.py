"""Fitted conditional MSC outputs, distinct from fixed MSC and D+L schemas."""

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict

import pandas as pd

from nwkit import __version__
from nwkit.mul_msc import _time_scale, dated_text, read_msc_inputs
from nwkit.mul_msc_fit import FitSettings, fit_candidate
from nwkit.mul_msc_model import dated_candidates, validate_sampling
from nwkit.output_transaction import output_transaction


def settings_from_args(args):
    mode = args.msc_fit
    age_fit, ne_fit = mode in {"age", "joint"}, mode in {"ne", "joint"}
    if ne_fit and (
        args.species_time_unit != "generations"
        or args.effective_population_size is not None
    ):
        raise ValueError(
            "Ne fitting requires generations and --population-size-bounds, not fixed Ne/coalescent units."
        )
    if age_fit and args.hybridization_age is not None:
        raise ValueError(
            "Fitted age uses --hybridization-age-bounds; omit fixed --hybridization-age."
        )
    if args.species_time_unit is None:
        raise ValueError("MSC fitting requires explicit --species-time-unit.")
    scale = None
    if not ne_fit:
        scale = _time_scale(args, require_age=not age_fit)
    defaults = {
        "grid_points": 5,
        "fit_starts": 3,
        "maxiter": 200,
        "max_evaluations": 5000,
    }
    values = {
        ("starts" if key == "fit_starts" else key): getattr(args, "msc_" + key)
        if getattr(args, "msc_" + key) is not None
        else value
        for key, value in defaults.items()
    }
    return FitSettings(
        mode,
        tuple(args.hybridization_age_bounds)
        if args.hybridization_age_bounds is not None
        else None,
        tuple(args.population_size_bounds)
        if args.population_size_bounds is not None
        else None,
        None if age_fit else args.hybridization_age,
        scale,
        **values,
        max_states=args.max_coalescent_states or 100000,
        max_assignments=args.max_coalescent_assignments or 10000,
    )


def _fit_task(task):
    return fit_candidate(*task)


def mul_msc_fit_main(args):
    species, genes, parser, paths = read_msc_inputs(
        args, extra_outputs=("msc_profile_out",)
    )
    settings = settings_from_args(args)
    candidates, polyploid = dated_candidates(
        species,
        args.h1,
        args.h2,
        settings.fixed_age,
        max_candidates=args.max_candidates,
        age_bounds=settings.age_bounds,
    )
    validate_sampling(genes, parser, polyploid, species.leaf_names())
    tasks = [
        (candidate, genes, parser, settings)
        for candidate in candidates
        if candidate.status == "evaluated"
    ]
    if args.cpus == 1:
        results = list(map(_fit_task, tasks))
    else:
        with ProcessPoolExecutor(max_workers=min(args.cpus, len(tasks))) as pool:
            results = list(pool.map(_fit_task, tasks, chunksize=1))
    results.sort(key=lambda result: (-result["log_likelihood"], result["id"]))
    write_fit_results(args, candidates, polyploid, genes, settings, results, paths)


def write_fit_results(args, candidates, polyploid, genes, settings, results, paths):
    best = results[0]
    tied = [
        result["id"]
        for result in results
        if best["log_likelihood"] - result["log_likelihood"] <= 1e-7
    ]
    if args.tree_out is not None and (
        len(tied) != 1 or best["diagnostics"]["status"] != "locally-distinguishable"
    ):
        raise ValueError(
            "MSC fitted tree is unresolved: tied parent candidates or non-interior/locally unidentified parameters; omit --tree-out and inspect diagnostics."
        )
    by_id = {candidate.id: candidate for candidate in candidates}
    rows = []
    for result in results:
        candidate = by_id[result["id"]]
        point = result["diagnostics"]["status"] == "locally-distinguishable"
        rows.append(
            {
                "mul.tree": candidate.id,
                "h1.node": candidate.h1,
                "h2.node": candidate.h2,
                "log_likelihood": result["log_likelihood"],
                "delta_log_likelihood": best["log_likelihood"]
                - result["log_likelihood"],
                "hybridization_age": result["estimates"].get(
                    "hybridization_age", settings.fixed_age
                ),
                "effective_population_size": result["estimates"].get(
                    "effective_population_size", args.effective_population_size
                ),
                "dated.tree": dated_text(result["candidate"].tree, args)
                if point
                else None,
                "status": result["diagnostics"]["status"],
                "reason": result["diagnostics"]["meaning"],
            }
        )
    rows.extend(
        {
            "mul.tree": c.id,
            "h1.node": c.h1,
            "h2.node": c.h2,
            "log_likelihood": None,
            "delta_log_likelihood": None,
            "hybridization_age": None,
            "effective_population_size": None,
            "dated.tree": None,
            "status": "excluded",
            "reason": c.reason,
        }
        for c in candidates
        if c.status == "excluded"
    )
    details = [row for result in results for row in result["gene_rows"]]
    profiles = [row for result in results for row in result["profiles"]]
    metadata = {
        "schema_version": 1,
        "method": "conditional-direct-parent-MUL-MSC-bounded-fit-v1",
        "score_schema": "nwkit-mul-msc-fit-v1",
        "nwkit_version": __version__,
        "species_time_unit": args.species_time_unit,
        "species_ages": "fixed input",
        "num_gene_trees": len(genes),
        "polyploid_species": list(polyploid),
        "fit_settings": asdict(settings),
        "fitted_parameters": list(best["estimates"]),
        "estimated_parameters": [
            name for name, value in best["estimates"].items() if value is not None
        ],
        "best_hypotheses": tied,
        "reported_hypothesis": best["id"],
        "score_tie_tolerance": 1e-7,
        "scores": rows,
        "gene_likelihoods": details,
        "fits": [
            {
                key: value
                for key, value in result.items()
                if key not in {"candidate", "gene_rows", "profiles"}
            }
            for result in results
        ],
        "profiles": profiles,
        "limitations": [
            "Assumes one disomic allopolyploid event; not a WGD test or DL+ILS model.",
            "Direct second-parent attachment, not a separately estimated ghost-parent divergence/hybridization history.",
            "Local observed-pattern sensitivity rank and nuisance-refitted profile grids do not prove global identifiability or give confidence intervals.",
            "Parameter values are withheld for boundary/locally unidentified fits; numerical solutions are retained only for reproducibility.",
            "All candidates use identical input families, fixed species ages and explicit bounds; no candidate-dependent filtering.",
            "Missing copies are conditioned on; no SSD, gene loss, detection, sequence likelihood or gene-tree/root uncertainty model.",
        ],
    }
    profile_table = pd.DataFrame(
        [
            {
                **row,
                "nuisance_parameters": json.dumps(
                    row["nuisance_parameters"], sort_keys=True
                ),
                "optimizer_attempts": json.dumps(
                    row["optimizer_attempts"], sort_keys=True
                ),
            }
            for row in profiles
        ]
    )
    tables = {
        "outfile": pd.DataFrame(rows),
        "report": pd.DataFrame(details),
        "check_out": pd.DataFrame(details),
        "msc_profile_out": profile_table,
    }
    with output_transaction(paths.values()) as staged:
        for name, path in paths.items():
            if name in tables:
                staged.write_text(
                    path,
                    lambda handle, name=name: tables[name].to_csv(
                        handle, sep="\t", index=False
                    ),
                )
            elif name == "tree_out":
                staged.write_text(
                    path,
                    lambda handle: handle.write(
                        dated_text(best["candidate"].tree, args) + "\n"
                    ),
                )
            else:
                staged.write_text(
                    path,
                    lambda handle: json.dump(
                        metadata, handle, indent=2, allow_nan=False
                    ),
                )
    if args.outfile == "-":
        tables["outfile"].to_csv(sys.stdout, sep="\t", index=False)
    sys.stderr.write(
        f"Conditional fitted MUL MSC: best={best['id']}, status={best['diagnostics']['status']}, log_likelihood={best['log_likelihood']}. Not a WGD test; profiles are not confidence intervals.\n"
    )
