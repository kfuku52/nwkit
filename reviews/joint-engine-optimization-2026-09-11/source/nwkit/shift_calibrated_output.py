"""CLI exports for calibrated contrast selection (schema 7)."""

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from nwkit import __version__
from nwkit.output_transaction import output_transaction
from nwkit.shift_calibration import CalibratedSearch, validate_calibration_options


def calibrated_main(args, tree, ids, tokens, newick, data, outputs):
    validate_calibration_options(args.calibration_replicates, args.calibration_level)
    if args.fit_out or args.bootstrap:
        raise ValueError(
            "--fit-out and --bootstrap are IC backend options; calibrated selection uses --calibration-replicates and exports its fit in JSON."
        )
    if args.search_strategy not in {"auto", "exhaustive"}:
        raise ValueError(
            "Calibrated selection performs exhaustive search; lasso/ensemble require --selection ic."
        )
    variances = (
        data.observation_variance.to_numpy()
        if "observation_variance" in data
        else np.zeros(len(data))
    )
    search = CalibratedSearch(
        tree,
        args.max_shifts,
        args.convergence,
        variances,
        args.exhaustive_max_configurations,
    )
    fit = search.fit(
        data.value.to_numpy(),
        args.seed,
        args.calibration_replicates,
        args.calibration_level,
    )
    selected = fit["model"]["shift_branch_ids"]
    aliases = {
        b: ("baseline" if 0 in group else f"shift_{min(group)}")
        for group in fit["model"]["groups"]
        for b in group
    }
    regimes: dict[Any, str] = {}
    optima: dict[Any, float | None] = {}
    identifiable = fit["optimum_identifiable"]
    effects_by_id = dict(zip(selected, fit["optimum_effects"], strict=False))
    effect_rows = []
    for node in tree.traverse("preorder"):
        branch = ids[node]
        regimes[node] = aliases[branch] if branch in aliases else regimes[node.up]
        optima[node] = (
            (fit["intercept"] if node.is_root else optima[node.up])
            if identifiable
            else None
        )
        if branch in selected:
            if identifiable:
                optima[node] += effects_by_id[branch]
            index = selected.index(branch)
            effect_rows.append(
                {
                    "branch_id": branch,
                    "regime": regimes[node],
                    "parent_regime": regimes[node.up],
                    "mean_effect": fit["mean_effects"][index],
                    "optimum_effect": fit["optimum_effects"][index],
                    "optimum_identifiable": identifiable,
                }
            )
    regime_rows = [
        {
            "regime": regimes[n],
            "branch_id": ids[n],
            "optimum": optima[n],
            "optimum_identifiable": identifiable,
        }
        for n in tree.traverse("preorder")
        if n.is_root or ids[n] in selected
    ]
    tip_rows = []
    for index, node in enumerate(tree.leaves()):
        observed, predicted = float(data.value.iloc[index]), fit["predicted"][index]
        tip_rows.append(
            {
                "leaf_name": node.name,
                "branch_id": ids[node],
                "regime": regimes[node],
                "standard_error": float(np.sqrt(variances[index])),
                "observation_variance": float(variances[index]),
                "observed": observed,
                "predicted": predicted,
                "residual": observed - predicted,
                "optimum": optima[node],
                "optimum_identifiable": identifiable,
            }
        )
    result = (
        pd.DataFrame([{"branch_id": ids[n], "regime": regimes[n]} for n in ids])
        .sort_values("branch_id")
        .to_csv(sep="\t", index=False)
    )
    model = {
        "schema_version": 7,
        "selection": "calibrated",
        "backend": "nwkit",
        "nwkit_version": __version__,
        "trait": args.state_column,
        "root_model": args.root_model,
        "effective_root_model": "root_mean_removed_by_contrasts",
        "time_scale": "original",
        "tree_normalized": False,
        "max_shifts": args.max_shifts,
        "convergence_searched": args.convergence,
        "model_family": {
            "finite_alpha": "OU_with_regime_optima",
            "zero_alpha": "scaled_effect_drift_limit_not_finite_OU_optima",
            "infinite_alpha": "independent_tip_process_with_regime_means",
            "root_mean": "removed_by_fixed_orthonormal_contrasts",
        },
        "shift_branch_ids": selected,
        "parameters": {
            k: fit[k]
            for k in (
                "alpha",
                "alpha_height",
                "alpha_status",
                "alpha_limit_supported",
                "sigma2",
                "process_tip_variance",
                "intercept",
                "contrast_log_likelihood",
            )
        },
        "alpha_profile": fit["alpha_profile"],
        "alpha_profile_note": "1.920729 log-likelihood support cutoff is diagnostic, not a post-selection confidence interval. Null alpha_height with limit=infinity is the exact independent limit.",
        "calibration": {
            "replicates": args.calibration_replicates,
            "level": args.calibration_level,
            "seed": args.seed,
            "tests": fit["tests"],
            "method": "sequential_complete_search_parametric_bootstrap",
            "null_calibration": fit["null_calibration"],
            "null_hypothesis": "no_effective_shift",
            "later_stage_interpretation": "family_selection_not_individual_branch_or_convergence_tests",
            "likelihood_process_scale": "profile_ML",
            "plugin_generation_process_scale": "RSS_over_residual_degrees_of_freedom"
            if not search.known_error
            else "fitted_process_variance_with_fixed_observation_variances",
            "replicate_failure_policy": "fail_entire_run",
            "note": "Without measurement error, the first test covers every fitted alpha grid point; early acceptance reports p_value=1 as an explicit upper bound. Later stages and known-error fits use plug-in calibration. This is not a uniform guarantee over continuous nuisance parameters. Every replicate repeats the candidate and nuisance-grid search.",
        },
        "search": {
            "strategy": "exhaustive",
            "evaluated_candidates": fit["candidate_count"],
            "likelihood": "Gaussian_orthonormal_contrast_density",
            "profile": "finite_grid_including_exact_alpha_limits",
            "continuous_global_optimum_certified": False,
            "process_variance_grid": search.variance_grid.tolist()
            if search.known_error
            else None,
        },
        "standard_error_column": args.standard_error_column,
        "effective_observation_error": "known_independent_variances"
        if search.known_error
        else "none",
        "optimum_convention": "root_mean_equals_baseline_optimum; suppressed if alpha limit is supported",
        "optimum_identifiability_rule": "diagnostic_only: finite_grid_maximum, positive_process_variance, and neither_alpha_limit_in_profile_support_set",
        "shift_effects": effect_rows,
        "regime_parameters": regime_rows,
        "tip_predictions": tip_rows,
        "tip_tokens": {token: node.name for node, token in tokens.items()},
        "analysis_tree_tokens": newick,
        "analysis_input_sha256": hashlib.sha256(
            (newick + "\n" + data.to_csv(sep="\t", index=False)).encode()
        ).hexdigest(),
        "branches": [
            {
                "branch_id": ids[n],
                "parent": -1 if n.is_root else ids[n.up],
                "name": n.name,
                "dist": None if n.is_root else n.dist,
                "regime": regimes[n],
            }
            for n in ids
        ],
        "asr_note": "Regime map only. ASR refits parameters; alpha-limit and selection uncertainty are not propagated. At the Brownian limit, group labels describe scaled-effect constraints, not identified finite OU optima.",
    }
    tables = [
        (
            args.effects_out,
            effect_rows,
            [
                "branch_id",
                "regime",
                "parent_regime",
                "mean_effect",
                "optimum_effect",
                "optimum_identifiable",
            ],
        ),
        (
            args.regime_parameters_out,
            regime_rows,
            ["regime", "branch_id", "optimum", "optimum_identifiable"],
        ),
        (
            args.tip_summary_out,
            tip_rows,
            [
                "leaf_name",
                "branch_id",
                "regime",
                "standard_error",
                "observation_variance",
                "observed",
                "predicted",
                "residual",
                "optimum",
                "optimum_identifiable",
            ],
        ),
    ]
    with output_transaction(outputs) as staged:
        Path(staged[args.model_out]).write_text(
            json.dumps(model, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        for path, rows, columns in tables:
            if path:
                pd.DataFrame(rows, columns=columns).to_csv(
                    staged[path], sep="\t", index=False, na_rep="NA"
                )
        if args.outfile != "-":
            Path(staged[args.outfile]).write_text(result, encoding="utf-8")
    if args.outfile == "-":
        sys.stdout.write(result)
