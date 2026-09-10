"""Original-unit, long-form native shift outputs and explicit alpha limits."""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit import __version__
from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_input import (
    native_fit_arguments,
    read_native_data,
    read_native_layout,
)
from nwkit.shift_native_provenance import native_provenance, read_native_resume


def _finite_optima_supported(record, fit):
    if (
        not math.isfinite(fit.alpha_height)
        or fit.alpha_height == 0
        or fit.process_variance == 0
        or not record["covariance_components_identifiable"]
    ):
        return False
    diagnostics = record["optimizer"]
    if not record["identifiability_diagnostic"]["finite_alpha_supported"]:
        return False
    if diagnostics["nuisance_variance_at_numerical_bound"]:
        return False
    if not diagnostics["alpha_estimated"]:
        return True
    lower, upper = diagnostics["alpha_height_bounds"]
    if fit.alpha_height <= lower * (1 + 1e-5) or fit.alpha_height >= upper * (1 - 1e-5):
        return False
    limits = [
        r["log_likelihood"]
        for r in diagnostics["alpha_candidates"]
        if r["success"] and r["alpha_mode"] in {"0.0", "inf"}
    ]
    return not limits or max(limits) < fit.log_likelihood - 1.920729410347062


def _remaining_times(tree):
    remaining = np.zeros(len(tree.branch_ids))
    for i in tree.compiled.postorder:
        if i:
            parent = tree.compiled.parents[i]
            remaining[parent] = max(remaining[parent], tree.times[i] + remaining[i])
    return remaining


def native_tables(data, result, regime_names=None):
    tree, layout = data.tree, result["layout"]
    tip_names = tree.leaf_names
    names = regime_names or [
        "baseline" if 0 in group else f"shift_{min(group)}" for group in layout.groups
    ]
    node_groups = layout.node_groups(tree)
    remaining = _remaining_times(tree)
    branches = [
        {
            "branch_id": branch,
            "parent": -1 if i == 0 else tree.branch_ids[tree.compiled.parents[i]],
            "name": str(tree.compiled.nodes[i].name or ""),
            "dist": None if i == 0 else float(tree.compiled.nodes[i].dist),
            "regime": names[node_groups[i]],
        }
        for i, branch in enumerate(tree.branch_ids)
    ]
    effects, regimes, tips = [], [], []
    for trait, (record, fit) in enumerate(
        zip(result["traits"], result["fits"], strict=True)
    ):
        identified = _finite_optima_supported(record, fit)
        record["optimum_identifiable"] = identified
        record["optimum_identifiability_rule"] = (
            "conditional_fixed_layout_diagnostic: local_covariance_rank_support, finite_alpha, positive_process_variance, interior_bounds, and neither_evaluated_alpha_limit_within_1.920729_loglik_units"
        )
        coefficients = np.asarray(record["scaled_regime_coefficients"])
        offsets = coefficients.copy()
        offsets[0] = 0
        optima = (
            coefficients[0] + offsets / -np.expm1(-fit.alpha_height)
            if identified
            else [None] * len(names)
        )
        for group, name in enumerate(names):
            regimes.append(
                {
                    "trait": data.trait_names[trait],
                    "regime": name,
                    "branch_id": min(layout.groups[group]),
                    "optimum": optima[group],
                    "scaled_optimum_offset": float(offsets[group]),
                    "optimum_identifiable": identified,
                }
            )
        for i, branch in enumerate(tree.branch_ids):
            if branch not in layout.shifts:
                continue
            parent = tree.compiled.parents[i]
            group, inherited = node_groups[i], node_groups[parent]
            age = remaining[parent]
            weight = (
                age
                if fit.alpha_height == 0
                else (
                    1.0
                    if math.isinf(fit.alpha_height)
                    else -math.expm1(-fit.alpha_height * age)
                    / -math.expm1(-fit.alpha_height)
                )
            )
            delta = offsets[group] - offsets[inherited]
            effects.append(
                {
                    "trait": data.trait_names[trait],
                    "branch_id": branch,
                    "regime": names[group],
                    "parent_regime": names[inherited],
                    "mean_effect": float(weight * delta),
                    "optimum_effect": float(optima[group] - optima[inherited])
                    if identified
                    else None,
                    "optimum_identifiable": identified,
                }
            )
        for row, i in enumerate(tree.compiled.leaf_indices):
            observed = (
                data.centers[trait] + data.scales[trait] * data.values[row, trait]
            )
            predicted = record["predicted"][row]
            tips.append(
                {
                    "trait": data.trait_names[trait],
                    "leaf_name": tip_names[row],
                    "branch_id": tree.branch_ids[i],
                    "regime": names[node_groups[i]],
                    "observed": float(observed) if math.isfinite(observed) else None,
                    "predicted": predicted,
                    "residual": float(observed - predicted)
                    if math.isfinite(observed)
                    else None,
                    "standard_error": float(
                        data.scales[trait] * np.sqrt(data.variances[row, trait])
                    ),
                    "optimum": optima[node_groups[i]],
                    "optimum_identifiable": identified,
                }
            )
    return branches, effects, regimes, tips


def _model(data, result, tables):
    branches, effects, regimes, tips = tables
    return {
        "schema_version": 8,
        "backend": "nwkit",
        "selection": "native",
        "inference_role": "fixed_layout_parameter_estimation",
        "nwkit_version": __version__,
        "trait_covariance": "diagonal",
        "trait_names": list(data.trait_names),
        "root_model": result["root_model"],
        "likelihood": "ordinary_Gaussian_ML_profiled_means",
        "time_scale": "original",
        "optimum_convention": "root_mean_equals_baseline_optimum; scaled_effect_drift_limit_at_alpha_zero",
        "shift_branch_ids": list(result["layout"].shifts),
        "groups": [list(g) for g in result["layout"].groups],
        "traits": result["traits"],
        "log_likelihood": result["log_likelihood"],
        "branches": branches,
        "shift_effects": effects,
        "regime_parameters": regimes,
        "tip_predictions": tips,
        "selection_calibration": None,
        "asr_note": "The regime map supplies assignments only. ASR refits and does not inherit these parameters, alpha limits, or selection uncertainty.",
    }


def write_native_outputs(args, model):
    tables = [
        (
            args.effects_out,
            model["shift_effects"],
            [
                "trait",
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
            model["regime_parameters"],
            [
                "trait",
                "regime",
                "branch_id",
                "optimum",
                "scaled_optimum_offset",
                "optimum_identifiable",
            ],
        ),
        (
            args.tip_summary_out,
            model["tip_predictions"],
            [
                "trait",
                "leaf_name",
                "branch_id",
                "regime",
                "observed",
                "predicted",
                "residual",
                "standard_error",
                "optimum",
                "optimum_identifiable",
            ],
        ),
    ]
    outputs = [args.model_out, *[path for path, _, _ in tables if path]]
    if args.outfile != "-":
        outputs.append(args.outfile)
    map_text = (
        pd.DataFrame(model["branches"])[["branch_id", "regime"]]
        .sort_values("branch_id")
        .to_csv(sep="\t", index=False)
    )
    text = json.dumps(model, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    with output_transaction(outputs) as staged:
        Path(staged[args.model_out]).write_text(text, encoding="utf-8")
        for path, rows, columns in tables:
            if path:
                pd.DataFrame(rows, columns=columns).to_csv(
                    staged[path], sep="\t", index=False, na_rep="NA"
                )
        if args.outfile != "-":
            Path(staged[args.outfile]).write_text(map_text, encoding="utf-8")
    if args.outfile == "-":
        sys.stdout.write(map_text)


def native_main(args):
    paths = [
        p
        for p in (
            args.model_out,
            args.effects_out,
            args.regime_parameters_out,
            args.tip_summary_out,
        )
        if p
    ]
    if "-" in paths:
        raise ValueError("Auxiliary native shift outputs require file paths.")
    if args.outfile != "-":
        paths.append(args.outfile)
    validate_output_targets(paths)
    validate_outputs_do_not_replace_inputs(
        [
            (name, path)
            for name, path in [
                ("--infile", args.infile),
                ("--trait", args.trait),
                ("--regime-map", args.regime_map),
                ("--resume-model", args.resume_model),
            ]
            if path
        ],
        [(path, path) for path in paths],
    )
    if args.regime_map and (args.bootstrap or args.convergence):
        raise ValueError(
            "A fixed native regime map is supplied; selection/bootstrap options cannot alter it."
        )
    data = read_native_data(args)
    arguments = native_fit_arguments(args, data)
    layout, names = (
        read_native_layout(args.regime_map, data.tree)
        if args.regime_map
        else (None, None)
    )
    provenance = native_provenance(data, args, layout, names)
    if args.resume_model:
        model = read_native_resume(args.resume_model, provenance)
        write_native_outputs(args, model)
        return
    metadata = {}
    if layout is None:
        from nwkit.shift_native_selection import select_native

        result, metadata = select_native(data, args, arguments)
    else:
        result = fit_native_layout(data, layout, **arguments)
    model = _model(data, result, native_tables(data, result, names))
    model.update(metadata)
    model.update(provenance)
    model["completion_status"] = "complete"
    if (
        native_provenance(data, args, layout, names)["implementation_sha256"]
        != provenance["implementation_sha256"]
    ):
        raise ValueError(
            "NWKIT implementation changed during native inference; rerun with a stable installation."
        )
    write_native_outputs(args, model)
