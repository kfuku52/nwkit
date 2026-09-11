"""Validate shift inputs and map external fit clades to NWKIT branch IDs."""

import hashlib
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

from nwkit import __version__
from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.rooting_state import require_rooted
from nwkit.shift_backend import run_backend
from nwkit.shift_bootstrap import collect_bootstrap
from nwkit.shift_convergence import collect_convergence, shared_regime_rows
from nwkit.shift_math import ancestral_effect_scales, observation_variances
from nwkit.shift_results import (
    align_regime_optima,
    collect_effects,
    collect_search,
    collect_tips,
    summarize_effects,
)
from nwkit.util import (
    assign_branch_ids,
    read_tip_table,
    read_tree,
    validate_unique_named_leaves,
)


def _validate_options(args):
    if getattr(args, "max_shifts", None) == "auto":
        raise ValueError("--max-shifts auto requires --selection native.")
    if getattr(args, "global_null_gate", False) and (
        args.selection != "native" or args.criterion != "AIC" or args.regime_map
    ):
        raise ValueError(
            "--global-null-gate requires native AIC search without --regime-map."
        )
    if args.search_strategy == "native-path" and args.selection != "native":
        raise ValueError("--search-strategy native-path requires --selection native.")
    if args.selection == "ic" and args.criterion is None:
        args.criterion = "pBIC"
    if args.selection == "calibrated" and args.criterion is not None:
        raise ValueError(
            "--criterion requires --selection ic or native; calibrated selection uses bootstrap tests."
        )
    if (
        args.selection == "ic"
        and args.convergence
        and args.criterion not in {"AIC", "AICc", "BIC", "pBIC"}
    ):
        raise ValueError(
            "--convergence supports only AIC, AICc, BIC and pBIC criteria."
        )
    if not 0 <= args.bootstrap <= 2147483647:
        raise ValueError("--bootstrap must be between 0 and 2147483647.")
    if not 0 <= args.bootstrap_seed <= 2147483647:
        raise ValueError("--bootstrap-seed must be between 0 and 2147483647.")
    if args.max_shifts < 0 or not 1 <= args.exhaustive_max_configurations <= 2147483647:
        raise ValueError(
            "--max-shifts must be non-negative and --exhaustive-max-configurations between 1 and 2147483647."
        )
    if not 0 <= args.seed <= 2147483647:
        raise ValueError("--seed must be between 0 and 2147483647.")


def _prepare(args):
    _validate_options(args)
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    require_rooted(tree, "Shift inference requires a rooted tree.")
    validate_unique_named_leaves(tree, "--infile")
    ids = assign_branch_ids(tree)
    distances = {tree: 0.0}
    for node in tree.traverse("preorder"):
        if not node.is_leaf and len(node.children) != 2:
            raise ValueError(
                "Initial shift inference requires a strictly bifurcating tree."
            )
        if not node.is_root:
            if node.dist is None or not math.isfinite(node.dist) or node.dist <= 0:
                raise ValueError(
                    "Shift inference requires finite, positive non-root branch lengths."
                )
            distances[node] = distances[node.up] + node.dist
    leaves = list(tree.leaves())
    heights = [distances[node] for node in leaves]
    if len(leaves) < 3 or not all(math.isfinite(x) for x in heights):
        raise ValueError(
            "Shift inference requires at least three tips and finite tree height."
        )
    if max(heights) - min(heights) > max(heights) * 1e-8:
        raise ValueError(
            "Shift inference requires an ultrametric tree; no automatic repair is performed."
        )
    if args.max_shifts >= len(leaves) - 1:
        raise ValueError(
            "--max-shifts must be smaller than the number of tips minus one."
        )
    if args.state_column == "leaf_name":
        raise ValueError("--state-column must differ from leaf_name.")
    if args.standard_error_column == "":
        raise ValueError("--standard-error-column requires a nonempty column name.")
    if args.standard_error_column in {"leaf_name", args.state_column}:
        raise ValueError(
            "--standard-error-column must differ from leaf_name and --state-column."
        )
    table, _, _ = read_tip_table(
        args.trait,
        tree_leaf_names=[n.name for n in leaves],
        required_columns=[args.state_column]
        + ([args.standard_error_column] if args.standard_error_column else []),
        unmatched="error",
    )
    values = pd.to_numeric(
        table.set_index("leaf_name")[args.state_column], errors="raise"
    )
    if values.isna().any() or not all(math.isfinite(float(x)) for x in values):
        raise ValueError(
            "Initial shift inference requires a finite observation for every tip."
        )
    if values.nunique() < 2:
        raise ValueError("Shift inference requires a non-invariant trait.")
    tokens = {node: f"t{i}" for i, node in enumerate(leaves)}
    clades: dict[Any, list[str]] = {}
    newicks: dict[Any, str] = {}
    for node in tree.traverse("postorder"):
        clades[node] = (
            [tokens[node]]
            if node.is_leaf
            else sorted(t for child in node.children for t in clades[child])
        )
        label = (
            tokens[node]
            if node.is_leaf
            else "(" + ",".join(newicks[c] for c in node.children) + ")"
        )
        newicks[node] = label + (";" if node.is_root else f":{node.dist:.17g}")
    data = pd.DataFrame(
        {
            "leaf_name": [tokens[n] for n in leaves],
            "value": [float(values[n.name]) for n in leaves],
        }
    )
    if args.standard_error_column:
        errors = pd.to_numeric(
            table.set_index("leaf_name")[args.standard_error_column], errors="raise"
        )
        if errors.isna().any():
            raise ValueError("A finite standard error is required for every tip.")
        data["standard_error"] = [float(errors[n.name]) for n in leaves]
        data["observation_variance"] = observation_variances(data.standard_error)
    mapping = {"/".join(clades[node]): ids[node] for node in ids}
    return tree, ids, tokens, mapping, newicks[tree], data


def _configuration(value, mapping):
    keys = value.split(";") if value else []
    if len(keys) != len(set(keys)) or any(
        key not in mapping or mapping[key] == 0 for key in keys
    ):
        raise ValueError("kfl1ou returned an invalid shift clade configuration.")
    return sorted(mapping[key] for key in keys)


def _metrics(path):
    model = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if len(model) != 1:
        raise ValueError("Expected one kfl1ou model.")
    metrics = {
        key: float(model.iloc[0][key])
        for key in ("alpha", "sigma2", "intercept", "log_likelihood", "score")
    }
    if (
        not all(math.isfinite(x) for x in metrics.values())
        or metrics["alpha"] < 0
        or metrics["sigma2"] <= 0
    ):
        raise ValueError("kfl1ou returned invalid model parameters or scores.")
    return metrics, str(model.iloc[0]["backend_version"])


def _collect(directory, mapping, convergence=False):
    shifts = pd.read_csv(directory / "shifts.tsv", sep="\t", keep_default_na=False)
    if list(shifts.columns) != ["clade", "mean_effect", "optimum_effect"]:
        raise ValueError("Invalid kfl1ou shifts table.")
    selected = _configuration(";".join(shifts.clade), mapping)
    metrics, version = _metrics(directory / "model.tsv")
    search_metrics = (
        _metrics(directory / "unconstrained-model.tsv")[0] if convergence else metrics
    )
    candidates = []
    for row in pd.read_csv(
        directory / "candidates.tsv", sep="\t", keep_default_na=False
    ).itertuples():
        score = float(row.score)
        candidates.append(
            {
                "score": score if math.isfinite(score) else None,
                "shift_branch_ids": _configuration(row.clades, mapping),
            }
        )
    if not any(
        c["shift_branch_ids"] == selected
        and c["score"] is not None
        and math.isclose(
            c["score"], search_metrics["score"], rel_tol=1e-10, abs_tol=1e-10
        )
        for c in candidates
    ):
        raise ValueError("Selected kfl1ou model is missing from its candidate profile.")
    return selected, metrics, version, candidates, search_metrics


def shift_main(args):
    tables_out = {
        "effects": args.effects_out,
        "regimes": args.regime_parameters_out,
        "tips": args.tip_summary_out,
    }
    outputs = (
        [args.model_out]
        + ([args.fit_out] if args.fit_out else [])
        + [p for p in tables_out.values() if p]
    )
    if any(path == "-" for path in outputs):
        raise ValueError("Auxiliary shift outputs require file paths.")
    if args.outfile != "-":
        outputs.append(args.outfile)
    targets = validate_output_targets(outputs)
    if any(not Path(target).parent.is_dir() for target in targets.values()):
        raise ValueError("Output parent directories must exist before inference.")
    validate_outputs_do_not_replace_inputs(
        [
            ("--infile", args.infile),
            ("--trait", args.trait),
            ("--rscript", shutil.which(args.rscript)),
        ],
        [(path, path) for path in outputs],
    )
    tree, ids, tokens, mapping, newick, data = _prepare(args)
    if args.selection == "calibrated":
        from nwkit.shift_calibrated_output import calibrated_main

        return calibrated_main(args, tree, ids, tokens, newick, data, outputs)
    with tempfile.TemporaryDirectory(prefix="nwkit-shift-") as temporary:
        directory = Path(temporary)
        (directory / "tree.nwk").write_text(newick, encoding="utf-8")
        trait_text = data.to_csv(sep="\t", index=False)
        (directory / "trait.tsv").write_text(trait_text, encoding="utf-8")
        execution = run_backend(directory, args)
        selected, metrics, version, candidates, search_metrics = _collect(
            directory, mapping, args.convergence
        )
        convergence, aliases = collect_convergence(
            directory, mapping, selected, args.convergence
        )
        if args.convergence and metrics["alpha"] <= 0:
            raise ValueError("Convergence requires identifiable OU optima (alpha > 0).")
        regimes: dict[Any, str] = {}
        for node in tree.traverse("preorder"):
            regimes[node] = (
                "baseline"
                if node.is_root
                else (aliases[ids[node]] if ids[node] in selected else regimes[node.up])
            )
        rows = [{"branch_id": ids[n], "regime": regimes[n]} for n in ids]
        result = (
            pd.DataFrame(rows).sort_values("branch_id").to_csv(sep="\t", index=False)
        )
        if metrics["alpha"] == 0 and args.root_model == "OUrandomRoot":
            raise ValueError(
                "Stationary/random-root OU is undefined at alpha=0; an explicit fixed-root model is required."
            )
        effects = collect_effects(directory, mapping, selected)
        shift_rows, regime_rows, optima = summarize_effects(
            tree, ids, regimes, metrics, effects
        )
        if args.convergence:
            scales = ancestral_effect_scales(
                tree, ids, metrics["intercept"], effects, "optimum_effect"
            )
            regime_rows = shared_regime_rows(
                regime_rows,
                {ids[node]: scale for node, scale in scales.items()},
            )
        tip_rows = collect_tips(
            directory, tokens, ids, regimes, data, metrics, effects, optima
        )
        align_regime_optima(regime_rows, tip_rows)
        search = collect_search(directory, search_metrics, tree)
        tables = {"effects": shift_rows, "regimes": regime_rows, "tips": tip_rows}
        model = {
            "schema_version": 5,
            "selection": "ic",
            "convergence": convergence,
            "unconstrained_parameters": search_metrics if args.convergence else None,
            "search_stage": "unconstrained_shift_discovery",
            "bootstrap": collect_bootstrap(directory, args, tree, ids, mapping),
            "standard_error_column": args.standard_error_column,
            "effective_observation_error": "known_independent_variances"
            if "observation_variance" in data and data.observation_variance.gt(0).any()
            else "none",
            "observation_error": "known_independent_variances"
            if args.standard_error_column
            else "none",
            "optimum_convention": "root_mean_equals_baseline_optimum",
            "effective_root_model": args.root_model,
            "shift_effects": shift_rows,
            "regime_parameters": regime_rows,
            "tip_predictions": tip_rows,
            "search": search,
            "nwkit_version": __version__,
            "backend": "kfl1ou",
            "backend_version": version,
            "trait": args.state_column,
            "root_model": args.root_model,
            "criterion": args.criterion,
            "search_strategy": args.search_strategy,
            "max_shifts": args.max_shifts,
            "exhaustive_max_configurations": args.exhaustive_max_configurations,
            "seed": args.seed,
            "time_scale": "original",
            "tree_normalized": False,
            "convergence_searched": args.convergence,
            "shift_branch_ids": selected,
            "parameters": metrics,
            "candidates": candidates,
            "tip_tokens": {token: node.name for node, token in tokens.items()},
            "analysis_tree_tokens": newick,
            "analysis_input_sha256": hashlib.sha256(
                (newick + "\n" + trait_text).encode()
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
            "search_diagnostics_r": (directory / "diagnostics.txt").read_text(
                encoding="utf-8"
            ),
            "execution": execution,
            "asr_note": "Regime map only. ASR refits a model; match root treatment, observation-error inputs and units before comparing. Uncertainty conditional on selected regimes excludes selection uncertainty.",
        }
        with output_transaction(outputs) as staged:
            Path(staged[args.model_out]).write_text(
                json.dumps(model, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            for role, path in tables_out.items():
                if path:
                    columns = {
                        "effects": [
                            "branch_id",
                            "regime",
                            "parent_regime",
                            "mean_effect",
                            "optimum_effect",
                            "optimum_identifiable",
                        ],
                        "regimes": [
                            "regime",
                            "branch_id",
                            "optimum",
                            "optimum_identifiable",
                        ],
                        "tips": [
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
                    }[role]
                    pd.DataFrame(tables[role], columns=columns).to_csv(
                        staged[path], sep="\t", index=False, na_rep="NA"
                    )
            if args.fit_out:
                shutil.copyfile(directory / "fit.rds", staged[args.fit_out])
            if args.outfile != "-":
                Path(staged[args.outfile]).write_text(result, encoding="utf-8")
        if execution.get("stderr"):
            sys.stderr.write(execution["stderr"])
        bootstrap = model["bootstrap"]
        if bootstrap is not None and bootstrap["failed"]:
            sys.stderr.write(
                f"Warning: {bootstrap['failed']}/{bootstrap['attempted']} bootstrap refits failed; "
                "support frequencies use successful refits only. See model JSON for details.\n"
            )
        if args.outfile == "-":
            sys.stdout.write(result)
