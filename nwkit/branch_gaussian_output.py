"""ASR output publication and prior-only runs for fixed branch processes."""

import io
import json
import sys
from contextlib import redirect_stdout
from copy import copy
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit import __version__
from nwkit.branch_gaussian import build_branch_gaussian_process
from nwkit.branch_gaussian_asr import (
    branch_fit,
    branch_model_table,
    load_branch_assignment,
)
from nwkit.branch_gaussian_options import branch_root
from nwkit.gaussian_inference import GaussianMarginal, simulate_gaussian_process
from nwkit.output_transaction import output_transaction
from nwkit.util import assign_branch_ids, get_node_class

_MAX_SIMULATION_VALUES = 2_000_000

_MODEL_COLUMNS = (
    "branch_id",
    "model",
    "sigma2",
    "alpha",
    "theta",
    "jump_mean",
    "jump_variance",
)


def run_branch_transaction(args, handler):
    from nwkit.asr import _validate_asr_output_paths

    _validate_asr_output_paths(args)
    paths = {
        name: path
        for name, path in vars(args).items()
        if (name == "outfile" or name.endswith("_out")) and path not in (None, "", "-")
    }
    staged_args = copy(args)
    staged_args._branch_output_staged = True
    if getattr(args, "figure_out", None):
        staged_args._branch_figure_format = Path(args.figure_out).suffix.lower()[1:]
    captured = io.StringIO()
    with output_transaction(paths.values()) as staged:
        for name, path in paths.items():
            setattr(staged_args, name, staged[path])
        with redirect_stdout(captured):
            result = handler(staged_args)
    sys.stdout.write(captured.getvalue())
    return result


def _node_rows(tree):
    ids = assign_branch_ids(tree)
    return [
        {
            "branch_id": identifier,
            "parent": -1 if node.is_root else ids[node.up],
            "node_class": get_node_class(node),
            "name": str(node.name or ""),
            "length": None if node.is_root else float(node.dist),
        }
        for node, identifier in ids.items()
    ]


def write_branch_outputs(tree, observed, errors, fit, args, settings):
    rows = fit.branch_assignment.model_rows()
    if getattr(args, "branch_models_out", None):
        pd.DataFrame(rows, columns=_MODEL_COLUMNS).to_csv(
            args.branch_models_out, sep="\t", index=False, float_format="%.17g"
        )
    if not getattr(args, "process_out", None):
        return
    ids = assign_branch_ids(tree)
    metadata = {
        "schema_version": 2,
        "nwkit_version": __version__,
        "command": "asr",
        "model": settings.model,
        "output": settings.output,
        "parameters": "fixed; no parameter estimation or model search"
        if fit.estimation is None
        else "estimated diffusion parameters; fixed assignments, jumps and root",
        "estimation": fit.estimation,
        "root": asdict(fit.process.root),
        "nodes": _node_rows(tree),
        "branch_models": rows,
        "branch_regimes": None
        if fit.branch_assignment.regime_by_branch_id is None
        else [
            {"branch_id": key, "regime": value}
            for key, value in sorted(fit.branch_assignment.regime_by_branch_id.items())
        ],
        "transitions": [
            {"branch_id": ids[node], **asdict(transition)}
            for node, transition in fit.process.transitions.items()
        ],
        "jump_position": "branch_end_after_diffusion",
        "observations": [
            {"leaf_name": name, "value": value, "standard_error": errors.get(name, 0.0)}
            for name, value in observed.items()
        ],
        "state_column": args.state_column,
        "standard_error_column": getattr(args, "standard_error_column", None),
        "ci_level": settings.ci_level,
        "target": getattr(args, "target", "all"),
        "missing_values": getattr(args, "missing_values", None),
        "unmatched": getattr(args, "unmatched", "warn"),
        "replicate_log_constant": getattr(args, "_replicate_log_constant", 0.0),
        "num_observed_positions": fit.num_observed_positions,
        "seed": getattr(args, "seed", None),
        "prior_samples": (getattr(args, "prior_samples", None) or 1000)
        if settings.output == "prior-samples"
        else None,
        "prior_root_value": fit.prior_root_value,
        "figure_simulations": getattr(args, "figure_simulations", 0),
        "figure_simulation_mode": getattr(args, "figure_simulation_mode", None)
        or "unconditional",
        "figure_simulation_steps": getattr(args, "figure_simulation_steps", None)
        or 200,
        "log_likelihood": fit.log_likelihood,
        "num_observed": fit.num_observed,
        "likelihood_rank": fit.num_effective_observations,
        "uncertainty": "Conditional on the input tree, fixed branch assignments and supplied/fitted parameters; parameter and assignment uncertainty is excluded. Prior draws are latent states without observation noise. Flat-root likelihoods use an improper constant-density root and are not comparable to proper-root likelihoods.",
    }
    Path(args.process_out).write_text(
        json.dumps(metadata, ensure_ascii=False, allow_nan=False, indent=2) + "\n",
        encoding="utf-8",
    )


def _prior_marginals(process, root_value):
    root = process.root
    mean = root_value if root.mode == "flat" else root.mean
    variance = 0.0 if root.mode == "flat" else root.variance
    result = {process.tree: GaussianMarginal(mean, variance)}
    for node in process.tree.traverse("preorder"):
        if node.is_root:
            continue
        transition, parent = process.transitions[node], result[node.up]
        mean = transition.slope * parent.mean + transition.intercept
        variance = transition.slope**2 * parent.variance + transition.variance
        if not np.isfinite([mean, variance]).all():
            raise ValueError(
                "Prior node moments exceed floating-point range; rescale trait units."
            )
        result[node] = GaussianMarginal(mean, variance)
    return result


def run_branch_prior(tree, args):
    from nwkit.asr import _write_table
    from nwkit.asr_figure import write_continuous_asr_figure
    from nwkit.asr_input import AsrSettings, asr_trait_columns, effective_asr_args

    if getattr(args, "trait_type", "auto") == "discrete":
        raise ValueError("BRANCH-GAUSSIAN prior sampling requires continuous traits.")
    settings = AsrSettings.from_args(args, "continuous")
    args = effective_asr_args(args, settings)
    args.state_column = asr_trait_columns(args.state_column or "trait", settings.model)[
        0
    ]
    assignment = load_branch_assignment(tree, args)
    process = build_branch_gaussian_process(
        tree, assignment.models_by_branch_id, root=branch_root(args)
    )
    count = getattr(args, "prior_samples", None) or 1000
    if (len(process.transitions) + 1) * count > _MAX_SIMULATION_VALUES:
        raise ValueError(
            "Prior simulation output exceeds the 2 million node-value limit; reduce --prior-samples."
        )
    root_value = getattr(args, "prior_root_value", None)
    fit = branch_fit(
        process, assignment, summary_kind="prior", prior_root_value=root_value
    )
    samples = simulate_gaussian_process(
        process,
        num_samples=count,
        seed=getattr(args, "seed", None),
        root_values=root_value,
    )
    ids = assign_branch_ids(tree)
    records = [
        {
            "branch_id": ids[node],
            "parent": -1 if node.is_root else ids[node.up],
            "node_class": get_node_class(node),
            "name": str(node.name or ""),
        }
        for node in samples.nodes
    ]
    table = pd.concat([pd.DataFrame(records)] * count, ignore_index=True)
    table.insert(
        0, "simulation", np.repeat(np.arange(1, count + 1), len(samples.nodes))
    )
    table["value"] = samples.values.ravel()
    _write_table(table, args.outfile)
    write_branch_outputs(tree, {}, {}, fit, args, settings)
    if getattr(args, "model_out", None):
        _write_table(branch_model_table(fit, args, settings.ci_level), args.model_out)
    if getattr(args, "figure_out", None):
        marginals = _prior_marginals(process, root_value)
        write_continuous_asr_figure(
            tree,
            {},
            {},
            marginals,
            [args.state_column],
            args,
            settings,
            fit,
            fit.display_assignment,
        )
