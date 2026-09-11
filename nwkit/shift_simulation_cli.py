"""Reproducible simulation from explicit parameters or a native shift model."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.rooting_state import require_rooted
from nwkit.shift_joint_model import integral_decay, original_diffusion_coordinate
from nwkit.shift_native_model import ShiftLayout, ShiftTree
from nwkit.shift_simulation import ShiftSimulation, simulate_shift
from nwkit.util import read_tree, validate_unique_named_leaves


def _alpha(value, p):
    if isinstance(value, (str, int, float)) and not isinstance(value, bool):
        value = [value]
    if not isinstance(value, list) or any(isinstance(v, bool) for v in value):
        raise ValueError(
            "Alpha requires a number or numeric list; use the string 'inf' for the independent limit."
        )
    try:
        alpha = np.broadcast_to(np.asarray([float(x) for x in value]), (p,)).copy()
    except (ValueError, TypeError) as exc:
        raise ValueError("Alpha must have one value or one per trait.") from exc
    if np.isnan(alpha).any() or np.any(alpha < 0):
        raise ValueError("Alpha must be nonnegative.")
    return alpha


def _covariance_coordinate(parameters, tree, alpha, p, root):
    supplied = [
        key
        for key in ("diffusion_covariance", "process_tip_covariance")
        if parameters.get(key) is not None
    ]
    if len(supplied) != 1:
        raise ValueError(
            "Supply exactly one of diffusion_covariance and process_tip_covariance."
        )
    covariance = np.asarray(parameters[supplied[0]], float)
    if covariance.shape != (p, p) or not np.isfinite(covariance).all():
        raise ValueError("Generating covariance must be a finite trait-square matrix.")
    if np.isinf(alpha).all():
        if supplied[0] != "process_tip_covariance":
            raise ValueError(
                "The independent limit requires process_tip_covariance, not diffusion."
            )
        return covariance
    if np.isinf(alpha).any():
        if supplied[0] == "process_tip_covariance" and not np.any(
            covariance - np.diag(np.diag(covariance))
        ):
            return covariance
        raise ValueError(
            "Simulation does not support mixed finite/infinite alpha with joint covariance."
        )
    if root == "OUrandomRoot" and (alpha == 0).any():
        raise ValueError("Stationary-root OU is undefined at alpha zero.")
    if supplied[0] == "process_tip_covariance":
        sums = alpha[:, None] + alpha[None]
        normalized_diffusion = (
            covariance * sums
            if root == "OUrandomRoot"
            else covariance / integral_decay(sums, 1.0)
        )
        covariance = normalized_diffusion / tree.height
    return original_diffusion_coordinate(tree, alpha, covariance, root)


def explicit_simulation(tree, parameters):
    allowed = {
        "schema",
        "trait_names",
        "alpha",
        "root_model",
        "shift_branch_ids",
        "groups",
        "regime_optima",
        "scaled_regime_coefficients",
        "diffusion_covariance",
        "process_tip_covariance",
        "sampling_standard_errors",
        "measurement_covariance",
        "missing",
    }
    unknown = set(parameters) - allowed
    if unknown:
        raise ValueError("Unknown generating parameters: " + ", ".join(sorted(unknown)))
    names = tuple(parameters["trait_names"])
    p, n = len(names), len(tree.leaf_names)
    original_alpha = _alpha(parameters["alpha"], p)
    with np.errstate(over="ignore", under="ignore"):
        alpha = original_alpha * tree.height
    if np.any(np.isfinite(original_alpha) & ~np.isfinite(alpha)) or np.any(
        (original_alpha > 0) & (alpha == 0)
    ):
        raise ValueError(
            "Alpha times tree height is not representable; rescale time units."
        )
    root = parameters.get("root_model", "OUfixedRoot")
    layout = ShiftLayout.build(
        tree, parameters.get("shift_branch_ids", ()), parameters.get("groups")
    )
    mean_keys = [
        k for k in ("regime_optima", "scaled_regime_coefficients") if k in parameters
    ]
    if len(mean_keys) != 1:
        raise ValueError(
            "Supply regime_optima or scaled_regime_coefficients, but not both."
        )
    coefficients = np.asarray(parameters[mean_keys[0]], float).copy()
    if coefficients.shape != (len(layout.groups), p):
        raise ValueError(
            "Generating means must have one row per canonical regime group and one column per trait."
        )
    if mean_keys[0] == "regime_optima":
        if (alpha == 0).any():
            raise ValueError(
                "The Brownian drift limit requires scaled_regime_coefficients."
            )
        coefficients[1:] = (coefficients[1:] - coefficients[0]) * -np.expm1(
            -alpha[None]
        )
    errors = np.broadcast_to(
        np.asarray(parameters.get("sampling_standard_errors", 0.0), float), (n, p)
    ).copy()
    if not np.isfinite(errors).all() or (errors < 0).any():
        raise ValueError("Sampling standard errors must be finite and nonnegative.")
    with np.errstate(over="ignore", under="ignore"):
        variances = errors**2
    if not np.isfinite(variances).all() or np.any((errors > 0) & (variances == 0)):
        raise ValueError("Squared standard errors are not representable.")
    missing = np.asarray(parameters.get("missing", np.zeros((n, p), bool)))
    if missing.dtype != bool:
        raise ValueError("Missing mask must contain JSON booleans.")
    spec = ShiftSimulation(
        tree,
        layout,
        names,
        alpha,
        _covariance_coordinate(parameters, tree, alpha, p, root),
        coefficients,
        root,
        variances,
        np.asarray(parameters.get("measurement_covariance", np.zeros((p, p))), float),
        missing,
    )
    spec.validate()
    return spec


def model_simulation(tree, model):
    if (
        model.get("selection") != "native"
        or model.get("completion_status") != "complete"
    ):
        raise ValueError("Simulation requires a completed native model.")
    names = tuple(model["trait_names"])
    p, n = len(names), len(tree.leaf_names)
    by_branch = {int(row["branch_id"]): row for row in model["branches"]}
    if set(by_branch) != set(tree.branch_ids):
        raise ValueError("Model and simulation tree branch IDs differ.")
    for i, b in enumerate(tree.branch_ids):
        row = by_branch[b]
        expected_parent = -1 if i == 0 else tree.branch_ids[tree.compiled.parents[i]]
        if row["parent"] != expected_parent or str(row["name"]) != str(
            tree.compiled.nodes[i].name or ""
        ):
            raise ValueError("Model and simulation tree topology/names differ.")
        if i and not np.isclose(
            float(row["dist"]), float(tree.compiled.nodes[i].dist), rtol=1e-12, atol=0
        ):
            raise ValueError("Model and simulation tree branch lengths differ.")
    traits = model["traits"]
    if [r["trait"] for r in traits] != list(names):
        raise ValueError("Model trait records are not in declared trait order.")
    alpha = []
    for record in traits:
        status = record["alpha_status"]
        if status == "independent_limit":
            alpha.append("inf")
        elif record.get("alpha") is not None:
            alpha.append(record["alpha"])
        else:
            raise ValueError(
                "Model does not export identified alpha; supply explicit generating parameters."
            )
    parameters = {
        "trait_names": list(names),
        "alpha": alpha,
        "root_model": model["root_model"],
        "shift_branch_ids": model["shift_branch_ids"],
        "groups": model["groups"],
        "scaled_regime_coefficients": np.column_stack(
            [r["scaled_regime_coefficients"] for r in traits]
        ).tolist(),
    }
    if "joint_covariance" in model:
        joint = model["joint_covariance"]
        key = (
            "process_tip_covariance"
            if joint["diffusion_covariance"] is None
            else "diffusion_covariance"
        )
        parameters[key] = joint[key]
        parameters["measurement_covariance"] = joint["measurement_covariance"]
    else:
        if any(
            r.get("process_tip_variance") is None
            or r.get("measurement_variance") is None
            for r in traits
        ):
            raise ValueError(
                "Model does not export a variance decomposition; supply explicit generating parameters."
            )
        parameters["process_tip_covariance"] = np.diag(
            [r["process_tip_variance"] for r in traits]
        ).tolist()
        parameters["measurement_covariance"] = np.diag(
            [r["measurement_variance"] for r in traits]
        ).tolist()
    errors = np.zeros((n, p))
    missing = np.zeros((n, p), bool)
    observed = set()
    tip_index = {name: i for i, name in enumerate(tree.leaf_names)}
    trait_index = {name: j for j, name in enumerate(names)}
    for row in model["tip_predictions"]:
        tip_key = (row["leaf_name"], row["trait"])
        if tip_key in observed:
            raise ValueError("Duplicate model tip/trait record.")
        observed.add(tip_key)
        i, j = tip_index[tip_key[0]], trait_index[tip_key[1]]
        errors[i, j] = row["standard_error"]
        missing[i, j] = row["observed"] is None
    if len(observed) != n * p:
        raise ValueError("Model is missing tip/trait records.")
    parameters["sampling_standard_errors"] = errors.tolist()
    parameters["missing"] = missing.tolist()
    return explicit_simulation(tree, parameters)


def shift_simulate_main(args):
    paths = [args.truth_out] + ([args.latent_out] if args.latent_out else [])
    if "-" in paths:
        raise ValueError("Simulation truth/latent outputs require file paths.")
    if args.outfile != "-":
        paths.append(args.outfile)
    validate_output_targets(paths)
    source = args.parameters or args.model_in
    validate_outputs_do_not_replace_inputs(
        [("tree", args.infile), ("generating model", source)], [(p, p) for p in paths]
    )
    if args.replicates < 1 or args.seed < 0 or args.memory_mb < 1:
        raise ValueError("Replicates/memory must be positive and seed nonnegative.")
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    require_rooted(tree, "Shift simulation requires a rooted tree.")
    validate_unique_named_leaves(tree, "--infile")
    prepared = ShiftTree.build(tree)
    payload = json.loads(Path(source).read_text(encoding="utf-8"))
    spec = (
        explicit_simulation(prepared, payload)
        if args.parameters
        else model_simulation(prepared, payload)
    )
    n, p = len(prepared.leaf_names), len(spec.trait_names)
    estimated = args.replicates * (len(prepared.branch_ids) + n) * p * 8 * 4
    if estimated > args.memory_mb * 1024**2:
        raise ValueError(
            f"Simulation allocation estimate {estimated} bytes exceeds --memory-mb."
        )
    reserved = {"replicate", "leaf_name", "branch_id", "node_name"}
    error_columns = [f"se_{name}" for name in spec.trait_names]
    if reserved.intersection(spec.trait_names) or set(error_columns).intersection(
        spec.trait_names
    ):
        raise ValueError(
            "Simulation trait names collide with output metadata or SE columns."
        )
    values, latent = simulate_shift(spec, args.replicates, seed=args.seed)
    table = pd.DataFrame(values.reshape(-1, p), columns=spec.trait_names)
    table.insert(0, "leaf_name", np.tile(prepared.leaf_names, args.replicates))
    table.insert(0, "replicate", np.repeat(np.arange(args.replicates), n))
    for j, name in enumerate(error_columns):
        table[name] = np.tile(np.sqrt(spec.sampling_variances[:, j]), args.replicates)
    labels = spec.layout.node_groups(prepared)
    truth = {
        "schema": "nwkit-shift-simulation-v1",
        "seed": args.seed,
        "replicates": args.replicates,
        "conditioning": "unconditional_generating_process",
        "trait_names": list(spec.trait_names),
        "tip_order": list(prepared.leaf_names),
        "root_model": spec.root_model,
        "alpha": [
            float(a / prepared.height) if np.isfinite(a) else "inf"
            for a in spec.alpha_height
        ],
        "covariance_coordinate_original_units": spec.covariance_coordinate.tolist(),
        "covariance_coordinate_definition": "variance-scaled diffusion; diagonal equals marginal tip variance",
        "scaled_regime_coefficients": spec.coefficients.tolist(),
        "shift_branch_ids": list(spec.layout.shifts),
        "groups": [list(g) for g in spec.layout.groups],
        "sampling_variances": spec.sampling_variances.tolist(),
        "measurement_covariance": spec.measurement_covariance.tolist(),
        "missing": spec.missing.tolist(),
        "branches": [
            {
                "branch_id": b,
                "parent": -1
                if i == 0
                else prepared.branch_ids[prepared.compiled.parents[i]],
                "name": str(prepared.compiled.nodes[i].name or ""),
                "dist": None if i == 0 else float(prepared.compiled.nodes[i].dist),
                "regime_group": int(labels[i]),
            }
            for i, b in enumerate(prepared.branch_ids)
        ],
    }
    from nwkit.shift_joint_model import joint_covariance_geometry

    _, _, _, tip, diffusion = joint_covariance_geometry(
        prepared, spec.alpha_height, spec.covariance_coordinate, spec.root_model
    )
    truth["process_tip_covariance"] = tip.tolist()
    truth["diffusion_covariance"] = (
        None if diffusion is None else (diffusion / prepared.height).tolist()
    )
    with output_transaction(paths) as staged:
        Path(staged[args.truth_out]).write_text(
            json.dumps(truth, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        if args.outfile != "-":
            table.to_csv(staged[args.outfile], sep="\t", index=False, na_rep="NA")
        if args.latent_out:
            nodes = len(prepared.branch_ids)
            latent_table = pd.DataFrame(latent.reshape(-1, p), columns=spec.trait_names)
            latent_table.insert(
                0,
                "node_name",
                np.tile(
                    [str(node.name or "") for node in prepared.compiled.nodes],
                    args.replicates,
                ),
            )
            latent_table.insert(
                0, "branch_id", np.tile(prepared.branch_ids, args.replicates)
            )
            latent_table.insert(
                0, "replicate", np.repeat(np.arange(args.replicates), nodes)
            )
            latent_table.to_csv(staged[args.latent_out], sep="\t", index=False)
    if args.outfile == "-":
        print(table.to_csv(sep="\t", index=False, na_rep="NA"), end="")
