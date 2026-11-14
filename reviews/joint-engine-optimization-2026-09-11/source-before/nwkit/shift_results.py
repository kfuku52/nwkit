"""Validated, label-restored shift effects, regime optima and tip predictions."""

import math
from typing import Any

import pandas as pd

from nwkit.shift_math import ancestral_effect_scales, close_in_units, remaining_heights


def _number(value, *, nullable=False):
    if nullable and str(value) == "NA":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Nonfinite value in kfl1ou result.")
    return result


def _table(path, columns):
    table = pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str)
    if list(table.columns) != columns:
        raise ValueError(f"Invalid kfl1ou {path.name} columns.")
    return table


def collect_effects(directory, mapping, selected):
    table = _table(directory / "shifts.tsv", ["clade", "mean_effect", "optimum_effect"])
    effects = {}
    for row in table.itertuples():
        branch = mapping.get(row.clade)
        if branch not in selected or branch in effects:
            raise ValueError("Invalid or duplicated effect clade.")
        effects[branch] = {
            "mean_effect": _number(row.mean_effect),
            "optimum_effect": _number(row.optimum_effect, nullable=True),
        }
    if set(effects) != set(selected):
        raise ValueError("Missing selected shift effects.")
    return effects


def summarize_effects(tree, ids, regimes, metrics, effects):
    alpha = metrics["alpha"]
    heights = remaining_heights(tree)
    baseline = metrics["intercept"] if alpha > 0 else None
    optima = {tree: baseline}
    shifts = []
    regime_rows = [
        {
            "regime": "baseline",
            "branch_id": 0,
            "optimum": baseline,
            "optimum_identifiable": baseline is not None,
        }
    ]
    for node in tree.traverse("preorder"):
        if node.is_root:
            continue
        optima[node] = optima[node.up]
        if ids[node] not in effects:
            continue
        effect = effects[ids[node]]
        delta = effect["optimum_effect"]
        expected = (
            None
            if alpha == 0
            else effect["mean_effect"] / -math.expm1(-alpha * heights[node.up])
        )
        if (expected is None) != (delta is None) or (
            expected is not None and not close_in_units(expected, delta)
        ):
            raise ValueError("kfl1ou mean and optimum shift effects are inconsistent.")
        optimum = (
            None
            if delta is None or optima[node.up] is None
            else optima[node.up] + delta
        )
        optima[node] = optimum
        shifts.append(
            {
                "branch_id": ids[node],
                "regime": regimes[node],
                "parent_regime": regimes[node.up],
                **effect,
                "optimum_identifiable": delta is not None,
            }
        )
        regime_rows.append(
            {
                "regime": regimes[node],
                "branch_id": ids[node],
                "optimum": optimum,
                "optimum_identifiable": optimum is not None,
            }
        )
    return sorted(shifts, key=lambda row: row["branch_id"]), regime_rows, optima


def collect_tips(directory, tokens, ids, regimes, data, metrics, effects, optima):
    table = _table(
        directory / "tips.tsv",
        ["token", "observed", "predicted", "residual", "optimum"],
    )
    by_token = {token: node for node, token in tokens.items()}
    if table.token.duplicated().any() or set(table.token) != set(by_token):
        raise ValueError("kfl1ou tip identities differ from the analysis tree.")
    observed = dict(zip(data.leaf_name, data.value, strict=True))
    errors = (
        dict(zip(data.leaf_name, data.standard_error, strict=True))
        if "standard_error" in data
        else {}
    )
    variances = (
        dict(zip(data.leaf_name, data.observation_variance, strict=True))
        if "observation_variance" in data
        else {}
    )
    expected_mean: dict[Any, float] = {}
    tree = next(iter(tokens)).root
    mean_scales = ancestral_effect_scales(
        tree, ids, metrics["intercept"], effects, "mean_effect"
    )
    optimum_scales = ancestral_effect_scales(
        tree, ids, metrics["intercept"], effects, "optimum_effect"
    )
    for node in tree.traverse("preorder"):
        expected_mean[node] = (
            metrics["intercept"]
            if node.is_root
            else expected_mean[node.up]
            + effects.get(ids[node], {}).get("mean_effect", 0.0)
        )
    rows = []
    for row in table.itertuples():
        node = by_token[row.token]
        values = {
            key: _number(getattr(row, key))
            for key in ("observed", "predicted", "residual")
        }
        if not close_in_units(values["observed"], observed[row.token], rtol=1e-12):
            raise ValueError("kfl1ou observations differ from the input trait values.")
        if not close_in_units(
            values["observed"] - values["predicted"],
            values["residual"],
            operands=(values["observed"], values["predicted"]),
        ) or not close_in_units(
            values["predicted"], expected_mean[node], operands=(mean_scales[node],)
        ):
            raise ValueError(
                "kfl1ou predictions, effects and residuals are inconsistent."
            )
        backend_optimum = _number(row.optimum, nullable=True)
        optimum = optima[node]
        # The BM intercept remains defined, but finite OU optima are not.
        if metrics["alpha"] > 0 and (
            (optimum is None) != (backend_optimum is None)
            or (
                optimum is not None
                and not close_in_units(
                    optimum, backend_optimum, operands=(optimum_scales[node],)
                )
            )
        ):
            raise ValueError("kfl1ou tip and regime optima are inconsistent.")
        rows.append(
            {
                "leaf_name": node.name,
                "branch_id": ids[node],
                "regime": regimes[node],
                "standard_error": errors.get(row.token, 0.0),
                "observation_variance": variances.get(row.token, 0.0),
                **values,
                "optimum": backend_optimum if metrics["alpha"] > 0 else None,
                "optimum_identifiable": optimum is not None,
            }
        )
    return sorted(rows, key=lambda row: row["branch_id"])


def align_regime_optima(regime_rows, tip_rows):
    """Keep directly fitted optima after checking reconstructed shift effects."""
    canonical = {
        row["regime"]: row["optimum"] for row in regime_rows if row["branch_id"] == 0
    }
    for row in tip_rows:
        regime, value = row["regime"], row["optimum"]
        if regime in canonical:
            expected = canonical[regime]
            if (expected is None) != (value is None) or (
                value is not None and not close_in_units(value, expected)
            ):
                raise ValueError("kfl1ou returned different optima for one regime.")
        else:
            canonical[regime] = value
    for row in [*regime_rows, *tip_rows]:
        if row["regime"] in canonical:
            row["optimum"] = canonical[row["regime"]]


def collect_search(directory, metrics, tree):
    columns = [
        "strategy",
        "configuration_space_size",
        "evaluated_configurations",
        "coverage",
        "globally_optimal",
        "ensemble_attempted",
        "ensemble_successful",
        "ensemble_failed",
        "alpha_lower",
        "alpha_upper",
    ]
    table = _table(directory / "search.tsv", columns)
    if len(table) != 1:
        raise ValueError("Expected one search summary.")
    raw = table.iloc[0]
    result = {
        key: _number(raw[key], nullable=True)
        for key in columns
        if key not in {"strategy", "globally_optimal"}
    }
    if raw.globally_optimal not in {"TRUE", "FALSE", "NA"}:
        raise ValueError("Invalid global-optimum diagnostic.")
    result.update(
        strategy=raw.strategy,
        globally_optimal={"TRUE": True, "FALSE": False, "NA": None}[
            raw.globally_optimal
        ],
    )
    _validate_search(result)
    height = remaining_heights(tree)[tree]
    for side in ("lower", "upper"):
        bound = result[f"alpha_{side}"]
        result[f"alpha_at_{side}_bound"] = (
            None
            if bound is None
            else math.isclose(
                metrics["alpha"], bound, rel_tol=1e-6, abs_tol=1e-8 / height
            )
        )
    return result


def _validate_search(result):
    if result["strategy"] not in {"none", "lasso", "ensemble", "exhaustive"}:
        raise ValueError("Unknown effective search strategy.")
    for key in (
        "configuration_space_size",
        "evaluated_configurations",
        "ensemble_attempted",
        "ensemble_successful",
        "ensemble_failed",
    ):
        value = result[key]
        if value is not None and (value < 0 or value != int(value)):
            raise ValueError("Search counts must be non-negative integers.")
    coverage = result["coverage"]
    if coverage is not None and not 0 <= coverage <= 1:
        raise ValueError("Search coverage must be between zero and one.")
    attempted, successful, failed = (
        result[k]
        for k in ("ensemble_attempted", "ensemble_successful", "ensemble_failed")
    )
    if (
        all(v is not None for v in (attempted, successful, failed))
        and attempted != successful + failed
    ):
        raise ValueError("Ensemble success/failure counts do not match attempts.")
    low, high = result["alpha_lower"], result["alpha_upper"]
    if any(v is not None and v < 0 for v in (low, high)) or (
        low is not None and high is not None and low > high
    ):
        raise ValueError("Invalid reported alpha bounds.")
