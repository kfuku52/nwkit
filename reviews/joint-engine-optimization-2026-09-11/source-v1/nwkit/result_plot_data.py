"""Read and validate saved reconciliation/dating results without fitting models."""

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit.clade_index import CladeIndex
from nwkit.draw_layouts import get_rectangular_coordinates
from nwkit.reconcile import _validate_rooted_binary_tree
from nwkit.util import read_tree

RADTE_PLOT_SUFFIXES = {
    "tree": ".dated.nwk",
    "nodes": ".nodes.tsv",
    "species": ".species.tsv",
    "manifest": ".manifest.json",
}
EVENT_TYPES = {"leaf", "speciation", "duplication", "transfer", "unresolved"}


@dataclass
class ResultPlotData:
    gene: object
    species: object
    gene_index: CladeIndex
    species_index: CladeIndex
    events: dict
    species_rows: dict
    manifest: dict
    dated: bool = False


def radte_plot_paths(prefix):
    if not prefix or str(prefix).strip() in {"", "-"}:
        raise ValueError("--radte-prefix requires a filesystem prefix.")
    return {key: str(prefix) + suffix for key, suffix in RADTE_PLOT_SUFFIXES.items()}


def radte_plot_protected_paths(prefix):
    from nwkit.radte import radte_paths

    return radte_paths(str(prefix))


def read_result_table(path):
    # Preserve literal identifiers such as NA and leading zeroes.
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _indexed_rows(table, key, required, label):
    missing = set(required) - set(table.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(sorted(missing))}.")
    if table[key].isna().any() or table[key].astype(str).eq("").any():
        raise ValueError(f"{label} contains an empty {key}.")
    if table[key].duplicated().any():
        raise ValueError(f"{label} contains duplicate {key} values.")
    return table.set_index(key, drop=False).to_dict("index")


def reconciliation_plot_data(gene, species, events):
    _validate_rooted_binary_tree(gene, "--infile")
    _validate_rooted_binary_tree(species, "--species-tree")
    gi, si = CladeIndex(gene), CladeIndex(species)
    rows = _indexed_rows(
        events,
        "gene_clade_id",
        {"gene_clade_id", "parent_gene_clade_id", "event_type", "species_event_id"},
        "Reconciliation table",
    )
    if set(rows) != {gi.clade_id_for_node(n) for n in gene.traverse()}:
        raise ValueError(
            "Reconciliation table must cover exactly the supplied gene tree."
        )
    species_ids = {si.clade_id_for_node(n) for n in species.traverse()}
    for node in gene.traverse():
        row = rows[gi.clade_id_for_node(node)]
        parent = "" if node is gene else gi.clade_id_for_node(node.up)
        if row["parent_gene_clade_id"] != parent:
            raise ValueError("Reconciliation parent clades do not match gene topology.")
        event = str(row["event_type"])
        if event not in EVENT_TYPES or ((event == "leaf") != node.is_leaf):
            raise ValueError(f"Invalid event type for gene clade: {event!r}.")
        sid = row["species_event_id"]
        if sid and sid not in species_ids:
            raise ValueError(
                "Reconciliation mapping does not match the supplied species tree."
            )
    return ResultPlotData(gene, species, gi, si, rows, {}, {})


def _number(value, label, *, optional=False):
    if optional and (pd.isna(value) or str(value) in {"", "NA", "NaN", "nan"}):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric.") from exc
    if not np.isfinite(result):
        raise ValueError(f"{label} must be finite.")
    return result


def _validate_age_rows(rows, label, *, allow_unestimated=False):
    for row in rows.values():
        for key in ("estimated_age", "age_min", "age_max"):
            row[key] = _number(
                row.get(key),
                f"{label} {key}",
                optional=allow_unestimated
                and key == "estimated_age"
                and row.get("estimation_status") == "not-sampled",
            )
        if (
            row["estimated_age"] is not None and row["estimated_age"] < 0
        ) or not 0 <= row["age_min"] <= row["age_max"]:
            raise ValueError(
                f"{label} ages and bounds must be nonnegative and ordered."
            )
        if label in {"Node", "Species"}:
            for key in ("interval_lower", "interval_upper"):
                row[key] = _number(row.get(key), key, optional=True)
            lower, upper = row["interval_lower"], row["interval_upper"]
            if (lower is None) != (upper is None):
                raise ValueError("Age intervals require both endpoints or neither.")
            if lower is not None and not 0 <= lower <= upper:
                raise ValueError("Age intervals must be nonnegative and ordered.")


def _validate_species_evidence(rows):
    for row in rows.values():
        for key in (
            "input_interval_lower",
            "input_interval_upper",
            "input_interval_level",
        ):
            row[key] = _number(row.get(key), key, optional=True)
        lo, hi, level = (
            row[k]
            for k in (
                "input_interval_lower",
                "input_interval_upper",
                "input_interval_level",
            )
        )
        if lo is None and hi is None and level is None:
            continue
        if (
            lo is None
            or hi is None
            or level is None
            or not 0 <= lo <= hi
            or not 0 < level < 1
        ):
            raise ValueError(
                "Species input intervals require ordered nonnegative endpoints and a probability level."
            )
        if row.get("input_interval_kind") not in {
            "confidence",
            "credible",
            "hpd",
            "percentile",
        } or not row.get("input_interval_source"):
            raise ValueError("Species input intervals require their kind and source.")


def dating_plot_data(gene, species, nodes, species_table, manifest):
    data = reconciliation_plot_data(gene, species, nodes)
    required = {
        "estimated_age",
        "age_min",
        "age_max",
        "shared_age_id",
        "interval_lower",
        "interval_upper",
        "interval_status",
    }
    if not required.issubset(nodes.columns):
        raise ValueError(
            "RADTE node table is missing age, interval, or shared-event columns."
        )
    species_rows = _indexed_rows(
        species_table,
        "species_event_id",
        {"species_event_id", "estimated_age", "age_min", "age_max", "age"},
        "Species table",
    )
    if set(species_rows) != {
        data.species_index.clade_id_for_node(n) for n in species.traverse()
    }:
        raise ValueError(
            "Species result table must cover exactly the supplied species tree."
        )
    _validate_age_rows(data.events, "Node")
    _validate_age_rows(
        species_rows,
        "Species",
        allow_unestimated=str(manifest.get("calibration_policy", "")).startswith(
            "PAML"
        ),
    )
    _validate_species_evidence(species_rows)
    scale = max(1.0, *(r["estimated_age"] for r in data.events.values()))
    tolerance = 1e-7 * scale
    if not str(manifest.get("calibration_policy", "hard")).startswith("PAML"):
        for row in [*data.events.values(), *species_rows.values()]:
            if (
                not row["age_min"] - tolerance
                <= row["estimated_age"]
                <= row["age_max"] + tolerance
            ):
                raise ValueError(
                    "Estimated age lies outside a saved hard calibration bound."
                )
    for node in gene.traverse():
        row = data.events[data.gene_index.clade_id_for_node(node)]
        if node.is_leaf and abs(row["estimated_age"]) > tolerance:
            raise ValueError("Dated gene tips must have age zero.")
        if node is not gene:
            parent = data.events[data.gene_index.clade_id_for_node(node.up)]
            length = _number(node.dist, "Dated branch length")
            if (
                length < 0
                or abs(parent["estimated_age"] - row["estimated_age"] - length)
                > tolerance
            ):
                raise ValueError(
                    "Dated tree branch lengths disagree with saved node ages."
                )
        if row["event_type"] == "speciation":
            sid = row["species_event_id"]
            if (
                not sid
                or species_rows[sid]["estimated_age"] is None
                or abs(row["estimated_age"] - species_rows[sid]["estimated_age"])
                > tolerance
            ):
                raise ValueError(
                    "Shared speciation ages disagree with the saved species ages."
                )
            if row.get("shared_age_id") != "S:" + sid:
                raise ValueError(
                    "Shared speciation group does not match its species event."
                )
    for node in species.traverse():
        if node is not species and _number(node.dist, "Species branch length") <= 0:
            raise ValueError("Species tree requires positive branch lengths.")
    depth, _, leaves = get_rectangular_coordinates(species)
    height = max(depth[n] for n in leaves)
    for node in species.traverse():
        row = species_rows[data.species_index.clade_id_for_node(node)]
        original = _number(row["age"], "Original species age")
        if abs(height - depth[node] - original) > tolerance:
            raise ValueError(
                "Species tree ages disagree with the saved calibration input."
            )
        if (
            node.is_leaf
            and row["estimated_age"] is not None
            and abs(row["estimated_age"]) > tolerance
        ):
            raise ValueError("Estimated species tip ages must be zero.")
        if node is not species:
            parent = species_rows[data.species_index.clade_id_for_node(node.up)]
            if (
                row["estimated_age"] is not None
                and parent["estimated_age"] is not None
                and row["estimated_age"] > parent["estimated_age"] + tolerance
            ):
                raise ValueError("Estimated species ages violate parent-child order.")
    data.species_rows = species_rows
    data.manifest = manifest
    data.dated = True
    return data


def load_dating_plot_data(
    prefix, species_path, *, species_rooted="auto", species_format="auto"
):
    paths = radte_plot_paths(prefix)
    manifest_bytes = Path(paths["manifest"]).read_bytes()
    manifest = json.loads(manifest_bytes)
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema") != "nwkit-radte-run-v1"
        or manifest.get("status") != "complete"
    ):
        raise ValueError("RADTE plotting requires a complete supported run manifest.")
    hashes = manifest.get("output_sha256", {})
    contents = {}
    for key, path in paths.items():
        if key == "manifest":
            continue
        contents[key] = Path(path).read_bytes()
        if hashlib.sha256(contents[key]).hexdigest() != hashes.get(key):
            raise ValueError(
                f"RADTE {key} output does not match its manifest; do not mix or edit result bundles."
            )
    if Path(paths["manifest"]).read_bytes() != manifest_bytes:
        raise ValueError(
            "RADTE results changed during reading; retry after the run finishes."
        )
    from io import StringIO

    gene = read_tree(contents["tree"].decode(), "auto", True)
    species = read_tree(species_path, species_format, True, rooted=species_rooted)
    frames = {
        key: read_result_table(StringIO(contents[key].decode()))
        for key in ("nodes", "species")
    }
    return dating_plot_data(gene, species, frames["nodes"], frames["species"], manifest)
