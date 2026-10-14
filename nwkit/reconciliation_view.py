"""Display state and arbitrary file-backed comparisons, with hash verification."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

from nwkit.reconciliation_document import _fields, _string, load_json, read_document


def validate_view(view, document=None):
    _fields(
        view,
        {"schema_version", "kind", "layout"},
        {"leaf_order", "collapsed_clades", "event_colors", "font", "output"},
        "view",
    )
    if (
        view["schema_version"] != "0.1.0"
        or view["kind"] != "reconciliation_view"
        or view["layout"] not in {"standalone", "topology_embedded", "paired"}
    ):
        raise ValueError("Unsupported reconciliation view version/layout.")
    for key in ("leaf_order", "event_colors"):
        if not isinstance(view.get(key, {}), dict):
            raise ValueError(f"View {key} must be an object.")
    if not isinstance(view.get("collapsed_clades", []), list):
        raise ValueError("Collapsed clades must be an array.")
    for tid, order in view.get("leaf_order", {}).items():
        if (
            not isinstance(order, list)
            or any(not isinstance(x, str) or not x for x in order)
            or len(set(order)) != len(order)
        ):
            raise ValueError("View leaf order requires unique leaf IDs.")
        if document:
            if tid not in document.trees or set(order) != {
                nid for nid, n in document.trees[tid].nodes.items() if not n["children"]
            }:
                raise ValueError(
                    "View leaf order must list every leaf ID in a known tree."
                )
    for item in view.get("collapsed_clades", []):
        _fields(item, {"tree", "node"}, set(), "collapsed clade")
        if document and (
            item["tree"] not in document.trees
            or item["node"] not in document.trees[item["tree"]].nodes
        ):
            raise ValueError("Collapsed clade is not in the scientific document.")
    if document and view["layout"] == "topology_embedded":
        upper = {r["upper_tree"] for r in document.record["reconciliations"]}
        if any(item["tree"] in upper for item in view.get("collapsed_clades", [])):
            raise ValueError(
                "Upper context collapse requires a standalone or paired view."
            )
    import re

    for key, color in view.get("event_colors", {}).items():
        _string(key, "event color name")
        if not isinstance(color, str) or not re.fullmatch(r"#[0-9a-fA-F]{6}", color):
            raise ValueError("Event colors require six-digit hex values.")
    if "font" in view:
        _fields(view["font"], {"family", "size"}, set(), "view font")
        _string(view["font"]["family"], "font family")
        size = view["font"]["size"]
        if (
            isinstance(size, bool)
            or not isinstance(size, (int, float))
            or not math.isfinite(size)
            or not 0 < size <= 72
        ):
            raise ValueError("Invalid view font size.")
    if "output" in view:
        _fields(
            view["output"], {"format"}, {"pdf_pages", "pdf_page_size"}, "view output"
        )
        if (
            view["output"]["format"] not in {"pdf", "svg", "png", "html"}
            or view["output"].get("pdf_pages", "auto")
            not in {"auto", "single", "tiled"}
            or view["output"].get("pdf_page_size", "A4") not in {"A4", "A3"}
        ):
            raise ValueError("Invalid view output options.")
    return view


def read_view(path, document=None):
    return validate_view(load_json(Path(path).read_text()), document)


def read_comparison(path):
    path = Path(path)
    record = load_json(path.read_text())
    _fields(
        record,
        {"schema_version", "kind", "mode", "panels", "shared_context", "view"},
        set(),
        "comparison",
    )
    if (
        record["schema_version"] != "0.1.0"
        or record["kind"] != "reconciliation_comparison"
        or record["mode"] not in {"general", "inferred_trees", "fixed_gene_tree"}
    ):
        raise ValueError("Unsupported comparison version/mode.")
    validate_view(record["view"])
    if not isinstance(record["panels"], list) or not record["panels"]:
        raise ValueError("Comparison requires at least one panel.")
    _fields(
        record["shared_context"],
        {"upper_tree_sha256", "gene_tree_sha256", "mapping_sha256"},
        set(),
        "shared context",
    )
    documents, paths, seen = {}, [], set()
    import re

    for value in record["shared_context"].values():
        if value is not None and (
            not isinstance(value, str) or not re.fullmatch("[a-f0-9]{64}", value)
        ):
            raise ValueError("Shared context hashes must be SHA-256 or null.")
    for panel in record["panels"]:
        _fields(
            panel,
            {
                "id",
                "label",
                "document",
                "document_sha256",
                "inference_status",
                "history_status",
                "inference_runtime",
                "history_runtime",
            },
            {"annotations"},
            "panel",
        )
        _string(panel["id"], "panel ID")
        _string(panel["label"], "panel label")
        if panel["id"] in seen:
            raise ValueError("Duplicate comparison panel ID.")
        seen.add(panel["id"])
        if panel["inference_status"] not in {
            "success",
            "partial",
            "failed",
            "timeout",
            "missing",
            "not_applicable",
        } or panel["history_status"] not in {
            "available",
            "partial",
            "failed",
            "timeout",
            "missing",
            "unsupported",
        }:
            raise ValueError("Unknown comparison inference/history status.")
        for kind in ("inference_runtime", "history_runtime"):
            runtime = panel[kind]
            _fields(runtime, {"seconds", "measured_at", "conditions"}, set(), "runtime")
            value = runtime["seconds"]
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError("Runtime must be nonnegative or null.")
            if runtime["measured_at"] is not None:
                _string(runtime["measured_at"], "runtime date")
            if not isinstance(runtime["conditions"], dict):
                raise ValueError("Runtime conditions require an object.")
        available = panel["history_status"] in {"available", "partial"}
        if available and (not panel["document"] or not panel["document_sha256"]):
            raise ValueError("Available history requires a document and hash.")
        if panel["document"] is not None:
            _string(panel["document"], "panel document")
            source = (path.parent / panel["document"]).resolve()
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            if digest != panel["document_sha256"]:
                raise ValueError(f"Comparison document hash mismatch: {panel['id']}.")
            document = read_document(source)
            for key, expected in record["shared_context"].items():
                if (
                    expected is not None
                    and document.record["provenance"].get(key) != expected
                ):
                    raise ValueError(
                        f"Panel does not verify shared {key}: {panel['id']}."
                    )
            documents[panel["id"]] = document
            paths.append(str(source))
    if (
        record["mode"] == "fixed_gene_tree"
        and record["shared_context"]["gene_tree_sha256"] is None
    ):
        raise ValueError(
            "Fixed gene tree comparison requires a verified shared gene tree hash."
        )
    if (
        record["mode"] == "fixed_gene_tree"
        and len({_fixed_gene_signature(doc) for doc in documents.values()}) > 1
    ):
        raise ValueError(
            "Fixed gene tree panels have different actual trees/lengths/units."
        )
    validate_comparison_view(record["view"], documents)
    return record, documents, paths


def view_for_document(view, document):
    """Apply global display selections only where that stable tree/clade exists."""
    result = dict(view)
    result["leaf_order"] = {
        tid: order
        for tid, order in view.get("leaf_order", {}).items()
        if tid in document.trees
    }
    result["collapsed_clades"] = [
        item
        for item in view.get("collapsed_clades", [])
        if item["tree"] in document.trees
        and item["node"] in document.trees[item["tree"]].nodes
    ]
    return result


def validate_comparison_view(view, documents):
    validate_view(view)
    tree_ids = {tid for document in documents.values() for tid in document.trees}
    if set(view.get("leaf_order", {})) - tree_ids:
        raise ValueError("Comparison leaf order references an unknown tree.")
    for item in view.get("collapsed_clades", []):
        if not any(
            item["tree"] in doc.trees and item["node"] in doc.trees[item["tree"]].nodes
            for doc in documents.values()
        ):
            raise ValueError("Comparison collapsed clade is absent from every source.")
    for document in documents.values():
        validate_view(view_for_document(view, document), document)
    return view


def _fixed_gene_signature(document):
    candidates = [
        tree
        for tree in document.trees.values()
        if tree.id == "G" or tree.record.get("role") == "observed"
    ]
    if not candidates and len(document.trees) == 1:
        candidates = list(document.trees.values())
    if len(candidates) != 1:
        raise ValueError("Fixed gene tree panels must identify one G/observed tree.")
    tree = candidates[0]
    signatures: dict[str, tuple[str, ...]] = {}
    nodes = []
    for nid in reversed(tree.order):
        node = tree.nodes[nid]
        clade = (
            tuple(sorted(label for c in node["children"] for label in signatures[c]))
            if node["children"]
            else (node["label"],)
        )
        signatures[nid] = clade
        nodes.append((clade, node["branch_length"]))
    return tree.record["branch_length_unit"], tuple(sorted(nodes))
