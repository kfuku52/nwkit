import csv
import math
import os
import stat
import sys
from copy import copy
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from typing import Any

from nwkit.clade_mapping import canonical_split
from nwkit.draw import _draw_tree
from nwkit.output_transaction import output_transaction
from nwkit.root import (
    DEFAULT_TAXONOMY_SOURCE_CHAIN,
    SUPPORTED_TAXONOMY_SOURCES,
    _parse_taxonomy_sources,
    mad_rooting,
    midpoint_rooting,
    mv_rooting,
    outgroup_rooting,
    reconciliation_rooting,
    taxonomy_rooting,
    transfer_root_with_taxon_mode,
)
from nwkit.root_evaluation import (
    RootingCandidate,
    RootingEvaluation,
    candidate_from_rooted_tree,
    format_taxa,
    root_split_id,
)
from nwkit.rooting_state import require_rooted
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    get_subtree_leaf_name_sets,
    read_tree,
    remove_singleton,
    validate_distinct_output_paths,
    validate_outputs_do_not_replace_inputs,
    validate_unique_named_leaves,
)

CORE_METHODS = ("midpoint", "mad", "mv")
CONTEXT_METHODS = ("outgroup", "transfer", "reconciliation")
TAXONOMY_METHOD_PREFIX = "taxonomy:"
SUPPORTED_METHODS = CORE_METHODS + CONTEXT_METHODS + ("taxonomy",)
DEFAULT_FIGURE_WIDTH = 7.2
MAXIMUM_AUTO_FIGURE_WIDTH = 200.0
AUTO_WIDTH_INCHES_PER_TIP_FONT_POINT = 0.015
AUTO_TIP_LABEL_LIMIT = 200

ROOT_COMPARE_COLUMNS = (
    "method",
    "status",
    "message",
    "selection_basis",
    "source",
    "tie_index",
    "num_best",
    "is_canonical",
    "root_split_id",
    "num_equivalent_splits",
    "equivalent_root_split_ids",
    "side_a_taxa",
    "side_b_taxa",
    "num_side_a_taxa",
    "num_side_b_taxa",
    "position_kind",
    "position_fraction_from_side_a",
    "edge_length",
    "score_name",
    "score",
    "score_unit",
    "num_evaluated_edges",
    "tie_rule",
    "duplications",
    "losses",
    "duplication_cost",
    "loss_cost",
)

_METHOD_LABELS = {
    "midpoint": "Midpoint",
    "mad": "MAD",
    "mv": "MV",
    "outgroup": "Outgroup",
    "transfer": "Root transfer",
    "reconciliation": "Reconciliation",
    "taxonomy:ncbi": "Taxonomy: NCBI",
    "taxonomy:opentree": "Taxonomy: OpenTree",
    "taxonomy:timetree": "Taxonomy: TimeTree",
}

_MARKER_STYLES = {
    "midpoint": ("#0072B2", "o", 7.0),
    "mad": ("#D55E00", "D", 6.5),
    "mv": ("#009E73", "s", 6.5),
    "outgroup": ("#CC79A7", "^", 7.0),
    "transfer": ("#E69F00", "v", 7.0),
    "reconciliation": ("#56B4E9", "P", 7.5),
    "taxonomy:ncbi": ("#000000", "X", 7.0),
    "taxonomy:opentree": ("#F0E442", "*", 9.5),
    "taxonomy:timetree": ("#777777", "h", 7.0),
}


def _method_tokens(value, option_name, *, allow_empty=False):
    if value is None:
        return []
    if value == "":
        if allow_empty:
            return []
        raise ValueError("'{}' must name at least one method.".format(option_name))
    raw_values = value if isinstance(value, (list, tuple)) else [value]
    tokens = []
    seen = set()
    for raw_value in raw_values:
        for raw_token in str(raw_value).split(","):
            token = raw_token.strip().lower()
            if token == "":
                continue
            if token in seen:
                continue
            seen.add(token)
            tokens.append(token)
    if not tokens and not allow_empty:
        raise ValueError("'{}' must name at least one method.".format(option_name))
    return tokens


def _validate_method_token(token, option_name, allow_all):
    taxonomy_methods = {
        "{}{}".format(TAXONOMY_METHOD_PREFIX, source)
        for source in SUPPORTED_TAXONOMY_SOURCES
    }
    supported = set(SUPPORTED_METHODS) | taxonomy_methods
    if allow_all:
        supported.add("all")
    if token not in supported:
        raise ValueError(
            "Unknown method in '{}': {}. Supported values are: {}.".format(
                option_name,
                token,
                ", ".join(sorted(supported)),
            )
        )


def _taxonomy_names_are_available(tree, args):
    parser = get_species_parser(args=args)
    parsed = [parser.parse(leaf.name) for leaf in tree.leaves()]
    if not parsed or any(
        record.species_label is None or record.taxonomy_query is None
        for record in parsed
    ):
        return False
    species_labels = {record.species_label for record in parsed}
    taxonomy_queries = {record.taxonomy_query for record in parsed}
    return len(species_labels) >= 2 and len(taxonomy_queries) >= 2


def _automatic_taxonomy_methods(tree, args, taxonomy_sources):
    names_are_available = _taxonomy_names_are_available(tree, args)
    taxids_are_available = getattr(args, "taxid_tsv", None) not in (None, "")
    methods: list[str] = []
    for source in taxonomy_sources:
        if source == "ncbi" and (names_are_available or taxids_are_available):
            methods.append("taxonomy:ncbi")
        elif source != "ncbi" and names_are_available:
            methods.append("{}{}".format(TAXONOMY_METHOD_PREFIX, source))
    return methods


def _expanded_explicit_methods(tokens, taxonomy_sources):
    methods: list[str] = []
    for token in tokens:
        if token == "taxonomy":
            methods.extend(
                "{}{}".format(TAXONOMY_METHOD_PREFIX, source)
                for source in taxonomy_sources
            )
        else:
            methods.append(token)
    return list(dict.fromkeys(methods))


def _method_is_excluded(method, excluded):
    if "all" in excluded:
        return True
    if method in excluded:
        return True
    return "taxonomy" in excluded and method.startswith(TAXONOMY_METHOD_PREFIX)


def resolve_root_compare_methods(tree, args):
    requested = _method_tokens(getattr(args, "methods", "all"), "--methods")
    for token in requested:
        _validate_method_token(token, "--methods", allow_all=True)
    if "all" in requested and requested != ["all"]:
        raise ValueError("'--methods all' cannot be combined with named methods.")

    excluded = _method_tokens(
        getattr(args, "exclude_methods", None),
        "--exclude-methods",
        allow_empty=True,
    )
    for token in excluded:
        _validate_method_token(token, "--exclude-methods", allow_all=True)

    automatic = requested == ["all"]
    automatic_taxonomy_disabled = automatic and bool(
        {"all", "taxonomy"}.intersection(excluded)
    )
    needs_configured_taxonomy_sources = "taxonomy" in requested or (
        automatic and not automatic_taxonomy_disabled
    )
    taxonomy_sources = (
        _parse_taxonomy_sources(
            getattr(args, "taxonomy_source", DEFAULT_TAXONOMY_SOURCE_CHAIN)
        )
        if needs_configured_taxonomy_sources
        else []
    )
    if automatic:
        methods = list(CORE_METHODS)
        if getattr(args, "outgroup", None) not in (None, ""):
            methods.append("outgroup")
        if getattr(args, "infile2", None) not in (None, ""):
            methods.append("transfer")
        if getattr(args, "species_tree", None) not in (None, ""):
            methods.append("reconciliation")
        if any(
            not _method_is_excluded(
                "{}{}".format(TAXONOMY_METHOD_PREFIX, source),
                excluded,
            )
            for source in taxonomy_sources
        ):
            methods.extend(_automatic_taxonomy_methods(tree, args, taxonomy_sources))
    else:
        methods = _expanded_explicit_methods(requested, taxonomy_sources)

    methods = [
        method for method in methods if not _method_is_excluded(method, excluded)
    ]
    if not methods:
        raise ValueError("No rooting methods remain after applying --exclude-methods.")
    return methods, automatic


def _edge_selection_evaluation(
    rooted_tree,
    *,
    method,
    selection_basis,
    source=None,
    position_kind="edge_unspecified",
):
    return RootingEvaluation(
        method=method,
        selection_basis=selection_basis,
        score_name="",
        score_unit="",
        candidates=(
            candidate_from_rooted_tree(
                rooted_tree,
                position_kind=position_kind,
            ),
        ),
        tie_rule="not_applicable",
        source=source,
    )


def _read_transfer_tree(args, cache):
    if "transfer_tree" not in cache:
        source = read_tree(
            args.infile2,
            getattr(args, "format2", "auto"),
            args.quoted_node_names,
            rooted=getattr(args, "infile2_rooted", "auto"),
        )
        source = remove_singleton(
            source,
            verbose=False,
            preserve_branch_length=True,
        )
        require_rooted(
            source,
            "'--infile2' must be rooted for root comparison.",
            "--infile2-rooted",
        )
        if len(list(source.leaves())) > 1 and len(source.get_children()) != 2:
            raise ValueError(
                "'--infile2' root must have exactly two children for root comparison."
            )
        cache["transfer_tree"] = source
    return cache["transfer_tree"]


def _finite_nonnegative_child_lengths(tree):
    children = tree.get_children()
    if len(children) != 2 or any(child.dist is None for child in children):
        return None
    try:
        lengths = [float(child.dist) for child in children]
    except (TypeError, ValueError, OverflowError):
        return None
    if any(not math.isfinite(length) or length < 0.0 for length in lengths):
        return None
    scale = max(lengths, default=0.0)
    if scale <= 0.0:
        return None
    scaled_total = math.fsum(length / scale for length in lengths)
    if not math.isfinite(scaled_total) or scaled_total <= 0.0:
        return None
    return children, lengths, scale, scaled_total


def _leaf_taxa(tree):
    return frozenset(map(str, tree.leaf_names()))


def _root_ratio_by_projected_side(tree, shared_taxa):
    values = _finite_nonnegative_child_lengths(tree)
    if values is None:
        return None
    children, lengths, scale, scaled_total = values
    ratios = {}
    shared_taxa = frozenset(str(taxon) for taxon in shared_taxa)
    for child, length in zip(children, lengths, strict=True):
        side = frozenset(
            str(name) for name in child.leaf_names() if str(name) in shared_taxa
        )
        if not side or side == shared_taxa or side in ratios:
            return None
        ratios[side] = (length / scale) / scaled_total
    return ratios


def _transferred_root_position_is_exact(rooted, source):
    shared_taxa = _leaf_taxa(rooted) & _leaf_taxa(source)
    source_ratios = _root_ratio_by_projected_side(source, shared_taxa)
    target_ratios = _root_ratio_by_projected_side(rooted, shared_taxa)
    if source_ratios is None or target_ratios is None:
        return False
    if set(source_ratios) != set(target_ratios):
        return False
    return all(
        math.isclose(
            target_ratios[side],
            source_ratios[side],
            rel_tol=1e-12,
            abs_tol=1e-15,
        )
        for side in source_ratios
    )


def _argument_source(args, name):
    value = getattr(args, name, None)
    return None if value in (None, "") else str(value)


def _read_species_tree(args, cache):
    if "species_tree" not in cache:
        cache["species_tree"] = read_tree(
            args.species_tree,
            getattr(args, "species_tree_format", "auto"),
            args.quoted_node_names,
            rooted=getattr(args, "species_tree_rooted", "auto"),
        )
    return cache["species_tree"]


def _evaluate_method(tree, method, args, cache):
    if method == "midpoint":
        _, evaluation = midpoint_rooting(tree, _return_evaluation=True)
        return evaluation
    if method == "mad":
        _, evaluation = mad_rooting(tree, _return_evaluation=True)
        return evaluation
    if method == "mv":
        _, evaluation = mv_rooting(tree, _return_evaluation=True)
        return evaluation
    if method == "outgroup":
        rooted = outgroup_rooting(tree, getattr(args, "outgroup", None))
        return _edge_selection_evaluation(
            rooted,
            method=method,
            selection_basis="specified_outgroup",
        )
    if method == "transfer":
        source = _read_transfer_tree(args, cache)
        rooted = transfer_root_with_taxon_mode(
            tree_to=tree,
            tree_from=source,
            taxon_mode=getattr(args, "taxon_mode", "exact"),
            verbose=False,
        )
        return _edge_selection_evaluation(
            rooted,
            method=method,
            selection_basis="reference_root",
            source=_argument_source(args, "infile2"),
            position_kind=(
                "exact"
                if _transferred_root_position_is_exact(rooted, source)
                else "edge_unspecified"
            ),
        )
    if method == "reconciliation":
        from nwkit.reconcile import _parsed_species_labels

        species_tree = _read_species_tree(args, cache)
        _, evaluation = reconciliation_rooting(
            tree,
            species_tree,
            _parsed_species_labels(tree, args),
            duplication_cost=getattr(args, "duplication_cost", 1.0),
            loss_cost=getattr(args, "loss_cost", 1.0),
            _return_evaluation=True,
        )
        return replace(evaluation, source=_argument_source(args, "species_tree"))
    if method.startswith(TAXONOMY_METHOD_PREFIX):
        source = method.removeprefix(TAXONOMY_METHOD_PREFIX)
        rooted = taxonomy_rooting(
            tree=tree,
            taxonomy_source=source,
            taxid_tsv=getattr(args, "taxid_tsv", None),
            rank=getattr(args, "rank", "no"),
            verbose=False,
            args=args,
        )
        return _edge_selection_evaluation(
            rooted,
            method=method,
            selection_basis="taxonomy_root",
            source=source,
        )
    raise ValueError("Unknown rooting method: {}".format(method))


def _split_sort_key(split):
    split = canonical_split(split[0], split[1])
    return tuple(sorted(split[0])), tuple(sorted(split[1]))


def _physical_edge_endpoints(tree):
    all_taxa = _leaf_taxa(tree)
    taxon_sets = get_subtree_leaf_name_sets(tree)
    nodes = list(tree.traverse())
    component_parent = {id(node): id(node) for node in nodes}

    def find_component(node_id):
        while component_parent[node_id] != node_id:
            component_parent[node_id] = component_parent[component_parent[node_id]]
            node_id = component_parent[node_id]
        return node_id

    def union_components(first, second):
        first_root = find_component(id(first))
        second_root = find_component(id(second))
        if first_root != second_root:
            component_parent[max(first_root, second_root)] = min(
                first_root,
                second_root,
            )

    def is_zero_length(value):
        if value is None:
            return True
        try:
            numeric = float(value)
        except (TypeError, ValueError, OverflowError):
            return False
        return math.isfinite(numeric) and numeric == 0.0

    root_children = tree.get_children()
    suppressed_root = tree if len(root_children) == 2 else None
    for node in nodes:
        if node.is_root:
            continue
        if suppressed_root is not None and node.up is suppressed_root:
            continue
        if is_zero_length(node.dist):
            union_components(node, node.up)
    if suppressed_root is not None and all(
        is_zero_length(child.dist) for child in root_children
    ):
        union_components(root_children[0], root_children[1])

    nodes_by_split: dict[Any, list[Any]] = {}
    for node in nodes:
        if node.is_root:
            continue
        side = frozenset(str(name) for name in taxon_sets[node])
        split = canonical_split(side, all_taxa - side)
        nodes_by_split.setdefault(split, []).append(node)

    endpoints = {}
    root_child_ids = {id(child) for child in root_children}
    for split, nodes in nodes_by_split.items():
        if len(nodes) == 1:
            node = nodes[0]
            side = frozenset(str(name) for name in taxon_sets[node])
            raw_endpoints = (node, node.up) if side == split[0] else (node.up, node)
            endpoints[split] = tuple(
                find_component(id(endpoint)) for endpoint in raw_endpoints
            )
            continue
        if len(nodes) == 2 and {id(node) for node in nodes} == root_child_ids:
            side_a = next(
                node
                for node in nodes
                if frozenset(str(name) for name in taxon_sets[node]) == split[0]
            )
            side_b = nodes[0] if nodes[1] is side_a else nodes[1]
            endpoints[split] = (
                find_component(id(side_a)),
                find_component(id(side_b)),
            )
    return endpoints


def _candidate_equivalent_splits(candidate):
    return candidate.equivalent_splits or (candidate.split,)


def _normalize_physical_root_positions(tree, evaluation, endpoints_by_split=None):
    if not evaluation.candidates:
        return evaluation
    if endpoints_by_split is None:
        endpoints_by_split = _physical_edge_endpoints(tree)
    node_groups: dict[int, tuple[Any, list[RootingCandidate]]] = {}
    retained = []
    for candidate in evaluation.candidates:
        fraction = candidate.position_fraction_from_side_a
        endpoints = endpoints_by_split.get(candidate.split)
        endpoint_index = None
        if (
            candidate.position_kind in {"exact", "node"}
            and fraction is not None
            and math.isfinite(float(fraction))
            and endpoints is not None
        ):
            if endpoints[0] == endpoints[1]:
                endpoint_index = 0
            elif math.isclose(float(fraction), 0.0, rel_tol=0.0, abs_tol=1e-12):
                endpoint_index = 0
            elif math.isclose(float(fraction), 1.0, rel_tol=0.0, abs_tol=1e-12):
                endpoint_index = 1
        if endpoint_index is None:
            retained.append(
                replace(
                    candidate,
                    equivalent_splits=tuple(
                        sorted(
                            set(_candidate_equivalent_splits(candidate)),
                            key=_split_sort_key,
                        )
                    ),
                )
            )
            continue
        endpoint = endpoints[endpoint_index]
        group = node_groups.setdefault(endpoint, (endpoint, []))[1]
        group.append(candidate)

    for endpoint, candidates in node_groups.values():
        representative = min(candidates, key=lambda item: _split_sort_key(item.split))
        representative_endpoints = endpoints_by_split[representative.split]
        fraction = 0.0 if representative_endpoints[0] == endpoint else 1.0
        equivalent_splits = {
            split
            for candidate in candidates
            for split in _candidate_equivalent_splits(candidate)
        }
        retained.append(
            replace(
                representative,
                position_kind="node",
                position_fraction_from_side_a=fraction,
                equivalent_splits=tuple(sorted(equivalent_splits, key=_split_sort_key)),
            )
        )

    retained.sort(
        key=lambda candidate: (
            _split_sort_key(candidate.split),
            (
                0.5
                if candidate.position_fraction_from_side_a is None
                else float(candidate.position_fraction_from_side_a)
            ),
        )
    )
    return replace(evaluation, candidates=tuple(retained))


def evaluate_root_compare_methods(tree, methods, args, automatic):
    evaluations = []
    cache: dict[str, Any] = {}
    endpoints_by_split = None
    for method in methods:
        try:
            evaluation = _evaluate_method(tree, method, args, cache)
        except Exception as exc:
            if not automatic:
                raise ValueError(
                    "Root comparison method '{}' failed: {}".format(method, exc)
                ) from exc
            if method.startswith(TAXONOMY_METHOD_PREFIX):
                source = method.removeprefix(TAXONOMY_METHOD_PREFIX)
            elif method == "transfer":
                source = _argument_source(args, "infile2")
            elif method == "reconciliation":
                source = _argument_source(args, "species_tree")
            else:
                source = None
            evaluation = RootingEvaluation.failed(method, exc, source=source)
            sys.stderr.write(
                "Warning: root comparison method '{}' failed: {}\n".format(
                    method,
                    exc,
                )
            )
        if evaluation.status == "ok":
            if endpoints_by_split is None:
                endpoints_by_split = _physical_edge_endpoints(tree)
            evaluation = _normalize_physical_root_positions(
                tree,
                evaluation,
                endpoints_by_split=endpoints_by_split,
            )
        evaluations.append(evaluation)
    return evaluations


def _format_value(value):
    if value in (None, ""):
        return ""
    if isinstance(value, Fraction):
        with localcontext() as context:
            context.prec = 15
            decimal_value = Decimal(value.numerator) / Decimal(value.denominator)
        return format(decimal_value.normalize(), ".15g")
    if isinstance(value, float):
        return format(value, ".15g") if math.isfinite(value) else ""
    return str(value)


def _single_line_message(message):
    return " ".join(str(message).split())


def root_compare_rows(evaluations):
    rows = []
    for evaluation in evaluations:
        if evaluation.status != "ok" or not evaluation.candidates:
            rows.append(
                {column: "" for column in ROOT_COMPARE_COLUMNS}
                | {
                    "method": evaluation.method,
                    "status": evaluation.status,
                    "message": _single_line_message(evaluation.message),
                    "source": evaluation.source or "",
                    "num_best": "0",
                }
            )
            continue
        num_best = len(evaluation.candidates)
        for tie_index, candidate in enumerate(evaluation.candidates, start=1):
            metrics = candidate.metrics
            equivalent_splits = _candidate_equivalent_splits(candidate)
            rows.append(
                {
                    "method": evaluation.method,
                    "status": evaluation.status,
                    "message": _single_line_message(evaluation.message),
                    "selection_basis": evaluation.selection_basis,
                    "source": evaluation.source or "",
                    "tie_index": str(tie_index),
                    "num_best": str(num_best),
                    "is_canonical": "yes" if tie_index == 1 else "no",
                    "root_split_id": root_split_id(candidate.split),
                    "num_equivalent_splits": str(len(equivalent_splits)),
                    "equivalent_root_split_ids": ",".join(
                        root_split_id(split) for split in equivalent_splits
                    ),
                    "side_a_taxa": format_taxa(candidate.split[0]),
                    "side_b_taxa": format_taxa(candidate.split[1]),
                    "num_side_a_taxa": str(len(candidate.split[0])),
                    "num_side_b_taxa": str(len(candidate.split[1])),
                    "position_kind": candidate.position_kind,
                    "position_fraction_from_side_a": _format_value(
                        candidate.position_fraction_from_side_a
                    ),
                    "edge_length": _format_value(candidate.edge_length),
                    "score_name": evaluation.score_name,
                    "score": _format_value(candidate.score),
                    "score_unit": evaluation.score_unit,
                    "num_evaluated_edges": _format_value(evaluation.evaluated_edges),
                    "tie_rule": evaluation.tie_rule,
                    "duplications": _format_value(metrics.get("duplications")),
                    "losses": _format_value(metrics.get("losses")),
                    "duplication_cost": _format_value(metrics.get("duplication_cost")),
                    "loss_cost": _format_value(metrics.get("loss_cost")),
                }
            )
    return rows


def _write_root_compare_rows(outfile, rows):
    if outfile == "-":
        writer = csv.DictWriter(
            sys.stdout,
            fieldnames=ROOT_COMPARE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        return
    with open(outfile, mode="w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ROOT_COMPARE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _marker_label(evaluation):
    label = _METHOD_LABELS.get(evaluation.method, evaluation.method)
    first = evaluation.candidates[0]
    details = []
    if evaluation.score_name and first.score is not None:
        score_label = {
            "ancestor_deviation": "score",
            "root_to_tip_variance": "variance",
            "weighted_duplication_loss_cost": "cost",
        }.get(evaluation.score_name, evaluation.score_name)
        try:
            numeric_score = float(first.score)
        except (TypeError, ValueError, OverflowError):
            displayed_score = _format_value(first.score)
        else:
            displayed_score = (
                format(numeric_score, ".6g")
                if math.isfinite(numeric_score)
                else _format_value(first.score)
            )
        details.append("{}={}".format(score_label, displayed_score))
    if len(evaluation.candidates) > 1:
        details.append("{} ties".format(len(evaluation.candidates)))
    equivalent_split_count = sum(
        len(_candidate_equivalent_splits(candidate))
        for candidate in evaluation.candidates
    )
    if equivalent_split_count > len(evaluation.candidates):
        details.append("{} equivalent splits".format(equivalent_split_count))
    if any(
        candidate.position_kind not in {"exact", "node"}
        or candidate.position_fraction_from_side_a is None
        for candidate in evaluation.candidates
    ):
        details.append("edge only")
    if details:
        return "{} ({})".format(label, "; ".join(details))
    return label


def _rooting_markers(evaluations):
    markers = []
    for evaluation in evaluations:
        if evaluation.status != "ok" or not evaluation.candidates:
            continue
        color, marker, size = _MARKER_STYLES[evaluation.method]
        label = _marker_label(evaluation)
        for candidate in evaluation.candidates:
            exact = (
                candidate.position_kind in {"exact", "node"}
                and candidate.position_fraction_from_side_a is not None
            )
            markers.append(
                {
                    "split": candidate.split,
                    "position_fraction_from_side_a": (
                        candidate.position_fraction_from_side_a if exact else None
                    ),
                    "color": color,
                    "marker": marker,
                    "filled": exact,
                    "size": size,
                    "label": label,
                }
            )
    return markers


def _resolved_unrooted_method(tree, requested):
    requested = str(requested or "auto").strip().lower()
    if requested != "auto":
        return requested
    displayed_nodes = sum(1 for _ in tree.traverse())
    if len(tree.get_children()) == 2:
        displayed_nodes -= 1
    return "equal-daylight" if displayed_nodes <= 2000 else "equal-angle"


def _resolved_tip_labels(tree, requested):
    if isinstance(requested, bool):
        return requested
    requested = str(requested or "auto").strip().lower()
    if requested == "yes":
        return True
    if requested == "no":
        return False
    return sum(1 for _ in tree.leaves()) <= AUTO_TIP_LABEL_LIMIT


def _resolved_figure_width(tree, requested, font_size, tip_labels):
    if requested is not None:
        return float(requested)
    if not tip_labels:
        return DEFAULT_FIGURE_WIDTH
    tip_count = sum(1 for _ in tree.leaves())
    density_width = tip_count * float(font_size) * AUTO_WIDTH_INCHES_PER_TIP_FONT_POINT
    return min(
        max(DEFAULT_FIGURE_WIDTH, density_width),
        MAXIMUM_AUTO_FIGURE_WIDTH,
    )


def _draw_root_comparison(tree, evaluations, args):
    markers = _rooting_markers(evaluations)
    if not markers:
        raise ValueError("No successful rooting position is available to draw.")
    requested_tip_labels = getattr(args, "tip_labels", "auto")
    tip_labels = _resolved_tip_labels(tree, requested_tip_labels)
    if str(requested_tip_labels or "auto").strip().lower() == "auto" and not tip_labels:
        sys.stderr.write(
            "Tip labels were omitted because the drawing exceeds {} tips; "
            "use '--tip-labels yes' to force them.\n".format(AUTO_TIP_LABEL_LIMIT)
        )
    _draw_tree(
        tree=tree,
        outfile=args.figure_out,
        image_format="pdf",
        node_type_by_node={},
        leaf_label_color_by_leaf={},
        group_color_by_name={},
        support_labels=False,
        figure_width=_resolved_figure_width(
            tree,
            getattr(args, "figure_width", None),
            getattr(args, "font_size", 8.0),
            tip_labels,
        ),
        figure_height=getattr(args, "figure_height", None),
        font_size=getattr(args, "font_size", 8.0),
        layout="unrooted",
        tip_label_position="branch-end",
        tip_labels=tip_labels,
        tip_label_wrap="auto",
        tip_spacing="label-aware",
        unrooted_method=_resolved_unrooted_method(
            tree,
            getattr(args, "unrooted_method", "auto"),
        ),
        root_marker="none",
        legend=True,
        legend_position="right",
        collision_policy="resolve",
        layout_report=getattr(args, "layout_report", None),
        branch_markers=markers,
        transactional_output=not getattr(args, "_nwkit_outputs_staged", False),
    )


def _validate_root_compare_paths(args):
    if args.figure_out in (None, "", "-"):
        raise ValueError("'--figure-out' must be a PDF file path, not STDOUT.")
    if os.path.splitext(str(args.figure_out))[1].lower() != ".pdf":
        raise ValueError("'--figure-out' must use the .pdf extension.")
    for option_name in ("figure_width", "figure_height", "font_size"):
        value = getattr(args, option_name, None)
        if value is not None and float(value) <= 0:
            raise ValueError(
                "'--{}' must be positive.".format(option_name.replace("_", "-"))
            )
    outputs = [
        ("--outfile", args.outfile),
        ("--figure-out", args.figure_out),
    ]
    layout_report = getattr(args, "layout_report", None)
    if layout_report not in (None, ""):
        if layout_report == "-":
            raise ValueError("'--layout-report' must be a file path, not STDOUT.")
        outputs.append(("--layout-report", layout_report))
    for option_name, path in outputs:
        if path in (None, "", "-"):
            continue
        real_path = os.path.realpath(os.fspath(path))
        if os.path.isdir(real_path):
            raise ValueError(
                "'{}' must be a file path, not a directory.".format(option_name)
            )
        if os.path.exists(real_path) and not stat.S_ISREG(os.stat(real_path).st_mode):
            raise ValueError("'{}' must identify a regular file.".format(option_name))
    validate_distinct_output_paths(outputs)
    validate_outputs_do_not_replace_inputs(
        [
            ("--infile", args.infile),
            ("--infile2", getattr(args, "infile2", None)),
            ("--species-tree", getattr(args, "species_tree", None)),
            ("--species-map-tsv", getattr(args, "species_map_tsv", None)),
            ("--taxid-tsv", getattr(args, "taxid_tsv", None)),
        ],
        outputs,
        label="Root comparison output",
    )


def _write_failed_evaluations(args, rows):
    if args.outfile == "-":
        _write_root_compare_rows("-", rows)
        return
    with output_transaction([args.outfile]) as staged:
        _write_root_compare_rows(staged[args.outfile], rows)


def _write_successful_outputs(tree, evaluations, args, rows):
    outputs = [args.figure_out]
    if args.outfile != "-":
        outputs.append(args.outfile)
    layout_report = getattr(args, "layout_report", None)
    if layout_report not in (None, ""):
        outputs.append(layout_report)
    write_stdout = (
        (lambda: _write_root_compare_rows("-", rows)) if args.outfile == "-" else None
    )
    with output_transaction(outputs, after_install=write_stdout) as staged:
        drawing_args = copy(args)
        drawing_args.figure_out = staged[args.figure_out]
        drawing_args._nwkit_outputs_staged = True
        if layout_report not in (None, ""):
            drawing_args.layout_report = staged[layout_report]
        _draw_root_comparison(tree, evaluations, drawing_args)
        if args.outfile != "-":
            _write_root_compare_rows(staged[args.outfile], rows)


def root_compare_main(args):
    _validate_root_compare_paths(args)
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    validate_unique_named_leaves(
        tree,
        option_name="--infile",
        context=" for root comparison",
    )
    if len(list(tree.leaves())) < 2:
        raise ValueError("Root comparison requires at least two tree tips.")
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
    methods, automatic = resolve_root_compare_methods(tree, args)
    evaluations = evaluate_root_compare_methods(tree, methods, args, automatic)
    rows = root_compare_rows(evaluations)
    successful = any(
        evaluation.status == "ok" and evaluation.candidates
        for evaluation in evaluations
    )
    if not successful:
        _write_failed_evaluations(args, rows)
        raise ValueError(
            "Every selected rooting method failed; the failure summary TSV was written, "
            "and the PDF was left unchanged."
        )
    _write_successful_outputs(tree, evaluations, args, rows)
    return evaluations
