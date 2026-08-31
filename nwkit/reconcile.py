import sys
from typing import Any

import pandas as pd

from nwkit.clade_index import CladeIndex, LcaIndex
from nwkit.reconciliation_properties import (
    COLLAPSED_EVENT_BOUNDARY_PROP,
    property_is_true,
)
from nwkit.rooting_state import require_rooted, rooting_option_for_input
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    assign_branch_ids,
    get_node_class,
    read_tree,
    validate_distinct_output_paths,
    validate_outputs_do_not_replace_inputs,
    validate_unique_named_leaves,
)

RECONCILIATION_COLUMNS = [
    "tree_id",
    "gene_branch_id",
    "gene_clade_id",
    "parent_gene_branch_id",
    "parent_gene_clade_id",
    "node_class",
    "gene_name",
    "lineage_id",
    "lineage_clade_id",
    "event_type",
    "event_source",
    "event_status",
    "collapsed_event_boundary",
    "transfer_source_species",
    "transfer_destination_species",
    "species_branch_id",
    "species_event_id",
    "species_node_class",
    "species_name",
    "mapping_status",
    "eligible",
    "coverage_status",
    "reason",
    "contrast_numerator_gene_branch_id",
    "contrast_numerator_gene_clade_id",
    "contrast_denominator_gene_branch_id",
    "contrast_denominator_gene_clade_id",
    "species_numerator_branch_id",
    "species_numerator_event_id",
    "species_denominator_branch_id",
    "species_denominator_event_id",
    "descendant_taxa",
    "num_taxa",
    "descendant_species_taxa",
    "num_descendant_species_taxa",
    "species_event_taxa",
    "num_species_event_taxa",
    "numerator_observed_species_taxa",
    "num_numerator_observed_species_taxa",
    "numerator_species_clade_taxa",
    "num_numerator_species_clade_taxa",
    "numerator_species_coverage",
    "denominator_observed_species_taxa",
    "num_denominator_observed_species_taxa",
    "denominator_species_clade_taxa",
    "num_denominator_species_clade_taxa",
    "denominator_species_coverage",
]


def _validate_rooted_binary_tree(tree, option_name):
    validate_unique_named_leaves(tree, option_name=option_name)
    if len(list(tree.leaves())) < 2:
        raise ValueError("'{}' must contain at least two tips.".format(option_name))
    require_rooted(
        tree,
        "'{}' must be rooted.".format(option_name),
        rooting_option_for_input(option_name),
    )
    multifurcations = [
        node for node in tree.traverse() if not node.is_leaf and len(node.children) != 2
    ]
    if multifurcations:
        raise ValueError(
            "'{}' must be strictly bifurcating; found {} non-binary internal node(s).".format(
                option_name, len(multifurcations)
            )
        )


def _parse_nhx_boolean(value):
    normalized = str(value).strip().upper()
    if normalized in {"Y", "YES", "TRUE", "1"}:
        return True
    if normalized in {"N", "NO", "FALSE", "0"}:
        return False
    return None


def _parse_nhx_transfer(value):
    parts = str(value).strip().split("@")
    is_transfer = _parse_nhx_boolean(parts[0])
    if is_transfer is None:
        return None, "", ""
    if not is_transfer:
        return (False, "", "") if len(parts) == 1 else (None, "", "")
    if len(parts) == 1:
        return True, "", ""
    if len(parts) == 3 and parts[1] != "" and parts[2] != "":
        return True, parts[1], parts[2]
    return None, "", ""


def _nhx_event(node):
    transfer_value = node.props.get("H")
    if transfer_value is not None:
        is_transfer, transfer_source, transfer_destination = _parse_nhx_transfer(
            transfer_value
        )
        if is_transfer is None:
            return "unresolved", "unresolved", "invalid_nhx_H", "", ""
        if is_transfer:
            return (
                "transfer",
                "resolved",
                "",
                transfer_source,
                transfer_destination,
            )
    duplication_value = node.props.get("D")
    if duplication_value is None:
        return "unresolved", "unresolved", "missing_nhx_D", "", ""
    is_duplication = _parse_nhx_boolean(duplication_value)
    if is_duplication is None:
        return "unresolved", "unresolved", "invalid_nhx_D", "", ""
    if is_duplication:
        return "duplication", "resolved", "", "", ""
    return "speciation", "resolved", "", "", ""


def _mapped_species_nodes(
    gene_tree,
    species_label_by_gene_leaf,
    species_leaf_by_name,
    species_clades,
    species_lca,
    event_source,
    species_node_by_name,
    ambiguous_species_names,
):
    species_node_by_gene_node: dict[Any, Any | None] = dict()
    species_mask_by_gene_node = dict()
    mapping_reason_by_node = dict()
    for node in gene_tree.traverse(strategy="postorder"):
        if node.is_leaf:
            species_label = species_label_by_gene_leaf.get(str(node.name))
            if species_label is None:
                species_node_by_gene_node[node] = None
                species_mask_by_gene_node[node] = 0
                mapping_reason_by_node[node] = "species_label_unresolved"
            elif species_label not in species_leaf_by_name:
                species_node_by_gene_node[node] = None
                species_mask_by_gene_node[node] = 0
                mapping_reason_by_node[node] = "species_not_in_species_tree"
            else:
                species_node = species_leaf_by_name[species_label]
                species_node_by_gene_node[node] = species_node
                species_mask_by_gene_node[node] = species_clades.mask_by_node[
                    species_node
                ]
                mapping_reason_by_node[node] = ""
            continue
        child_nodes = [species_node_by_gene_node[child] for child in node.children]
        species_mask_by_gene_node[node] = 0
        for child in node.children:
            species_mask_by_gene_node[node] |= species_mask_by_gene_node[child]
        if any(child_node is None for child_node in child_nodes):
            species_node_by_gene_node[node] = None
            mapping_reason_by_node[node] = "descendant_mapping_unresolved"
        else:
            species_node_by_gene_node[node] = species_lca.common_ancestor(
                child_nodes[0], child_nodes[1]
            )
            mapping_reason_by_node[node] = ""
    if event_source == "nhx":
        for node in gene_tree.traverse():
            annotated_species = node.props.get("S")
            if annotated_species in (None, ""):
                continue
            annotated_species = str(annotated_species)
            if annotated_species in ambiguous_species_names:
                species_node_by_gene_node[node] = None
                mapping_reason_by_node[node] = "ambiguous_nhx_S"
                continue
            mapped_node = species_node_by_name.get(annotated_species)
            if mapped_node is None:
                species_node_by_gene_node[node] = None
                mapping_reason_by_node[node] = "nhx_S_not_in_species_tree"
                continue
            if node.is_leaf:
                parsed_species = species_label_by_gene_leaf.get(str(node.name))
                if parsed_species != annotated_species:
                    species_node_by_gene_node[node] = None
                    mapping_reason_by_node[node] = "nhx_S_conflicts_with_species_map"
                    continue
            species_node_by_gene_node[node] = mapped_node
            mapping_reason_by_node[node] = ""
    return (
        species_node_by_gene_node,
        species_mask_by_gene_node,
        mapping_reason_by_node,
    )


def _event_annotations(
    gene_tree,
    event_source,
    species_node_by_gene_node,
    species_mask_by_gene_node,
    species_node_by_name,
    ambiguous_species_names,
):
    event_type_by_node = dict()
    event_status_by_node = dict()
    event_reason_by_node = dict()
    transfer_by_node = dict()
    for node in gene_tree.traverse(strategy="postorder"):
        if node.is_leaf:
            event_type_by_node[node] = "leaf"
            event_status_by_node[node] = "not-applicable"
            event_reason_by_node[node] = ""
            transfer_by_node[node] = ("", "")
            continue
        mapped_node = species_node_by_gene_node[node]
        transfer_source = transfer_destination = ""
        if event_source == "nhx":
            (
                event_type,
                event_status,
                event_reason,
                transfer_source,
                transfer_destination,
            ) = _nhx_event(node)
            if event_type == "transfer" and transfer_source != "":
                endpoints = (transfer_source, transfer_destination)
                if transfer_source == transfer_destination:
                    event_type, event_status = "unresolved", "unresolved"
                    event_reason = "invalid_nhx_H_endpoints"
                    transfer_source = transfer_destination = ""
                elif any(name in ambiguous_species_names for name in endpoints):
                    event_type, event_status = "unresolved", "unresolved"
                    event_reason = "ambiguous_nhx_H_endpoint"
                    transfer_source = transfer_destination = ""
                elif any(name not in species_node_by_name for name in endpoints):
                    event_type, event_status = "unresolved", "unresolved"
                    event_reason = "nhx_H_endpoint_not_in_species_tree"
                    transfer_source = transfer_destination = ""
                elif node.props.get("S") not in (None, "", transfer_source):
                    event_type, event_status = "unresolved", "unresolved"
                    event_reason = "nhx_H_source_conflicts_with_S"
                    transfer_source = transfer_destination = ""
        elif mapped_node is None:
            event_type, event_status = "unresolved", "unresolved"
            event_reason = "mapping_unresolved"
        elif event_source == "species-overlap":
            child1, child2 = node.children
            event_type = (
                "duplication"
                if species_mask_by_gene_node[child1] & species_mask_by_gene_node[child2]
                else "speciation"
            )
            event_status, event_reason = "resolved", ""
        elif event_source == "lca":
            child_species_nodes = [
                species_node_by_gene_node[child] for child in node.children
            ]
            event_type = (
                "duplication"
                if any(child is mapped_node for child in child_species_nodes)
                else "speciation"
            )
            event_status, event_reason = "resolved", ""
        else:
            raise ValueError("Unsupported event source: {}".format(event_source))
        event_type_by_node[node] = event_type
        event_status_by_node[node] = event_status
        event_reason_by_node[node] = event_reason
        transfer_by_node[node] = (transfer_source, transfer_destination)
    return (
        event_type_by_node,
        event_status_by_node,
        event_reason_by_node,
        transfer_by_node,
    )


def lca_duplication_loss_contribution(
    mapped_species_node,
    child_species_nodes,
    species_depth_by_node,
):
    """Return the standard LCA-reconciliation D/L contribution of one event.

    Losses are counted on species-tree edges.  A duplication contributes the
    complete mapped distance to each child lineage, whereas a speciation
    consumes the first edge on each child path and therefore contributes
    ``distance - 1``.
    """
    child_species_nodes = tuple(child_species_nodes)
    parent_depth = species_depth_by_node[mapped_species_node]
    child_depths = tuple(species_depth_by_node[child] for child in child_species_nodes)
    is_duplication = any(child is mapped_species_node for child in child_species_nodes)
    speciation_offset = 0 if is_duplication else 1
    losses = sum(
        child_depth - parent_depth - speciation_offset for child_depth in child_depths
    )
    return int(is_duplication), losses


def _is_subset(mask, container_mask):
    return mask != 0 and mask & ~container_mask == 0


def _speciation_orientation(
    gene_node,
    mapped_species_node,
    species_mask_by_gene_node,
    species_clades,
):
    if mapped_species_node.is_leaf:
        return None, "mapped_to_species_leaf"
    species_children = sorted(
        mapped_species_node.children,
        key=species_clades.names_for_node,
    )
    numerator_species, denominator_species = species_children
    numerator_mask = species_clades.mask_by_node[numerator_species]
    denominator_mask = species_clades.mask_by_node[denominator_species]
    gene_child1, gene_child2 = gene_node.children
    child1_mask = species_mask_by_gene_node[gene_child1]
    child2_mask = species_mask_by_gene_node[gene_child2]
    if _is_subset(child1_mask, numerator_mask) and _is_subset(
        child2_mask, denominator_mask
    ):
        return (
            gene_child1,
            gene_child2,
            numerator_species,
            denominator_species,
        ), ""
    if _is_subset(child2_mask, numerator_mask) and _is_subset(
        child1_mask, denominator_mask
    ):
        return (
            gene_child2,
            gene_child1,
            numerator_species,
            denominator_species,
        ), ""
    return None, "children_do_not_match_distinct_species_clades"


def _lineage_roots(gene_tree, event_type_by_node):
    lineage_root_by_node = {gene_tree: gene_tree}
    for node in gene_tree.traverse(strategy="preorder"):
        if property_is_true(node.props.get(COLLAPSED_EVENT_BOUNDARY_PROP, "N")):
            lineage_root_by_node[node] = node
        if node.is_leaf:
            continue
        starts_new_lineage = event_type_by_node[node] != "speciation"
        for child in node.children:
            lineage_root_by_node[child] = (
                child if starts_new_lineage else lineage_root_by_node[node]
            )
    return lineage_root_by_node


def _coverage_fields(orientation, species_mask_by_gene_node, species_clades):
    empty = {
        "coverage_status": "not-applicable",
        "numerator_observed_species_taxa": "",
        "num_numerator_observed_species_taxa": "",
        "numerator_species_clade_taxa": "",
        "num_numerator_species_clade_taxa": "",
        "numerator_species_coverage": "",
        "denominator_observed_species_taxa": "",
        "num_denominator_observed_species_taxa": "",
        "denominator_species_clade_taxa": "",
        "num_denominator_species_clade_taxa": "",
        "denominator_species_coverage": "",
    }
    if orientation is None:
        return empty
    numerator_gene, denominator_gene, numerator_species, denominator_species = (
        orientation
    )
    numerator_observed = species_mask_by_gene_node[numerator_gene]
    denominator_observed = species_mask_by_gene_node[denominator_gene]
    numerator_full = species_clades.mask_by_node[numerator_species]
    denominator_full = species_clades.mask_by_node[denominator_species]
    numerator_observed_count = species_clades.count_for_mask(numerator_observed)
    denominator_observed_count = species_clades.count_for_mask(denominator_observed)
    numerator_full_count = species_clades.count_for_mask(numerator_full)
    denominator_full_count = species_clades.count_for_mask(denominator_full)
    complete = (
        numerator_observed == numerator_full
        and denominator_observed == denominator_full
    )
    return {
        "coverage_status": "complete" if complete else "partial",
        "numerator_observed_species_taxa": species_clades.csv_for_mask(
            numerator_observed
        ),
        "num_numerator_observed_species_taxa": numerator_observed_count,
        "numerator_species_clade_taxa": species_clades.csv_for_mask(numerator_full),
        "num_numerator_species_clade_taxa": numerator_full_count,
        "numerator_species_coverage": numerator_observed_count / numerator_full_count,
        "denominator_observed_species_taxa": species_clades.csv_for_mask(
            denominator_observed
        ),
        "num_denominator_observed_species_taxa": denominator_observed_count,
        "denominator_species_clade_taxa": species_clades.csv_for_mask(denominator_full),
        "num_denominator_species_clade_taxa": denominator_full_count,
        "denominator_species_coverage": denominator_observed_count
        / denominator_full_count,
    }


def _normalized_species_mapping(gene_tree, species_label_by_gene_leaf):
    if not isinstance(species_label_by_gene_leaf, dict):
        raise ValueError("Gene-tip species labels must be supplied as a dictionary.")
    gene_leaf_names = {str(leaf.name) for leaf in gene_tree.leaves()}
    mapping_names = {str(name) for name in species_label_by_gene_leaf}
    if len(mapping_names) != len(species_label_by_gene_leaf):
        raise ValueError(
            "Gene-tip species mapping keys must remain unique when converted to strings."
        )
    missing_mapping = sorted(gene_leaf_names - mapping_names)
    extra_mapping = sorted(mapping_names - gene_leaf_names)
    if missing_mapping or extra_mapping:
        raise ValueError(
            "Gene-tip species mapping must exactly cover gene-tree tips "
            "(missing={}; extra={}).".format(
                ",".join(missing_mapping), ",".join(extra_mapping)
            )
        )
    normalized: dict[str, str | None] = {}
    for name, species_label in species_label_by_gene_leaf.items():
        normalized_name = str(name)
        if species_label is None:
            normalized[normalized_name] = None
            continue
        normalized_label = str(species_label).strip()
        if normalized_label == "":
            raise ValueError(
                "Gene-tip species mapping values must be non-empty strings or None."
            )
        normalized[normalized_name] = normalized_label
    return normalized


def _named_species_nodes(species_tree):
    species_nodes_by_name = {}
    ambiguous_species_names = set()
    for node in species_tree.traverse():
        if node.name in (None, ""):
            continue
        name = str(node.name)
        if name in species_nodes_by_name:
            ambiguous_species_names.add(name)
        else:
            species_nodes_by_name[name] = node
    return species_nodes_by_name, ambiguous_species_names


def _reconciliation_state(
    node,
    mapped_species_node,
    mapping_reason,
    event_type,
    event_status,
    event_reason,
    species_mask_by_gene_node,
    species_clades,
):
    mapping_status = "mapped" if mapped_species_node is not None else "unmapped"
    orientation = None
    reason = "leaf" if node.is_leaf else ""
    if mapping_status != "mapped":
        reason = mapping_reason
    elif event_status != "resolved" and not node.is_leaf:
        reason = event_reason
    elif event_type != "speciation" and not node.is_leaf:
        reason = "event_not_speciation"
    elif not node.is_leaf:
        orientation, reason = _speciation_orientation(
            node,
            mapped_species_node,
            species_mask_by_gene_node,
            species_clades,
        )
    return mapping_status, orientation, reason


def _orientation_nodes(orientation):
    if orientation is None:
        return None, None, None, None
    return orientation


def _optional_node_value(node, value):
    return "" if node is None else value


def _parent_branch_id(node, branch_ids):
    return -1 if node.is_root else branch_ids[node.up]


def _parent_clade_id(node, clades):
    return "" if node.is_root else clades.clade_id_for_node(node.up)


def _optional_clade_id(node, clades):
    return "" if node is None else clades.clade_id_for_node(node)


def _optional_node_class(node):
    return "" if node is None else get_node_class(node)


def _node_name(node):
    return "" if node is None or node.name in (None, "") else str(node.name)


def _yes_no(value):
    return "yes" if value else "no"


def _reconciliation_row(
    node,
    tree_id,
    mapped_species_node,
    mapping_status,
    orientation,
    reason,
    lineage_root,
    transfer,
    event_source,
    event_type,
    event_status,
    gene_branch_ids,
    species_branch_ids,
    gene_clades,
    species_clades,
    species_mask_by_gene_node,
):
    numerator_gene, denominator_gene, numerator_species, denominator_species = (
        _orientation_nodes(orientation)
    )
    transfer_source, transfer_destination = transfer
    species_event_mask = (
        0
        if mapped_species_node is None
        else species_clades.mask_by_node[mapped_species_node]
    )
    row = {
        "tree_id": tree_id,
        "gene_branch_id": gene_branch_ids[node],
        "gene_clade_id": gene_clades.clade_id_for_node(node),
        "parent_gene_branch_id": _parent_branch_id(node, gene_branch_ids),
        "parent_gene_clade_id": _parent_clade_id(node, gene_clades),
        "node_class": get_node_class(node),
        "gene_name": _node_name(node),
        "lineage_id": gene_branch_ids[lineage_root],
        "lineage_clade_id": gene_clades.clade_id_for_node(lineage_root),
        "event_type": event_type,
        "event_source": event_source,
        "event_status": event_status,
        "collapsed_event_boundary": _yes_no(
            property_is_true(node.props.get(COLLAPSED_EVENT_BOUNDARY_PROP, "N"))
        ),
        "transfer_source_species": transfer_source,
        "transfer_destination_species": transfer_destination,
        "species_branch_id": _optional_node_value(
            mapped_species_node, species_branch_ids.get(mapped_species_node)
        ),
        "species_event_id": _optional_clade_id(mapped_species_node, species_clades),
        "species_node_class": _optional_node_class(mapped_species_node),
        "species_name": _node_name(mapped_species_node),
        "mapping_status": mapping_status,
        "eligible": _yes_no(orientation is not None),
        "reason": reason,
        "contrast_numerator_gene_branch_id": _optional_node_value(
            numerator_gene, gene_branch_ids.get(numerator_gene)
        ),
        "contrast_numerator_gene_clade_id": _optional_clade_id(
            numerator_gene, gene_clades
        ),
        "contrast_denominator_gene_branch_id": _optional_node_value(
            denominator_gene, gene_branch_ids.get(denominator_gene)
        ),
        "contrast_denominator_gene_clade_id": _optional_clade_id(
            denominator_gene, gene_clades
        ),
        "species_numerator_branch_id": _optional_node_value(
            numerator_species, species_branch_ids.get(numerator_species)
        ),
        "species_numerator_event_id": _optional_clade_id(
            numerator_species, species_clades
        ),
        "species_denominator_branch_id": _optional_node_value(
            denominator_species, species_branch_ids.get(denominator_species)
        ),
        "species_denominator_event_id": _optional_clade_id(
            denominator_species, species_clades
        ),
        "descendant_taxa": gene_clades.csv_for_node(node),
        "num_taxa": gene_clades.count_for_node(node),
        "descendant_species_taxa": species_clades.csv_for_mask(
            species_mask_by_gene_node[node]
        ),
        "num_descendant_species_taxa": species_clades.count_for_mask(
            species_mask_by_gene_node[node]
        ),
        "species_event_taxa": species_clades.csv_for_mask(species_event_mask),
        "num_species_event_taxa": species_clades.count_for_mask(species_event_mask),
    }
    row.update(_coverage_fields(orientation, species_mask_by_gene_node, species_clades))
    return row


def build_reconciliation_table(
    gene_tree,
    species_tree,
    species_label_by_gene_leaf,
    event_source="lca",
    tree_id="",
):
    _validate_rooted_binary_tree(gene_tree, "--infile")
    _validate_rooted_binary_tree(species_tree, "--species-tree")
    if event_source not in {"lca", "nhx", "species-overlap"}:
        raise ValueError("Unsupported event source: {}".format(event_source))
    species_label_by_gene_leaf = _normalized_species_mapping(
        gene_tree, species_label_by_gene_leaf
    )
    gene_branch_ids = assign_branch_ids(gene_tree)
    species_branch_ids = assign_branch_ids(species_tree)
    gene_clades = CladeIndex(gene_tree)
    species_clades = CladeIndex(species_tree)
    species_lca = LcaIndex(species_tree)
    species_leaf_by_name = {str(leaf.name): leaf for leaf in species_tree.leaves()}
    species_nodes_by_name, ambiguous_species_names = _named_species_nodes(species_tree)
    (
        species_node_by_gene_node,
        species_mask_by_gene_node,
        mapping_reason_by_node,
    ) = _mapped_species_nodes(
        gene_tree,
        species_label_by_gene_leaf,
        species_leaf_by_name,
        species_clades,
        species_lca,
        event_source,
        species_nodes_by_name,
        ambiguous_species_names,
    )
    (
        event_type_by_node,
        event_status_by_node,
        event_reason_by_node,
        transfer_by_node,
    ) = _event_annotations(
        gene_tree,
        event_source,
        species_node_by_gene_node,
        species_mask_by_gene_node,
        species_nodes_by_name,
        ambiguous_species_names,
    )
    lineage_root_by_node = _lineage_roots(gene_tree, event_type_by_node)
    rows = list()
    for node in gene_tree.traverse():
        mapped_species_node = species_node_by_gene_node[node]
        mapping_status, orientation, reason = _reconciliation_state(
            node,
            mapped_species_node,
            mapping_reason_by_node[node],
            event_type_by_node[node],
            event_status_by_node[node],
            event_reason_by_node[node],
            species_mask_by_gene_node,
            species_clades,
        )
        rows.append(
            _reconciliation_row(
                node,
                tree_id,
                mapped_species_node,
                mapping_status,
                orientation,
                reason,
                lineage_root_by_node[node],
                transfer_by_node[node],
                event_source,
                event_type_by_node[node],
                event_status_by_node[node],
                gene_branch_ids,
                species_branch_ids,
                gene_clades,
                species_clades,
                species_mask_by_gene_node,
            )
        )
    return pd.DataFrame(rows, columns=RECONCILIATION_COLUMNS)


def _parsed_species_labels(gene_tree, args):
    parser = get_species_parser(args=args)
    return {
        str(leaf.name): parser.parse(leaf.name).species_label
        for leaf in gene_tree.leaves()
    }


def _report_unmatched_species(species_label_by_gene_leaf, species_tree, policy):
    species_tree_labels = {str(leaf.name) for leaf in species_tree.leaves()}
    unresolved = sorted(
        leaf_name
        for leaf_name, species_label in species_label_by_gene_leaf.items()
        if species_label is None
    )
    absent = sorted(
        "{}={}".format(leaf_name, species_label)
        for leaf_name, species_label in species_label_by_gene_leaf.items()
        if species_label is not None and species_label not in species_tree_labels
    )
    if policy == "error" and (unresolved or absent):
        raise ValueError(
            "Gene tips could not be mapped to --species-tree "
            "(unresolved={}; absent={}).".format(",".join(unresolved), ",".join(absent))
        )
    if policy == "warn":
        if unresolved:
            sys.stderr.write(
                "Warning: species labels could not be parsed for gene tips: {}\n".format(
                    " ".join(unresolved)
                )
            )
        if absent:
            sys.stderr.write(
                "Warning: parsed species were absent from --species-tree: {}\n".format(
                    " ".join(absent)
                )
            )


def reconcile_main(args):
    outputs = [("--outfile", args.outfile)]
    validate_distinct_output_paths(outputs)
    validate_outputs_do_not_replace_inputs(
        [
            ("--infile", args.infile),
            ("--species-tree", args.species_tree),
            ("--species-map-tsv", getattr(args, "species_map_tsv", None)),
        ],
        outputs,
        label="Reconciliation output",
    )
    gene_tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    species_tree = read_tree(
        args.species_tree,
        args.species_tree_format,
        args.quoted_node_names,
        rooted=getattr(args, "species_tree_rooted", "auto"),
    )
    species_label_by_gene_leaf = _parsed_species_labels(gene_tree, args)
    _report_unmatched_species(
        species_label_by_gene_leaf, species_tree, policy=args.unmatched
    )
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        species_label_by_gene_leaf,
        event_source=args.event_source,
        tree_id=getattr(args, "tree_id", ""),
    )
    if args.outfile == "-":
        print(table.to_csv(sep="\t", index=False), end="")
    else:
        from nwkit.regression_pipeline import _write_dataframes_transactionally

        _write_dataframes_transactionally([(args.outfile, table)])
