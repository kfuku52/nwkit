"""Reconciliation adapters and shared-event chronology for RADTE.

All input routes produce one event per gene clade. Species event IDs, rather
than traversal numbers, identify shared ages. No reconciliation executable runs.
"""

from dataclasses import dataclass
from graphlib import CycleError, TopologicalSorter

import numpy as np
import pandas as pd

from nwkit.clade_index import CladeIndex
from nwkit.reconcile import (
    _parsed_species_labels,
    _report_unmatched_species,
    _validate_rooted_binary_tree,
    build_reconciliation_table,
)
from nwkit.util import read_tree


@dataclass
class Chronology:
    gene: object
    species: object
    nodes: list
    edges: list
    groups: list[str]
    group_by_node: np.ndarray
    parent: np.ndarray
    child: np.ndarray
    constraint_parent: np.ndarray
    constraint_child: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    initial: np.ndarray
    scale: float
    events: pd.DataFrame
    species_table: pd.DataFrame
    min_duration: float

    def durations(self, ages):
        return ages[self.parent] - ages[self.child]


def _read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _unique_names(tree):
    result = {}
    for node in tree.traverse():
        name = str(node.name or "")
        if name:
            if name in result:
                raise ValueError(f"Ambiguous node label: {name}")
            result[name] = node
    return result


def species_ages(tree, bounds_path=None):
    _validate_rooted_binary_tree(tree, "--species-tree")
    depth = {tree: 0.0}
    for node in tree.traverse():
        if node is tree:
            continue
        if node.dist is None or not np.isfinite(node.dist) or node.dist <= 0:
            raise ValueError("Species tree requires finite, positive branch lengths.")
        depth[node] = depth[node.up] + node.dist
    height = max(depth.values())
    if max(abs(depth[n] - height) for n in tree.leaves()) > 1e-7 * height:
        raise ValueError(
            "Species tree must be ultrametric; ages are never silently repaired."
        )
    index = CladeIndex(tree)
    rows = []
    for node in tree.traverse():
        age = 0.0 if node.is_leaf else height - depth[node]
        rows.append(
            dict(
                species_event_id=index.clade_id_for_node(node),
                node=str(node.name or ""),
                age=age,
                age_min=age,
                age_max=age,
            )
        )
    frame = pd.DataFrame(rows)
    if bounds_path:
        bounds = _read_tsv(bounds_path)
        if not {"node", "age_min", "age_max"}.issubset(bounds.columns):
            raise ValueError("Species bounds require node, age_min, age_max columns.")
        names = _unique_names(tree)
        if bounds.node.duplicated().any():
            raise ValueError("Duplicate species bounds node.")
        for record in bounds.to_dict("records"):
            if record["node"] not in names:
                raise ValueError(f"Unknown species bounds node: {record['node']}")
            lo, hi = float(record["age_min"]), float(record["age_max"])
            if not np.isfinite([lo, hi]).all() or not 0 <= lo <= hi:
                raise ValueError("Species bounds must be finite with 0 <= min <= max.")
            node = names[record["node"]]
            if node.is_leaf and (lo != 0 or hi != 0):
                raise ValueError(
                    "RADTE currently requires contemporaneous tips at age zero."
                )
            frame.loc[
                frame.species_event_id == index.clade_id_for_node(node),
                ["age_min", "age_max"],
            ] = [lo, hi]
    return frame, index, height


def _notung_events(path, gene, species, table):
    names, species_names = _unique_names(gene), _unique_names(species)
    index, sp_index = CladeIndex(gene), CladeIndex(species)
    records = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].upper() != "#D":
                continue
            if len(fields) >= 2 and fields[1].upper() == "DUPLICATION":
                continue
            if len(fields) >= 3 and [f.upper() for f in fields[1:3]] == [
                "GENE",
                "NODE",
            ]:
                continue
            if len(fields) != 4:
                raise ValueError(
                    "Malformed Notung #D record; expected node, lower, upper."
                )
            _, name, lower, upper = fields
            if name not in names or names[name].is_leaf:
                raise ValueError(f"Unknown/noninternal Notung duplication: {name}")
            if lower not in species_names:
                raise ValueError(f"Unknown Notung lower species node: {lower}")
            low_node = species_names[lower]
            if upper not in species_names:
                if low_node is not species or upper.upper() not in {
                    "-",
                    "NA",
                    "N/A",
                    "NONE",
                    "NULL",
                    "NIL",
                }:
                    raise ValueError(f"Unknown Notung upper species node: {upper}")
                up_id = ""
            else:
                up_node = species_names[upper]
                if up_node not in set(low_node.ancestors()):
                    raise ValueError(
                        "Notung upper node must be an ancestor of its lower node."
                    )
                up_id = sp_index.clade_id_for_node(up_node)
            records.append(
                (
                    index.clade_id_for_node(names[name]),
                    sp_index.clade_id_for_node(low_node),
                    up_id,
                )
            )
    if len({row[0] for row in records}) != len(records):
        raise ValueError("Duplicate Notung duplication record.")
    table = table.copy()
    # The Notung #D list is authoritative, including an empty list.
    table.loc[table.event_type != "leaf", "event_type"] = "speciation"
    table["upper_species_event_id"] = ""
    for gene_id, low_id, up_id in records:
        mask = table.gene_clade_id == gene_id
        table.loc[
            mask, ["event_type", "species_event_id", "upper_species_event_id"]
        ] = ["duplication", low_id, up_id]
    table["event_source"] = "notung"
    return table


def remap_generax_species(gene, species, source_path):
    """Translate GeneRax S labels by topology, never borrow its branch ages."""
    source = read_tree(source_path, 1, True, rooted="auto")
    _validate_rooted_binary_tree(source, "--reconciliation-species-tree")
    source_names = _unique_names(source)
    _unique_names(species)
    old, new = CladeIndex(source), CladeIndex(species)
    target = {new.clade_id_for_node(n): n for n in species.traverse()}
    if {old.clade_id_for_node(n) for n in source.traverse()} != set(target):
        raise ValueError(
            "Reconciliation and dated species trees must have identical rooted clades."
        )
    for node in gene.traverse():
        label = str(node.props.get("S", ""))
        if label not in source_names:
            raise ValueError(f"Unknown GeneRax species annotation: {label}")
        sid = old.clade_id_for_node(source_names[label])
        destination = target[sid]
        if not destination.name:
            destination.name = sid
        node.props["S"] = str(destination.name)
    _unique_names(species)


def read_inputs(args):
    if bool(args.generax_nhx) == bool(args.gene_tree):
        raise ValueError("Provide exactly one of --generax-nhx or --gene-tree.")
    modes = [
        bool(args.notung_parsable),
        bool(args.reconciliation),
        bool(args.reconcile),
    ]
    if args.generax_nhx and any(modes):
        raise ValueError(
            "GeneRax NHX cannot be combined with another reconciliation source."
        )
    if args.gene_tree and sum(modes) != 1:
        raise ValueError(
            "--gene-tree requires exactly one of --notung-parsable, --reconciliation, --reconcile lca."
        )
    gene = read_tree(
        args.generax_nhx or args.gene_tree,
        args.gene_tree_format,
        True,
        rooted=args.input_rooted,
    )
    species = read_tree(
        args.species_tree,
        args.species_tree_format,
        True,
        rooted=args.species_tree_rooted,
    )
    _validate_rooted_binary_tree(gene, "--gene-tree")
    _validate_rooted_binary_tree(species, "--species-tree")
    source_species = getattr(args, "reconciliation_species_tree", None)
    if source_species:
        if not args.generax_nhx:
            raise ValueError("--reconciliation-species-tree requires --generax-nhx.")
        remap_generax_species(gene, species, source_species)
    if args.generax_nhx:
        # Never replace an absent or invalid GeneRax S annotation by LCA.
        if any(str(node.props.get("S", "")) == "" for node in gene.traverse()):
            raise ValueError("GeneRax requires an S annotation on every node.")
        for node in gene.traverse():
            duplication = str(node.props.get("D", "N")).upper()
            if duplication in {"Y", "YES", "TRUE", "T", "1"}:
                node.props["D"] = "Y"
            elif duplication in {"N", "NO", "FALSE", "F", "0", ""}:
                node.props["D"] = "N"
            else:
                raise ValueError(
                    f"Unsupported GeneRax duplication annotation: {duplication}"
                )
            transfer = str(node.props.get("H", "N")).upper()
            if transfer not in {"N", "NO", "FALSE", "F", "0", ""}:
                raise ValueError(
                    "RADTE does not date transfer or unresolved GeneRax H events."
                )
            node.props["H"] = "N"
        labels = {str(n.name): str(n.props["S"]) for n in gene.leaves()}
    else:
        labels = _parsed_species_labels(gene, args)
    _report_unmatched_species(labels, species, "error")
    if args.reconciliation:
        table = _read_tsv(args.reconciliation)
    else:
        table = build_reconciliation_table(
            gene, species, labels, event_source="nhx" if args.generax_nhx else "lca"
        )
    if args.notung_parsable:
        table = _notung_events(args.notung_parsable, gene, species, table)
    if {"gene_clade_id", "species_event_id"}.issubset(table.columns):
        gene_index, species_index = CladeIndex(gene), CladeIndex(species)
        placements = table.set_index("gene_clade_id")["species_event_id"].to_dict()
        species_leaves = {str(n.name): n for n in species.leaves()}
        for tip in gene.leaves():
            if placements.get(
                gene_index.clade_id_for_node(tip)
            ) != species_index.clade_id_for_node(species_leaves[labels[str(tip.name)]]):
                raise ValueError(
                    "Reconciliation tip placements disagree with the gene-to-species mapping."
                )
    chronology = build_chronology(
        gene, species, table, args.species_node_bounds_tsv, args.max_age
    )
    from nwkit.radte_species import attach_species_intervals

    attach_species_intervals(
        chronology, getattr(args, "species_node_intervals_tsv", None)
    )
    return chronology


def _event_bounds(node, gid, rec, sp_by_id, bounds, sp_index, max_age):
    extra_constraints = []
    sid, event = rec["species_event_id"], rec["event_type"]
    if sid not in sp_by_id or event not in {"leaf", "speciation", "duplication"}:
        raise ValueError(
            "RADTE requires resolved duplication/loss events; transfer/unresolved events are unsupported."
        )
    if rec.get("collapsed_event_boundary", "no") not in {"", "no"}:
        raise ValueError("Date the complete gene tree before pruning collapsed events.")
    if node.is_leaf != (event == "leaf"):
        raise ValueError("Reconciliation leaf/internal event mismatch.")
    if node.is_leaf and not sp_by_id[sid].is_leaf:
        raise ValueError("A gene tip must map to a species tip.")
    if not node.is_leaf and event == "speciation" and sp_by_id[sid].is_leaf:
        raise ValueError("A speciation event cannot map to a species tip.")
    key = "D:" + gid if event == "duplication" else "S:" + sid
    lower, upper = bounds[sid]["age_min"], bounds[sid]["age_max"]
    upper_id = sid
    if event == "duplication":
        parent = sp_by_id[sid].up
        upper_id = rec.get("upper_species_event_id", "")
        if not upper_id and parent is not None:
            upper_id = sp_index.clade_id_for_node(parent)
        if upper_id:
            if upper_id not in bounds or sp_by_id[upper_id] not in set(
                sp_by_id[sid].ancestors()
            ):
                raise ValueError("Invalid duplication upper species bound.")
            upper = bounds[upper_id]["age_max"]
        else:
            if max_age is None or not np.isfinite(max_age) or max_age <= lower:
                raise ValueError(
                    "Duplication above the species root requires --max-age greater than its lower bound."
                )
            upper = max_age
        extra_constraints.append((key, "S:" + sid))
        if upper_id:
            extra_constraints.append(("S:" + upper_id, key))
    if node.is_leaf:
        lower = upper = 0.0
    return sid, event, key, lower, upper, upper_id, extra_constraints


def _validate_species_placements(gene, by_id, gene_index, sp_by_id, sp_index):
    observed = {}
    for node in gene.traverse(strategy="postorder"):
        record = by_id[gene_index.clade_id_for_node(node)]
        sid = record["species_event_id"]
        if sid not in sp_by_id:
            raise ValueError("Unknown species event in reconciliation.")
        mapped_mask = sp_index.mask_by_node[sp_by_id[sid]]
        if node.is_leaf:
            observed[node] = mapped_mask
            continue
        first, second = node.children
        mask = observed[first] | observed[second]
        observed[node] = mask
        if mask & mapped_mask != mask:
            raise ValueError(
                "Gene descendants fall outside the reconciled species lineage."
            )
        if record["event_type"] == "speciation" and not sp_by_id[sid].is_leaf:
            daughters = [sp_index.mask_by_node[n] for n in sp_by_id[sid].children]
            orientations = [(daughters[0], daughters[1]), (daughters[1], daughters[0])]
            if not any(
                observed[first] & a == observed[first]
                and observed[second] & b == observed[second]
                for a, b in orientations
            ):
                raise ValueError(
                    "Shared speciation mapping has overlapping descendant lineages or a zero-duration path."
                )


def build_chronology(gene, species, events, bounds_path=None, max_age=None):
    _validate_rooted_binary_tree(gene, "--gene-tree")
    sp_table, sp_index, scale = species_ages(species, bounds_path)
    gene_index = CladeIndex(gene)
    nodes = list(gene.traverse())
    required = {
        "gene_clade_id",
        "parent_gene_clade_id",
        "species_event_id",
        "event_type",
    }
    if not required.issubset(events.columns):
        raise ValueError(f"Reconciliation requires columns: {sorted(required)}")
    if events.gene_clade_id.duplicated().any():
        raise ValueError("Duplicate gene clade in reconciliation.")
    by_id = events.set_index("gene_clade_id").to_dict("index")
    if set(by_id) != {gene_index.clade_id_for_node(n) for n in nodes}:
        raise ValueError("Reconciliation must cover exactly the supplied gene tree.")
    sp_by_id = {sp_index.clade_id_for_node(n): n for n in species.traverse()}
    _validate_species_placements(gene, by_id, gene_index, sp_by_id, sp_index)
    bounds = sp_table.set_index("species_event_id").to_dict("index")
    groups = ["S:" + sid for sid in bounds]
    node_keys, rows = [], []
    lo_by_key = {"S:" + sid: rec["age_min"] / scale for sid, rec in bounds.items()}
    hi_by_key = {"S:" + sid: rec["age_max"] / scale for sid, rec in bounds.items()}
    extra_constraints = [
        ("S:" + sp_index.clade_id_for_node(n.up), "S:" + sid)
        for sid, n in sp_by_id.items()
        if n.up is not None
    ]
    for node in nodes:
        gid = gene_index.clade_id_for_node(node)
        rec = by_id[gid]
        expected_parent = "" if node is gene else gene_index.clade_id_for_node(node.up)
        if rec["parent_gene_clade_id"] != expected_parent:
            raise ValueError("Reconciliation parent clades do not match gene topology.")
        sid, event, key, lower, upper, upper_id, constraints = _event_bounds(
            node, gid, rec, sp_by_id, bounds, sp_index, max_age
        )
        extra_constraints.extend(constraints)
        if key not in lo_by_key:
            groups.append(key)
            lo_by_key[key], hi_by_key[key] = lower / scale, upper / scale
        node_keys.append(key)
        rows.append(
            dict(
                gene_clade_id=gid,
                gene_name=str(node.name or ""),
                parent_gene_clade_id=expected_parent,
                event_type=event,
                event_source=rec.get("event_source", "precomputed"),
                species_event_id=sid,
                species_name=bounds[sid]["node"],
                upper_species_event_id=upper_id,
                shared_age_id=key,
                age_min=lower,
                age_max=upper,
            )
        )
    group_index = {key: i for i, key in enumerate(groups)}
    group_by_node = np.array([group_index[k] for k in node_keys])
    node_index = {n: i for i, n in enumerate(nodes)}
    edges = [n for n in nodes if n is not gene]
    parent = np.array([group_by_node[node_index[n.up]] for n in edges])
    child = np.array([group_by_node[node_index[n]] for n in edges])
    all_constraints = set(zip(parent.tolist(), child.tolist(), strict=True))
    all_constraints.update(
        (group_index[p], group_index[c]) for p, c in extra_constraints
    )
    cp, cc = np.array(sorted(all_constraints), dtype=int).T
    lower = np.array([lo_by_key[k] for k in groups])
    upper = np.array([hi_by_key[k] for k in groups])
    margin = 1e-10
    predecessors: dict[int, set[int]] = {i: set() for i in range(len(groups))}
    for p, c in zip(cp, cc, strict=True):
        predecessors[int(p)].add(int(c))
    try:
        order = list(TopologicalSorter(predecessors).static_order())
    except CycleError as exc:
        raise ValueError(
            "Shared speciation ages create a cycle or zero-duration path."
        ) from exc
    effective_lo, effective_hi = lower.copy(), upper.copy()
    for p in order:
        for c in predecessors[p]:
            effective_lo[p] = max(effective_lo[p], effective_lo[c] + margin)
    for p in reversed(order):
        for c in predecessors[p]:
            effective_hi[c] = min(effective_hi[c], effective_hi[p] - margin)
    if np.any(effective_lo > effective_hi):
        raise ValueError(
            "Infeasible shared ages or calibration bounds with positive durations."
        )
    initial = effective_lo.copy()
    for p in reversed(order):
        initial[p] = (effective_lo[p] + effective_hi[p]) / 2
        for c in predecessors[p]:
            effective_hi[c] = min(effective_hi[c], initial[p] - margin)
    for node in edges:
        if node.dist is None or not np.isfinite(node.dist) or node.dist <= 0:
            raise ValueError(
                "Gene tree requires finite, positive substitution branch lengths."
            )
    return Chronology(
        gene,
        species,
        nodes,
        edges,
        groups,
        group_by_node,
        parent,
        child,
        cp,
        cc,
        lower,
        upper,
        initial,
        scale,
        pd.DataFrame(rows),
        sp_table,
        margin,
    )
