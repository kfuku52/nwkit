"""Conservative fossil placement and auditable MCMCtree constraints."""

import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict

from nwkit.angiocal import fossil_id_text, positive_number
from nwkit.time_tree import parse_mcmctree_calibration
from nwkit.util import (
    get_ete_ncbitaxa,
    get_species_group_records,
    get_subtree_leaf_name_sets,
    warn_cleanup_failure,
)

REPORT_COLUMNS = (
    "dataset_version",
    "source",
    "source_sha256",
    "source_row",
    "fossil_id",
    "fossil_taxon",
    "clade",
    "placement",
    "node_calibrated",
    "minimum_age_ma",
    "safe_minimum_age",
    "age_quality_score",
    "node_assignment_score",
    "reconciliation_score",
    "relationship_reference",
    "age_reference",
    "status",
    "reason",
    "mapping_method",
    "node_clade_id",
    "node_tips",
    "selected_minimum_age_ma",
    "time_unit_ma",
    "constraint",
)


def validate_angiocal_tree(tree):
    if any(not node.is_leaf and len(node.children) != 2 for node in tree.traverse()):
        raise ValueError("AngioCal PAML output requires a fully bifurcating tree.")
    names = list(tree.leaf_names())
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("AngioCal requires nonempty, unique leaf names.")
    if any(
        len(name) > 100
        or not name.isascii()
        or not name.isprintable()
        or name.isdecimal()
        or any(char.isspace() or char in "()[]':;,#\"" for char in name)
        for name in names
    ):
        raise ValueError(
            "AngioCal PAML output requires unquoted ASCII tip identifiers of at most 100 characters without "
            "whitespace, control characters, Newick delimiters, quotes or '#', and not consisting solely of digits. "
            "Rename tips and alignment records together."
        )
    _normalize_existing_calibrations(tree)


def _normalize_existing_calibrations(tree):
    for node in tree.traverse():
        if node.is_leaf:
            continue
        text = str(node.name or "").strip().strip("'").strip('"')
        calibration = parse_mcmctree_calibration(text)
        if calibration is not None:
            _validate_fossil_prior(calibration)
            # NWKIT accepts lowercase constructors; stock PAML only recognizes
            # uppercase L/U/B. Keep all numeric text and prior parameters.
            raw = calibration["raw"]
            node.name = raw[:1].upper() + raw[1:]
        elif (
            text.startswith(("@", ">", "<"))
            or "#" in text
            or re.search(r"\b(?:L|U|B|G|SN|ST|S2N)\s*[({]", text, re.IGNORECASE)
        ):
            raise ValueError(
                f"Unsupported existing PAML calibration: {text!r}. AngioCal retains fully specified "
                "L(age, offset, scale, tail), U(age, tail), B(lower, upper, lower_tail, upper_tail) "
                "or legacy >/< bounds. Expand shorthand or remove unsupported priors before importing."
            )


def read_calibration_map(path, records, tree, args):
    """Each pair asserts the final calibrated divergence, including for stem fossils."""
    if not path:
        return {}
    valid_ids = {record.fossil_id for record in records}
    leaf_names = set(tree.leaf_names())
    _, groups, _ = get_species_group_records(tree, args=args)
    handle = sys.stdin if path == "-" else open(path, encoding="utf-8-sig", newline="")
    mapped = {}
    try:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        if (
            not all(headers)
            or len(set(headers)) != len(headers)
            or not {
                "fossil_id",
                "left_species",
                "right_species",
            }.issubset(headers)
        ):
            raise ValueError(
                "Calibration map needs unique fossil_id/left_species/right_species columns."
            )
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed calibration mapping row.")
            fossil_id = fossil_id_text(row["fossil_id"])
            if fossil_id not in valid_ids or fossil_id in mapped:
                raise ValueError(
                    f"Unknown or duplicate calibration-map fossil ID: {fossil_id}."
                )
            sides = []
            for key in ("left_species", "right_species"):
                label = row[key].strip()
                tips = {label} if label in leaf_names else set(groups.get(label, []))
                if not tips:
                    raise ValueError(f"Calibration-map species not found: {label}.")
                sides.append(tips)
            target = tree.common_ancestor(sorted(sides[0] | sides[1]))
            if sides[0] & sides[1] or target.is_leaf:
                raise ValueError(
                    "Calibration-map anchors must identify different lineages."
                )
            species = {
                label
                for label, tips in groups.items()
                if set(tips) & (sides[0] | sides[1])
            }
            if len(species) < 2:
                raise ValueError(
                    "Calibration-map anchors must identify different species."
                )
            mapped[fossil_id] = target
    finally:
        if handle is not sys.stdin:
            handle.close()
    return mapped


def _named_target(tree, record):
    # Internal labels explicitly assert a biological clade; ordinary tip MRCA
    # alone cannot establish the crown split when basal lineages are absent.
    final_labels = {f"{record.placement} {record.clade}", record.node_calibrated} - {""}
    internal = [node for node in tree.traverse() if not node.is_leaf]
    direct = [node for node in internal if node.name in final_labels]
    if len(direct) == 1:
        return direct[0], "", "node_label"
    if direct:
        return None, "ambiguous_node_label", "node_label"
    crowns = [node for node in internal if node.name == record.clade]
    if len(crowns) != 1:
        return None, "ambiguous_node_label" if crowns else "", "node_label"
    crown = crowns[0]
    if record.placement == "crown":
        return crown, "", "node_label"
    if crown.is_root:
        return None, "stem_outgroup_missing", "node_label"
    return crown.up, "", "node_label"


def _taxonomy_context(tree, ncbi, args):
    leaf_to_species, groups, queries = get_species_group_records(tree, args=args)
    names = ncbi.get_name_translator(sorted(set(queries.values())))
    lineages = {}
    for label, tips in groups.items():
        ids = names.get(queries[label], [])
        if len(ids) != 1:
            return leaf_to_species, {}, "unresolved_tip_taxonomy"
        lineage = set(ncbi.get_lineage(ids[0]))
        for tip in tips:
            lineages[tip] = lineage
    return leaf_to_species, lineages, ""


def _taxonomy_tip_clades(lineages):
    """Project represented NCBI clades onto the input tip set."""
    taxon_tips = defaultdict(set)
    for tip, lineage in lineages.items():
        for taxid in lineage:
            taxon_tips[taxid].add(tip)
    total = len(lineages)
    return {frozenset(tips) for tips in taxon_tips.values() if 1 < len(tips) < total}


def _stem_conflicts_with_taxonomy(target, leaf_sets, taxon_clades, cache):
    # A projected taxon cannot partly overlap either stem branch or their
    # parent. A conflict here can move a stem fossil onto a false ancestor.
    for node in (target, *target.children):
        if node not in cache:
            tips = leaf_sets[node]
            cache[node] = any(
                bool(tips & clade) and not (tips <= clade or clade <= tips)
                for clade in taxon_clades
            )
        if cache[node]:
            return True
    return False


def _taxonomy_stem_target(
    tree, record, ncbi, context, leaf_sets, taxon_clades, conflict_cache
):
    _, lineages, error = context
    if error:
        return None, error, "taxonomy_stem"
    components = [component.strip() for component in record.clade.split("+")]
    names = ncbi.get_name_translator(components)
    if any(len(names.get(component, [])) != 1 for component in components):
        return None, "unresolved_clade_taxonomy", "taxonomy_stem"
    ids = {names[component][0] for component in components}
    tips = {tip for tip, lineage in lineages.items() if ids & lineage}
    if not tips:
        return None, "clade_not_sampled", "taxonomy_stem"
    crown = tree.common_ancestor(sorted(tips))
    if leaf_sets[crown] != tips:
        return None, "nonmonophyletic_clade", "taxonomy_stem"
    if crown.is_root:
        return None, "stem_outgroup_missing", "taxonomy_stem"
    if _stem_conflicts_with_taxonomy(crown.up, leaf_sets, taxon_clades, conflict_cache):
        return None, "conflicting_tree_taxonomy", "taxonomy_stem"
    # The nearest sampled outside lineage may diverge earlier than the true
    # stem when its sister is missing. A minimum remains valid on that ancestor.
    return crown.up, "sampled_stem_ancestor", "taxonomy_stem"


def place_fossils(tree, dataset, args):
    leaf_to_species, _, _ = get_species_group_records(tree, args=args)
    explicit = read_calibration_map(
        getattr(args, "calibration_map_tsv", None), dataset.records, tree, args
    )
    leaf_sets = get_subtree_leaf_name_sets(tree)
    placements = {}
    pending = []
    for record in dataset.records:
        if record.fossil_id in explicit:
            placements[record.fossil_id] = (
                explicit[record.fossil_id],
                "",
                "explicit_map",
            )
            continue
        target, reason, method = _named_target(tree, record)
        if target is not None and record.placement == "crown":
            if len({leaf_to_species[tip] for tip in leaf_sets[target]}) < 2:
                target, reason = None, "crown_not_represented"
        if target is not None or reason:
            placements[record.fossil_id] = (target, reason, method)
        elif record.placement == "crown":
            placements[record.fossil_id] = (None, "crown_requires_anchors", "")
        elif not getattr(args, "angiocal_taxonomy", True):
            placements[record.fossil_id] = (None, "taxonomy_disabled", "")
        else:
            pending.append(record)
    if pending:
        ncbi = get_ete_ncbitaxa(args=args)
        try:
            context = _taxonomy_context(tree, ncbi, args)
            taxon_clades = _taxonomy_tip_clades(context[1])
            conflict_cache: dict[object, bool] = {}
            for record in pending:
                placements[record.fossil_id] = _taxonomy_stem_target(
                    tree,
                    record,
                    ncbi,
                    context,
                    leaf_sets,
                    taxon_clades,
                    conflict_cache,
                )
        finally:
            db = getattr(ncbi, "db", None)
            if db is not None:
                try:
                    db.close()
                except Exception as exc:
                    warn_cleanup_failure("NCBI taxonomy database handle", exc)
    return placements, leaf_sets


def _lower_constraint(age, args):
    from nwkit.mcmctree import _finite_number, _number_text, _tail_probability

    age = positive_number(age, "Converted AngioCal minimum age")
    offset = _finite_number(args.lower_offset, "'--lower-offset'", minimum=0)
    scale = positive_number(args.lower_scale, "'--lower-scale'")
    parameters = [
        repr(float(age)),
        _number_text(offset),
        _number_text(scale),
        _tail_probability(args, "lower"),
    ]
    return "L(" + ", ".join(parameters) + ")"


def _combine_existing(node, minimum, args):
    from nwkit.mcmctree import _number_text, _tail_probability

    lower_text = _lower_constraint(minimum, args)
    existing = parse_mcmctree_calibration(node.name)
    if existing is None:
        return lower_text
    _validate_fossil_prior(existing)
    if existing.get("upper", float("inf")) < minimum:
        raise ValueError(
            "Existing calibration upper/point age is younger than an AngioCal minimum."
        )
    if existing.get("lower", -1) >= minimum:
        return existing["raw"]
    if "lower" in existing:
        raise ValueError(
            "An existing weaker lower calibration has a separate prior. "
            "Remove that calibration before importing AngioCal."
        )
    if existing["upper"] == minimum:
        raise ValueError(
            "An existing upper-only age must exceed the AngioCal minimum to form a bounded prior."
        )
    return (
        "B("
        + ", ".join(
            [
                repr(float(minimum)),
                repr(float(existing["upper"])),
                _tail_probability(args, "lower"),
                _number_text(
                    existing.get("upper_tail", float(_tail_probability(args, "upper")))
                ),
            ]
        )
        + ")"
    )


def _validate_fossil_prior(calibration):
    if calibration.get("type") == "point":
        raise ValueError(
            "AngioCal cannot retain '@age' annotations: PAML does not treat them as "
            "fossil priors. Remove them or provide an L/U/B fossil prior."
        )
    if calibration.get("type") == "lower":
        age = calibration["lower"]
        location = age * (1.0 + calibration.get("offset", 0.1))
        scale = age * calibration.get("scale", 1.0)
        if age <= 0 or scale <= 0 or not all(map(math.isfinite, (location, scale))):
            raise ValueError(
                "A PAML lower fossil prior's age and derived location/scale must be finite and positive."
            )
    if calibration.get("type") == "bounded":
        if calibration["lower"] >= calibration["upper"]:
            raise ValueError(
                "A PAML bounded fossil prior needs a strictly positive age range."
            )
        if (
            1.0
            - calibration.get("lower_tail", 0.025)
            - calibration.get("upper_tail", 0.025)
            <= 0
        ):
            raise ValueError(
                "Bounded fossil prior tail probabilities must sum to less than 1."
            )


def validate_calibration_order(tree):
    """Check nominal age bounds throughout the tree without changing priors."""
    descendant_minima: dict[object, float] = {}
    for node in tree.traverse(strategy="postorder"):
        child_minimum = max(
            (descendant_minima[child] for child in node.children), default=0.0
        )
        calibration = (
            (parse_mcmctree_calibration(node.name) or {}) if not node.is_leaf else {}
        )
        _validate_fossil_prior(calibration)
        lower = max(child_minimum, calibration.get("lower", 0.0))
        upper = calibration.get("upper", float("inf"))
        if lower > upper or (node.children and child_minimum >= upper):
            raise ValueError(
                "Nominal calibration bounds conflict between an ancestor and its descendants."
            )
        descendant_minima[node] = lower


def _report_row(record, dataset, placement, leaf_sets, unit):
    node, reason, method = placement
    row = {key: getattr(record, key, "") for key in REPORT_COLUMNS}
    row.update(
        dataset_version=dataset.version,
        source=dataset.source,
        source_sha256=dataset.sha256,
        status="candidate" if node is not None else "skipped",
        reason=reason,
        mapping_method=method,
        time_unit_ma=unit,
    )
    if node is not None:
        tips = json.dumps(
            sorted(leaf_sets[node]), ensure_ascii=False, separators=(",", ":")
        )
        row.update(
            node_tips=tips, node_clade_id=hashlib.sha256(tips.encode()).hexdigest()
        )
    return row


def add_angiocal_constraints(tree, dataset, args):
    unit = positive_number(getattr(args, "time_unit_ma", 1), "'--time-unit-ma'")
    # Validate prior parameters even if every fossil is skipped.
    _lower_constraint(1, args)
    placements, leaf_sets = place_fossils(tree, dataset, args)
    rows = []
    by_node = defaultdict(list)
    tree_size = len(leaf_sets[tree])
    for record in dataset.records:
        placement = placements[record.fossil_id]
        row = _report_row(record, dataset, placement, leaf_sets, unit)
        rows.append(row)
        node = placement[0]
        if node is None:
            continue
        if not node.is_root and len(leaf_sets[node]) < args.min_clade_prop * tree_size:
            row.update(status="skipped", reason="min_clade_prop")
            continue
        by_node[node].append((record, row))
    for node, entries in by_node.items():
        minimum_ma = max(record.minimum_age_ma for record, _ in entries)
        constraint = _combine_existing(node, minimum_ma / unit, args)
        node.name = "'" + constraint + "'"
        for record, row in entries:
            row.update(
                status="applied"
                if record.minimum_age_ma == minimum_ma
                else "supporting",
                selected_minimum_age_ma=minimum_ma,
                constraint=constraint,
            )
    validate_calibration_order(tree)
    counts = Counter(row["status"] for row in rows)
    sys.stderr.write(
        f"AngioCal {dataset.version}: {len(by_node)} calibrated node(s), "
        f"{counts['skipped']} skipped fossil(s).\n"
    )
    return rows, len(by_node)


def write_angiocal_report(handle, rows):
    writer = csv.DictWriter(
        handle, fieldnames=REPORT_COLUMNS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
