"""Species-tree Ks correction, multiple-outgroup diagnostics and outputs."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit import __version__
from nwkit.clade_index import CladeIndex
from nwkit.gene_family_input import read_gene_family_table
from nwkit.ksrate_model import (
    KsObservation,
    KsTrio,
    PairKs,
    bootstrap_corrections,
    grouped_corrections,
    simultaneous_median_corrections,
    trio_values,
)
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.rooting_state import require_rooted
from nwkit.util import (
    assign_branch_ids,
    read_tree,
    validate_outputs_do_not_replace_inputs,
    validate_unique_named_leaves,
)


def make_trios(tree, focal_names=None, outgroup_policy="nearest", max_trios=100000):
    require_rooted(tree, "Ks correction requires a rooted species tree.")
    validate_unique_named_leaves(tree, "--infile")
    if any(len(node.children) == 1 for node in tree.traverse()):
        raise ValueError("Ks correction does not accept unary species-tree nodes.")
    clades = CladeIndex(tree)
    branches = assign_branch_ids(tree)
    leaves = {leaf.name: leaf for leaf in tree.leaves()}
    if len(leaves) < 2:
        raise ValueError("Ks correction needs at least two species-tree tips.")
    focal_names = sorted(leaves) if focal_names is None else sorted(focal_names)
    if len(set(focal_names)) != len(focal_names) or any(
        name not in leaves for name in focal_names
    ):
        raise ValueError("Focal species must be unique species-tree tips.")
    if outgroup_policy not in {"nearest", "all"}:
        raise ValueError("Outgroup policy must be nearest or all.")
    trios: list[KsTrio] = []
    events = []
    for focal in focal_names:
        child = leaves[focal]
        while child.up is not None:
            node = child.up
            sisters = sorted(
                set(clades.names_for_node(node)) - set(clades.names_for_node(child))
            )
            if node.is_root:
                outgroups = []
            else:
                comparison = node.up if outgroup_policy == "nearest" else tree
                outgroups = sorted(
                    set(clades.names_for_node(comparison))
                    - set(clades.names_for_node(node))
                )
            event_id = clades.clade_id_for_node(node)
            indices = []
            for sister in sisters:
                for outgroup in outgroups:
                    indices.append(len(trios))
                    trios.append(
                        KsTrio(focal, sister, outgroup, event_id, branches[node])
                    )
                    if len(trios) > max_trios:
                        raise ValueError(
                            "Too many Ks trios; restrict --focals or use nearest outgroups."
                        )
            events.append(
                {
                    "focal": focal,
                    "species_event_id": event_id,
                    "branch_id": branches[node],
                    "descendant_taxa": clades.csv_for_node(node),
                    "trio_indices": indices,
                }
            )
            child = node
    return trios, events


def _paths(args):
    paths = {}
    for role in ("outfile", "trios_out", "model_out"):
        path = getattr(args, role)
        if path == "" or (path == "-" and role != "outfile"):
            raise ValueError(
                f"--{role.replace('_', '-')} requires a nonempty file path."
            )
        if path not in (None, "-"):
            paths[role] = path
    validate_outputs_do_not_replace_inputs(
        [("--infile", args.infile), ("--ks-tsv", args.ks_tsv)],
        [(role, path) for role, path in paths.items()],
    )
    validate_output_targets(paths.values())
    return paths


def ksrate_main(args):
    paths = _paths(args)
    if (
        args.bootstrap < 0
        or args.seed < 0
        or not np.isfinite(args.ci_level)
        or not 0 < args.ci_level < 1
    ):
        raise ValueError("Bootstrap/seed must be nonnegative and ci-level in (0, 1).")
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    trios, events = make_trios(
        tree,
        None if args.focals is None else args.focals.split(","),
        args.outgroup_policy,
    )
    _, rows = read_gene_family_table(
        args.ks_tsv, ["species_a", "species_b", "family_id", "ks"]
    )
    observations = [
        KsObservation(
            row["species_a"], row["species_b"], row["family_id"], float(row["ks"])
        )
        for row in rows
    ]
    pairs = PairKs(observations, [leaf.name for leaf in tree.leaves()])
    estimates = pairs.estimates()
    values = trio_values(trios, estimates)
    corrections = grouped_corrections(trios, values)
    bootstrap = bootstrap_corrections(pairs, trios, args.bootstrap, args.seed)
    simultaneous, num_simultaneous_pairs = (
        simultaneous_median_corrections(pairs, trios, args.ci_level)
        if args.ci_method == "pair-median-bonferroni"
        else ({}, 0)
    )
    trio_rows = []
    for trio, value in zip(trios, values, strict=True):
        trio_rows.append(
            {
                "focal": trio.focal,
                "sister": trio.sister,
                "outgroup": trio.outgroup,
                "species_event_id": trio.species_event_id,
                "branch_id": trio.branch_id,
                "ks_focal_sister": estimates.get(
                    tuple(sorted((trio.focal, trio.sister))), np.nan
                ),
                "ks_focal_outgroup": estimates.get(
                    tuple(sorted((trio.focal, trio.outgroup))), np.nan
                ),
                "ks_sister_outgroup": estimates.get(
                    tuple(sorted((trio.sister, trio.outgroup))), np.nan
                ),
                "num_families_focal_sister": pairs.count(trio.focal, trio.sister),
                "num_families_focal_outgroup": pairs.count(trio.focal, trio.outgroup),
                "num_families_sister_outgroup": pairs.count(trio.sister, trio.outgroup),
                "raw_corrected_ks": value,
                "status": "missing_comparison"
                if not np.isfinite(value)
                else "negative_correction"
                if value < 0
                else "ok",
            }
        )
    output = []
    previous_by_focal: dict[str, float] = {}
    for event in events:
        key = event["focal"], event["species_event_id"]
        selected = values[event["trio_indices"]]
        complete = selected[np.isfinite(selected)]
        raw = corrections.get(key, np.nan)
        status = "ok"
        if not event["trio_indices"]:
            status = "no_external_outgroup"
        elif not complete.size:
            status = "missing_comparisons"
        elif raw < 0:
            status = "negative_correction"
        elif np.any(complete < 0):
            status = "inconsistent_trios"
        elif complete.size != len(selected):
            status = "partial_comparisons"
        corrected = raw if np.isfinite(raw) and raw >= 0 else np.nan
        previous = previous_by_focal.get(event["focal"])
        monotone = (
            "not_comparable"
            if previous is None or not np.isfinite(corrected)
            else "yes"
            if corrected >= previous
            else "no"
        )
        if np.isfinite(corrected):
            previous_by_focal[event["focal"]] = corrected
        draws = bootstrap.get(key, np.full(args.bootstrap, np.nan))
        estimable = int(np.isfinite(draws).sum())
        bootstrap_status = (
            "not_run"
            if args.bootstrap == 0
            else "ok"
            if estimable == args.bootstrap
            else "unavailable_missing_bootstrap_comparisons"
        )
        bootstrap_low = bootstrap_high = np.nan
        if args.bootstrap and estimable == args.bootstrap:
            bootstrap_low, bootstrap_high = np.quantile(
                draws, [(1 - args.ci_level) / 2, (1 + args.ci_level) / 2]
            )
        if args.ci_method == "pair-median-bonferroni":
            low, high = simultaneous.get(key, (np.nan, np.nan))
            interval_status = (
                "ok"
                if np.isfinite(low) and np.isfinite(high)
                else "unavailable_insufficient_pair_families"
                if complete.size
                else "unavailable_missing_comparisons"
            )
        else:
            low, high = bootstrap_low, bootstrap_high
            interval_status = bootstrap_status
        output.append(
            {
                **{
                    name: event[name]
                    for name in (
                        "focal",
                        "species_event_id",
                        "branch_id",
                        "descendant_taxa",
                    )
                },
                "raw_corrected_ks": raw,
                "corrected_ks": corrected,
                "ci_lower": low,
                "ci_upper": high,
                "ci_level": args.ci_level,
                "ci_method": args.ci_method
                if args.ci_method == "pair-median-bonferroni" or args.bootstrap
                else "not-run",
                "num_simultaneous_pairs": num_simultaneous_pairs,
                "bootstrap_ci_lower": bootstrap_low,
                "bootstrap_ci_upper": bootstrap_high,
                "bootstrap_interval_status": bootstrap_status,
                "num_bootstrap": args.bootstrap,
                "num_bootstrap_estimable": estimable,
                "interval_status": interval_status,
                "num_trios": len(selected),
                "num_complete_trios": int(complete.size),
                "num_negative_trios": int((complete < 0).sum()),
                "monotone_from_younger_node": monotone,
                "status": status,
            }
        )
    serialized = {
        "outfile": pd.DataFrame(output).to_csv(
            sep="\t", index=False, float_format="%.12g", na_rep="NA"
        ),
        "trios_out": pd.DataFrame(
            trio_rows,
            columns=[
                "focal",
                "sister",
                "outgroup",
                "species_event_id",
                "branch_id",
                "ks_focal_sister",
                "ks_focal_outgroup",
                "ks_sister_outgroup",
                "num_families_focal_sister",
                "num_families_focal_outgroup",
                "num_families_sister_outgroup",
                "raw_corrected_ks",
                "status",
            ],
        ).to_csv(sep="\t", index=False, float_format="%.12g", na_rep="NA"),
    }
    metadata = {
        "schema_version": 1,
        "nwkit_version": __version__,
        "method": "focal-Ks-trio-median-family-bootstrap-v1",
        "formula": "Ks(F,S)+Ks(F,O)-Ks(S,O)",
        "interpretation": "twice the focal lineage synonymous distance from the focal/sister ancestor; not absolute age",
        "pair_estimator": "median, one representative observation per unordered species-pair/family",
        "trio_estimator": "median across complete sister/outgroup combinations; negative raw values retained for diagnostics",
        "outgroup_policy": args.outgroup_policy,
        "num_families": len(pairs.families),
        "num_pairs": len(pairs.pairs),
        "num_trios": len(trios),
        "bootstrap": args.bootstrap,
        "seed": args.seed,
        "ci_level": args.ci_level,
        "ci_method": args.ci_method,
        "num_simultaneous_pairs": num_simultaneous_pairs,
        "bootstrap_unit": "global shared family IDs; a family has one jointly resampled multiplicity across all species pairs",
        "uncertainty": "Intervals are conditional on the input tree, family membership and selected pairs; not posterior probability, a genome-event age interval, or a guarantee under synonymous saturation. Family-bootstrap percentile intervals are approximate and can under-cover at finite sample sizes. Pair-median Bonferroni intervals simultaneously bound compared population medians and propagate through trio/node medians; their finite-sample guarantee assumes iid independent families within each pair, not independence across species pairs.",
        "limitations": [
            "Root nodes have no external outgroup and remain unresolved.",
            "Nonmonotone/negative corrections are diagnosed, not clipped or projected.",
            "Between-pair distances are medians, not a fitted additive substitution tree; estimator choice can affect corrections.",
            "Missing bootstrap comparisons make the entire interval unavailable; incomplete draws are not silently dropped.",
            "Simultaneous pair-median bounds can be wide or unavailable when too few independent families give finite binomial order-statistic bounds; correlated/nonrepresentative families violate their guarantee.",
            "Input KS method, genetic code, alignment/orthology quality, saturation and gene conversion must be checked by the producer.",
        ],
    }
    serialized["model_out"] = (
        json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    with output_transaction(paths.values()) as staged:
        for role, path in paths.items():
            Path(staged[path]).write_text(serialized[role], encoding="utf-8")
    if args.outfile == "-":
        sys.stdout.write(serialized["outfile"])
    print(
        f"Ks correction: {len(pairs.families)} families, {len(trios)} trios; root nodes unresolved without external outgroups.",
        file=sys.stderr,
    )
