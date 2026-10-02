"""Plug-in fixed-gene-tree likelihoods using an inspectable native count fit."""

import json
import sys
from pathlib import Path

import pandas as pd

from nwkit.clade_index import CladeIndex
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.rooting_state import require_rooted
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    assign_branch_ids,
    read_tree,
    validate_outputs_do_not_replace_inputs,
    validate_unique_named_leaves,
)
from nwkit.wgd_count import _count_tree
from nwkit.wgd_count_model import MultiplicationEvent
from nwkit.wgd_tree_model import GeneTopology, TopologyLikelihood


def gene_topology(tree, parser):
    require_rooted(tree, "Fixed-topology likelihood requires a rooted gene tree.")
    validate_unique_named_leaves(tree, "--infile")
    nodes = tuple(tree.traverse("postorder"))
    indices = {node: i for i, node in enumerate(nodes)}
    clades = CladeIndex(tree)
    design = GeneTopology(
        tuple(tuple(indices[child] for child in node.children) for node in nodes),
        tuple(
            parser.parse(node.name).species_label or "" if node.is_leaf else ""
            for node in nodes
        ),
        tuple(clades.clade_id_for_node(node) for node in nodes),
    )
    return design, nodes


def model_mapping(metadata, design):
    expected = metadata.get("species_tree")
    if expected is None:
        raise ValueError(
            "Count model lacks full species-tree identity; regenerate the native fit."
        )
    source_ids = metadata["species_event_ids"]
    if len(source_ids) != len(set(source_ids)) or set(source_ids) != set(
        design.clade_ids
    ):
        raise ValueError("Count model and full species tree have different clades.")
    mapping = {i: design.clade_ids.index(clade) for i, clade in enumerate(source_ids)}
    if list(design.tip_names) != metadata["tip_names"]:
        raise ValueError("Count model and species tree have different tips.")
    groups = [0] * len(design.parents)
    for old, new in mapping.items():
        parent = expected["parents"][old]
        if (-1 if parent == -1 else mapping[parent]) != design.parents[new] or float(
            expected["lengths"][old]
        ) != design.lengths[new]:
            raise ValueError(
                "Count model and species tree have different branch lengths or parents."
            )
        groups[new] = metadata["branch_rate_groups"][old]
    return mapping, tuple(groups)


def wgd_tree_main(args):
    paths = {
        role: getattr(args, role)
        for role in ("outfile", "model_out")
        if getattr(args, role) not in (None, "-")
    }
    if any(value == "" for value in paths.values()) or args.model_out == "-":
        raise ValueError(
            "Likelihood output paths must be nonempty; only the primary table accepts stdout."
        )
    validate_outputs_do_not_replace_inputs(
        [
            ("--infile", args.infile),
            ("--species-tree", args.species_tree),
            ("--count-model", args.count_model),
            ("--species-map-tsv", args.species_map_tsv),
        ],
        list(paths.items()),
    )
    validate_output_targets(paths.values())
    species = read_tree(
        args.species_tree,
        "auto",
        args.quoted_node_names,
        rooted=args.species_tree_rooted,
    )
    design, _, _ = _count_tree(species)
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    if any(
        str(node.props.get("H", "")).upper().split("@")[0] in {"Y", "YES", "1", "TRUE"}
        for node in tree.traverse()
    ):
        raise ValueError(
            "Fixed DL/WGD topology likelihood does not model annotated transfers."
        )
    gene, nodes = gene_topology(tree, get_species_parser(args=args))
    metadata = json.loads(
        sys.stdin.read()
        if args.count_model == "-"
        else Path(args.count_model).read_text()
    )
    mapping, groups = model_mapping(metadata, design)
    model = TopologyLikelihood(
        design,
        gene,
        detection=metadata["detection_probabilities"],
        branch_groups=groups,
        rate_scales=metadata["family_rate_scales"],
        ascertainment=args.ascertainment,
        max_tips=args.max_tips,
    )
    options = {"tolerance": args.tolerance, "origin_tolerance": args.origin_tolerance}
    background = metadata["background"]
    background_result = model.evaluate(
        background["rates"], background["root_mean"], **options
    )
    requested = None if args.event_ids is None else args.event_ids.split(",")
    if requested is not None and len(set(requested)) != len(requested):
        raise ValueError("Event identifiers must be unique.")
    candidates = metadata["candidates"]
    known = {
        design.clade_ids[mapping[item["event"]["event"]["node"]]] for item in candidates
    }
    if requested is not None and set(requested) - known:
        raise ValueError("Requested event is absent from the supplied count scan.")
    branches = assign_branch_ids(tree)
    by_id = {
        identifier: branches[node]
        for node, identifier in zip(nodes, gene.node_ids, strict=True)
    }
    rows, results = [], []
    for item in candidates:
        fit = item["event"]
        selected = fit["event"]
        node = mapping[selected["node"]]
        clade = design.clade_ids[node]
        if requested is not None and clade not in requested:
            continue
        event = MultiplicationEvent(
            node, selected["retention"], selected["fraction"], selected["multiplicity"]
        )
        result = model.evaluate(fit["rates"], fit["root_mean"], event, **options)
        results.append(
            {
                "species_event_id": clade,
                "gene_log_likelihood": result.log_likelihood,
                "background_gene_log_likelihood": background_result.log_likelihood,
                "log_likelihood_error": result.log_likelihood_error,
                "origin_probability_error": result.origin_probability_error,
                "arithmetic_precision_bits": result.precision_bits,
            }
        )
        for gene_id, probability in zip(
            result.origin_node_ids, result.origin_probabilities, strict=True
        ):
            rows.append(
                {
                    "tree_id": args.tree_id,
                    "species_event_id": clade,
                    "gene_clade_id": gene_id,
                    "gene_branch_id": by_id[gene_id],
                    "conditional_wgd_probability": probability,
                    "gene_log_likelihood": result.log_likelihood,
                    "background_gene_log_likelihood": background_result.log_likelihood,
                    "conditional_log_likelihood_difference": result.log_likelihood
                    - background_result.log_likelihood,
                    "parameter_source": "supplied_count_fit",
                    "probability_meaning": "latent_node_origin_given_topology_parameters_and_one_event",
                }
            )
    columns = [
        "tree_id",
        "species_event_id",
        "gene_clade_id",
        "gene_branch_id",
        "conditional_wgd_probability",
        "gene_log_likelihood",
        "background_gene_log_likelihood",
        "conditional_log_likelihood_difference",
        "parameter_source",
        "probability_meaning",
    ]
    table = pd.DataFrame(rows, columns=columns).to_csv(
        sep="\t", index=False, float_format="%.12g", na_rep="NA"
    )
    report = {
        "schema_version": 1,
        "method": "native-fixed-colored-topology-DL-WGD-v1",
        "ascertainment": args.ascertainment,
        "root_topology_prior": "unit Yule stem rate log(root_mean); not multiplied by family-rate categories",
        "background_gene_log_likelihood": background_result.log_likelihood,
        "candidates": results,
        "numerics": "Exact positive polynomial branch flow in log space; double/extended-precision agreement checked against likelihood and assignment tolerances. Extended precision is platform-dependent; this is not a certified roundoff bound.",
        "limitations": [
            "One fixed rooted species-colored topology; observed gene branch lengths and sequence likelihood are unused.",
            "Parameters are supplied count-fit estimates, not jointly fitted to the gene tree; this is not an independent significance test or Bayes factor.",
            "Assignments condition on topology, parameters, root prior, sampling, and one supplied doubling. They are not probabilities that a WGD occurred or parameter/tree posteriors.",
            "Later family origins, HGT, ILS, allopolyploid donor trees and triplication polytomy resolution are not modeled.",
        ],
    }
    with output_transaction(paths.values()) as staged:
        for role, path in paths.items():
            Path(staged[path]).write_text(
                table
                if role == "outfile"
                else json.dumps(report, indent=2, allow_nan=False) + "\n"
            )
    if args.outfile == "-":
        sys.stdout.write(table)
