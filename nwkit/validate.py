import math
from typing import Any

import pandas as pd

from nwkit.clade_mapping import projected_root_split
from nwkit.util import (
    compute_node_ages,
    get_species_group_records,
    inspect_tree_text,
    is_rooted,
    normalize_phylogenetic_tree_text,
    read_input_text,
    read_tree,
    split_newick_stream,
    support_is_missing,
)

VALIDATE_COLUMNS = (
    "tree_id",
    "status",
    "parse_ok",
    "parse_error",
    "input_format",
    "format_ambiguous",
    "has_quoted_node_names",
    "has_quoted_internal_node_names",
    "num_leaves",
    "num_nodes",
    "num_singleton_nodes",
    "num_multifurcation_nodes",
    "num_zero_branch_nodes",
    "num_negative_branch_nodes",
    "num_non_finite_branch_nodes",
    "num_missing_support_internal_nodes",
    "num_non_finite_support_internal_nodes",
    "is_rooted",
    "is_ultrametric",
    "leaf_names_unique",
    "num_duplicate_leaf_names",
    "num_empty_leaf_names",
    "leaf_set_matches_first",
    "rooting_matches_first",
    "species_parseable",
    "species_groups_monophyletic",
    "issues",
)


def _get_duplicate_leaf_names(tree):
    leaf_names = list(tree.leaf_names())
    seen = set()
    duplicates = list()
    for leaf_name in leaf_names:
        if leaf_name in seen:
            duplicates.append(str(leaf_name))
        else:
            seen.add(leaf_name)
    return sorted(set(duplicates))


def _get_empty_leaf_names(tree):
    return [
        str(leaf.name)
        for leaf in tree.leaves()
        if (leaf.name is None) or (str(leaf.name).strip() == "")
    ]


def _collect_tree_metrics(tree):
    num_nodes = 0
    num_singleton_nodes = 0
    num_multifurcation_nodes = 0
    num_zero_branch_nodes = 0
    num_negative_branch_nodes = 0
    num_non_finite_branch_nodes = 0
    num_missing_support_internal_nodes = 0
    num_non_finite_support_internal_nodes = 0
    for node in tree.traverse():
        num_nodes += 1
        if not node.is_leaf:
            num_children = len(node.get_children())
            if num_children == 1:
                num_singleton_nodes += 1
            elif num_children > 2:
                num_multifurcation_nodes += 1
            if (not node.is_root) and support_is_missing(node.support):
                num_missing_support_internal_nodes += 1
            if node.support is not None:
                try:
                    if not math.isfinite(float(node.support)):
                        num_non_finite_support_internal_nodes += 1
                except (TypeError, ValueError):
                    num_non_finite_support_internal_nodes += 1
        if node.dist is not None:
            distance = float(node.dist)
            if not math.isfinite(distance):
                num_non_finite_branch_nodes += 1
            if math.isfinite(distance) and distance < 0.0:
                num_negative_branch_nodes += 1
            if not node.is_root:
                if distance == 0.0:
                    num_zero_branch_nodes += 1
    try:
        compute_node_ages(tree)
        is_ultrametric = True
    except ValueError:
        is_ultrametric = False
    return {
        "num_leaves": len(list(tree.leaves())),
        "num_nodes": num_nodes,
        "num_singleton_nodes": num_singleton_nodes,
        "num_multifurcation_nodes": num_multifurcation_nodes,
        "num_zero_branch_nodes": num_zero_branch_nodes,
        "num_negative_branch_nodes": num_negative_branch_nodes,
        "num_non_finite_branch_nodes": num_non_finite_branch_nodes,
        "num_missing_support_internal_nodes": num_missing_support_internal_nodes,
        "num_non_finite_support_internal_nodes": num_non_finite_support_internal_nodes,
        "is_rooted": is_rooted(tree),
        "is_ultrametric": is_ultrametric,
    }


def _validate_species_labels(tree, args):
    try:
        get_species_group_records(
            tree,
            option_name="--infile",
            context=" for 'validate'",
            args=args,
        )
    except ValueError as exc:
        message = str(exc)
        if "could not be parsed" in message:
            return False, False
        if "not monophyletic" in message:
            return True, False
        raise
    return True, True


def _build_parse_error_row(tree_id, inspection, issues):
    return {
        "tree_id": tree_id,
        "status": "invalid",
        "parse_ok": False,
        "parse_error": inspection["parse_error"],
        "input_format": inspection["input_format"],
        "format_ambiguous": inspection["format_ambiguous"],
        "has_quoted_node_names": inspection["has_quoted_node_names"],
        "has_quoted_internal_node_names": inspection["has_quoted_internal_node_names"],
        "num_leaves": "",
        "num_nodes": "",
        "num_singleton_nodes": "",
        "num_multifurcation_nodes": "",
        "num_zero_branch_nodes": "",
        "num_negative_branch_nodes": "",
        "num_non_finite_branch_nodes": "",
        "num_missing_support_internal_nodes": "",
        "num_non_finite_support_internal_nodes": "",
        "is_rooted": "",
        "is_ultrametric": "",
        "leaf_names_unique": "",
        "num_duplicate_leaf_names": "",
        "num_empty_leaf_names": "",
        "leaf_set_matches_first": "",
        "rooting_matches_first": "",
        "species_parseable": "",
        "species_groups_monophyletic": "",
        "issues": ",".join(issues),
    }


def _validation_options(args):
    return {
        "require_rooted": getattr(args, "require_rooted", False),
        "require_ultrametric": getattr(args, "require_ultrametric", False),
        "require_same_leaf_set": getattr(args, "require_same_leaf_set", True),
        "require_same_rooting": getattr(args, "require_same_rooting", False),
        "require_binary": getattr(args, "require_binary", False),
        "require_all_support": getattr(args, "require_all_support", False),
        "require_unambiguous_format": getattr(
            args, "require_unambiguous_format", False
        ),
        "require_unquoted_names": getattr(args, "require_unquoted_names", False),
        "check_species": getattr(args, "check_species", False),
        "fail_on_issue": getattr(args, "fail_on_issue", False),
    }


def _inspection_issues(inspection, options):
    issues = []
    if inspection["format_ambiguous"] and options["require_unambiguous_format"]:
        issues.append("format_ambiguous")
    if inspection["has_quoted_node_names"] and options["require_unquoted_names"]:
        issues.append("quoted_node_names")
    return issues


def _metric_issues(metrics, duplicate_leaf_names, empty_leaf_names, options):
    checks = (
        (duplicate_leaf_names, "duplicate_leaf_names"),
        (empty_leaf_names, "empty_leaf_names"),
        (metrics["num_negative_branch_nodes"] > 0, "negative_branch_length"),
        (metrics["num_non_finite_branch_nodes"] > 0, "non_finite_branch_length"),
        (metrics["num_non_finite_support_internal_nodes"] > 0, "non_finite_support"),
        (options["require_rooted"] and not metrics["is_rooted"], "not_rooted"),
        (
            options["require_ultrametric"] and not metrics["is_ultrametric"],
            "not_ultrametric",
        ),
        (
            options["require_binary"]
            and (
                metrics["num_singleton_nodes"] > 0
                or metrics["num_multifurcation_nodes"] > 0
            ),
            "not_binary",
        ),
        (
            options["require_all_support"]
            and metrics["num_missing_support_internal_nodes"] > 0,
            "missing_support",
        ),
    )
    return [issue for condition, issue in checks if condition]


def _comparison_values(tree_id, tree, metrics, references, options, issues):
    leaf_name_set = set(tree.leaf_names())
    if tree_id == 1:
        references.update(
            {
                "leaf_set": leaf_name_set,
                "rooted_state": metrics["is_rooted"],
                "root_split": (
                    projected_root_split(tree, leaf_name_set)
                    if metrics["is_rooted"]
                    else None
                ),
                "first_tree_parsed": True,
            }
        )
        return True, True
    if not references["first_tree_parsed"]:
        return "", ""
    leaf_set_matches = leaf_name_set == references["leaf_set"]
    if options["require_same_leaf_set"] and not leaf_set_matches:
        issues.append("leaf_set_mismatch")
    rooting_matches: bool | str
    if not leaf_set_matches:
        rooting_matches = ""
    elif metrics["is_rooted"] != references["rooted_state"]:
        rooting_matches = False
    elif not references["rooted_state"]:
        rooting_matches = True
    else:
        rooting_matches = (
            projected_root_split(tree, leaf_name_set) == references["root_split"]
        )
    if options["require_same_rooting"] and rooting_matches is False:
        issues.append("rooting_mismatch")
    return leaf_set_matches, rooting_matches


def _species_values(tree, args, enabled, issues):
    if not enabled:
        return "", ""
    parseable, monophyletic = _validate_species_labels(tree, args)
    if not parseable:
        issues.append("species_parse_failed")
    elif not monophyletic:
        issues.append("species_not_monophyletic")
    return parseable, monophyletic


def _build_valid_tree_row(tree_id, tree, inspection, issues, references, options, args):
    metrics = _collect_tree_metrics(tree)
    duplicate_leaf_names = _get_duplicate_leaf_names(tree)
    empty_leaf_names = _get_empty_leaf_names(tree)
    issues.extend(
        _metric_issues(metrics, duplicate_leaf_names, empty_leaf_names, options)
    )
    leaf_set_matches_first, rooting_matches_first = _comparison_values(
        tree_id, tree, metrics, references, options, issues
    )
    species_parseable, species_groups_monophyletic = _species_values(
        tree, args, options["check_species"], issues
    )
    status = "ok" if not issues else "invalid"
    return {
        "tree_id": tree_id,
        "status": status,
        "parse_ok": True,
        "parse_error": "",
        "input_format": inspection["input_format"],
        "format_ambiguous": inspection["format_ambiguous"],
        "has_quoted_node_names": inspection["has_quoted_node_names"],
        "has_quoted_internal_node_names": inspection["has_quoted_internal_node_names"],
        **metrics,
        "leaf_names_unique": not duplicate_leaf_names,
        "num_duplicate_leaf_names": len(duplicate_leaf_names),
        "num_empty_leaf_names": len(empty_leaf_names),
        "leaf_set_matches_first": leaf_set_matches_first,
        "rooting_matches_first": rooting_matches_first,
        "species_parseable": species_parseable,
        "species_groups_monophyletic": species_groups_monophyletic,
        "issues": ",".join(issues),
    }


def validate_main(args):
    options = _validation_options(args)
    raw_text = normalize_phylogenetic_tree_text(
        read_input_text(args.infile),
        collection=True,
    )
    tree_strings = split_newick_stream(raw_text)
    if len(tree_strings) == 0:
        raise Exception("Failed to parse the input trees.")
    rows = list()
    invalid_tree_ids = list()
    references: dict[str, Any] = {
        "leaf_set": None,
        "rooted_state": None,
        "root_split": None,
        "first_tree_parsed": None,
    }
    for tree_id, tree_string in enumerate(tree_strings, start=1):
        inspection = inspect_tree_text(
            newick_text=tree_string,
            format=args.format,
            quoted_node_names=args.quoted_node_names,
        )
        issues = _inspection_issues(inspection, options)
        if not inspection["parse_ok"]:
            issues.append("parse_error")
            if tree_id == 1:
                references["first_tree_parsed"] = False
            row = _build_parse_error_row(tree_id, inspection, issues)
            rows.append(row)
            invalid_tree_ids.append(tree_id)
            continue
        tree = read_tree(
            tree_string,
            args.format,
            args.quoted_node_names,
            quiet=True,
            allow_non_finite=True,
        )
        row = _build_valid_tree_row(
            tree_id, tree, inspection, issues, references, options, args
        )
        if row["status"] != "ok":
            invalid_tree_ids.append(tree_id)
        rows.append(row)
    out = pd.DataFrame(rows, columns=VALIDATE_COLUMNS)
    if args.outfile == "-":
        print(out.to_csv(sep="\t", index=False), end="")
    else:
        out.to_csv(args.outfile, sep="\t", index=False)
    if options["fail_on_issue"] and invalid_tree_ids:
        raise ValueError(
            "Validation failed for tree_id(s): {}".format(
                ", ".join(str(tree_id) for tree_id in invalid_tree_ids)
            )
        )
