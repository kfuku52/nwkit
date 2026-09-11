import re
import sys
from collections import Counter, defaultdict
from importlib import resources
from typing import Any

import ete4
import numpy as np
import pandas as pd

from nwkit.rooting_state import set_rooting_info
from nwkit.util import (
    extract_taxonomy_query,
    get_ete_ncbitaxa,
    label2sciname,
    read_item_per_line_file,
    read_tip_table,
    read_tree,
    remove_singleton,
    validate_unique_named_leaves,
    warn_cleanup_failure,
    write_tree,
)


def _close_ncbi_db(ncbi):
    if ncbi is None:
        return
    db = getattr(ncbi, "db", None)
    if db is not None:
        try:
            db.close()
        except Exception as exc:
            warn_cleanup_failure("NCBI taxonomy database handle", exc)


def read_taxid_tsv(taxid_tsv):
    taxid_df, _, _ = read_tip_table(
        taxid_tsv,
        option_name="--taxid-tsv",
        required_columns=("taxid",),
        unmatched="ignore",
    )
    if taxid_df.empty:
        raise ValueError("--taxid-tsv is empty.")
    if taxid_df["taxid"].isna().any():
        raise ValueError('--taxid-tsv contains missing values in the "taxid" column.')
    taxid_numeric = pd.to_numeric(taxid_df["taxid"], errors="coerce")
    if taxid_numeric.isna().any():
        raise ValueError(
            '--taxid-tsv contains non-numeric values in the "taxid" column.'
        )
    if (taxid_numeric % 1 != 0).any():
        raise ValueError(
            '--taxid-tsv contains non-integer values in the "taxid" column.'
        )
    taxid_df["taxid"] = taxid_numeric.astype(int)
    return taxid_df


def check_input_file(args):
    labels = None
    taxid_df = None
    if (args.species_list is None) and (args.taxid_tsv is None):
        raise ValueError("Either --species-list or --taxid-tsv must be specified.")
    if (args.species_list is not None) and (args.taxid_tsv is not None):
        raise ValueError("Only one of --species-list or --taxid-tsv can be specified.")
    if (args.backbone != "ncbi") and (args.taxid_tsv is not None):
        raise ValueError(
            "--taxid-tsv is currently compatible only with --backbone ncbi."
        )
    if args.species_list is not None:
        labels = read_item_per_line_file(args.species_list)
        if len(labels) == 0:
            raise ValueError("--species-list is empty.")
        if len(labels) != len(set(labels)):
            raise ValueError("Duplicate entries in --species-list are not supported.")
    if args.taxid_tsv is not None:
        taxid_df = read_taxid_tsv(args.taxid_tsv)
    return labels, taxid_df


def initialize_tree(tree):
    for leaf in tree.leaves():
        leaf.add_props(has_taxon=False, taxon_names=[])
    return tree


def name_to_taxid(sp, ncbi):
    name2taxid = ncbi.get_name_translator(
        [
            sp,
        ]
    )
    if len(name2taxid) == 0:
        genus = re.sub(" .*", "", sp)
        name2taxid = ncbi.get_name_translator(
            [
                genus,
            ]
        )
    return name2taxid


def limit_lineage_to_rank(lineage, ncbi, rank):
    if rank == "no":
        return lineage
    ranks = ncbi.get_rank(lineage)
    sorted_ranks = [ranks[taxid] for taxid in lineage]
    lineage_ge_rank = list()
    for taxid, my_rank in zip(lineage, sorted_ranks, strict=True):
        lineage_ge_rank.append(int(taxid))
        if my_rank == rank:
            break
    return lineage_ge_rank


def get_lineage(sp, ncbi, rank):
    name2taxid = name_to_taxid(sp, ncbi)
    if len(name2taxid) == 0:
        txt = "Genus-level match was not found in the NCBI database: Excluded from the output: {}\n"
        sys.stderr.write(txt.format(sp))
        return []
    key = sp if sp in name2taxid.keys() else re.sub(" .*", "", sp)
    lineage = ncbi.get_lineage(name2taxid[key][0])
    return limit_lineage_to_rank(lineage, ncbi, rank)


def get_lineages(labels, rank, args=None, ncbi=None):
    splist = [
        extract_taxonomy_query(label, args=args) or str(label) for label in labels
    ]
    owns_ncbi = ncbi is None
    if ncbi is None:
        ncbi = get_ete_ncbitaxa(args=args)
    try:
        lineages = dict()
        for sp, label in zip(splist, labels, strict=True):
            lineages[label] = get_lineage(sp, ncbi, rank)
        return lineages
    finally:
        if owns_ncbi:
            _close_ncbi_db(ncbi)


def get_lineage_from_taxid(taxid, ncbi, rank):
    lineage = ncbi.get_lineage(taxid)
    return limit_lineage_to_rank(lineage, ncbi, rank)


def get_lineages_from_taxid(taxid_df, rank, args=None, ncbi=None):
    owns_ncbi = ncbi is None
    if ncbi is None:
        ncbi = get_ete_ncbitaxa(args=args)
    try:
        taxid_numeric = pd.to_numeric(taxid_df["taxid"], errors="coerce")
        if taxid_numeric.isna().any():
            raise ValueError(
                '--taxid-tsv contains non-numeric values in the "taxid" column.'
            )
        if (taxid_numeric % 1 != 0).any():
            raise ValueError(
                '--taxid-tsv contains non-integer values in the "taxid" column.'
            )
        lineages = dict()
        for label, taxid in zip(
            taxid_df["leaf_name"], taxid_numeric.astype(int), strict=True
        ):
            lineages[label] = get_lineage_from_taxid(taxid, ncbi, rank)
        return lineages
    finally:
        if owns_ncbi:
            _close_ncbi_db(ncbi)


def match_taxa(tree, labels, backbone_method, args=None):
    splist = [
        extract_taxonomy_query(label, args=args) or str(label) for label in labels
    ]
    ncbi = get_ete_ncbitaxa(args=args) if backbone_method.startswith("ncbi") else None
    try:
        leaf_names = [ln.replace("_", " ") for ln in tree.leaf_names()]
        leaf_name_set = set(leaf_names)
        leaf_by_name = {leaf.name: leaf for leaf in tree.leaves()}
        for sp, label in zip(splist, labels, strict=True):
            ancestor: list[Any]
            if backbone_method.startswith("ncbi"):
                assert ncbi is not None
                lineage = get_lineage(sp, ncbi, rank="no")
                ancestor_names = ncbi.get_taxid_translator(lineage)
                ancestor = list(
                    leaf_name_set.intersection(set(ancestor_names.values()))
                )
                if len(ancestor) > 1:
                    txt = "Multiple hits. Excluded from the output. Taxon in the list = {}, Taxa in the tree = {}\n"
                    sys.stderr.write(txt.format(sp, ",".join(list(ancestor))))
                    continue
                elif len(ancestor) == 0:
                    txt = "No hit. Excluded from the output. Taxon = {}\n"
                    sys.stderr.write(txt.format(sp, ",".join(list(ancestor))))
                    continue
            elif backbone_method == "user":
                ancestor = [sp]
            leaf = leaf_by_name.get(ancestor[0])
            if leaf is not None:
                leaf.props["has_taxon"] = True
                leaf.props["taxon_names"].append(label)
        return tree
    finally:
        _close_ncbi_db(ncbi)


def delete_nomatch_leaves(tree):
    for leaf in list(tree.leaves()):
        if not leaf.props.get("has_taxon"):
            sys.stderr.write(
                "Deleting taxon in the tree with no match: {}\n".format(leaf.name)
            )
            leaf.delete()
    return tree


def polytomize_one2many_matches(tree):
    for leaf in tree.leaves():
        taxon_names = leaf.props.get("taxon_names", [])
        if len(taxon_names) == 0:
            raise ValueError(
                "Leaf should have at least one match: {}".format(leaf.name)
            )
        if len(taxon_names) == 1:
            leaf.name = taxon_names[0]
        else:
            for i in range(len(taxon_names)):
                leaf.add_child(name=taxon_names[i])
            leaf.name = ""
    return tree


def get_taxid_counts(lineages):
    taxid_counter: Counter[Any] = Counter()
    for lineage in lineages.values():
        taxid_counter.update(lineage)
    uniq_taxids = list(taxid_counter.keys())
    counts = [taxid_counter[taxid] for taxid in uniq_taxids]
    taxid_counts = np.array([uniq_taxids, counts]).T
    count_order = np.argsort(taxid_counts[:, 1])
    taxid_counts = taxid_counts[count_order, :]
    return taxid_counts


def get_mrca_taxid(multi_counts, args=None):
    ncbi = get_ete_ncbitaxa(args=args)
    try:
        max_count = multi_counts[:, 1].max()
        is_max_count = multi_counts[:, 1] == max_count
        max_taxids = multi_counts[is_max_count, 0]
        mrca_taxid = 1
        max_ancestor_num = 0
        for mt in max_taxids:
            ancestor_num = len(ncbi.get_lineage(mt))
            if ancestor_num > max_ancestor_num:
                mrca_taxid = mt
                max_ancestor_num = ancestor_num
        return mrca_taxid
    finally:
        _close_ncbi_db(ncbi)


def get_max_ancestor_overlap_node(tree, ancestors):
    ancestors = set(ancestors)
    max_num_overlap = 0
    return_node = None
    for node in tree.traverse():
        overlap_taxids = set(node.props.get("ancestors", [])).intersection(ancestors)
        num_overlap = len(overlap_taxids)
        if num_overlap > max_num_overlap:
            max_num_overlap = num_overlap
            return_node = node
    return return_node


def populate_leaves(new_clade, taxid, lineages):
    for sp, lineage in lineages.items():
        if taxid not in lineage:
            continue
        new_leaf = ete4.Tree({"name": sp})
        new_leaf.add_props(ancestors=lineage)
        new_clade.add_child(new_leaf)
    return new_clade


def _add_new_clade_with_leaf_set_cache(
    clades, clade_leaf_sets, new_clade, new_leaf_set
):
    flag_append = True
    delete_index = list()
    reference_new_leaf_set = set(new_leaf_set)
    leaf_name_to_nodes = None
    if len(reference_new_leaf_set) > 0:
        for j, (clade, clade_leaf_set) in enumerate(
            zip(clades, clade_leaf_sets, strict=True)
        ):
            if reference_new_leaf_set <= clade_leaf_set:
                flag_append = False
                continue
            if clade_leaf_set <= reference_new_leaf_set:
                if leaf_name_to_nodes is None:
                    leaf_name_to_nodes = defaultdict(list)
                    for leaf in new_clade.leaves():
                        leaf_name_to_nodes[leaf.name].append(leaf)
                for leaf_name in clade_leaf_set:
                    nodes = leaf_name_to_nodes.get(leaf_name, [])
                    if len(nodes) == 0:
                        continue
                    for nnode in nodes:
                        nnode.delete()
                    nodes.clear()
                new_clade.add_child(clade)
                for leaf in clade.leaves():
                    leaf_name_to_nodes[leaf.name].append(leaf)
                delete_index.append(j)
    for di in sorted(delete_index, reverse=True):
        del clades[di]
        del clade_leaf_sets[di]
    if flag_append:
        new_clade_leaves = list(new_clade.leaf_names())
        if len(new_clade_leaves) != len(set(new_clade_leaves)):
            txt = "Redundant leaves in the new clade: {}".format(
                ", ".join(new_clade_leaves)
            )
            raise Exception(txt)
        clades.append(new_clade)
        clade_leaf_sets.append(reference_new_leaf_set)
    return clades, clade_leaf_sets


def add_new_clade(clades, new_clade):
    new_leaf_set = set(new_clade.leaf_names())
    clade_leaf_sets = [set(clade.leaf_names()) for clade in clades]
    clades, _ = _add_new_clade_with_leaf_set_cache(
        clades=clades,
        clade_leaf_sets=clade_leaf_sets,
        new_clade=new_clade,
        new_leaf_set=new_leaf_set,
    )
    return clades


def get_polytomy_index(tree):
    num_polytomy = 0
    for node in tree.traverse():
        if node.is_leaf:
            continue
        num_polytomy += len(node.get_children()) - 2
    return num_polytomy


def taxid2tree(lineages, taxid_counts, args=None, ncbi=None):
    if len(lineages) == 0:
        raise ValueError("No valid taxa found to build a constraint tree.")
    if len(lineages) == 1:
        sp = next(iter(lineages.keys()))
        tree = ete4.Tree({"name": sp})
        tree.add_props(ancestors=lineages[sp])
        return tree
    is_multiple = taxid_counts[:, 1] > 1
    multi_counts = taxid_counts[is_multiple, :]
    if multi_counts.shape[0] == 0:
        raise ValueError("No shared taxonomic ranks were found across input taxa.")
    shared_taxids = {int(taxid) for taxid, count in multi_counts if int(count) > 1}
    node_by_taxid: dict[int, Any] = {}
    parent_taxid_by_taxid: dict[int, int | None] = {}
    root_taxids: set[int] = set()
    for species_name, raw_lineage in lineages.items():
        lineage = [int(taxid) for taxid in raw_lineage]
        shared_lineage = [taxid for taxid in lineage if taxid in shared_taxids]
        if not shared_lineage:
            continue
        parent_taxid = None
        for lineage_index, taxid in enumerate(lineage):
            if taxid not in shared_taxids:
                continue
            node = node_by_taxid.get(taxid)
            if node is None:
                node = ete4.Tree()
                node.add_props(ancestors=lineage[: lineage_index + 1])
                node_by_taxid[taxid] = node
                parent_taxid_by_taxid[taxid] = parent_taxid
                if parent_taxid is None:
                    root_taxids.add(taxid)
                else:
                    node_by_taxid[parent_taxid].add_child(node)
            elif parent_taxid_by_taxid[taxid] != parent_taxid:
                raise ValueError(
                    "Taxonomy lineages do not form a consistent hierarchy."
                )
            parent_taxid = taxid
        leaf = ete4.Tree({"name": species_name})
        leaf.add_props(ancestors=lineage)
        assert parent_taxid is not None
        node_by_taxid[parent_taxid].add_child(leaf)
    if len(root_taxids) != 1:
        raise ValueError("Failed to merge clades into a single tree.")
    tree = node_by_taxid[next(iter(root_taxids))]
    # Successive ranks can describe the same species set. The former clade
    # merger retained the first such rank instead of emitting a unary chain.
    for node in list(tree.traverse(strategy="postorder")):
        children = node.get_children()
        if len(children) != 1 or children[0].is_leaf:
            continue
        redundant_child = children[0]
        node.remove_child(redundant_child)
        for grandchild in list(redundant_child.get_children()):
            node.add_child(grandchild)
    return tree


def tree_sciname2label(tree, labels, args=None):
    leaf_name_to_nodes = defaultdict(list)
    for leaf in tree.leaves():
        leaf_name_to_nodes[leaf.name].append(leaf)
    for label in labels:
        label_sciname = label2sciname(labels=label, args=args)
        candidate_nodes = leaf_name_to_nodes.get(label_sciname, [])
        if len(candidate_nodes) > 0:
            candidate_nodes.pop(0).name = label
        else:
            sys.stderr.write("Original label not found: {}\n".format(label))
    return tree


def collapse_genes(tree):
    for leaf in tree.leaves():
        leaf_name = leaf.name or ""
        splitted = leaf_name.split("_")
        if len(splitted) < 2:
            raise ValueError(
                "Leaf name does not match the expected GENUS_SPECIES[_...] format for --collapse: {}".format(
                    leaf_name
                )
            )
        genus_species = splitted[0] + "_" + splitted[1]
        leaf.name = genus_species
    leaf_name_counts = Counter(tree.leaf_names())
    for leaf in list(tree.leaves()):
        if leaf_name_counts[leaf.name] > 1:
            leaf.delete(prevent_nondicotomic=False, preserve_branch_length=True)
            leaf_name_counts[leaf.name] -= 1
    return tree


def _prepare_constraint_output(tree, backbone):
    if not backbone.endswith("user"):
        set_rooting_info(tree, True, source="taxonomy")
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=False)
    for node in tree.traverse():
        node.name = (node.name or "").replace(" ", "_")
    return tree


def constrain_main(args):
    labels, taxid_df = check_input_file(args)
    if args.backbone == "ncbi":
        ncbi = get_ete_ncbitaxa(args=args)
        try:
            if args.taxid_tsv is not None:
                lineages = get_lineages_from_taxid(
                    taxid_df,
                    rank=args.rank,
                    args=args,
                    ncbi=ncbi,
                )
            else:
                lineages = get_lineages(
                    labels=labels,
                    rank=args.rank,
                    args=args,
                    ncbi=ncbi,
                )
            taxid_counts = get_taxid_counts(lineages)
            tree = taxid2tree(
                lineages,
                taxid_counts,
                args=args,
                ncbi=ncbi,
            )
        finally:
            _close_ncbi_db(ncbi)
    else:
        if args.backbone.endswith("user"):
            tree = read_tree(
                args.infile,
                args.format,
                args.quoted_node_names,
                rooted=getattr(args, "input_rooted", "auto"),
            )
            for node in tree.traverse():
                node.name = (node.name or "").replace("_", " ")
            validate_unique_named_leaves(
                tree, option_name="--infile", context=" for '--backbone user'"
            )
        elif args.backbone == "ncbi_apgiv":
            file_path = "data_tree/apgiv.nwk"
            package_name = __package__ or "nwkit"
            nwk_string = (
                resources.files(package_name)
                .joinpath(file_path)
                .read_text(encoding="utf-8")
            )
            tree = ete4.Tree(nwk_string, parser=0)
        tree = initialize_tree(tree)
        tree = match_taxa(tree, labels, args.backbone, args=args)
        if not any(leaf.props.get("has_taxon") for leaf in tree.leaves()):
            raise ValueError(
                "No taxa from --species-list matched the backbone tree labels."
            )
        tree = delete_nomatch_leaves(tree)
        tree = polytomize_one2many_matches(tree)
    if args.collapse:
        tree = collapse_genes(tree)
    tree = _prepare_constraint_output(tree, args.backbone)
    write_tree(tree, args, format=9)
