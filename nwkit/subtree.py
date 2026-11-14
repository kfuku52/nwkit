import sys

from nwkit.util import (
    annotate_duplication_confidence_and_count_species,
    annotate_scientific_names,
    read_tree,
    validate_unique_named_leaves,
    write_tree,
)


def subtree_main(args):
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    validate_unique_named_leaves(
        tree,
        option_name="--infile",
        context=" for 'subtree'",
    )
    leaf_nodes = list(tree.leaves())
    leaf_positions = {leaf.name: index for index, leaf in enumerate(leaf_nodes)}
    if args.leaves is not None:
        seed_leaf_names = [
            name.strip() for name in args.leaves.split(",") if name.strip() != ""
        ]
    else:
        seed_leaf_names = [args.left_leaf, args.right_leaf]
        seed_leaf_names = [
            name.strip()
            for name in seed_leaf_names
            if name is not None and name.strip() != ""
        ]
    if len(seed_leaf_names) == 0:
        raise ValueError("Specify either '--leaves' or '--left-leaf/--right-leaf'.")
    first = len(leaf_nodes)
    last = -1
    for leaf_name in seed_leaf_names:
        position = leaf_positions.get(leaf_name)
        if position is None:
            raise Exception(
                "Specified leaf not found in the tree: {}".format(leaf_name)
            )
        first = min(first, position)
        last = max(last, position)
    # All descendant leaves of a clade are contiguous in DFS order.
    mrca = (
        leaf_nodes[first]
        if first == last
        else tree.common_ancestor([leaf_nodes[first], leaf_nodes[last]])
    )
    if args.orthogroup:
        tree = annotate_scientific_names(tree, args=args)
        species_counts = annotate_duplication_confidence_and_count_species(tree)
        num_all_species = species_counts[tree]
        sys.stderr.write(
            "Number of species in input tree: {:,}\n".format(num_all_species)
        )
        num_species_last_updated = species_counts[mrca]
        node_species_last_updated = mrca
        for current_node in mrca.ancestors():
            if current_node.is_root:
                break
            if num_species_last_updated == num_all_species:
                break
            current_num_species = species_counts[current_node]
            is_species_num_increased = current_num_species > num_species_last_updated
            is_sufficiently_small_dup_conf_score = (
                current_node.props.get("dup_conf_score")
                <= args.dup_conf_score_threshold
            )
            if is_species_num_increased and is_sufficiently_small_dup_conf_score:
                node_species_last_updated = current_node
                num_species_last_updated = current_num_species
            # txt = 'i: {}, dup_conf_score: {}, num_species: {}'
            # print(txt.format(i, current_node.props.get('dup_conf_score'), current_num_species))
        subtree = node_species_last_updated
        num_species_subtree = species_counts[subtree]
        sys.stderr.write(
            "Number of species in output orthogroup tree: {:,}\n".format(
                num_species_subtree
            )
        )
    else:
        subtree = mrca
    write_tree(subtree, args, format=args.outformat)
