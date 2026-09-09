"""Known chronologies for the RADTE benchmark (independent of dating code)."""

from ete4 import Tree


def renamed_copy(source, index):
    tree = source.copy()
    for node in tree.leaves():
        node.name += f"_{index}"
    return tree


def gene_chronogram(species, scenario):
    if scenario == "nested":
        left, right = species.children
        if left.is_leaf:
            raise ValueError(
                "Nested-duplication simulation requires at least four species"
            )
        gene = Tree()
        gene.name = species.name
        duplicate = gene.add_child(name="D", dist=2.5)
        for copy in [1, 2]:
            subtree = renamed_copy(left, copy)
            subtree.dist = 2.5
            duplicate.add_child(subtree)
        other = renamed_copy(right, 1)
        gene.add_child(other)
        truth = 7.5
    else:
        gene = Tree()
        gene.name = "D"
        for copy in [1, 2]:
            subtree = renamed_copy(species, copy)
            subtree.dist = 10
            gene.add_child(subtree)
        truth = 20.0
    if scenario == "loss":
        # Keep both species-root lineages represented in each copy. Pruning
        # collapses lost paths while preserving elapsed time on the survivors.
        tips = list(gene.leaves())
        keep = [n.name for i, n in enumerate(tips) if i % 3 != 1]
        for subtree in gene.children:
            keep.extend(next(child.leaves()).name for child in subtree.children)
        gene.prune(set(keep), preserve_branch_length=True)
    return gene, truth


def calibration_intervals(species, fractional_width):
    if fractional_width == 0:
        return None
    rows = ["node\tage_min\tage_max\n"]
    depth = {species: 0.0}
    for node in species.traverse():
        if node is species:
            continue
        depth[node] = depth[node.up] + node.dist
        if not node.is_leaf:
            age = 10 - depth[node]
            rows.append(
                f"{node.name}\t{age * (1 - fractional_width)}\t{age * (1 + fractional_width)}\n"
            )
    # The species root stays fixed, anchoring absolute time scale.
    return "".join(rows)
