"""Auditable node mappings for every co-optimal D+L hypothesis."""

import csv
import hashlib
import json
from functools import lru_cache

from nwkit.clade_index import CladeIndex
from nwkit.mul_reconcile_model import Reconciliation
from nwkit.rooting_state import require_rooted
from nwkit.util import validate_unique_named_leaves

NODE_FIELDS = (
    "gene.tree",
    "gene_topology_id",
    "species_topology_id",
    "gene_clade_id",
    "gene_node",
    "gene_label",
    "is_leaf",
    "mul.tree",
    "h1.node",
    "h2.node",
    "hypothesis.kind",
    "mapping.id",
    "optimal.mappings",
    "total.score",
    "mul_node",
    "mul_label",
    "mul_descendant_tips",
    "mapped_species",
    "duplication",
    "child_edge_losses",
    "root_losses",
)


def _rooted_topology_index(tree):
    validate_unique_named_leaves(tree, "Node diagnostic tree")
    require_rooted(tree, "Node diagnostic tree must be rooted.")
    index = CladeIndex(tree)
    digest = hashlib.sha256(b"nwkit-rooted-topology-v1\0")
    for clade in sorted(index.clade_id_for_node(n) for n in tree.traverse()):
        digest.update(clade.encode("ascii") + b"\0")
    return "rooted-topology-sha256:" + digest.hexdigest(), index


def rooted_topology_id(tree):
    """Identify labeled rooted topology, ignoring lengths and child order."""
    identity, _ = _rooted_topology_index(tree)
    return identity


def _mapped_taxa_lookup(clades):
    @lru_cache(maxsize=64)
    def lookup(mapped):
        return (
            json.dumps(clades.names_for_node(mapped), ensure_ascii=True),
            json.dumps(
                sorted({n.props.get("mul_species", n.name) for n in mapped.leaves()}),
                ensure_ascii=True,
            ),
        )

    return lookup


def write_node_diagnostics(
    handle, candidates, genes, parser, species, *, max_state_pairs, max_maps
):
    """Stream complete assignments; limits fail the enclosing transaction."""
    writer = csv.DictWriter(handle, fieldnames=NODE_FIELDS, delimiter="\t")
    writer.writeheader()
    species_id = rooted_topology_id(species)
    for candidate in candidates:
        # Occurrence clades cannot be compared using duplicate literal labels.
        validate_unique_named_leaves(candidate.tree, "MUL node diagnostic tree")
        mapped_taxa = _mapped_taxa_lookup(CladeIndex(candidate.tree))
        species_index = None
        for number, gene in enumerate(genes, 1):
            gene_id, gene_clades = _rooted_topology_index(gene)
            result = Reconciliation(
                gene,
                candidate.tree,
                parser,
                max_state_pairs,
                _species_index=species_index,
            )
            species_index = result._species_index
            for mapping_id, (_, _, mapping) in enumerate(result.mappings(max_maps), 1):
                for node, row in zip(result.nodes, mapping, strict=True):
                    descendants, species_labels = mapped_taxa(
                        result.species[row["mul_node"]]
                    )
                    writer.writerow(
                        {
                            **row,
                            "gene.tree": number,
                            "gene_topology_id": gene_id,
                            "species_topology_id": species_id,
                            "gene_clade_id": gene_clades.clade_id_for_node(node),
                            "is_leaf": int(node.is_leaf),
                            "mul.tree": candidate.id,
                            "h1.node": candidate.h1,
                            "h2.node": candidate.h2,
                            "hypothesis.kind": candidate.kind,
                            "mapping.id": mapping_id,
                            "optimal.mappings": result.num_maps,
                            "total.score": result.score,
                            "mul_descendant_tips": descendants,
                            "mapped_species": species_labels,
                        }
                    )
