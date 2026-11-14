"""Restricted one-event, direct-second-parent MUL/MSC candidate model."""

import math
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from itertools import permutations, product
from typing import Any

from scipy.special import logsumexp

from nwkit.mul_coalescent import CoalescentTopology
from nwkit.mul_reconcile_model import (
    hypotheses,
    label_nodes,
    select_nodes,
    validate_binary,
)
from nwkit.util import compute_node_ages


@dataclass(frozen=True)
class MscHypothesis:
    id: int
    tree: Any
    h1: str
    h2: str
    status: str
    reason: str = ""
    age_bounds: tuple[float, float] | None = None


def _attachment_interval(ages, polyploid, parent, bounds):
    lower_edges = (ages[polyploid], ages[parent])
    upper_edges = (
        ages[polyploid.up],
        ages[parent.up] if parent.up is not None else math.inf,
    )
    lower, upper = max(bounds[0], *lower_edges), min(bounds[1], *upper_edges)
    if lower in lower_edges:
        lower = math.nextafter(lower, upper)
    if upper in upper_edges:
        upper = math.nextafter(upper, lower)
    return (lower, upper) if lower < upper else None


def _dated_input(species, h1, hybridization_age, age_bounds):
    """Validate the chronogram and fixed H1 before generating donor candidates."""
    validate_binary(species, "Species tree")
    if not h1:
        raise ValueError("MSC requires one fixed polyploid clade via --h1.")
    if age_bounds is None and (
        hybridization_age is None
        or not math.isfinite(hybridization_age)
        or hybridization_age <= 0
    ):
        raise ValueError("--hybridization-age must be finite and positive.")
    if age_bounds is not None and (
        len(age_bounds) != 2
        or any(not math.isfinite(value) or value <= 0 for value in age_bounds)
        or age_bounds[0] >= age_bounds[1]
    ):
        raise ValueError("Hybridization age bounds require finite 0 < lower < upper.")
    original = species.copy()
    lengths = [node.dist for node in original.traverse() if not node.is_root]
    if any(
        length is None or not math.isfinite(length) or length < 0 for length in lengths
    ):
        raise ValueError("MSC requires all finite, nonnegative species branch lengths.")
    scale = max(lengths)
    ages = compute_node_ages(original, tolerance=scale * 1e-10)
    if any(not math.isfinite(age) for age in ages.values()):
        raise ValueError("MSC species node ages must be finite.")
    label_nodes(original)
    selected = select_nodes(original, h1)
    if len(selected) != 1 or selected[0].is_root:
        raise ValueError("MSC requires a single non-root --h1 clade.")
    polyploid = selected[0]
    if age_bounds is None and not (
        ages[polyploid] < hybridization_age < ages[polyploid.up]
    ):
        raise ValueError("--hybridization-age must lie strictly inside the H1 stem.")
    for node in original.traverse():
        node.props.pop("msc_attachment", None)
        node.add_prop("msc_age", ages[node])
    return original, ages, polyploid


def dated_candidates(
    species,
    h1,
    h2,
    hybridization_age,
    *,
    max_candidates=10000,
    age_bounds=None,
    require_evaluated=True,
):
    """Unfold a direct H2 donor edge; H1's stem represents the first parent."""
    original, ages, polyploid = _dated_input(species, h1, hybridization_age, age_bounds)
    candidates = hypotheses(original, h1, h2, max_candidates=max_candidates)
    source = {node.name: node for node in original.traverse()}
    result = []
    for candidate in candidates:
        reason = ""
        interval = None
        attachment_age = hybridization_age
        if candidate.kind == "no-polyploidy":
            reason = "No-polyploidy DL+ILS comparison is outside conditional MSC; use the separate locus-mc model."
        elif candidate.kind == "autopolyploid":
            reason = "Autopolyploid inheritance is outside the disomic MSC prototype."
        else:
            parent = source[candidate.h2]
            upper = ages[parent.up] if parent.up is not None else math.inf
            if age_bounds is not None:
                interval = _attachment_interval(ages, polyploid, parent, age_bounds)
                if interval is not None:
                    lower, upper = interval
                    attachment_age = lower + (upper - lower) / 2
                else:
                    reason = "No H1/H2 branch overlap within the supplied age bounds."
            elif not ages[parent] < hybridization_age < upper:
                reason = (
                    "The H2 branch does not exist at the supplied hybridization age."
                )
        if not reason:
            for node in candidate.tree.traverse():
                if "msc_age" not in node.props:
                    node.add_prop("msc_age", attachment_age)
                    node.add_prop("msc_attachment", True)
            for node in candidate.tree.traverse():
                node.dist = (
                    0.0
                    if node.is_root
                    else node.up.props["msc_age"] - node.props["msc_age"]
                )
        result.append(
            MscHypothesis(
                candidate.id,
                candidate.tree,
                candidate.h1,
                candidate.h2,
                "excluded" if reason else "evaluated",
                reason,
                interval,
            )
        )
    if require_evaluated and not any(
        candidate.status == "evaluated" for candidate in result
    ):
        raise ValueError("No temporally compatible allopolyploid H2 candidates remain.")
    return tuple(result), tuple(sorted(polyploid.leaf_names()))


def retime_candidate(candidate, age):
    """Change the direct donor attachment without moving supplied species ages."""
    if candidate.age_bounds is not None and not (
        candidate.age_bounds[0] <= age <= candidate.age_bounds[1]
    ):
        raise ValueError("MSC attachment age is outside candidate bounds.")
    tree = candidate.tree.copy()
    for node in tree.traverse():
        if node.props.get("msc_attachment"):
            node.props["msc_age"] = age
    for node in tree.traverse():
        node.dist = (
            0.0 if node.is_root else node.up.props["msc_age"] - node.props["msc_age"]
        )
        if not math.isfinite(node.dist) or node.dist < 0:
            raise ValueError(
                "MSC attachment creates invalid population branch lengths."
            )
    return replace(candidate, tree=tree)


def copy_assignments(gene, species, parser, *, max_assignments=10000):
    """Uniform prior on injective, species-preserving homoeolog assignments."""
    by_species = defaultdict(list)
    occurrences = defaultdict(list)
    for node in gene.leaves():
        by_species[parser.parse(node.name).species_label].append(node.name)
    for node in species.leaves():
        occurrences[node.props.get("mul_species", node.name)].append(node)
    groups = []
    count = 1
    for name, genes in by_species.items():
        tips = occurrences.get(name, ())
        if not tips or len(genes) > len(tips):
            raise ValueError(
                f"MSC input has unmatched or excess copies for species {name!r}; "
                "at most one sample per diploid/subgenome, no SSD or allelic replicates."
            )
        ways = math.perm(len(tips), len(genes))
        count *= ways
        if count > max_assignments:
            raise ValueError("MSC assignments exceed --max-coalescent-assignments.")
        groups.append(
            tuple(
                tuple(zip(genes, assignment, strict=True))
                for assignment in permutations(tips, len(genes))
            )
        )
    return count, (
        dict(pair for group in choice for pair in group) for choice in product(*groups)
    )


def validate_sampling(genes, parser, polyploid_species, species_names):
    """Reject ambiguous sampling rather than silently pruning or calling copies alleles."""
    polyploid = set(polyploid_species)
    known = set(species_names)
    for gene in genes:
        validate_binary(gene, "Gene tree")
        counts = Counter(
            parser.parse(node.name).species_label for node in gene.leaves()
        )
        for name, count in counts.items():
            if name not in known or count > (2 if name in polyploid else 1):
                raise ValueError(
                    f"MSC unmatched/excess copies for species {name!r}: {count}; "
                    "requires <=1 diploid copy or <=2 distinct homoeologs."
                )


def score_gene(
    gene, candidate, parser, *, time_scale=1.0, max_states=100000, max_assignments=10000
):
    count, assignments = copy_assignments(
        gene, candidate.tree, parser, max_assignments=max_assignments
    )
    solver = CoalescentTopology(gene, max_states=max_states)
    likelihood = float(
        logsumexp(
            [
                solver.log_probability(
                    candidate.tree, assignment, time_scale=time_scale
                )
                for assignment in assignments
            ]
        )
    ) - math.log(count)
    return likelihood, count, solver.work
