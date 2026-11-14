"""Independent CTMC/dense-covariance generator for RADTE stress validation.

Does not call native rate precision, transition or sequence simulation code.
Shared deterministic scenario definitions describe the chronology only.
"""

import json

import numpy as np
from ete4 import Tree
from radte_benchmark_cases import calibration_intervals, gene_chronogram
from scipy.linalg import expm


def species_text(names, shape, age=10.0):
    if len(names) == 1:
        return names[0], 0.0
    split = len(names) // 2 if shape == "balanced" else len(names) - 1
    children = [
        species_text(part, shape, age / 2) for part in (names[:split], names[split:])
    ]
    return "(" + ",".join(
        f"{text}:{age - child_age}" for text, child_age in children
    ) + ")", age


def simulate(
    directory,
    count,
    sites,
    seed,
    sd,
    *,
    scenario="root",
    shape="balanced",
    rho=0.0,
    width=0.0,
    model="jc69",
    kappa=2.0,
    copy_ratio=1.0,
):
    directory.mkdir(parents=True, exist_ok=False)
    text, _ = species_text([f"S{i}" for i in range(count)], shape)
    species = Tree(text + ";", parser=1)
    for i, node in enumerate(species.traverse()):
        if not node.is_leaf:
            node.name = f"N{i}"
    species.write(
        outfile=str(directory / "species.nwk"), parser=1, format_root_node=True
    )
    gene, truth = gene_chronogram(species, scenario)
    edges = [node for node in gene.traverse() if node is not gene]
    ancestry = {}
    for node in edges:
        path = []
        current = node
        while current is not None:
            path.append(current)
            current = current.up
        ancestry[node] = path
    covariance = np.empty((len(edges), len(edges)))
    for i, a in enumerate(edges):
        for j, b in enumerate(edges):
            common = next(node for node in ancestry[a] if node in ancestry[b])
            distance = ancestry[a].index(common) + ancestry[b].index(common)
            covariance[i, j] = rho**distance
    rng = np.random.default_rng(seed)
    deviations = rng.multivariate_normal(np.zeros(len(edges)), sd**2 * covariance)
    realized = []
    for node, deviation in zip(edges, deviations, strict=True):
        shift = (
            copy_ratio if all(tip.name.endswith("_2") for tip in node.leaves()) else 1.0
        )
        node.dist *= 0.01 * shift * np.exp(deviation)
        realized.append(
            dict(
                edge=len(realized), name=node.name, log_rate_deviation=float(deviation)
            )
        )
    gene.write(outfile=str(directory / "gene.nwk"), parser=1, format_root_node=True)
    (directory / "mapping.tsv").write_text(
        "leaf_name\tspecies_label\n"
        + "".join(
            f"{node.name}\t{node.name.rsplit('_', 1)[0]}\n" for node in gene.leaves()
        )
    )
    bounds = calibration_intervals(species, width)
    if bounds is not None:
        (directory / "bounds.tsv").write_text(bounds)
    q = np.ones((4, 4))
    if model == "hky":
        q[0, 2] = q[2, 0] = q[1, 3] = q[3, 1] = kappa
    np.fill_diagonal(q, 0)
    np.fill_diagonal(q, -q.sum(axis=1))
    q /= -np.trace(q) / 4
    states = {gene: rng.integers(4, size=sites)}
    for node in edges:
        p = expm(q * node.dist)
        cumulative = np.cumsum(p[states[node.up]], axis=1)
        states[node] = (rng.random(sites)[:, None] > cumulative).sum(axis=1)
    alphabet = np.array(list("ACGT"))
    (directory / "alignment.fasta").write_text(
        "".join(
            f">{node.name}\n{''.join(alphabet[states[node]])}\n"
            for node in gene.leaves()
        )
    )
    (directory / "truth.json").write_text(
        json.dumps(
            dict(
                duplication_age=truth,
                target_ages={"D": truth},
                scenario=scenario,
                rate_deviations=realized,
                generator="independent-dense-covariance-expm",
                seed=seed,
            )
        )
    )
