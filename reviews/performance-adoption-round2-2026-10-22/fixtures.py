"""Deterministic fixtures copied from the independently reviewed candidate study."""

import numpy as np

from nwkit import asr
from nwkit import full_ou as module
from nwkit.evolution import build_evolutionary_process
from nwkit.multivariate_gaussian_asr import _prepare_observations
from nwkit.util import read_tree


def signal_fixture(n, *, shape="balanced", seed=71, truth=0.65, nonultrametric=False):
    rng = np.random.default_rng(seed)
    covariance = np.zeros((n, n))
    if shape == "balanced":
        size = 2
        while size < n:
            for start in range(0, n, size):
                covariance[start : start + size, start : start + size] += 0.25
            size *= 2
    elif shape == "pectinate":
        for start in range(1, n - 1):
            covariance[start:, start:] += 1 / n
    covariance[np.diag_indices(n)] += rng.lognormal(0, 0.7, n) if nonultrametric else 1
    diagonal = np.diag(np.diag(covariance))
    generating = diagonal + truth * (covariance - diagonal)
    values = 1.7 + np.linalg.cholesky(generating) @ rng.standard_normal(n)
    return covariance, values, np.zeros(n)


def gaussian_fixture(size, topology="balanced", model="brownian", parameter=None):
    leaves = [f"T{i}:{0.2 + (i % 11) / 20:.12g}" for i in range(size)]
    if topology == "star":
        newick = "(" + ",".join(leaves) + ")R;"
    elif topology == "pectinate":
        newick = f"({leaves[0]},{leaves[1]})I0:0.2"
        for index, leaf in enumerate(leaves[2:], 1):
            newick = f"({newick},{leaf})I{index}:0.2"
        newick = newick.rsplit(":", 1)[0] + ";"
    else:
        current = leaves
        counter = 0
        while len(current) > 1:
            next_nodes = []
            for index in range(0, len(current), 2):
                if index + 1 == len(current):
                    next_nodes.append(current[index])
                else:
                    next_nodes.append(
                        f"({current[index]},{current[index + 1]})I{counter}:0.3"
                    )
                    counter += 1
            current = next_nodes
        newick = current[0].rsplit(":", 1)[0] + ";"
    tree = read_tree(newick, "1", True, quiet=True, rooted="yes")
    process = build_evolutionary_process(tree, model=model, parameter=parameter)
    return tree, process, tuple(tree.leaves())


def count_fixture(counts, *, detection=None, scales=None, ascertainment="observed"):
    import numpy as np

    from nwkit.wgd_count_model import CountLikelihood, CountTree

    tree = CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 0.4, 0.25, 0.45, 0.7),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )
    return CountLikelihood(
        tree,
        np.asarray(counts, dtype=float),
        detection=detection,
        rate_scales=scales,
        ascertainment=ascertainment,
    )


def rare_counts(families, maximum=84, *, missing=True):
    import numpy as np

    counts = np.ones((families, 3))
    counts[:, 0] = maximum - np.arange(families) % min(16, maximum - 1)
    if missing:
        counts[1::4, 1] = np.nan
        counts[2::4, 2] = np.nan
    return counts


def pca_fixture(n, p, *, seed=7125, nonultrametric=False):
    covariance = np.eye(n) * 0.75
    width = 2
    while width < n:
        for start in range(0, n, width):
            covariance[start : start + width, start : start + width] += 0.25
        width *= 2
    if nonultrametric:
        covariance[np.diag_indices(n)] += np.linspace(0.01, 0.9, n)
    rng = np.random.default_rng(seed + n + p)
    values = np.linalg.cholesky(covariance) @ rng.normal(size=(n, p))
    values += np.linspace(-1, 2, p)
    return covariance, values


def map_fixture(tips, states=3, *, unique=False, regime=False, projection=False):
    lengths = iter(np.linspace(0.15, 0.8, 4 * tips))
    leaves = [f"T{i}:{next(lengths) if unique else 0.4:.16g}" for i in range(tips)]
    while len(leaves) > 1:
        leaves = [
            "("
            + ",".join(leaves[i : i + 2])
            + f"):{next(lengths) if unique else 0.4:.16g}"
            if i + 1 < len(leaves)
            else leaves[i]
            for i in range(0, len(leaves), 2)
        ]
    tree = read_tree(leaves[0] + ";", "1", True, quiet=True, rooted="yes")
    labels = tuple(str(i) for i in range(states))
    observed = {
        f"T{i}": None if i % 7 == 0 else labels[i % states] for i in range(tips)
    }
    evidence = {
        name: np.ones(states)
        if value is None
        else np.array([float(value == x) for x in labels])
        for name, value in observed.items()
    }
    _, fit = asr.compute_mk_marginals(tree, labels, observed, evidence, rate=0.45)
    if regime:
        matrices = [fit["rate_matrix"].copy(), fit["rate_matrix"].copy() * 1.7]
        by_node = {
            node: matrices[i % 2]
            for i, node in enumerate(tree.traverse("preorder"))
            if not node.is_root
        }
        inside, _, terms, transitions = asr._compute_log_inside_likelihoods(
            tree, evidence, by_node
        )
        outside = asr._compute_log_outside_likelihoods(
            tree, terms, transitions, fit["root_prior"]
        )
        fit.update(
            rate_matrix_by_node=by_node,
            log_inside=inside,
            transition_matrices=transitions,
            posterior_by_node={
                node: asr._probabilities_from_logs(inside[node] + outside[node])
                for node in tree.traverse()
            },
        )
    if projection:
        fit["state_projection"] = np.arange(states) % 2
        labels = ("0", "1")
    return tree, labels, fit


def er_fixture(tips, states, *, unique=False):
    rng = np.random.default_rng(9090 + tips + states)
    count = 0

    def length():
        return float(rng.uniform(0.2, 0.8)) if unique else 0.4

    parts = [f"T{i}:{length():.17g}" for i in range(tips)]
    while len(parts) > 1:
        following = []
        for start in range(0, len(parts), 2):
            if start + 1 == len(parts):
                following.append(parts[start])
            else:
                following.append(
                    f"({parts[start]},{parts[start + 1]})I{count}:{length():.17g}"
                )
                count += 1
        parts = following
    tree = read_tree(parts[0] + ";", "1", True, quiet=True, rooted="yes")
    labels = tuple(f"s{i}" for i in range(states))
    observed, likelihood = {}, {}
    for i, leaf in enumerate(tree.leaves()):
        row = np.zeros(states)
        row[i % states] = 1.0
        value = labels[i % states]
        if i % 11 == 0:
            row[:] = 1.0
            value = None
        elif i % 7 == 0:
            row[(i + 1) % states] = 1.0
            value = (labels[i % states], labels[(i + 1) % states])
        observed[leaf.name] = value
        likelihood[leaf.name] = row
    return tree, labels, observed, likelihood


def ou_parameters(dimension, seed=1010, *, initial=False):
    values, _ = module.full_ou_initial(dimension)
    values = np.asarray(values)
    if not initial:
        values = values + np.random.default_rng(seed).normal(0, 0.12, len(values))
    return values


def ou_design(dimension, tips, *, equal=False, gap=None, limit=4096):
    if gap is not None:
        lengths = [1 + gap * i for i in range(tips)]
    elif equal:
        lengths = [1.0] * tips
    else:
        lengths = [0.25 + (i % 8) * 0.125 for i in range(tips)]
    leaves = [f"T{i}:{length:.16g}" for i, length in enumerate(lengths)]
    tree = read_tree(
        "(" + ",".join(leaves) + ")R;", "1", True, quiet=True, rooted="yes"
    )
    values = {
        f"T{i}": list(np.random.default_rng(710 + i).normal(size=dimension))
        for i in range(tips)
    }
    data = _prepare_observations(
        tree, values, tuple(f"trait{i}" for i in range(dimension)), None
    )
    return module.observation_moment_design(data, limit=limit)
