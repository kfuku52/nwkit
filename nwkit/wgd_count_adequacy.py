"""Descriptive checks of the fitted count background, not WGD tests."""

import numpy as np


def _subtree_totals(values, node_count, tip_nodes, parents, *, dtype=None):
    totals = np.zeros(
        (node_count, len(values)), dtype=values.dtype if dtype is None else dtype
    )
    totals[tip_nodes] = values.T
    for node in range(node_count - 1, 0, -1):
        totals[parents[node]] += totals[node]
    return totals


def _statistics(model):
    tree = model.tree
    observed = ~np.isnan(model.counts)
    tip_labels = [
        (i, name) for i, name in enumerate(tree.tip_names) if observed[:, i].any()
    ]
    node_count = len(tree.parents)
    tip_nodes = np.array(tree.tip_nodes)
    observed_totals = _subtree_totals(
        observed, node_count, tip_nodes, tree.parents, dtype=np.intp
    )
    contrasts = []
    count_dtype = np.int32 if len(tip_nodes) <= np.iinfo(np.int32).max else np.intp
    for node in sorted(range(1, node_count), key=tree.branch_ids.__getitem__):
        n_inside = observed_totals[node]
        n_outside = observed_totals[0] - n_inside
        paired = np.flatnonzero((n_inside > 0) & (n_outside > 0))
        if len(paired):
            contrasts.append(
                (
                    tree.branch_ids[node],
                    node,
                    paired,
                    n_inside[paired].astype(count_dtype, copy=False),
                    n_outside[paired].astype(count_dtype, copy=False),
                )
            )
    columns = [i for i, _ in tip_labels]
    n_observed = observed[:, columns].sum(axis=0)
    del observed_totals
    tip_names = [name for _, name in tip_labels]
    labels = {
        "tip_mean": tip_names,
        "tip_zero_rate": tip_names,
        "tip_multicopy_rate": tip_names,
        "mean_copy_contrast": [branch for branch, *_ in contrasts],
    }

    def evaluate(counts):
        tips = counts[:, columns]
        yield from np.nansum(tips, axis=0) / n_observed
        yield from (tips == 0).sum(axis=0) / n_observed
        yield from (tips >= 2).sum(axis=0) / n_observed
        # Counts are nonnegative integers. Below 2**53 their sums and the
        # complement subtraction are exact for every subtree and family.
        # For larger counts, sum the original tip selections independently.
        integer_sums_are_exact = np.nanmax(counts) <= 2.0**53 / len(tip_nodes)
        if integer_sums_are_exact:
            totals = _subtree_totals(
                np.where(observed, counts, 0.0), node_count, tip_nodes, tree.parents
            )
        else:
            clades: list[set[int]] = [set() for _ in tree.parents]
            for column, node in enumerate(tree.tip_nodes):
                clades[node].add(column)
            for node in range(node_count - 1, 0, -1):
                clades[tree.parents[node]].update(clades[node])
        for _, node, paired, ni, no in contrasts:
            if integer_sums_are_exact:
                inside = totals[node, paired]
                outside = totals[0, paired] - inside
            else:
                inside = np.nansum(counts[np.ix_(paired, sorted(clades[node]))], axis=1)
                outside = np.nansum(
                    counts[np.ix_(paired, sorted(clades[0] - clades[node]))], axis=1
                )
            yield np.mean(inside / ni - outside / no)

    return evaluate, labels


def count_background_adequacy(model, fit, draws, seed):
    """Compare observed summaries with simulations at a fixed fitted background."""
    if isinstance(draws, bool) or not isinstance(draws, int) or draws < 4:
        raise ValueError(
            "Count model adequacy requires at least four simulation draws."
        )
    rng = np.random.default_rng(np.random.SeedSequence([seed, 0x57474441]))
    evaluate, labels = _statistics(model)
    observed = np.array(list(evaluate(model.counts)))
    simulated = np.array(
        [
            list(evaluate(model.simulate(fit.rates, fit.root_mean, rng)))
            for _ in range(draws)
        ]
    )
    lower, median, upper = np.quantile(simulated, [0.025, 0.5, 0.975], axis=0)
    training, evaluation = np.split(simulated, [draws // 2])
    center, spread = training.mean(axis=0), training.std(axis=0, ddof=1)
    # Exact constants occasionally have a nonzero rounded standard deviation.
    varying = np.any(training != training[0], axis=0) & (spread > 0)
    statistics = {}
    offset = 0
    for name, targets in labels.items():
        indices = np.arange(offset, offset + len(targets))
        offset += len(targets)
        target_key = "branch_id" if name == "mean_copy_contrast" else "tip_name"
        rows = [
            {
                target_key: target,
                "observed": float(observed[i]),
                "simulated_median": float(median[i]),
                "simulated_lower": float(lower[i]),
                "simulated_upper": float(upper[i]),
                "used_in_joint_discrepancy": bool(varying[i]),
            }
            for i, target in zip(indices, targets, strict=True)
        ]
        indices = indices[varying[indices]]
        joint = None
        if len(indices):
            observed_z = (observed[indices] - center[indices]) / spread[indices]
            simulated_z = (evaluation[:, indices] - center[indices]) / spread[indices]
            observed_rms = float(np.sqrt(np.mean(observed_z**2)))
            simulated_rms = np.sqrt(np.mean(simulated_z**2, axis=1))
            hits = int(np.count_nonzero(simulated_rms >= observed_rms))
            joint = {
                "observed_standardized_rms": observed_rms,
                "simulated_95th_percentile": float(np.quantile(simulated_rms, 0.95)),
                "num_evaluation_exceedances": hits,
                "simulated_upper_tail_fraction": (1 + hits) / (1 + len(evaluation)),
            }
        statistics[name] = {"summaries": rows, "joint_discrepancy": joint}
    return {
        "method": "fitted-background-simulation-discrepancy",
        "seed": seed,
        "num_draws": draws,
        "num_training_draws": len(training),
        "num_evaluation_draws": len(evaluation),
        "statistics": statistics,
        "interpretation": "Descriptive comparison at the fitted background, with original missingness and ascertainment. No refitting, formal goodness-of-fit test, or change to genome-event P-values.",
    }
