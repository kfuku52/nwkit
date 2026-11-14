"""Descriptive checks of the fitted count background, not WGD tests."""

import numpy as np


def _statistics(model):
    tree = model.tree
    observed = ~np.isnan(model.counts)
    tip_labels = [
        (i, name) for i, name in enumerate(tree.tip_names) if observed[:, i].any()
    ]
    clades: list[set[int]] = [set() for _ in tree.parents]
    for column, node in enumerate(tree.tip_nodes):
        clades[node].add(column)
    for node in range(len(tree.parents) - 1, 0, -1):
        clades[tree.parents[node]].update(clades[node])
    contrasts = []
    for node in sorted(range(1, len(clades)), key=tree.branch_ids.__getitem__):
        inside = sorted(clades[node])
        outside = sorted(clades[0] - clades[node])
        n_inside = observed[:, inside].sum(axis=1)
        n_outside = observed[:, outside].sum(axis=1)
        paired = np.flatnonzero((n_inside > 0) & (n_outside > 0))
        if len(paired):
            contrasts.append(
                (
                    tree.branch_ids[node],
                    np.ix_(paired, inside),
                    np.ix_(paired, outside),
                    n_inside[paired],
                    n_outside[paired],
                )
            )
    columns = [i for i, _ in tip_labels]
    n_observed = observed[:, columns].sum(axis=0)
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
        for _, inside, outside, ni, no in contrasts:
            yield np.mean(
                np.nansum(counts[inside], axis=1) / ni
                - np.nansum(counts[outside], axis=1) / no
            )

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
