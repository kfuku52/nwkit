"""Parametric rate and sequence draws for conditional RADTE validation.

Rates are redrawn on the rooted genealogy before simulating alignment sites.
This differs from resampling columns of one realized family's alignment.
"""

import numpy as np

from nwkit.radte_sequence import SequenceLikelihood


def simulate_rate_lengths(chronology, ages, mean, sd, rho, rng):
    """Return substitution lengths; mean is in normalized chronology units."""
    if not np.isfinite([mean, sd, rho]).all() or sd < 0 or not 0 <= rho < 1:
        raise ValueError("Finite mean, nonnegative SD and 0 <= rho < 1 required.")
    deviations = {chronology.gene: rng.normal(0, sd)}
    for node in chronology.gene.traverse():
        if node is not chronology.gene:
            deviations[node] = rho * deviations[node.up] + rng.normal(
                0, sd * np.sqrt(1 - rho**2)
            )
    return chronology.durations(ages) * np.exp(
        mean + np.array([deviations[node] for node in chronology.edges])
    )


def simulate_sequence_likelihood(exact, lengths, rng, *, chronology=None):
    """Simulate the fitted CTMC and preserve the observed missing-data mask.

    Gamma category is drawn once per site, shared by all branches at that site.
    The reconstructed likelihood re-estimates empirical frequencies, as the usual
    sequence input route does. Other nuisance refits use the saved fit_settings.
    """
    if not isinstance(exact, SequenceLikelihood):
        raise ValueError(
            "Parametric sequence simulation requires the native alignment."
        )
    full_mask = (1 << len(exact.pi)) - 1
    ambiguous = (exact.raw_matrix & (exact.raw_matrix - 1)) != 0
    if np.any(ambiguous & (exact.raw_matrix != full_mask)) or np.any(
        exact.raw_matrix == 0
    ):
        raise ValueError(
            "Parametric simulation of informative ambiguity requires an observation model; only complete states and fully missing sites are supported."
        )
    lengths = np.asarray(lengths)
    if (
        lengths.shape != (len(exact.edges),)
        or not np.isfinite(lengths).all()
        or np.any(lengths < 0)
    ):
        raise ValueError("One finite nonnegative length per rooted edge is required.")
    sites = exact.raw_matrix.shape[1]
    categories = rng.integers(len(exact.rates), size=sites)
    states = {exact.chronology.gene: rng.choice(len(exact.pi), sites, p=exact.pi)}
    for node in exact.preorder:
        if node is exact.chronology.gene:
            continue
        parent = states[node.up]
        child = np.empty(sites, dtype=int)
        for category, rate in enumerate(exact.rates):
            selected = np.flatnonzero(categories == category)
            transition, _ = exact.transition(lengths[exact.edge_id[node]], rate)
            transition /= transition.sum(axis=1, keepdims=True)
            cumulative = np.cumsum(transition[parent[selected]], axis=1)
            child[selected] = np.sum(
                rng.random(len(selected))[:, None] > cumulative, axis=1
            )
        states[node] = child
    matrix = np.array(
        [
            np.left_shift(np.uint64(1), states[node].astype(np.uint64))
            for node in exact.chronology.gene.leaves()
        ],
        dtype=exact.raw_matrix.dtype,
    )
    matrix[ambiguous] = exact.raw_matrix[ambiguous]
    result = SequenceLikelihood(
        chronology or exact.chronology,
        None,
        model=exact.model,
        kappa=exact.kappa,
        gamma_shape=exact.gamma_shape,
        gamma_categories=exact.gamma_categories,
        matrix=matrix,
        exchangeabilities=exact.exchangeabilities,
        omega=exact.omega,
        codon_frequencies=exact.codon_frequencies,
        genetic_code=exact.genetic_code,
    )
    result.fit_settings = exact.fit_settings.copy()
    return result
