"""Unclocked nuisance-model estimation before conditional sequence dating."""

import numpy as np
from scipy.optimize import minimize

from nwkit.fasta import parse_fasta
from nwkit.radte_sequence import DNA_CODES


def default_sequence_model(path):
    with open(path, encoding="utf-8") as handle:
        records = parse_fasta(handle)
    letters = set(
        "".join("".join(record.raw.splitlines()[1:]) for record in records).upper()
    ) - {" ", "\t"}
    return "gtr" if letters.issubset(DNA_CODES) else "lg"


def unrooted_lengths(likelihood, lengths):
    roots = {likelihood.edge_id[n] for n in likelihood.chronology.gene.children}
    mapping, next_id = [], 1
    for index in range(len(lengths)):
        if index in roots:
            mapping.append(0)
        else:
            mapping.append(next_id)
            next_id += 1
    mapping = np.asarray(mapping)
    combined = np.bincount(mapping, weights=lengths)
    fractions = lengths / combined[mapping]
    return mapping, combined, fractions


def fit_sequence_model(
    likelihood,
    lengths,
    *,
    fit_kappa=False,
    fit_gamma=False,
    fit_gtr=False,
    maxiter=1000,
):
    likelihood.fit_settings = dict(
        fit_kappa=fit_kappa, fit_gamma=fit_gamma, fit_gtr=fit_gtr
    )
    settings: list[tuple[str, int | None, float, float, float]] = []
    if fit_kappa:
        settings.append(
            ("kappa", None, np.log(likelihood.kappa), np.log(0.05), np.log(100))
        )
    if fit_gtr:
        likelihood.exchangeabilities /= likelihood.exchangeabilities[-1]
        for j in range(5):
            settings.append(
                (
                    "exchangeabilities",
                    j,
                    np.log(likelihood.exchangeabilities[j]),
                    np.log(0.01),
                    np.log(100),
                )
            )
    if fit_gamma and likelihood.gamma_categories > 1:
        settings.append(
            (
                "gamma_shape",
                None,
                np.log(likelihood.gamma_shape),
                np.log(0.02),
                np.log(50),
            )
        )
    mapping, combined, fractions = unrooted_lengths(likelihood, np.asarray(lengths))
    count = len(combined)
    initial = np.concatenate([np.log(combined), [s[2] for s in settings]])
    bounds = [(-18.0, 4.0)] * count + [(s[3], s[4]) for s in settings]

    def assign(parameters):
        if not settings:
            return
        for value, (attribute, index, _, _, _) in zip(
            parameters, settings, strict=True
        ):
            if index is None:
                setattr(likelihood, attribute, float(np.exp(value)))
            else:
                getattr(likelihood, attribute)[index] = np.exp(value)
        likelihood.update_parameters()

    def objective(x):
        assign(x[count:])
        branch_lengths = np.exp(x[:count][mapping]) * fractions
        nll, branch_gradient = likelihood.value_gradient(branch_lengths)
        gradient = np.zeros(len(x))
        gradient[:count] = np.bincount(
            mapping, weights=branch_gradient * branch_lengths
        )
        for j in range(len(settings)):
            plus, minus = x[count:].copy(), x[count:].copy()
            plus[j] += 1e-4
            minus[j] -= 1e-4
            assign(plus)
            above = likelihood.value_gradient(branch_lengths)[0]
            assign(minus)
            below = likelihood.value_gradient(branch_lengths)[0]
            gradient[count + j] = (above - below) / 2e-4
        assign(x[count:])
        return nll, gradient

    result = minimize(
        objective,
        initial,
        jac=True,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": maxiter, "ftol": 1e-11, "gtol": 1e-4},
    )
    if not result.success or not np.isfinite(result.fun):
        raise ValueError(
            "Unclocked sequence fit did not converge: " + str(result.message)
        )
    assign(result.x[count:])
    likelihood.initial_lengths = np.exp(result.x[:count][mapping]) * fractions
    active = [
        settings[j][0]
        for j, value in enumerate(result.x[count:])
        if min(abs(value - settings[j][3]), abs(value - settings[j][4])) < 1e-4
    ]
    return dict(
        status="estimated-unclocked-conditional-model"
        if settings
        else "estimated-unclocked-lengths-fixed-substitution-model",
        objective=float(result.fun),
        iterations=int(result.nit),
        boundary_parameters=sorted(set(active)),
        boundary_unrooted_branches=[
            j
            for j, value in enumerate(result.x[:count])
            if min(abs(value - bounds[j][0]), abs(value - bounds[j][1])) < 1e-4
        ],
    )
