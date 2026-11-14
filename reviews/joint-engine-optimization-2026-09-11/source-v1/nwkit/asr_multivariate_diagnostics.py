"""Joint samples, predictive discrepancies and bootstrap for vector ASR."""

import numpy as np
import pandas as pd

from nwkit.util import assign_branch_ids, get_node_class
from nwkit.vector_gaussian import condition_vector_tree
from nwkit.vector_processes import vector_brownian_process, vector_ou_process


def fitted_vector_process(tree, model, fit):
    dimension = len(fit.trait_names)
    if model == "MV-BM":
        return vector_brownian_process(tree, fit.sigma)
    if model == "MV-OU":
        attraction = np.eye(dimension) * fit.alpha
        diffusion = 2 * fit.alpha * fit.sigma
    elif model == "MV-OU-DIAG":
        attraction = np.diag(fit.alpha_by_trait)
        diffusion = fit.diffusion_sigma
    elif model == "MV-OU-FULL":
        attraction = fit.attraction_matrix
        diffusion = fit.diffusion_sigma
    else:
        raise ValueError(f"Unsupported multivariate diagnostic model: {model}.")
    return vector_ou_process(tree, attraction, diffusion, fit.theta)


def vector_error_covariances(observed, errors, fit=None):
    if getattr(fit, "measurement_covariances", None) is not None:
        return fit.measurement_covariances
    if errors is None:
        return None
    return {
        name: np.diag(
            [
                0.0 if value is None else float(error) ** 2
                for value, error in zip(values, errors[name], strict=True)
            ]
        )
        for name, values in observed.items()
        if values is not None
    }


def write_multivariate_samples(tree, observed, errors, args, settings, fit):
    path = getattr(args, "posterior_samples_out", None)
    if path in (None, ""):
        return
    from nwkit.asr_continuous_diagnostics import _positive_count

    count = _positive_count(
        getattr(args, "posterior_samples", None), "--posterior-samples", 1000
    )
    process = fitted_vector_process(tree, settings.model, fit)
    posterior = condition_vector_tree(
        process,
        observed,
        error_covariances=vector_error_covariances(observed, errors, fit),
    )
    samples = posterior.sample(count, seed=getattr(args, "seed", None))
    ids = assign_branch_ids(tree)
    rows = [
        {
            "sample": sample,
            "branch_id": ids[node],
            "parent": -1 if node.is_root else ids[node.up],
            "node_class": get_node_class(node),
            "name": "" if node.name in (None, "") else str(node.name),
            "trait": trait,
            "value": samples[sample, index, k],
        }
        for sample in range(count)
        for index, node in enumerate(posterior.nodes)
        for k, trait in enumerate(fit.trait_names)
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def vector_predictive_statistics(tree, observed, trait_names):
    from nwkit.asr_diagnostics import _summary_statistics
    from nwkit.phylogenetic_predictive import sister_clade_contrast

    result = {}
    for k, trait in enumerate(trait_names):
        values = {
            name: vector[k]
            for name, vector in observed.items()
            if vector is not None and vector[k] is not None
        }
        if len(values) < 2:
            raise ValueError(
                "Vector predictive checks require two observations per trait."
            )
        for statistic, value in _summary_statistics(list(values.values())).items():
            result[(trait, "", statistic)] = value
        result[(trait, "", "sister_clade_mean_squared_difference")] = (
            sister_clade_contrast(tree, values)
        )
        for j in range(k):
            pairs = np.array(
                [
                    [vector[k], vector[j]]
                    for vector in observed.values()
                    if vector is not None
                    and vector[k] is not None
                    and vector[j] is not None
                ],
                dtype=float,
            )
            if len(pairs) < 2:
                continue
            result[(trait, trait_names[j], "covariance")] = float(
                np.cov(pairs.T, ddof=1)[0, 1]
            )
    return result


def vector_posterior_predictive(
    process, observed, trait_names, *, errors=None, num_simulations=1000, seed=None
):
    from nwkit.phylogenetic_predictive import predictive_summary
    from nwkit.vector_simulation import simulated_vector_observations

    posterior = condition_vector_tree(process, observed, error_covariances=errors)
    root_seed, simulation_seed = np.random.SeedSequence(seed).spawn(2)
    roots = np.random.default_rng(root_seed).multivariate_normal(
        posterior.means[0],
        posterior.covariances[0],
        size=num_simulations,
        check_valid="raise",
    )
    replicates = simulated_vector_observations(
        process,
        observed,
        errors,
        num_simulations,
        int(simulation_seed.generate_state(1)[0]),
        root_values=roots,
    )
    statistics = vector_predictive_statistics(process.tree, observed, trait_names)
    replicated = [
        vector_predictive_statistics(process.tree, values, trait_names)
        for values in replicates
    ]
    return pd.DataFrame(
        [
            {
                "trait": key[0],
                "other_trait": key[1],
                **predictive_summary(key[2], value, [row[key] for row in replicated]),
            }
            for key, value in statistics.items()
        ]
    )


def refit_vector_model(tree, values, errors, args, settings, fit):
    if settings.model == "MV-BM":
        from nwkit.multivariate_asr import compute_mvbm_marginals

        return compute_mvbm_marginals(
            tree,
            values,
            fit.trait_names,
            standard_errors=errors,
            measurement_covariances=getattr(fit, "measurement_covariances", None),
            compute_posterior=False,
        )[1]
    from nwkit.multivariate_gaussian_asr import fit_dense_mvou, fit_dense_mvou_diag

    keywords = dict(
        standard_errors=errors,
        compute_posterior=False,
        measurement_covariances=getattr(fit, "measurement_covariances", None),
    )
    if settings.model == "MV-OU-FULL":
        from nwkit.full_ou_fit import fit_full_mvou

        return fit_full_mvou(
            tree,
            values,
            fit.trait_names,
            attraction=None if fit.attraction_estimated else fit.attraction_matrix,
            diffusion=None if fit.attraction_estimated else fit.diffusion_sigma,
            **keywords,
        )[1]
    if getattr(args, "alpha_bounds", None) not in (None, ""):
        from nwkit.ou_asr import parse_alpha_bounds

        keywords["alpha_bounds"] = parse_alpha_bounds(args.alpha_bounds, tree)
    if settings.model == "MV-OU":
        return fit_dense_mvou(
            tree,
            values,
            fit.trait_names,
            alpha=None if fit.alpha_estimated else fit.alpha,
            **keywords,
        )[1]
    return fit_dense_mvou_diag(
        tree,
        values,
        fit.trait_names,
        alpha_by_trait=None if fit.alpha_estimated else fit.alpha_by_trait,
        **keywords,
    )[1]


def vector_fit_row(fit):
    row = {"model": getattr(fit, "model", "MV-BM"), "fit_model_status": fit.fit_status}
    for name in (
        "sigma",
        "diffusion_sigma",
        "theta",
        "alpha_by_trait",
        "attraction_matrix",
    ):
        value = getattr(fit, name, None)
        if value is not None:
            for indices, scalar in np.ndenumerate(value):
                row[name + "_" + "_".join(map(str, indices))] = float(scalar)
    for name in ("alpha", "log_likelihood", "restricted_log_likelihood"):
        value = getattr(fit, name, None)
        if value is not None:
            row[name] = float(value)
    return row


def write_multivariate_diagnostics(tree, observed, errors, args, settings, fit):
    from nwkit.asr_continuous_diagnostics import _positive_count
    from nwkit.asr_diagnostics import parametric_bootstrap
    from nwkit.vector_simulation import simulated_vector_observations

    write_multivariate_samples(tree, observed, errors, args, settings, fit)
    predictive_path = getattr(args, "posterior_predictive_out", None)
    bootstrap_path = getattr(args, "bootstrap_out", None)
    if predictive_path in (None, "") and bootstrap_path in (None, ""):
        return
    process = fitted_vector_process(tree, settings.model, fit)
    covariance = vector_error_covariances(observed, errors, fit)
    predictive_seed, bootstrap_seed = np.random.SeedSequence(
        getattr(args, "seed", None)
    ).spawn(2)
    if predictive_path not in (None, ""):
        vector_posterior_predictive(
            process,
            observed,
            fit.trait_names,
            errors=covariance,
            num_simulations=_positive_count(
                getattr(args, "posterior_predictive_simulations", None),
                "--posterior-predictive-simulations",
                1000,
            ),
            seed=int(predictive_seed.generate_state(1)[0]),
        ).to_csv(predictive_path, sep="\t", index=False)
    if bootstrap_path not in (None, ""):
        root = (
            condition_vector_tree(
                process, observed, error_covariances=covariance
            ).means[0]
            if process.root_mean is None
            else None
        )
        rows, failures = parametric_bootstrap(
            lambda seed: simulated_vector_observations(
                process, observed, covariance, 1, seed, root_values=root
            )[0],
            lambda values: refit_vector_model(
                tree, values, errors, args, settings, fit
            ),
            vector_fit_row,
            num_simulations=_positive_count(
                getattr(args, "bootstrap_simulations", None),
                "--bootstrap-simulations",
                100,
            ),
            seed=int(bootstrap_seed.generate_state(1)[0]),
        )
        if len(failures):
            rows = pd.concat((rows, failures), ignore_index=True, sort=False)
        rows.sort_values("replicate").to_csv(bootstrap_path, sep="\t", index=False)
