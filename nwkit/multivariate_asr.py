"""Linear-time complete-case multivariate Brownian ancestral reconstruction."""

import math
from collections.abc import Mapping
from copy import copy
from dataclasses import dataclass

import numpy as np
import pandas as pd

from nwkit.continuous_asr import (
    GaussianMarginal,
    _finite_number,
    _tree_groups,
    compute_bm_marginals,
)
from nwkit.continuous_asr_io import (
    _asr_node_metadata,
    _asr_table_from_columns,
    _marginal_columns,
    _summary,
)
from nwkit.util import assign_branch_ids, write_tree

_LOG_2PI = math.log(2.0 * math.pi)
_LOG_2 = math.log(2.0)


@dataclass(frozen=True)
class MultivariateGaussianMarginal:
    mean: np.ndarray
    covariance: np.ndarray


@dataclass(frozen=True)
class MultivariateBrownianFit:
    trait_names: tuple[str, ...]
    sigma: np.ndarray
    sigma_rank: int
    sigma_estimated: bool
    restricted_log_likelihood: float | None
    num_observed: int
    num_effective_observations: int
    residual_df: int
    fit_status: str


def _validated_vectors(values_by_leaf, trait_names):
    dimension = len(trait_names)
    if dimension < 2:
        raise ValueError("MV-BM requires at least two comma-separated trait columns.")
    result: dict[str, np.ndarray | None] = {}
    partial: list[str] = []
    for name, raw in values_by_leaf.items():
        if raw is None:
            result[name] = None
            continue
        if len(raw) != dimension:
            raise ValueError(f"MV-BM trait dimension mismatch for tip '{name}'.")
        missing = [value is None for value in raw]
        if any(missing):
            if not all(missing):
                partial.append(name)
            result[name] = None
            continue
        result[name] = np.asarray(
            [
                _finite_number(value, f"Trait '{trait}' for '{name}'")
                for trait, value in zip(trait_names, raw, strict=True)
            ],
            dtype=float,
        )
    if partial:
        raise ValueError(
            "MV-BM currently requires complete trait vectors at observed tips; "
            "partially missing tip(s): " + ", ".join(sorted(partial))
        )
    if not any(value is not None for value in result.values()):
        raise ValueError("MV-BM requires at least one complete observed tip vector.")
    return result


def _scaled_vectors(vectors, trait_names):
    observed = [value for value in vectors.values() if value is not None]
    center = observed[0].copy()
    differences = [value - center for value in observed]
    if any(np.any(~np.isfinite(value)) for value in differences):
        raise ValueError(
            "MV-BM trait ranges exceed floating-point range; rescale units."
        )
    exponents = []
    for trait_index in range(len(trait_names)):
        size = max(abs(value[trait_index]) for value in differences)
        exponents.append(math.frexp(size)[1] - 1 if size > 0.0 else 0)
    scales = np.asarray([math.ldexp(1.0, exponent) for exponent in exponents])
    return {
        name: None if value is None else (value - center) / scales
        for name, value in vectors.items()
    }, np.asarray(exponents, dtype=int)


def _independent_contrasts(tree, vectors, *, tree_validated=False):
    node_groups, parents, lengths = _tree_groups(tree, validated=tree_validated)
    dimension = len(next(value for value in vectors.values() if value is not None))
    local = [None] * len(parents)
    num_observed = 0
    longest = max(lengths)
    time_exponent = math.frexp(longest)[1] - 1 if longest else 0
    scaled_lengths = [math.ldexp(length, -time_exponent) for length in lengths]
    if any(
        scaled_length == 0.0 < length
        for scaled_length, length in zip(scaled_lengths, lengths, strict=True)
    ):
        time_exponent = 0
    else:
        lengths = scaled_lengths
    for leaf in tree.leaves():
        value = vectors.get(str(leaf.name))
        if value is None:
            continue
        num_observed += 1
        group = node_groups[leaf]
        if local[group] is not None:
            if not np.array_equal(local[group], value):
                raise ValueError(
                    "Conflicting exact multivariate observations are connected by zero-length branches."
                )
            continue
        local[group] = value
    children: list[list[int]] = [[] for _ in parents]
    for child, parent in enumerate(parents):
        if parent >= 0:
            children[parent].append(child)
    inside_mean = list(local)
    inside_variance = [0.0 if value is not None else None for value in local]
    contrasts = []
    log_variances = []
    for group in range(len(parents) - 1, -1, -1):
        mean = inside_mean[group]
        variance = inside_variance[group]
        for child in children[group]:
            child_mean = inside_mean[child]
            if child_mean is None:
                continue
            child_variance = inside_variance[child] + lengths[child]
            if mean is None:
                mean, variance = child_mean, child_variance
                continue
            total = variance + child_variance
            if total == 0.0:
                if not np.array_equal(mean, child_mean):
                    raise ValueError(
                        "Conflicting exact multivariate observations have zero likelihood."
                    )
                continue
            if not math.isfinite(total) or total < 0.0:
                raise ValueError("MV-BM contrast variance is not positive and finite.")
            contrasts.append((mean - child_mean) / math.sqrt(total))
            log_variances.append(math.log(total))
            left_weight = child_variance / total
            right_weight = variance / total
            mean = left_weight * mean + right_weight * child_mean
            variance = min(variance, child_variance) * (
                max(variance, child_variance) / total
            )
        inside_mean[group], inside_variance[group] = mean, variance
    if inside_mean[0] is None:
        raise ValueError("MV-BM flat root prior requires observed tips.")
    if contrasts:
        matrix = np.vstack(contrasts)
    else:
        matrix = np.empty((0, dimension), dtype=float)
    return (
        matrix,
        math.fsum(log_variances),
        num_observed,
        sum(value is not None for value in local),
        time_exponent,
    )


def _measurement_errors_are_all_zero(
    tree, values_by_leaf, trait_names, standard_errors
):
    """Validate complete-case errors before deciding whether the fast path applies."""

    if not isinstance(standard_errors, Mapping):
        raise ValueError("MV-BM measurement errors must be supplied as a tip mapping.")
    dimension = len(trait_names)
    has_nonzero = False
    for leaf in tree.leaves():
        name = str(leaf.name)
        vector = values_by_leaf.get(name)
        if vector is None:
            continue
        try:
            vector_dimension = len(vector)
        except TypeError as exc:
            raise ValueError(
                f"MV-BM trait dimension mismatch for tip '{name}'."
            ) from exc
        if vector_dimension != dimension:
            raise ValueError(f"MV-BM trait dimension mismatch for tip '{name}'.")
        error_vector = standard_errors.get(name)
        if error_vector is None:
            raise ValueError(
                f"Measurement errors are required for observed tip '{name}'."
            )
        try:
            error_dimension = len(error_vector)
        except TypeError as exc:
            raise ValueError(
                f"Measurement-error dimension mismatch for tip '{name}'."
            ) from exc
        if error_dimension != dimension:
            raise ValueError(f"Measurement-error dimension mismatch for tip '{name}'.")
        for trait, raw_value, raw_error in zip(
            trait_names, vector, error_vector, strict=True
        ):
            if raw_value is None:
                continue
            error = _finite_number(
                raw_error,
                f"Standard error for trait '{trait}' at '{name}'",
                nonnegative=True,
            )
            has_nonzero = has_nonzero or error != 0.0
    return not has_nonzero


def compute_mvbm_marginals(
    tree,
    values_by_leaf,
    trait_names,
    *,
    standard_errors=None,
    _tree_validated=False,
    compute_posterior=True,
    _geometry_cache=None,
):
    """Estimate a trait covariance and return all-node MV-BM marginals.

    Complete exact observations use the linear-time path. Partial vectors or
    non-zero known measurement errors use the dense Gaussian implementation;
    complete all-zero error mappings retain the linear-time path after validation.
    """

    trait_names = tuple(str(value) for value in trait_names)
    partial = any(
        vector is not None
        and any(value is None for value in vector)
        and not all(value is None for value in vector)
        for vector in values_by_leaf.values()
    )
    if standard_errors is not None:
        if _measurement_errors_are_all_zero(
            tree, values_by_leaf, trait_names, standard_errors
        ):
            standard_errors = None
    if partial or standard_errors is not None:
        from nwkit.multivariate_gaussian_asr import fit_dense_mvbm

        return fit_dense_mvbm(
            tree,
            values_by_leaf,
            trait_names,
            standard_errors=standard_errors,
            compute_posterior=compute_posterior,
            _geometry_cache=_geometry_cache,
        )
    vectors = _validated_vectors(values_by_leaf, trait_names)
    scaled, exponents = _scaled_vectors(vectors, trait_names)
    contrasts, log_variance, num_observed, num_effective, time_exponent = (
        _independent_contrasts(tree, scaled, tree_validated=_tree_validated)
    )
    residual_df = len(contrasts)
    if residual_df == 0:
        raise ValueError(
            "MV-BM evolutionary covariance is not identifiable from one distinct observed position."
        )
    sigma_scaled = contrasts.T @ contrasts / residual_df
    sigma_scaled = (sigma_scaled + sigma_scaled.T) / 2.0
    sigma = np.ldexp(
        sigma_scaled, exponents[:, None] + exponents[None, :] - time_exponent
    )
    if np.any(~np.isfinite(sigma)):
        raise ValueError(
            "MV-BM covariance exceeds floating-point range; rescale units."
        )
    eigenvalues = np.linalg.eigvalsh(sigma_scaled)
    tolerance = (
        np.finfo(float).eps * float(np.max(eigenvalues)) * max(sigma_scaled.shape)
    )
    rank = int(np.sum(eigenvalues > tolerance))
    if rank < len(trait_names):
        log_likelihood = None
        fit_status = "singular_covariance"
    else:
        sign, log_determinant = np.linalg.slogdet(sigma_scaled)
        if sign <= 0.0 or not math.isfinite(float(log_determinant)):
            raise ValueError("MV-BM fitted covariance is not positive definite.")
        inverse = np.linalg.inv(sigma_scaled)
        quadratic = float(np.einsum("ni,ij,nj->", contrasts, inverse, contrasts))
        scaled_log_likelihood = -0.5 * (
            residual_df * len(trait_names) * _LOG_2PI
            + len(trait_names) * log_variance
            + residual_df * float(log_determinant)
            + quadratic
        )
        log_likelihood = scaled_log_likelihood - residual_df * math.fsum(
            int(value) * _LOG_2 for value in exponents
        )
        fit_status = "ok"
    posterior = (
        _mvbm_posterior(tree, vectors, trait_names, sigma) if compute_posterior else {}
    )
    fit = MultivariateBrownianFit(
        trait_names=trait_names,
        sigma=sigma,
        sigma_rank=rank,
        sigma_estimated=True,
        restricted_log_likelihood=log_likelihood,
        num_observed=num_observed,
        num_effective_observations=num_effective,
        residual_df=residual_df,
        fit_status=fit_status,
    )
    return posterior, fit


def _mvbm_posterior(tree, vectors, trait_names, sigma):
    scalar_results = []
    for trait_index in range(len(trait_names)):
        scalar_values = {
            name: None if value is None else float(value[trait_index])
            for name, value in vectors.items()
        }
        posterior, _ = compute_bm_marginals(
            tree, scalar_values, sigma2=1.0, _tree_validated=True
        )
        scalar_results.append(posterior)
    posterior = {}
    for node in tree.traverse():
        mean = np.asarray([result[node].mean for result in scalar_results], dtype=float)
        variance_factor = scalar_results[0][node].variance
        posterior[node] = MultivariateGaussianMarginal(
            mean=mean, covariance=variance_factor * sigma
        )
    return posterior


def _trait_id(trait):
    return str(trait).encode("utf-8").hex()


def multivariate_output_table(
    tree,
    selected_nodes,
    observed,
    posterior,
    trait_names,
    ci_level,
    *,
    errors=None,
):
    selected_nodes = tuple(selected_nodes)
    traits = np.asarray(trait_names, dtype=object)
    trait_count = len(traits)
    metadata = _asr_node_metadata(selected_nodes, assign_branch_ids(tree))
    columns = {
        **{name: np.repeat(values, trait_count) for name, values in metadata.items()},
        "trait": np.tile(traits, len(selected_nodes)),
        **_multivariate_observed_columns(selected_nodes, observed, errors, trait_count),
        **_multivariate_marginal_columns(
            selected_nodes, posterior, trait_count, ci_level
        ),
    }
    return _asr_table_from_columns(columns)


def _multivariate_observed_columns(nodes, observed, errors, trait_count):
    shape = (len(nodes), trait_count)
    value = np.full(shape, "", dtype=object)
    se = np.full(shape, "", dtype=object)
    imputed = np.zeros(shape, dtype=bool)
    for index, node in enumerate(nodes):
        if not node.is_leaf:
            continue
        name = str(node.name)
        vector = observed.get(name)
        error = None if errors is None else errors.get(name)
        for trait in range(trait_count):
            measurement = None if vector is None else vector[trait]
            imputed[index, trait] = measurement is None
            if measurement is not None:
                value[index, trait] = measurement
                deviation = None if error is None else error[trait]
                se[index, trait] = 0.0 if deviation is None else deviation
    return {
        "observed_value": _numeric_observations(value.ravel()),
        "observed_se": _numeric_observations(se.ravel()),
        "is_imputed": imputed.ravel(),
    }


def _numeric_observations(values):
    return (
        pd.Series(values.tolist()).to_numpy()
        if len(values) and not np.any(values == "")
        else values
    )


def _multivariate_marginal_columns(nodes, posterior, trait_count, ci_level):
    means = np.empty((len(nodes), trait_count), dtype=float)
    variances = np.empty_like(means)
    for index, node in enumerate(nodes):
        marginal = posterior[node]
        means[index] = marginal.mean[:trait_count]
        variances[index] = np.diag(marginal.covariance)[:trait_count]
    return _marginal_columns(means.ravel(), variances.ravel(), ci_level)


def multivariate_covariance_table(tree, selected_nodes, posterior, trait_names):
    ids = assign_branch_ids(tree)
    selected_nodes = tuple(selected_nodes)
    traits = np.asarray(trait_names, dtype=object)
    first, second = np.triu_indices(len(traits))
    pair_count = len(first)
    node_count = len(selected_nodes)
    metadata = _asr_node_metadata(selected_nodes, ids)
    covariance = np.empty((node_count, pair_count), dtype=float)
    correlation = np.empty((node_count, pair_count), dtype=object)
    for index, node in enumerate(selected_nodes):
        matrix = posterior[node].covariance
        covariance[index] = matrix[first, second]
        correlation[index] = _upper_triangle_correlations(
            matrix, covariance[index], first, second
        )
    covariance_column = covariance.ravel()
    correlation_column = correlation.ravel()
    if correlation_column.size and not np.any(correlation_column == ""):
        correlation_column = correlation_column.astype(float)
    return _asr_table_from_columns(
        {
            **{
                name: np.repeat(values, pair_count) for name, values in metadata.items()
            },
            "trait_1": np.tile(traits[first], node_count),
            "trait_2": np.tile(traits[second], node_count),
            "covariance": covariance_column,
            "correlation": correlation_column,
        }
    )


def _upper_triangle_correlations(matrix, covariance, first, second):
    sd = np.sqrt(np.maximum(np.diag(matrix), 0.0))
    denominator = sd[first] * sd[second]
    defined = denominator != 0.0
    correlation = np.full(len(first), "", dtype=object)
    correlation[defined] = np.clip(
        covariance[defined] / denominator[defined], -1.0, 1.0
    )
    return correlation


def multivariate_model_table(fit, args, ci_level):
    model = getattr(fit, "model", "MV-BM")
    is_ou = model in {"MV-OU", "MV-OU-DIAG"}
    row = {
        "trait_type": "continuous",
        "trait_type_requested": getattr(args, "trait_type", "auto"),
        "trait": ",".join(fit.trait_names),
        "model": model,
        "root_prior": "stationary" if is_ou else "flat",
        "sigma_estimated": fit.sigma_estimated,
        "estimation_method": "ML" if is_ou else "REML",
        "restricted_log_likelihood": fit.restricted_log_likelihood,
        "log_likelihood": getattr(fit, "log_likelihood", None),
        "likelihood_kind": (
            "multivariate_stationary_root_ml"
            if is_ou
            else "multivariate_flat_root_integrated"
        ),
        "num_observed": fit.num_observed,
        "num_effective_observations": fit.num_effective_observations,
        "residual_df": fit.residual_df,
        "sigma_rank": fit.sigma_rank,
        "fit_status": fit.fit_status,
        "standard_error_column": getattr(args, "standard_error_column", None) or "",
        "ci_level": ci_level,
        "interval_kind": "conditional_on_covariance",
        "parameter_uncertainty_included": False,
        "tree_uncertainty_included": False,
    }
    if is_ou:
        row["alpha"] = fit.alpha
        row["alpha_estimated"] = fit.alpha_estimated
        for index, trait in enumerate(fit.trait_names):
            row[f"theta_{_trait_id(trait)}"] = float(fit.theta[index])
            if model == "MV-OU-DIAG":
                row[f"alpha_{_trait_id(trait)}"] = float(fit.alpha_by_trait[index])
    row["sigma_interpretation"] = (
        "stationary_trait_covariance" if is_ou else "diffusion_covariance"
    )
    for first, first_trait in enumerate(fit.trait_names):
        for second, second_trait in enumerate(fit.trait_names):
            row[f"sigma_{_trait_id(first_trait)}_to_{_trait_id(second_trait)}"] = float(
                fit.sigma[first, second]
            )
            if model == "MV-OU-DIAG":
                row[
                    f"diffusion_sigma_{_trait_id(first_trait)}_to_"
                    f"{_trait_id(second_trait)}"
                ] = float(fit.diffusion_sigma[first, second])
    return pd.DataFrame([row])


def write_multivariate_tree(
    tree, observed, posterior, trait_names, args, ci_level, *, errors=None
):
    if getattr(args, "tree_out", None) in (None, ""):
        return
    mean_props = []
    summary_props = []
    observed_props = []
    for node in tree.traverse():
        marginal = posterior[node]
        node.add_props(
            asr_trait_type="continuous",
            asr_model=getattr(args, "model", "MV-BM"),
            asr_ci_level=ci_level,
            asr_interval_kind="conditional_on_covariance",
        )
        observed_vector = observed.get(str(node.name)) if node.is_leaf else None
        for trait_index, trait in enumerate(trait_names):
            trait_id = _trait_id(trait)
            scalar = GaussianMarginal(
                float(marginal.mean[trait_index]),
                float(marginal.covariance[trait_index, trait_index]),
            )
            summary = _summary(scalar, ci_level)
            mean_property = f"asr_mean_{trait_id}"
            mean_props.append(mean_property)
            node.add_props(**{mean_property: summary["mean"]})
            for key in ("variance", "sd", "ci_lower", "ci_upper"):
                prop = f"asr_{key}_{trait_id}"
                summary_props.append(prop)
                node.add_props(**{prop: summary[key]})
            if node.is_leaf:
                prop = f"asr_observed_value_{trait_id}"
                observed_props.append(prop)
                node.add_props(
                    **{
                        prop: ""
                        if observed_vector is None
                        or observed_vector[trait_index] is None
                        else float(observed_vector[trait_index])
                    }
                )
                if errors is not None:
                    error_prop = f"asr_observed_se_{trait_id}"
                    observed_props.append(error_prop)
                    error_vector = errors.get(str(node.name))
                    node.add_props(
                        **{
                            error_prop: ""
                            if error_vector is None or error_vector[trait_index] is None
                            else float(error_vector[trait_index])
                        }
                    )
        if node.is_leaf:
            node.add_props(
                asr_is_imputed="yes"
                if observed_vector is None
                or any(value is None for value in observed_vector)
                else "no"
            )
        for first in range(len(trait_names)):
            for second in range(first + 1, len(trait_names)):
                prop = (
                    f"asr_covariance_{_trait_id(trait_names[first])}_"
                    f"{_trait_id(trait_names[second])}"
                )
                summary_props.append(prop)
                node.add_props(**{prop: float(marginal.covariance[first, second])})
    props = ["asr_trait_type", "asr_model", *dict.fromkeys(mean_props)]
    if args.tree_annotation != "mean":
        props.extend(dict.fromkeys(summary_props))
        props.extend(["asr_ci_level", "asr_interval_kind"])
    if args.tree_annotation == "all":
        props.extend(dict.fromkeys(observed_props))
        props.append("asr_is_imputed")
    output_args = copy(args)
    output_args.outfile = args.tree_out
    write_tree(
        tree,
        output_args,
        format=getattr(args, "tree_outformat", "auto"),
        quiet=True,
        props=props,
    )
