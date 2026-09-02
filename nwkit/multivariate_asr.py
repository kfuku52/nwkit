"""Linear-time complete-case multivariate Brownian ancestral reconstruction."""

import math
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
from nwkit.continuous_asr_io import _summary
from nwkit.util import assign_branch_ids, get_node_class, write_tree

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
            variance = variance * child_variance / total
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
    )


def compute_mvbm_marginals(
    tree,
    values_by_leaf,
    trait_names,
    *,
    standard_errors=None,
    _tree_validated=False,
):
    """Estimate a trait covariance and return all-node MV-BM marginals.

    Observed tips must contain either every selected trait or none. Known
    measurement errors are intentionally excluded because arbitrary errors
    break the Kronecker structure used by the linear-time implementation.
    """

    trait_names = tuple(str(value) for value in trait_names)
    partial = any(
        vector is not None
        and any(value is None for value in vector)
        and not all(value is None for value in vector)
        for vector in values_by_leaf.values()
    )
    if standard_errors is not None:
        try:
            has_nonzero_error = any(
                error is not None and float(error) != 0.0
                for vector in standard_errors.values()
                if vector is not None
                for error in vector
            )
        except (TypeError, ValueError, OverflowError):
            # Let the dense implementation report the precise validation error.
            has_nonzero_error = True
        if not has_nonzero_error:
            standard_errors = None
    if partial or standard_errors is not None:
        from nwkit.multivariate_gaussian_asr import fit_dense_mvbm

        return fit_dense_mvbm(
            tree,
            values_by_leaf,
            trait_names,
            standard_errors=standard_errors,
        )
    vectors = _validated_vectors(values_by_leaf, trait_names)
    scaled, exponents = _scaled_vectors(vectors, trait_names)
    contrasts, log_variance, num_observed, num_effective = _independent_contrasts(
        tree, scaled, tree_validated=_tree_validated
    )
    residual_df = len(contrasts)
    if residual_df == 0:
        raise ValueError(
            "MV-BM evolutionary covariance is not identifiable from one distinct observed position."
        )
    sigma_scaled = contrasts.T @ contrasts / residual_df
    sigma_scaled = (sigma_scaled + sigma_scaled.T) / 2.0
    scales = np.asarray([math.ldexp(1.0, int(value)) for value in exponents])
    sigma = sigma_scaled * scales[:, None] * scales[None, :]
    if np.any(~np.isfinite(sigma)):
        raise ValueError(
            "MV-BM covariance exceeds floating-point range; rescale units."
        )
    eigenvalues = np.linalg.eigvalsh(sigma_scaled)
    tolerance = (
        np.finfo(float).eps
        * max(1.0, float(np.max(eigenvalues)))
        * max(sigma_scaled.shape)
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
    ids = assign_branch_ids(tree)
    rows = []
    for node in selected_nodes:
        observed_vector = observed.get(str(node.name)) if node.is_leaf else None
        marginal = posterior[node]
        for trait_index, trait in enumerate(trait_names):
            value = None if observed_vector is None else observed_vector[trait_index]
            error_vector = errors.get(str(node.name)) if errors is not None else None
            error = None if error_vector is None else error_vector[trait_index]
            scalar = GaussianMarginal(
                float(marginal.mean[trait_index]),
                float(marginal.covariance[trait_index, trait_index]),
            )
            rows.append(
                {
                    "branch_id": ids[node],
                    "parent": -1 if node.is_root else ids[node.up],
                    "node_class": get_node_class(node),
                    "name": "" if node.name in (None, "") else str(node.name),
                    "trait": trait,
                    "observed_value": "" if value is None else value,
                    "observed_se": ""
                    if value is None
                    else (0.0 if error is None else error),
                    "is_imputed": bool(node.is_leaf and value is None),
                    **_summary(scalar, ci_level),
                }
            )
    columns = [
        "branch_id",
        "parent",
        "node_class",
        "name",
        "trait",
        "observed_value",
        "observed_se",
        "is_imputed",
        "mean",
        "variance",
        "sd",
        "ci_lower",
        "ci_upper",
        "ci_level",
    ]
    return pd.DataFrame(rows, columns=columns)


def multivariate_covariance_table(tree, selected_nodes, posterior, trait_names):
    ids = assign_branch_ids(tree)
    rows = []
    for node in selected_nodes:
        covariance = posterior[node].covariance
        diagonal = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        for first in range(len(trait_names)):
            for second in range(first, len(trait_names)):
                denominator = diagonal[first] * diagonal[second]
                correlation: str | float
                if denominator == 0.0:
                    correlation = ""
                else:
                    raw_correlation = float(covariance[first, second] / denominator)
                    correlation = min(1.0, max(-1.0, raw_correlation))
                rows.append(
                    {
                        "branch_id": ids[node],
                        "parent": -1 if node.is_root else ids[node.up],
                        "node_class": get_node_class(node),
                        "name": "" if node.name in (None, "") else str(node.name),
                        "trait_1": trait_names[first],
                        "trait_2": trait_names[second],
                        "covariance": float(covariance[first, second]),
                        "correlation": correlation,
                    }
                )
    return pd.DataFrame(
        rows,
        columns=[
            "branch_id",
            "parent",
            "node_class",
            "name",
            "trait_1",
            "trait_2",
            "covariance",
            "correlation",
        ],
    )


def multivariate_model_table(fit, args, ci_level):
    model = getattr(fit, "model", "MV-BM")
    is_ou = model == "MV-OU"
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
    for first, first_trait in enumerate(fit.trait_names):
        for second, second_trait in enumerate(fit.trait_names):
            row[f"sigma_{_trait_id(first_trait)}_to_{_trait_id(second_trait)}"] = float(
                fit.sigma[first, second]
            )
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
