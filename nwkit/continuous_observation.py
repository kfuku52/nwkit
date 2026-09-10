"""Explicit correlated Gaussian measurement-error input contracts."""

import csv
import math

import numpy as np

from nwkit.vector_gaussian import _array, _inverse


def validate_measurement_covariances(observed, trait_names, covariances):
    """Validate full matrices; exact coordinates must have zero rows/columns."""
    dimension = len(trait_names)
    if set(covariances) - set(observed):
        raise ValueError("Measurement covariance contains unknown tips.")
    result = {}
    for name, vector in observed.items():
        if vector is None or all(value is None for value in vector):
            continue
        if name not in covariances:
            raise ValueError(f"Measurement covariance is required for tip '{name}'.")
        matrix = _array(
            covariances[name], (dimension, dimension), "Measurement covariance"
        )
        if not np.allclose(matrix, matrix.T, rtol=1e-12, atol=0):
            raise ValueError("Measurement covariance must be symmetric.")
        exact = np.flatnonzero(np.diag(matrix) == 0)
        if len(exact) and np.any(matrix[exact] != 0):
            raise ValueError(
                "Exact measurement coordinates require zero covariance rows."
            )
        noisy = np.flatnonzero(np.diag(matrix) != 0)
        if len(noisy):
            _inverse(matrix[np.ix_(noisy, noisy)], "Measurement covariance noisy block")
        result[name] = matrix.copy()
    return result


def read_measurement_covariances(path, observed, trait_names):
    """Read every ordered matrix entry from a long-format TSV."""
    names = tuple(trait_names)
    matrices: dict[str, np.ndarray] = {}
    seen = set()
    with open(path, newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["leaf_name", "trait", "other_trait", "covariance"]:
            raise ValueError(
                "Measurement covariance TSV requires leaf_name, trait, other_trait, covariance columns in that order."
            )
        for row in reader:
            name, first, second = row["leaf_name"], row["trait"], row["other_trait"]
            if name not in observed or first not in names or second not in names:
                raise ValueError(
                    "Measurement covariance contains an unknown tip or trait."
                )
            key = name, first, second
            if key in seen:
                raise ValueError("Duplicated measurement covariance entry.")
            seen.add(key)
            matrix = matrices.setdefault(
                name, np.full((len(names), len(names)), np.nan)
            )
            matrix[names.index(first), names.index(second)] = float(row["covariance"])
    return validate_measurement_covariances(observed, names, matrices)


def prepare_correlated_observations(
    tree, observed, trait_names, covariances, standard_errors=None
):
    from nwkit.multivariate_gaussian_asr import _prepare_observations
    from nwkit.vector_fit import normalized_vector_observations

    if standard_errors is not None:
        raise ValueError(
            "Full measurement covariance cannot be combined with standard errors."
        )
    covariances = validate_measurement_covariances(observed, trait_names, covariances)
    errors = {name: np.sqrt(np.diag(matrix)) for name, matrix in covariances.items()}
    for name, vector in observed.items():
        if vector is not None and all(value is None for value in vector):
            errors[name] = np.zeros(len(trait_names))
    data = _prepare_observations(tree, observed, trait_names, errors, dense_limit=False)
    values, normalized = normalized_vector_observations(data)
    for name in normalized:
        normalized[name] = covariances[name] / np.outer(data.scales, data.scales)
    return data, values, normalized, covariances


def summarize_independent_replicates(values, standard_errors):
    """Return mean, SE and log constant in product N(y_i|x,s_i²)=c N(mean|x,se²).

    SEs describe known total within-species observation variability. They are
    not estimated from these replicates. Noisy observations must have positive
    SEs, so the input always defines an ordinary density.
    """
    values = np.asarray(values, dtype=float)
    errors = np.asarray(standard_errors, dtype=float)
    if (
        values.ndim != 1
        or not len(values)
        or errors.shape != values.shape
        or not np.isfinite(values).all()
        or not np.isfinite(errors).all()
        or np.any(errors <= 0)
    ):
        raise ValueError(
            "Replicate observations need finite values and positive finite standard errors."
        )
    scale = float(np.min(errors))
    weights = (scale / errors) ** 2
    weight = float(np.sum(weights))
    center = float(values[0])
    mean = center + float(np.dot(weights, values - center)) / weight
    error = scale / math.sqrt(weight)
    residual = (values - mean) / errors
    log_constant = -0.5 * (
        (len(values) - 1) * math.log(2 * math.pi)
        + 2 * float(np.sum(np.log(errors)))
        - 2 * math.log(error)
        + float(residual @ residual)
    )
    if not all(math.isfinite(v) for v in (mean, error, log_constant)) or error <= 0:
        raise ValueError("Replicate summary exceeds floating-point range.")
    return mean, error, log_constant


def apply_replicate_observations(observed, errors, trait_names, args):
    """Replace supplied tip/trait cells by Gaussian replicate sufficient statistics."""
    path = getattr(args, "replicate_observations", None)
    if path in (None, ""):
        return observed, errors
    if getattr(args, "measurement_covariance", None) not in (None, ""):
        raise ValueError(
            "--replicate-observations cannot be combined with --measurement-covariance."
        )
    scalar = len(trait_names) == 1
    dimension = len(trait_names)
    values = {
        name: [vector]
        if scalar
        else list(vector)
        if vector is not None
        else [None] * dimension
        for name, vector in observed.items()
    }
    uncertainty = {
        name: [0.0] * dimension
        if errors is None or all(value is None for value in values[name])
        else [errors[name]]
        if scalar
        else list(errors[name])
        for name in observed
    }
    records: dict[tuple[str, str], list[tuple[float, float]]] = {}
    with open(path, newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["leaf_name", "trait", "value", "standard_error"]:
            raise ValueError(
                "Replicate TSV requires leaf_name, trait, value, standard_error columns in that order."
            )
        for row in reader:
            key = row["leaf_name"], row["trait"]
            if key[0] not in values or key[1] not in trait_names:
                raise ValueError(
                    "Replicate observations contain unknown tips or traits."
                )
            records.setdefault(key, []).append(
                (float(row["value"]), float(row["standard_error"]))
            )
    if not records:
        raise ValueError("Replicate observation table is empty.")
    constants = []
    for (name, trait), entries in records.items():
        mean, error, constant = summarize_independent_replicates(
            *zip(*entries, strict=True)
        )
        index = trait_names.index(trait)
        values[name][index], uncertainty[name][index] = mean, error
        constants.append(constant)
    args._replicate_log_constant = math.fsum(constants)
    if scalar:
        return {name: vector[0] for name, vector in values.items()}, {
            name: vector[0] for name, vector in uncertainty.items()
        }
    return values, uncertainty


def restore_replicate_likelihood(fit, args):
    from dataclasses import fields, replace

    constant = getattr(args, "_replicate_log_constant", 0.0)
    if not constant:
        return fit
    likelihood_fields = {item.name for item in fields(fit)} & {
        "log_likelihood",
        "restricted_log_likelihood",
    }
    changes = {
        name: getattr(fit, name) + constant
        for name in likelihood_fields
        if getattr(fit, name) is not None
    }
    return replace(fit, **changes)
