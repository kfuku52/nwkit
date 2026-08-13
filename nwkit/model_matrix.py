"""Typed predictor encoding shared by ordinary and reconciled models."""

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class CategoricalObservation:
    """Categorical replicate proportions and their biological sample size."""

    probabilities: Mapping[str, float]
    n_observations: int


@dataclass(frozen=True)
class ReplicatedObservation:
    """Independent observations sharing one phylogenetic tip effect."""

    values: tuple[float, ...]


@dataclass(frozen=True)
class PredictorTerm:
    """Metadata for one numeric design-matrix column."""

    name: str
    source: str
    kind: str
    level: str = ""
    reference: str = ""
    coding: str = ""


@dataclass(frozen=True)
class PredictorUncertainty:
    """Within-observation covariance among columns from one predictor."""

    source: str
    term_names: tuple[str, ...]
    covariance_by_observation: np.ndarray


@dataclass(frozen=True)
class EncodedPredictors:
    """Numeric predictor traits plus their source-level metadata."""

    values_by_trait: dict[str, dict[str, float]]
    terms: tuple[PredictorTerm, ...]
    groups: dict[str, tuple[str, ...]]
    uncertainties: tuple[PredictorUncertainty, ...]

    @property
    def term_names(self) -> list[str]:
        return [term.name for term in self.terms]

    @property
    def metadata_by_term(self) -> dict[str, PredictorTerm]:
        return {term.name: term for term in self.terms}


@dataclass(frozen=True)
class ResponseSpec:
    """Resolved likelihood and levels for one response trait."""

    name: str
    kind: str
    family: str
    levels: tuple[str, ...] = ()
    reference: str = ""


RESPONSE_FAMILIES = {
    "gaussian": ("continuous", "identity"),
    "binomial": ("categorical", "logit"),
    "multinomial": ("categorical", "reference-logit"),
    "ordinal": ("ordered", "cumulative-logit"),
    "poisson": ("count", "log"),
    "negative-binomial": ("count", "log"),
    "zero-inflated-poisson": ("count", "log/zero-logit"),
    "zero-inflated-negative-binomial": ("count", "log/zero-logit"),
    "hurdle-poisson": ("count", "log/zero-logit"),
    "hurdle-negative-binomial": ("count", "log/zero-logit"),
    "gamma": ("positive", "log"),
    "lognormal": ("positive", "log"),
    "beta": ("proportion", "logit"),
    "beta-binomial": ("proportion", "logit"),
    "censored-gaussian": ("continuous", "identity"),
}


def response_family_link(family: str) -> str:
    """Return the canonical link for one supported response family."""
    try:
        return RESPONSE_FAMILIES[family][1]
    except KeyError as exc:
        raise ValueError("Unsupported response family: {}.".format(family)) from exc


def _validate_auxiliary_applicability(specs, checks):
    for settings, allowed, option in checks:
        for response in settings:
            if response not in specs:
                raise ValueError(
                    "{} names unknown response '{}'.".format(option, response)
                )
            family = specs[response].family
            if family not in allowed:
                raise ValueError(
                    "{} does not apply to response '{}' with family '{}'.".format(
                        option, response, family
                    )
                )


def _validate_positive_response_settings(settings, option):
    for response, value in settings.items():
        try:
            numeric = float(value)  # type: ignore[arg-type]
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(
                "{} for '{}' must be positive.".format(option, response)
            ) from exc
        if not np.isfinite(numeric) or numeric <= 0.0:
            raise ValueError("{} for '{}' must be positive.".format(option, response))


def _validate_zero_probabilities(settings):
    for response, value in settings.items():
        message = (
            "--response-zero-probability for '{}' must lie strictly in (0, 1).".format(
                response
            )
        )
        try:
            numeric = float(value)  # type: ignore[arg-type]
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(message) from exc
        if not np.isfinite(numeric) or not 0.0 < numeric < 1.0:
            raise ValueError(message)


def _require_beta_binomial_trials(specs, trials):
    for response, spec in specs.items():
        if spec.family == "beta-binomial" and response not in trials:
            raise ValueError(
                "Beta-binomial response '{}' requires --response-trials.".format(
                    response
                )
            )


def validate_response_auxiliaries(
    specs: Mapping[str, ResponseSpec],
    *,
    offsets: Mapping[str, object] | None = None,
    trials: Mapping[str, object] | None = None,
    censor_lower: Mapping[str, object] | None = None,
    censor_upper: Mapping[str, object] | None = None,
    dispersions: Mapping[str, object] | None = None,
    zero_probabilities: Mapping[str, object] | None = None,
) -> None:
    """Reject response-family settings that would otherwise be silently ignored."""
    count_families = {
        "poisson",
        "negative-binomial",
        "zero-inflated-poisson",
        "zero-inflated-negative-binomial",
        "hurdle-poisson",
        "hurdle-negative-binomial",
    }
    dispersion_families = {
        "negative-binomial",
        "zero-inflated-negative-binomial",
        "hurdle-negative-binomial",
        "gamma",
        "lognormal",
        "beta",
        "beta-binomial",
        "censored-gaussian",
    }
    zero_families = {
        "zero-inflated-poisson",
        "zero-inflated-negative-binomial",
        "hurdle-poisson",
        "hurdle-negative-binomial",
    }
    offsets = {} if offsets is None else offsets
    trials = {} if trials is None else trials
    censor_lower = {} if censor_lower is None else censor_lower
    censor_upper = {} if censor_upper is None else censor_upper
    dispersions = {} if dispersions is None else dispersions
    zero_probabilities = {} if zero_probabilities is None else zero_probabilities
    _validate_auxiliary_applicability(
        specs,
        [
            (offsets, count_families, "--response-offset"),
            (trials, {"beta-binomial"}, "--response-trials"),
            (censor_lower, {"censored-gaussian"}, "--response-censor-lower"),
            (censor_upper, {"censored-gaussian"}, "--response-censor-upper"),
            (dispersions, dispersion_families, "--response-dispersion"),
            (zero_probabilities, zero_families, "--response-zero-probability"),
        ],
    )
    _validate_positive_response_settings(dispersions, "--response-dispersion")
    _validate_zero_probabilities(zero_probabilities)
    _require_beta_binomial_trials(specs, trials)


def _finite_number(value: object) -> float | None:
    if isinstance(value, (CategoricalObservation, ReplicatedObservation)):
        return None
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError, OverflowError):
        return None
    return number if np.isfinite(number) else None


def _missing_value(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, (CategoricalObservation, ReplicatedObservation, str)):
        return False
    try:
        return bool(np.isnan(value))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return False


def predictor_kind(
    name: str,
    values: Sequence[object],
    *,
    categorical: set[str],
    ordered_levels: Mapping[str, Sequence[str]],
) -> str:
    """Resolve one predictor as continuous, categorical, or ordered."""
    if name in ordered_levels:
        return "ordered"
    if name in categorical or any(
        isinstance(value, CategoricalObservation) for value in values
    ):
        return "categorical"
    return (
        "continuous"
        if all(_finite_number(value) is not None for value in values)
        else "categorical"
    )


def _probabilities(value: object, levels: Sequence[str]) -> np.ndarray:
    if isinstance(value, CategoricalObservation):
        probabilities = np.asarray(
            [float(value.probabilities.get(level, 0.0)) for level in levels],
            dtype=float,
        )
        if (
            not np.isfinite(probabilities).all()
            or np.any(probabilities < 0.0)
            or not np.isclose(probabilities.sum(), 1.0)
        ):
            raise ValueError(
                "Categorical-state probabilities must be non-negative and sum to one."
            )
        return probabilities
    label = str(value)
    if label not in levels:
        raise ValueError(
            "Categorical value '{}' is not an encoded level.".format(label)
        )
    probabilities = np.zeros(len(levels), dtype=float)
    probabilities[levels.index(label)] = 1.0
    return probabilities


def _observed_levels(values: Sequence[object]) -> list[str]:
    levels: set[str] = set()
    for value in values:
        if isinstance(value, CategoricalObservation):
            levels.update(str(level) for level in value.probabilities)
        else:
            levels.add(str(value))
    return sorted(levels)


def _validate_typed_names(
    selected: Sequence[str], configured: set[str], label: str
) -> None:
    unknown = configured - set(selected)
    if unknown:
        raise ValueError(
            "{} options name unknown traits: {}.".format(
                label, ", ".join(sorted(unknown))
            )
        )


def _values_for_samples(
    values_by_trait: Mapping[str, Mapping[str, object]],
    trait: str,
    sample_names: Sequence[str],
    role: str,
) -> list[object]:
    if trait not in values_by_trait:
        raise ValueError("{} trait '{}' is absent.".format(role, trait))
    values_by_sample = values_by_trait[trait]
    missing = sorted(set(sample_names) - set(values_by_sample))
    if missing:
        raise ValueError(
            "{} '{}' is missing samples: {}.".format(role, trait, ", ".join(missing))
        )
    return [values_by_sample[sample] for sample in sample_names]


def _ordered_response_spec(
    response: str, values: Sequence[object], ordered_levels
) -> ResponseSpec:
    observed = _observed_levels(values)
    levels = [str(level) for level in ordered_levels[response]]
    valid = (
        len(levels) >= 2
        and len(levels) == len(set(levels))
        and set(observed) == set(levels)
    )
    if not valid:
        raise ValueError(
            "Ordered response '{}' needs unique levels covering all observations.".format(
                response
            )
        )
    return ResponseSpec(response, "ordered", "ordinal", tuple(levels))


def _categorical_response_spec(
    response: str, values: Sequence[object], references
) -> ResponseSpec:
    observed = _observed_levels(values)
    if len(observed) < 2:
        raise ValueError(
            "Categorical response '{}' needs at least two levels.".format(response)
        )
    reference = str(references.get(response, observed[0]))
    if reference not in observed:
        raise ValueError(
            "Response reference '{}' is not observed for '{}'.".format(
                reference, response
            )
        )
    family = "binomial" if len(observed) == 2 else "multinomial"
    return ResponseSpec(response, "categorical", family, tuple(observed), reference)


def _explicit_response_spec(
    response, values, explicit_family, ordered_levels, references, allow_missing
):
    if explicit_family not in RESPONSE_FAMILIES:
        raise ValueError(
            "Unsupported response family for '{}': {}.".format(
                response, explicit_family
            )
        )
    if explicit_family == "ordinal":
        if response not in ordered_levels:
            raise ValueError(
                "Ordinal response '{}' requires --ordered-responses.".format(response)
            )
        return _ordered_response_spec(response, values, ordered_levels)
    if explicit_family in {"binomial", "multinomial"}:
        spec = _categorical_response_spec(response, values, references)
        if spec.family != explicit_family:
            raise ValueError(
                "Response '{}' has {} observed levels, incompatible with family "
                "'{}'.".format(response, len(spec.levels), explicit_family)
            )
        return spec
    kind, _link = RESPONSE_FAMILIES[explicit_family]
    invalid = any(
        _finite_number(value) is None
        and not isinstance(value, ReplicatedObservation)
        and not (allow_missing and _missing_value(value))
        for value in values
    )
    if invalid:
        raise ValueError(
            "Response '{}' must be numeric for family '{}'.".format(
                response, explicit_family
            )
        )
    return ResponseSpec(response, kind, explicit_family)


def _inferred_response_spec(
    response, values, categorical_set, ordered_levels, references, allow_missing
):
    if response in ordered_levels:
        return _ordered_response_spec(response, values, ordered_levels)
    evaluated_values = [
        value for value in values if not (allow_missing and _missing_value(value))
    ]
    if not evaluated_values:
        raise ValueError("Response '{}' has no observed values.".format(response))
    discrete = response in categorical_set or any(
        isinstance(value, CategoricalObservation) or _finite_number(value) is None
        for value in evaluated_values
    )
    if discrete:
        return _categorical_response_spec(response, values, references)
    return ResponseSpec(response, "continuous", "gaussian")


def _resolve_response_spec(
    response,
    values,
    categorical_set,
    ordered_levels,
    references,
    families,
    allow_missing,
) -> ResponseSpec:
    explicit_family = families.get(response)
    if explicit_family is not None:
        return _explicit_response_spec(
            response,
            values,
            explicit_family,
            ordered_levels,
            references,
            allow_missing,
        )
    return _inferred_response_spec(
        response,
        values,
        categorical_set,
        ordered_levels,
        references,
        allow_missing,
    )


def _validate_response_type_configuration(categorical_set, ordered_levels, families):
    conflicting_kinds = categorical_set & set(ordered_levels)
    if conflicting_kinds:
        raise ValueError(
            "Responses cannot be both unordered and ordered: {}.".format(
                ", ".join(sorted(conflicting_kinds))
            )
        )
    for response, family in families.items():
        if response in ordered_levels and family != "ordinal":
            raise ValueError(
                "Response '{}' is ordered but explicitly uses family '{}'.".format(
                    response, family
                )
            )
        if response in categorical_set and family not in {"binomial", "multinomial"}:
            raise ValueError(
                "Response '{}' is unordered categorical but explicitly uses family '{}'.".format(
                    response, family
                )
            )


def _validate_response_references(resolved, references):
    invalid_references = [
        response
        for response in references
        if resolved[response].family not in {"binomial", "multinomial"}
    ]
    if invalid_references:
        raise ValueError(
            "Response references apply only to binomial or multinomial responses: {}.".format(
                ", ".join(sorted(invalid_references))
            )
        )


def _configured_response_names(categorical_set, ordered_levels, references, families):
    return categorical_set | set(ordered_levels) | set(references) | set(families)


def resolve_response_specs(
    values_by_trait: Mapping[str, Mapping[str, object]],
    responses: Sequence[str],
    sample_names: Sequence[str],
    *,
    categorical: Sequence[str] = (),
    ordered_levels: Mapping[str, Sequence[str]] | None = None,
    references: Mapping[str, str] | None = None,
    families: Mapping[str, str] | None = None,
    allow_missing: bool = False,
) -> dict[str, ResponseSpec]:
    """Resolve Gaussian, binomial, multinomial, and ordinal responses."""
    categorical_set = set(categorical)
    ordered_levels = {} if ordered_levels is None else ordered_levels
    references = {} if references is None else references
    families = {} if families is None else families
    _validate_typed_names(
        responses,
        _configured_response_names(
            categorical_set, ordered_levels, references, families
        ),
        "Response type/reference",
    )
    _validate_response_type_configuration(categorical_set, ordered_levels, families)
    resolved = {
        response: _resolve_response_spec(
            response,
            _values_for_samples(values_by_trait, response, sample_names, "Response"),
            categorical_set,
            ordered_levels,
            references,
            families,
            allow_missing,
        )
        for response in responses
    }
    _validate_response_references(resolved, references)
    return resolved


def _ordered_contrast(levels: Sequence[str]) -> tuple[np.ndarray, list[str]]:
    scores = np.arange(len(levels), dtype=float)
    powers = np.column_stack([scores**degree for degree in range(len(levels))])
    orthogonal, triangular = np.linalg.qr(powers)
    signs = np.sign(np.diag(triangular))
    signs[signs == 0.0] = 1.0
    contrasts = orthogonal[:, 1:] * signs[1:]
    labels = ["linear", "quadratic"]
    names = [
        labels[index] if index < len(labels) else "degree{}".format(index + 1)
        for index in range(len(levels) - 1)
    ]
    return contrasts, names


def _factor_contrast(
    levels: list[str], reference: str, coding: str
) -> tuple[np.ndarray, list[str]]:
    if reference not in levels:
        raise ValueError(
            "Factor reference '{}' is not an observed level (levels: {}).".format(
                reference, ", ".join(levels)
            )
        )
    nonreference = [level for level in levels if level != reference]
    matrix = np.zeros((len(levels), len(nonreference)), dtype=float)
    for column, level in enumerate(nonreference):
        matrix[levels.index(level), column] = 1.0
    if coding == "sum":
        matrix[levels.index(reference), :] = -1.0
    elif coding != "treatment":
        raise ValueError("Unsupported factor coding: {}.".format(coding))
    return matrix, nonreference


def _encode_discrete_predictor(
    name: str,
    values_by_sample: Mapping[str, object],
    sample_names: Sequence[str],
    *,
    kind: str,
    ordered_levels: Mapping[str, Sequence[str]],
    factor_references: Mapping[str, str],
    factor_coding: str,
) -> tuple[
    list[PredictorTerm], dict[str, dict[str, float]], PredictorUncertainty | None
]:
    values = [values_by_sample[sample] for sample in sample_names]
    if kind == "ordered":
        levels = [str(level) for level in ordered_levels[name]]
        observed = set(_observed_levels(values))
        if (
            len(levels) < 2
            or len(levels) != len(set(levels))
            or not observed <= set(levels)
        ):
            raise ValueError(
                "Ordered predictor '{}' needs unique levels covering all observations.".format(
                    name
                )
            )
        contrast, column_labels = _ordered_contrast(levels)
        reference = ""
        coding = "polynomial"
    else:
        levels = _observed_levels(values)
        if len(levels) < 2:
            raise ValueError(
                "Categorical predictor '{}' needs at least two levels.".format(name)
            )
        reference = str(factor_references.get(name, levels[0]))
        contrast, column_labels = _factor_contrast(levels, reference, factor_coding)
        coding = factor_coding
    term_names = ["{}[{}]".format(name, label) for label in column_labels]
    terms = [
        PredictorTerm(
            name=term_name,
            source=name,
            kind=kind,
            level=column_labels[index],
            reference=reference,
            coding=coding,
        )
        for index, term_name in enumerate(term_names)
    ]
    encoded: dict[str, dict[str, float]] = {term_name: {} for term_name in term_names}
    covariance_by_observation = np.zeros(
        (len(sample_names), len(term_names), len(term_names)), dtype=float
    )
    uncertain = False
    for sample_index, sample in enumerate(sample_names):
        probabilities = _probabilities(values_by_sample[sample], levels)
        expected = probabilities @ contrast
        second_moment = contrast.T @ (probabilities[:, None] * contrast)
        covariance = second_moment - np.outer(expected, expected)
        observation = values_by_sample[sample]
        if isinstance(observation, CategoricalObservation):
            if observation.n_observations <= 0:
                raise ValueError(
                    "Categorical biological replicate count must be positive."
                )
            covariance /= float(observation.n_observations)
        covariance_by_observation[sample_index] = (covariance + covariance.T) / 2.0
        uncertain = uncertain or bool(np.max(np.abs(covariance)) > 1e-15)
        for column, term_name in enumerate(term_names):
            encoded[term_name][sample] = float(expected[column])
    uncertainty = (
        PredictorUncertainty(
            source=name,
            term_names=tuple(term_names),
            covariance_by_observation=covariance_by_observation,
        )
        if uncertain
        else None
    )
    return terms, encoded, uncertainty


def _encode_continuous_predictor(
    predictor: str, values: Sequence[object], sample_names: Sequence[str]
) -> tuple[list[PredictorTerm], dict[str, dict[str, float]], None]:
    numeric = [_finite_number(value) for value in values]
    if any(value is None for value in numeric):
        raise ValueError("Continuous predictor '{}' must be finite.".format(predictor))
    encoded: dict[str, dict[str, float]] = {predictor: {}}
    for sample, value in zip(sample_names, numeric, strict=True):
        if value is None:
            raise AssertionError("Finite predictor validation was inconsistent.")
        encoded[predictor][sample] = value
    return [PredictorTerm(predictor, predictor, "continuous")], encoded, None


def _encode_one_predictor(
    predictor,
    values_by_trait,
    sample_names,
    categorical_set,
    ordered_levels,
    factor_references,
    factor_coding,
):
    values = _values_for_samples(values_by_trait, predictor, sample_names, "Predictor")
    values_by_sample = values_by_trait[predictor]
    kind = predictor_kind(
        predictor,
        values,
        categorical=categorical_set,
        ordered_levels=ordered_levels,
    )
    if kind == "continuous":
        return _encode_continuous_predictor(predictor, values, sample_names)
    return _encode_discrete_predictor(
        predictor,
        values_by_sample,
        sample_names,
        kind=kind,
        ordered_levels=ordered_levels,
        factor_references=factor_references,
        factor_coding=factor_coding,
    )


def _validate_predictor_encoding_options(
    predictors,
    categorical_set,
    ordered_levels,
    factor_references,
    factor_coding,
):
    _validate_typed_names(
        predictors,
        categorical_set | set(ordered_levels) | set(factor_references),
        "Predictor type/reference",
    )
    if factor_coding not in {"treatment", "sum"}:
        raise ValueError("Unsupported factor coding: {}.".format(factor_coding))


def encode_predictors(
    values_by_trait: Mapping[str, Mapping[str, object]],
    predictors: Sequence[str],
    sample_names: Sequence[str],
    *,
    categorical: Sequence[str] = (),
    ordered_levels: Mapping[str, Sequence[str]] | None = None,
    factor_references: Mapping[str, str] | None = None,
    factor_coding: str = "treatment",
) -> EncodedPredictors:
    """Expand typed predictors into finite numeric design-matrix columns."""
    predictors = list(predictors)
    sample_names = [str(name) for name in sample_names]
    categorical_set = set(categorical)
    ordered_levels = {} if ordered_levels is None else ordered_levels
    factor_references = {} if factor_references is None else factor_references
    _validate_predictor_encoding_options(
        predictors,
        categorical_set,
        ordered_levels,
        factor_references,
        factor_coding,
    )
    encoded: dict[str, dict[str, float]] = {}
    terms: list[PredictorTerm] = []
    groups: dict[str, tuple[str, ...]] = {}
    uncertainties: list[PredictorUncertainty] = []
    for predictor in predictors:
        predictor_terms, predictor_values, uncertainty = _encode_one_predictor(
            predictor,
            values_by_trait,
            sample_names,
            categorical_set,
            ordered_levels,
            factor_references,
            factor_coding,
        )
        encoded.update(predictor_values)
        term_names = tuple(term.name for term in predictor_terms)
        if set(term_names) & set(encoded) - set(term_names):
            raise ValueError("Encoded predictor term names are not unique.")
        terms.extend(predictor_terms)
        groups[predictor] = term_names
        if uncertainty is not None:
            uncertainties.append(uncertainty)
    if len(terms) != len({term.name for term in terms}):
        raise ValueError("Encoded predictor term names are not unique.")
    return EncodedPredictors(
        values_by_trait=encoded,
        terms=tuple(terms),
        groups=groups,
        uncertainties=tuple(uncertainties),
    )


def parse_name_list(value: object) -> list[str]:
    """Parse an optional comma-separated list without accepting empty names."""
    if value in (None, ""):
        return []
    names = [item.strip() for item in str(value).split(",")]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise ValueError("A predictor-name list contains empty or duplicate names.")
    return names


def parse_key_values(value: object, option_name: str) -> dict[str, str]:
    """Parse comma-separated NAME=VALUE settings."""
    if value in (None, ""):
        return {}
    parsed: dict[str, str] = {}
    for item in str(value).split(","):
        name, separator, setting = item.strip().partition("=")
        if not separator or not name or not setting or name in parsed:
            raise ValueError(
                "'{}' requires unique NAME=VALUE items.".format(option_name)
            )
        parsed[name] = setting
    return parsed


def parse_numeric_key_values(
    value: object,
    option_name: str,
    *,
    lower: float | None = None,
    upper: float | None = None,
) -> dict[str, float]:
    """Parse finite NAME=FLOAT settings with optional strict bounds."""
    parsed = parse_key_values(value, option_name)
    numeric: dict[str, float] = {}
    for name, setting in parsed.items():
        try:
            number = float(setting)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(
                "'{}' value for '{}' must be finite numeric.".format(option_name, name)
            ) from exc
        if not np.isfinite(number):
            raise ValueError(
                "'{}' value for '{}' must be finite numeric.".format(option_name, name)
            )
        if lower is not None and number <= lower:
            raise ValueError(
                "'{}' value for '{}' must be greater than {}.".format(
                    option_name, name, lower
                )
            )
        if upper is not None and number >= upper:
            raise ValueError(
                "'{}' value for '{}' must be less than {}.".format(
                    option_name, name, upper
                )
            )
        numeric[name] = number
    return numeric


def parse_ordered_levels(
    value: object, option_name: str = "--ordered-predictors"
) -> dict[str, tuple[str, ...]]:
    """Parse NAME=LEVEL1|LEVEL2|... ordered-predictor declarations."""
    settings = parse_key_values(value, option_name)
    parsed = {}
    for name, raw_levels in settings.items():
        levels = tuple(level.strip() for level in raw_levels.split("|"))
        if (
            len(levels) < 2
            or any(not level for level in levels)
            or len(levels) != len(set(levels))
        ):
            raise ValueError(
                "Ordered predictor '{}' needs at least two unique levels.".format(name)
            )
        parsed[name] = levels
    return parsed
