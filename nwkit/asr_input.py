"""Trait-type inference and mode-specific ASR input contracts."""

import math
from copy import copy
from dataclasses import dataclass

from nwkit.asr_models import default_model, model_definition, model_names
from nwkit.util import is_missing_table_value, read_tip_table

_MODE_OPTIONS = {
    "discrete": {
        "output": ("probabilities", "map"),
        "tree_annotation": ("map", "state", "probability", "all"),
    },
    "continuous": {
        "output": ("summary",),
        "tree_annotation": ("summary", "mean", "all"),
    },
}
_DISCRETE_ONLY = (
    "states",
    "rate",
    "rate_bounds",
    "rate_matrix",
    "transition_graph",
    "ambiguous_separator",
    "stochastic_map_out",
    "n_sim",
    "threads",
    "seed",
)
_CONTINUOUS_ONLY = (
    "sigma2",
    "standard_error_column",
    "ci_level",
    "alpha",
    "alpha_bounds",
    "theta",
)


def _mode_switch_hint(trait_type):
    if trait_type == "continuous":
        return "Use --trait-type discrete for numeric categories."
    return "Use --trait-type continuous for a continuous numeric trait."


@dataclass(frozen=True)
class AsrSettings:
    trait_type: str
    model: str
    root_prior: str
    output: str
    tree_annotation: str
    ci_level: float | None

    @classmethod
    def from_args(cls, args, trait_type="discrete"):
        options = _resolve_mode_options(args, trait_type)
        _validate_mode_arguments(args, trait_type)
        _validate_model_arguments(args, trait_type, options["model"])
        ci_level = _confidence_level(args) if trait_type == "continuous" else None
        return cls(trait_type=trait_type, ci_level=ci_level, **options)


def _resolve_mode_options(args, trait_type):
    if trait_type not in _MODE_OPTIONS:
        raise ValueError(f"Unsupported resolved trait type: {trait_type}.")
    model = getattr(args, "model", None) or default_model(trait_type)
    definition = model_definition(model)
    if definition.trait_type != trait_type:
        raise ValueError(
            f"--model={model} is not supported for {trait_type} traits; choose "
            f"{', '.join(model_names(trait_type))}. {_mode_switch_hint(trait_type)}"
        )
    root_prior = getattr(args, "root_prior", None) or definition.default_root_prior
    if root_prior not in definition.root_priors:
        raise ValueError(
            f"--root-prior={root_prior} is not supported for {trait_type} traits "
            f"with --model={model}; "
            f"choose {', '.join(definition.root_priors)}."
        )
    resolved = {"model": model, "root_prior": root_prior}
    for name, choices in _MODE_OPTIONS[trait_type].items():
        value = getattr(args, name, None)
        value = choices[0] if value is None else value
        if value not in choices:
            option = "--" + name.replace("_", "-")
            raise ValueError(
                f"{option}={value} is not supported for {trait_type} traits; "
                f"choose {', '.join(choices)}. {_mode_switch_hint(trait_type)}"
            )
        resolved[name] = value
    return resolved


def _validate_model_arguments(args, trait_type, model):
    if trait_type == "discrete":
        rate_matrix = getattr(args, "rate_matrix", None)
        transition_graph = getattr(args, "transition_graph", None)
        if model == "CUSTOM" and rate_matrix in (None, ""):
            raise ValueError("--model CUSTOM requires --rate-matrix.")
        if model != "CUSTOM" and rate_matrix not in (None, ""):
            raise ValueError("--rate-matrix requires --model CUSTOM.")
        if model == "CUSTOM" and transition_graph not in (None, ""):
            raise ValueError(
                "--transition-graph cannot be combined with --model CUSTOM."
            )
        if model == "CUSTOM" and getattr(args, "rate", None) is not None:
            raise ValueError("--rate cannot be combined with --model CUSTOM.")
        if model == "CUSTOM" and getattr(args, "rate_bounds", None) is not None:
            raise ValueError("--rate-bounds cannot be combined with --model CUSTOM.")
        if transition_graph == "ordered" and getattr(args, "states", None) in (
            None,
            "",
        ):
            raise ValueError("--transition-graph ordered requires explicit --states.")
        return
    if model == "BM":
        supplied = [
            "--" + name.replace("_", "-")
            for name in ("alpha", "alpha_bounds", "theta")
            if getattr(args, name, None) is not None
        ]
        if supplied:
            raise ValueError(f"Options requiring --model OU: {', '.join(supplied)}.")
    elif (
        getattr(args, "alpha", None) is not None
        and getattr(args, "alpha_bounds", None) is not None
    ):
        raise ValueError("--alpha-bounds cannot be combined with fixed --alpha.")


def _validate_mode_arguments(args, trait_type):
    incompatible = _DISCRETE_ONLY if trait_type == "continuous" else _CONTINUOUS_ONLY
    supplied = [
        "--" + name.replace("_", "-")
        for name in incompatible
        if getattr(args, name, None) is not None
    ]
    if supplied:
        raise ValueError(
            f"Options not supported for {trait_type} traits: {', '.join(supplied)}. "
            f"{_mode_switch_hint(trait_type)}"
        )


def _confidence_level(args):
    value = getattr(args, "ci_level", None)
    value = 0.95 if value is None else value
    if isinstance(value, bool):
        raise ValueError(
            "--ci-level must be a finite number strictly between zero and one."
        )
    try:
        level = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            "--ci-level must be a finite number strictly between zero and one."
        ) from exc
    if not math.isfinite(level) or not 0.0 < level < 1.0:
        raise ValueError(
            "--ci-level must be a finite number strictly between zero and one."
        )
    return level


def effective_asr_args(args, settings):
    effective = copy(args)
    for name in ("model", "root_prior", "output", "tree_annotation", "ci_level"):
        setattr(effective, name, getattr(settings, name))
    if settings.trait_type == "discrete":
        for name, default in (
            ("ambiguous_separator", "|"),
            ("n_sim", 100),
            ("threads", 1),
        ):
            if getattr(effective, name, None) is None:
                setattr(effective, name, default)
    return effective


def read_asr_table(
    trait_path,
    state_column,
    tree_leaf_names,
    *,
    missing_values=None,
    unmatched="warn",
    standard_error_column=None,
):
    if trait_path in ("", None):
        raise ValueError("'--trait' is required.")
    if state_column in ("", None):
        raise ValueError("'--state-column' is required.")
    dataframe, _, _ = read_tip_table(
        trait_path,
        option_name="--trait",
        tree_leaf_names=tree_leaf_names,
        required_columns=(state_column,),
        unmatched=unmatched,
        missing_values=missing_values,
        preserve_columns=tuple(
            column
            for column in (state_column, standard_error_column)
            if column is not None
        ),
    )
    return dataframe[dataframe["leaf_name"].isin(set(tree_leaf_names))].copy()


def _numeric_value(value):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return False
    if not math.isfinite(number):
        raise ValueError(
            "Trait values that parse as numbers must be finite for --trait-type auto."
        )
    return True


def resolve_trait_type(requested, dataframe, state_column):
    if requested not in {"auto", "discrete", "continuous"}:
        raise ValueError(f"Unsupported --trait-type: {requested}.")
    if requested != "auto":
        return requested
    observed = [
        value
        for value in dataframe[state_column]
        if not is_missing_table_value(value, set())
    ]
    if not observed:
        raise ValueError(
            "Cannot infer --trait-type from all-missing values. Specify --trait-type "
            "discrete with --states for a discrete prior-only analysis; continuous ASR needs observations."
        )
    # Inspect every token before deciding: ``all(generator)`` would stop at the
    # first text category and could hide a later numeric-looking infinity.
    numeric = [_numeric_value(value) for value in observed]
    return "continuous" if all(numeric) else "discrete"


def _continuous_number(value, label, *, nonnegative=False):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be numeric and finite.") from exc
    if not math.isfinite(number) or (nonnegative and number < 0.0):
        qualifier = "non-negative and finite" if nonnegative else "finite"
        raise ValueError(f"{label} must be {qualifier}.")
    return number


def continuous_tip_values(
    dataframe, state_column, leaf_names, standard_error_column=None
):
    if (
        standard_error_column is not None
        and standard_error_column not in dataframe.columns
    ):
        raise ValueError(
            f"Missing required standard-error column: {standard_error_column}."
        )
    if standard_error_column == state_column:
        raise ValueError("--standard-error-column must differ from --state-column.")
    values = dict.fromkeys(leaf_names)
    errors = {}
    for _, row in dataframe.iterrows():
        name = str(row["leaf_name"])
        raw_value = row[state_column]
        if is_missing_table_value(raw_value, set()):
            continue
        values[name] = _continuous_number(raw_value, f"Trait value for '{name}'")
        if standard_error_column is not None:
            raw_error = row[standard_error_column]
            if is_missing_table_value(raw_error, set()):
                raise ValueError(
                    f"A standard error is required for observed tip '{name}'."
                )
            errors[name] = _continuous_number(
                raw_error, f"Standard error for '{name}'", nonnegative=True
            )
    return values, errors if standard_error_column is not None else None
