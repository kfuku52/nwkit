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
    "rate_design",
    "transition_graph",
    "hidden_categories",
    "mixture_model",
    "rate_mixture",
    "rate_categories",
    "thresholds",
    "liability_samples",
    "liability_burnin",
    "liability_thin",
    "liability_chains",
    "liability_out",
    "ambiguous_separator",
    "stochastic_map_out",
    "n_sim",
    "threads",
)
_CONTINUOUS_ONLY = (
    "sigma2",
    "evolution_parameter",
    "evolution_parameter_bounds",
    "standard_error_column",
    "ci_level",
    "alpha",
    "alpha_by_trait",
    "alpha_bounds",
    "theta",
    "eb_rate",
    "eb_rate_bounds",
    "drift",
    "covariance_out",
    "root_mean",
    "root_variance",
    "profile_ci_level",
    "posterior_samples_out",
    "posterior_samples",
    "posterior_predictive_out",
    "posterior_predictive_simulations",
    "bootstrap_out",
    "bootstrap_simulations",
)
_REGIME_MODELS = frozenset(
    {"MK-REGIME", "BMS", "BMS-DRIFT", "OUM", "OUMA", "OUMV", "OUMVA"}
)
_CONTINUOUS_REGIME_MODELS = frozenset(_REGIME_MODELS - {"MK-REGIME"})
_TRANSFORMED_MODELS = frozenset({"LAMBDA", "KAPPA", "DELTA", "EB", "ACDC"})
_OU_MODELS = frozenset({"OU", "MV-OU", "MV-OU-DIAG", "OUM", "OUMA", "OUMV", "OUMVA"})
_MULTIVARIATE_MODELS = frozenset({"MV-BM", "MV-OU", "MV-OU-DIAG"})
_PAGEL_MODELS = frozenset({"PAGEL-INDEPENDENT", "PAGEL-DEPENDENT"})


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
        _validate_model_arguments(
            args, trait_type, options["model"], options["root_prior"]
        )
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


def _validate_regime_arguments(args, model):
    regime_map = getattr(args, "regime_map", None)
    regime_model = getattr(args, "regime_model", None)
    regime_parameters = getattr(args, "regime_parameters", None)
    if model in _REGIME_MODELS and regime_map in (None, ""):
        raise ValueError(f"--model {model} requires --regime-map.")
    if model not in _REGIME_MODELS and regime_map not in (None, ""):
        raise ValueError("--regime-map requires a regime model.")
    if model != "MK-REGIME" and regime_model not in (None, ""):
        raise ValueError("--regime-model requires --model MK-REGIME.")
    if model not in _CONTINUOUS_REGIME_MODELS and regime_parameters not in (None, ""):
        raise ValueError("--regime-parameters requires a continuous regime model.")


def _validate_custom_and_graph_options(args, model):
    rate_matrix = getattr(args, "rate_matrix", None)
    rate_design = getattr(args, "rate_design", None)
    transition_graph = getattr(args, "transition_graph", None)
    if model == "CUSTOM" and rate_matrix in (None, ""):
        raise ValueError("--model CUSTOM requires --rate-matrix.")
    if model != "CUSTOM" and rate_matrix not in (None, ""):
        raise ValueError("--rate-matrix requires --model CUSTOM.")
    if model == "CUSTOM" and transition_graph not in (None, ""):
        raise ValueError("--transition-graph cannot be combined with --model CUSTOM.")
    if model == "CUSTOM" and getattr(args, "rate", None) is not None:
        raise ValueError("--rate cannot be combined with --model CUSTOM.")
    if model == "CUSTOM" and getattr(args, "rate_bounds", None) is not None:
        raise ValueError("--rate-bounds cannot be combined with --model CUSTOM.")
    if model == "MK-REGIME" and rate_matrix not in (None, ""):
        raise ValueError("--rate-matrix is not supported with --model MK-REGIME.")
    if model == "MK-DESIGN" and rate_design in (None, ""):
        raise ValueError("--model MK-DESIGN requires --rate-design.")
    if model != "MK-DESIGN" and rate_design not in (None, ""):
        raise ValueError("--rate-design requires --model MK-DESIGN.")
    if model == "MK-DESIGN" and transition_graph not in (None, ""):
        raise ValueError(
            "--transition-graph cannot be combined with --model MK-DESIGN; "
            "the rate design defines its allowed edges."
        )
    if model in _PAGEL_MODELS and transition_graph not in (None, ""):
        raise ValueError(
            "--transition-graph is not supported by Pagel models; their "
            "four-state transition graph is fixed."
        )
    if transition_graph == "ordered" and getattr(args, "states", None) in (None, ""):
        raise ValueError("--transition-graph ordered requires explicit --states.")
    hidden_categories = getattr(args, "hidden_categories", None)
    if model not in {"HRM", "COVARION"} and hidden_categories is not None:
        raise ValueError("--hidden-categories requires --model HRM or COVARION.")


def _validate_mixture_options(args, model):
    option_names = ("mixture_model", "rate_mixture", "rate_categories")
    supplied = [name for name in option_names if getattr(args, name, None) is not None]
    if model != "MK-MIXTURE" and supplied:
        raise ValueError(
            "--mixture-model, --rate-mixture, and --rate-categories require "
            "--model MK-MIXTURE."
        )
    if model == "MK-MIXTURE" and getattr(args, "stochastic_map_out", None) not in (
        None,
        "",
    ):
        raise ValueError(
            "--stochastic-map-out is not defined for a jointly fitted "
            "across-character rate mixture."
        )
    if model == "MK-MIXTURE" and getattr(args, "compare_models", None) not in (
        None,
        "",
    ):
        raise ValueError(
            "--compare-models is not defined for the multi-character "
            "--model MK-MIXTURE input contract."
        )


def _validate_threshold_options(args, model):
    option_names = (
        "thresholds",
        "liability_samples",
        "liability_burnin",
        "liability_thin",
        "liability_chains",
        "liability_out",
    )
    if model != "THRESHOLD" and any(
        getattr(args, name, None) is not None for name in option_names
    ):
        raise ValueError(
            "--thresholds and --liability-* options require --model THRESHOLD."
        )
    if model != "THRESHOLD":
        return
    if getattr(args, "states", None) in (None, ""):
        raise ValueError("--model THRESHOLD requires ordered categories in --states.")
    forbidden = [
        "--" + name.replace("_", "-")
        for name in (
            "rate",
            "rate_bounds",
            "rate_matrix",
            "transition_graph",
            "hidden_categories",
            "mixture_model",
            "rate_mixture",
            "rate_categories",
        )
        if getattr(args, name, None) is not None
    ]
    if forbidden:
        raise ValueError(
            "Options not defined for --model THRESHOLD: " + ", ".join(forbidden) + "."
        )
    if getattr(args, "stochastic_map_out", None) not in (None, ""):
        raise ValueError("--stochastic-map-out is not defined for --model THRESHOLD.")
    if getattr(args, "compare_models", None) not in (None, ""):
        raise ValueError("--compare-models is not defined for --model THRESHOLD.")


def _validate_discrete_model_arguments(args, model):
    _validate_custom_and_graph_options(args, model)
    _validate_mixture_options(args, model)
    _validate_threshold_options(args, model)
    if getattr(args, "stochastic_map_out", None) in (None, ""):
        supplied = [
            "--" + name.replace("_", "-")
            for name in ("n_sim", "threads")
            if getattr(args, name, None) is not None
        ]
        if supplied:
            raise ValueError(
                "These options require --stochastic-map-out: "
                + ", ".join(supplied)
                + "."
            )


def _validate_transformed_options(args, model):
    if model not in {"EB", "ACDC"} and any(
        getattr(args, name, None) is not None for name in ("eb_rate", "eb_rate_bounds")
    ):
        raise ValueError("--eb-rate and --eb-rate-bounds require --model EB or ACDC.")
    if model in {"EB", "ACDC"} and all(
        getattr(args, name, None) is not None for name in ("eb_rate", "eb_rate_bounds")
    ):
        raise ValueError("--eb-rate-bounds cannot be combined with fixed --eb-rate.")
    parameter = getattr(args, "evolution_parameter", None)
    bounds = getattr(args, "evolution_parameter_bounds", None)
    if model not in _TRANSFORMED_MODELS and any(
        value is not None for value in (parameter, bounds)
    ):
        raise ValueError(
            "--evolution-parameter and --evolution-parameter-bounds require "
            "--model LAMBDA, KAPPA, DELTA, EB, or ACDC."
        )
    if parameter is not None and bounds is not None:
        raise ValueError(
            "--evolution-parameter-bounds cannot be combined with fixed "
            "--evolution-parameter."
        )
    if (
        model in {"EB", "ACDC"}
        and parameter is not None
        and getattr(args, "eb_rate", None) is not None
    ):
        raise ValueError("--evolution-parameter cannot be combined with --eb-rate.")
    if (
        model in {"EB", "ACDC"}
        and bounds is not None
        and getattr(args, "eb_rate_bounds", None) is not None
    ):
        raise ValueError(
            "--evolution-parameter-bounds cannot be combined with --eb-rate-bounds."
        )
    _validate_profile_options(args, model, parameter)


def _validate_profile_options(args, model, evolution_parameter):
    level = getattr(args, "profile_ci_level", None)
    if level is None:
        return
    if model not in _TRANSFORMED_MODELS:
        raise ValueError(
            "--profile-ci-level requires --model LAMBDA, KAPPA, DELTA, EB, or ACDC."
        )
    if not 0.0 < float(level) < 1.0:
        raise ValueError("--profile-ci-level must be strictly between zero and one.")
    if evolution_parameter is not None or getattr(args, "eb_rate", None) is not None:
        raise ValueError(
            "--profile-ci-level requires the evolution parameter to be estimated."
        )


def _validate_ou_root_options(args, model, root_prior):
    root_mean = getattr(args, "root_mean", None)
    root_variance = getattr(args, "root_variance", None)
    if model != "OU" and any(value is not None for value in (root_mean, root_variance)):
        raise ValueError(
            "--root-mean and --root-variance currently require --model OU."
        )
    if model != "OU":
        return
    if root_prior == "stationary" and any(
        value is not None for value in (root_mean, root_variance)
    ):
        raise ValueError(
            "Stationary OU determines its root distribution; remove "
            "--root-mean and --root-variance."
        )
    if root_prior in {"fixed", "gaussian"} and root_mean is None:
        raise ValueError(f"--root-prior {root_prior} requires --root-mean.")
    if root_prior == "fixed" and root_variance is not None:
        raise ValueError("--root-prior fixed cannot take --root-variance.")
    if root_prior == "gaussian" and root_variance is None:
        raise ValueError("--root-prior gaussian requires --root-variance.")


def _validate_multivariate_options(args, model):
    diagnostic_options = (
        "posterior_samples_out",
        "posterior_samples",
        "posterior_predictive_out",
        "posterior_predictive_simulations",
        "bootstrap_out",
        "bootstrap_simulations",
    )
    if model in _MULTIVARIATE_MODELS:
        if getattr(args, "sigma2", None) is not None:
            raise ValueError("--sigma2 is not supported for multivariate ASR models.")
        if (
            model in {"MV-OU", "MV-OU-DIAG"}
            and getattr(args, "theta", None) is not None
        ):
            raise ValueError(
                "--theta is scalar; multivariate OU models estimate one optimum "
                "per trait."
            )
        if any(getattr(args, name, None) is not None for name in diagnostic_options):
            raise ValueError(
                "Simulation diagnostics currently support scalar continuous models, "
                "not MV-BM, MV-OU, or MV-OU-DIAG."
            )
        if getattr(args, "compare_models", None) not in (None, ""):
            raise ValueError(
                "--compare-models currently supports scalar continuous traits, "
                "not MV-BM, MV-OU, or MV-OU-DIAG."
            )
    elif getattr(args, "covariance_out", None) not in (None, ""):
        raise ValueError(
            "--covariance-out requires --model MV-BM, MV-OU, or MV-OU-DIAG."
        )


def _validate_diagnostic_counts(args):
    for count_name, output_name in (
        ("posterior_samples", "posterior_samples_out"),
        ("posterior_predictive_simulations", "posterior_predictive_out"),
        ("bootstrap_simulations", "bootstrap_out"),
    ):
        if getattr(args, count_name, None) is not None and getattr(
            args, output_name, None
        ) in (None, ""):
            raise ValueError(
                "--{} requires --{}.".format(
                    count_name.replace("_", "-"), output_name.replace("_", "-")
                )
            )


def _validate_continuous_model_arguments(args, model, root_prior):
    if model not in _OU_MODELS:
        supplied = [
            "--" + name.replace("_", "-")
            for name in ("alpha", "alpha_by_trait", "alpha_bounds", "theta")
            if getattr(args, name, None) is not None
        ]
        if supplied:
            raise ValueError(
                "Options requiring an OU model: " + ", ".join(supplied) + "."
            )
    else:
        alpha = getattr(args, "alpha", None)
        alpha_by_trait = getattr(args, "alpha_by_trait", None)
        alpha_bounds = getattr(args, "alpha_bounds", None)
        if alpha_by_trait is not None and model != "MV-OU-DIAG":
            raise ValueError("--alpha-by-trait requires --model MV-OU-DIAG.")
        if alpha is not None and alpha_by_trait is not None:
            raise ValueError("--alpha cannot be combined with --alpha-by-trait.")
        if alpha_bounds is not None and (
            alpha is not None or alpha_by_trait is not None
        ):
            raise ValueError(
                "--alpha-bounds cannot be combined with fixed --alpha or "
                "--alpha-by-trait."
            )
    _validate_transformed_options(args, model)
    if (
        model not in {"BM-DRIFT", "BMS-DRIFT"}
        and getattr(args, "drift", None) is not None
    ):
        raise ValueError("--drift requires --model BM-DRIFT or BMS-DRIFT.")
    _validate_ou_root_options(args, model, root_prior)
    _validate_multivariate_options(args, model)
    _validate_diagnostic_counts(args)


def _validate_model_arguments(args, trait_type, model, root_prior=None):
    _validate_regime_arguments(args, model)
    if trait_type == "discrete":
        _validate_discrete_model_arguments(args, model)
        return
    _validate_continuous_model_arguments(args, model, root_prior)


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
            ("hidden_categories", 2),
        ):
            if getattr(effective, name, None) is None:
                setattr(effective, name, default)
        if settings.model == "THRESHOLD":
            for name, default in (
                ("liability_samples", 1000),
                ("liability_burnin", 500),
                ("liability_thin", 1),
                ("liability_chains", 4),
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
    state_columns = (
        tuple(state_column)
        if isinstance(state_column, (list, tuple))
        else (state_column,)
    )
    standard_error_columns = (
        tuple(standard_error_column)
        if isinstance(standard_error_column, (list, tuple))
        else ()
        if standard_error_column is None
        else (standard_error_column,)
    )
    dataframe, _, _ = read_tip_table(
        trait_path,
        option_name="--trait",
        tree_leaf_names=tree_leaf_names,
        required_columns=state_columns,
        unmatched=unmatched,
        missing_values=missing_values,
        preserve_columns=tuple(
            column
            for column in (*state_columns, *standard_error_columns)
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
    state_columns = (
        tuple(state_column)
        if isinstance(state_column, (list, tuple))
        else (state_column,)
    )
    observed = [
        value
        for column in state_columns
        for value in dataframe[column]
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


def asr_trait_columns(state_column, model):
    if model == "MK-MIXTURE":
        columns = tuple(item.strip() for item in str(state_column).split(","))
        if len(columns) < 2 or any(item == "" for item in columns):
            raise ValueError(
                "--model MK-MIXTURE requires at least two comma-separated "
                "--state-column values."
            )
        if len(columns) != len(set(columns)):
            raise ValueError(
                "MK-MIXTURE --state-column contains duplicated character columns."
            )
        return columns
    if model in _PAGEL_MODELS:
        columns = tuple(item.strip() for item in str(state_column).split(","))
        if len(columns) != 2 or any(item == "" for item in columns):
            raise ValueError(
                "Pagel models require exactly two comma-separated "
                "--state-column values."
            )
        if columns[0] == columns[1]:
            raise ValueError(
                "Pagel --state-column values must name two distinct traits."
            )
        return columns
    if model not in _MULTIVARIATE_MODELS:
        return (state_column,)
    columns = tuple(item.strip() for item in str(state_column).split(","))
    if len(columns) < 2 or any(item == "" for item in columns):
        raise ValueError(
            "A multivariate ASR model requires at least two comma-separated "
            "--state-column values."
        )
    if len(columns) != len(set(columns)):
        raise ValueError(
            "Multivariate --state-column contains duplicated trait columns."
        )
    return columns


def asr_standard_error_columns(value, model, trait_columns):
    if value in (None, ""):
        return None
    if model not in _MULTIVARIATE_MODELS:
        return value
    columns = tuple(item.strip() for item in str(value).split(","))
    if len(columns) != len(trait_columns) or any(item == "" for item in columns):
        raise ValueError(
            "Multivariate --standard-error-column must list one column per trait."
        )
    if len(columns) != len(set(columns)):
        raise ValueError(
            "Multivariate --standard-error-column contains duplicated columns."
        )
    if set(columns) & set(trait_columns):
        raise ValueError(
            "Multivariate standard-error columns must differ from trait columns."
        )
    return columns


def continuous_tip_vectors(dataframe, state_columns, leaf_names):
    vectors = dict.fromkeys(leaf_names)
    for _, row in dataframe.iterrows():
        name = str(row["leaf_name"])
        vector: list[float | None] = []
        for column in state_columns:
            raw = row[column]
            vector.append(
                None
                if is_missing_table_value(raw, set())
                else _continuous_number(raw, f"Trait '{column}' for '{name}'")
            )
        vectors[name] = (
            None if all(value is None for value in vector) else tuple(vector)
        )
    return vectors


def continuous_tip_vector_errors(dataframe, error_columns, state_columns, leaf_names):
    if error_columns is None:
        return None
    missing = [column for column in error_columns if column not in dataframe.columns]
    if missing:
        raise ValueError(
            "Missing required standard-error column(s): " + ", ".join(missing)
        )
    errors = dict.fromkeys(leaf_names)
    for _, row in dataframe.iterrows():
        name = str(row["leaf_name"])
        vector: list[float | None] = []
        for trait_column, error_column in zip(
            state_columns, error_columns, strict=True
        ):
            trait_missing = is_missing_table_value(row[trait_column], set())
            error_missing = is_missing_table_value(row[error_column], set())
            if trait_missing:
                vector.append(None)
            elif error_missing:
                raise ValueError(
                    f"Missing standard error '{error_column}' for observed trait "
                    f"'{trait_column}' at '{name}'."
                )
            else:
                vector.append(
                    _continuous_number(
                        row[error_column],
                        f"Standard error '{error_column}' for '{name}'",
                        nonnegative=True,
                    )
                )
        errors[name] = tuple(vector)
    return errors


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
