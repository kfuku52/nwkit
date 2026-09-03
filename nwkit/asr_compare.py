"""Batch comparison of ancestral-state reconstruction models."""

import math
import os
import re
import sys
import time
from copy import copy
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from nwkit.asr_compare_figure import (
    _criterion_label,
    _font_family_for_text,
    draw_comparison_figure,
)
from nwkit.asr_comparison import (
    grouped_model_comparison_table,
    has_nonregular_variance_boundary,
    summarize_fit,
)
from nwkit.asr_input import (
    AsrSettings,
    continuous_tip_values,
    continuous_tip_vector_errors,
    continuous_tip_vectors,
    read_asr_table,
    resolve_trait_type,
)
from nwkit.asr_models import model_definition, model_names
from nwkit.discrete_asr_models import model_equivalence_family
from nwkit.output_transaction import output_transaction
from nwkit.util import (
    read_tree,
    validate_distinct_output_paths,
    validate_outputs_do_not_replace_inputs,
)

REGIME_MODELS = frozenset(
    {"MK-REGIME", "BMS", "BMS-DRIFT", "OUM", "OUMA", "OUMV", "OUMVA"}
)
MULTIVARIATE_MODELS = frozenset({"MV-BM", "MV-OU"})
ROOT_VARIANT_MODEL = "OU"
ROOT_VARIANTS = frozenset({"stationary", "fixed", "gaussian"})
MODEL_TOKEN = re.compile(r"^([A-Z][A-Z0-9-]*)(?:\[([a-z]+)\])?$", re.IGNORECASE)

ASR_COMPARE_COLUMNS = (
    "model_id",
    "model",
    "trait_type",
    "trait_type_requested",
    "trait_columns",
    "num_traits",
    "comparison_group",
    "likelihood_kind",
    "root_prior",
    "root_prior_source",
    "equivalent_to",
    "status",
    "rankable",
    "message",
    "fit_status",
    "optimizer_success",
    "optimizer_message",
    "sample_size",
    "num_parameters",
    "estimated_parameters",
    "fixed_parameters",
    "log_likelihood",
    "aic",
    "aicc",
    "bic",
    "delta_aic",
    "delta_aicc",
    "delta_bic",
    "aic_weight",
    "aicc_weight",
    "bic_weight",
    "num_comparable_models",
    "criterion",
    "criterion_value",
    "criterion_rank",
    "is_best",
    "shared_preparation_seconds",
    "elapsed_seconds",
)

INTEGER_COLUMNS = (
    "num_traits",
    "sample_size",
    "num_parameters",
    "num_comparable_models",
    "criterion_rank",
)

TRANSFORMED_MODELS = frozenset({"LAMBDA", "KAPPA", "DELTA", "EB", "ACDC"})
OU_MODELS = frozenset({"OU", "MV-OU", "OUM", "OUMA", "OUMV", "OUMVA"})
DISCRETE_RATE_MODELS = frozenset(
    {
        "ER",
        "SYM",
        "ARD",
        "F81",
        "GTR",
        "MK-REGIME",
        "HRM",
        "COVARION",
        "MK-MIXTURE",
    }
)
CONTINUOUS_SCALAR_MODELS = frozenset(model_names("continuous")) - MULTIVARIATE_MODELS

MODEL_OPTION_CONSUMERS = {
    "rate": DISCRETE_RATE_MODELS,
    "rate_bounds": DISCRETE_RATE_MODELS,
    "transition_graph": DISCRETE_RATE_MODELS,
    "regime_map": REGIME_MODELS,
    "regime_model": frozenset({"MK-REGIME"}),
    "regime_parameters": REGIME_MODELS - {"MK-REGIME"},
    "hidden_categories": frozenset({"HRM", "COVARION"}),
    "mixture_model": frozenset({"MK-MIXTURE"}),
    "rate_mixture": frozenset({"MK-MIXTURE"}),
    "rate_categories": frozenset({"MK-MIXTURE"}),
    "rate_matrix": frozenset({"CUSTOM"}),
    "sigma2": CONTINUOUS_SCALAR_MODELS,
    "evolution_parameter": TRANSFORMED_MODELS,
    "evolution_parameter_bounds": TRANSFORMED_MODELS,
    "alpha": OU_MODELS,
    "alpha_bounds": OU_MODELS,
    "theta": OU_MODELS - {"MV-OU"},
    "eb_rate": frozenset({"EB", "ACDC"}),
    "eb_rate_bounds": frozenset({"EB", "ACDC"}),
    "drift": frozenset({"BM-DRIFT", "BMS-DRIFT"}),
}

REGIME_PARAMETER_COLUMNS = {
    "BMS": ("sigma2",),
    "BMS-DRIFT": ("sigma2", "drift"),
    "OUM": ("theta",),
    "OUMA": ("theta", "alpha"),
    "OUMV": ("theta", "sigma2"),
    "OUMVA": ("theta", "alpha", "sigma2"),
}


@dataclass(frozen=True)
class ComparisonCandidate:
    model: str
    root_prior: str
    root_prior_source: str

    @property
    def model_id(self):
        if self.model == ROOT_VARIANT_MODEL:
            return f"{self.model}[{self.root_prior}]"
        return self.model


@dataclass
class ComparisonContext:
    tree: Any
    trait_df: pd.DataFrame
    trait_type: str
    trait_columns: tuple[str, ...]
    error_columns: tuple[str, ...] | None
    args: Any
    cache: dict[str, Any] = field(default_factory=dict)

    @property
    def leaf_names(self):
        return list(self.tree.leaf_names())


@dataclass(frozen=True)
class _CachedFailure:
    error: Exception


def _comma_tokens(value, option_name, *, allow_empty=False):
    if value in (None, ""):
        if allow_empty:
            return []
        raise ValueError(f"'{option_name}' must name at least one model.")
    tokens = [item.strip() for item in str(value).split(",")]
    if any(token == "" for token in tokens):
        raise ValueError(f"'{option_name}' contains an empty model name.")
    return tokens


def _parse_model_token(token, option_name):
    normalized = str(token).strip()
    match = MODEL_TOKEN.fullmatch(normalized)
    if match is None:
        raise ValueError(
            f"Invalid model in '{option_name}': {token}. "
            "OU root variants use OU[stationary], OU[fixed], or OU[gaussian]."
        )
    raw_model, raw_root_prior = match.groups()
    model = raw_model.upper()
    root_prior = None if raw_root_prior is None else raw_root_prior.lower()
    if model not in model_names():
        raise ValueError(
            f"Unknown model in '{option_name}': {model}. Supported models are: "
            + ", ".join(model_names())
            + "."
        )
    if root_prior is not None and (
        model != ROOT_VARIANT_MODEL or root_prior not in ROOT_VARIANTS
    ):
        raise ValueError(
            f"Unsupported model variant in '{option_name}': {token}. "
            "Only OU accepts [stationary], [fixed], or [gaussian]."
        )
    return model, root_prior


def _candidate(model, variant, args):
    definition = model_definition(model)
    requested_prior = getattr(args, "root_prior", None)
    if variant is not None and requested_prior not in (None, variant):
        raise ValueError(
            f"Model {model}[{variant}] conflicts with --root-prior {requested_prior}."
        )
    if variant is not None:
        prior, source = variant, "model-token"
    elif requested_prior is not None:
        prior, source = requested_prior, "--root-prior"
    else:
        prior, source = definition.default_root_prior, "model-default"
    return ComparisonCandidate(model, prior, source)


def _excluded(candidate, exclusions):
    return any(
        model == "all"
        or (model == candidate.model and variant is None)
        or (model == candidate.model and variant == candidate.root_prior)
        for model, variant in exclusions
    )


def _parse_exclusions(args):
    excluded_tokens = _comma_tokens(
        getattr(args, "exclude_models", None),
        "--exclude-models",
        allow_empty=True,
    )
    exclusions: list[tuple[str, str | None]] = []
    for token in excluded_tokens:
        if token.lower() == "all":
            exclusions.append(("all", None))
        else:
            exclusions.append(_parse_model_token(token, "--exclude-models"))
    return exclusions


def _resolve_requested_models(trait_type, requested, automatic):
    if automatic:
        return [(model, None) for model in model_names(trait_type)]
    parsed = [_parse_model_token(token, "--models") for token in requested]
    wrong_type = [
        model
        for model, _variant in parsed
        if model_definition(model).trait_type != trait_type
    ]
    if wrong_type:
        raise ValueError(
            "Models incompatible with the resolved {} trait: {}.".format(
                trait_type, ", ".join(wrong_type)
            )
        )
    return parsed


def _reject_duplicate_candidates(candidates):
    identifiers = [candidate.model_id for candidate in candidates]
    duplicated = sorted(
        identifier
        for identifier in set(identifiers)
        if identifiers.count(identifier) > 1
    )
    if duplicated:
        raise ValueError(
            "'--models' contains duplicated candidates: " + ", ".join(duplicated)
        )


def resolve_comparison_candidates(trait_type, args):
    requested = _comma_tokens(getattr(args, "models", "all"), "--models")
    automatic = len(requested) == 1 and requested[0].lower() == "all"
    if any(token.lower() == "all" for token in requested) and not automatic:
        raise ValueError("'--models all' cannot be combined with named models.")

    parsed = _resolve_requested_models(trait_type, requested, automatic)
    candidates = [_candidate(model, variant, args) for model, variant in parsed]
    _reject_duplicate_candidates(candidates)
    exclusions = _parse_exclusions(args)
    candidates = [
        candidate for candidate in candidates if not _excluded(candidate, exclusions)
    ]
    if not candidates:
        raise ValueError("No ASR models remain after applying --exclude-models.")
    return candidates, automatic


def _parse_columns(value, option_name):
    columns = tuple(item.strip() for item in str(value).split(","))
    if any(column == "" for column in columns):
        raise ValueError(f"'{option_name}' contains an empty column name.")
    if len(columns) != len(set(columns)):
        raise ValueError(f"'{option_name}' contains duplicated column names.")
    return columns


def _parse_error_columns(value, trait_columns):
    if value in (None, ""):
        return None
    columns = _parse_columns(value, "--standard-error-column")
    if len(columns) != len(trait_columns):
        raise ValueError(
            "--standard-error-column must list one column for every --state-column."
        )
    if set(columns) & set(trait_columns):
        raise ValueError(
            "Standard-error columns must differ from the corresponding trait columns."
        )
    return columns


def _option_supplied(args, name):
    return getattr(args, name, None) not in (None, "")


def _candidate_consumes_option(candidate, name):
    if name == "root_mean":
        return candidate.model == "OU" and candidate.root_prior in {"fixed", "gaussian"}
    if name == "root_variance":
        return candidate.model == "OU" and candidate.root_prior == "gaussian"
    return candidate.model in MODEL_OPTION_CONSUMERS.get(name, ())


def _candidate_args(args, candidate):
    candidate_args = copy(args)
    candidate_args.model = candidate.model
    candidate_args.root_prior = candidate.root_prior
    for name in MODEL_OPTION_CONSUMERS:
        if not _candidate_consumes_option(candidate, name):
            setattr(candidate_args, name, None)
    for name in ("root_mean", "root_variance"):
        if not _candidate_consumes_option(candidate, name):
            setattr(candidate_args, name, None)
    return candidate_args


def _validate_option_consumers(context, candidates, *, automatic):
    args = context.args
    for name in MODEL_OPTION_CONSUMERS:
        if not _option_supplied(args, name):
            continue
        consumers = [
            candidate.model_id
            for candidate in candidates
            if _candidate_consumes_option(candidate, name)
            and _candidate_inapplicability(context, candidate) is None
            and not _automatic_diagnostic_skip(candidate, automatic)
        ]
        if not consumers:
            option = "--" + name.replace("_", "-")
            raise ValueError(
                f"{option} is not used by any selected applicable ASR comparison model."
            )


def _validate_transform_values(context, candidates):
    from nwkit.evolution import validate_evolution_parameter
    from nwkit.transformed_continuous_asr import parse_parameter_bounds

    args = context.args
    fixed_options = (
        ("evolution_parameter", getattr(args, "evolution_parameter", None)),
        ("eb_rate", getattr(args, "eb_rate", None)),
    )
    for option_name, value in fixed_options:
        if value is None:
            continue
        for candidate in candidates:
            if not _candidate_consumes_option(candidate, option_name):
                continue
            try:
                validate_evolution_parameter(candidate.model.lower(), value)
            except ValueError as exc:
                option = "--" + option_name.replace("_", "-")
                raise ValueError(
                    f"{option}={value} is invalid for selected model "
                    f"{candidate.model_id}: {exc}"
                ) from exc

    bound_options = (
        (
            "evolution_parameter_bounds",
            getattr(args, "evolution_parameter_bounds", None),
        ),
        ("eb_rate_bounds", getattr(args, "eb_rate_bounds", None)),
    )
    for option_name, value in bound_options:
        if value in (None, ""):
            continue
        for candidate in candidates:
            if not _candidate_consumes_option(candidate, option_name):
                continue
            try:
                parse_parameter_bounds(
                    value,
                    context.tree,
                    candidate.model.lower(),
                )
            except ValueError as exc:
                option = "--" + option_name.replace("_", "-")
                raise ValueError(
                    f"{option}={value} is invalid for selected model "
                    f"{candidate.model_id}: {exc}"
                ) from exc


def _validate_root_options(context, candidates):
    args = context.args
    trait_type = context.trait_type
    root_prior = getattr(args, "root_prior", None)
    allowed_roots = (
        {"equal", "empirical", "stationary", "gaussian"}
        if trait_type == "discrete"
        else {"flat", "stationary", "fixed", "gaussian"}
    )
    if root_prior is not None and root_prior not in allowed_roots:
        raise ValueError(
            "--root-prior={} is not available for {} comparison; choose {}.".format(
                root_prior, trait_type, ", ".join(sorted(allowed_roots))
            )
        )
    root_mean = getattr(args, "root_mean", None)
    root_variance = getattr(args, "root_variance", None)
    candidate_roots = {
        candidate.root_prior
        for candidate in candidates
        if candidate.model == ROOT_VARIANT_MODEL
    }
    proper_roots = candidate_roots & {"fixed", "gaussian"}
    if root_mean is not None and not proper_roots:
        raise ValueError("--root-mean requires an OU[fixed] or OU[gaussian] candidate.")
    if root_variance is not None and "gaussian" not in candidate_roots:
        raise ValueError("--root-variance requires an OU[gaussian] candidate.")
    if proper_roots and root_mean is None:
        raise ValueError("OU[fixed] and OU[gaussian] candidates require --root-mean.")
    if "gaussian" in candidate_roots and (
        root_variance is None or float(root_variance) <= 0.0
    ):
        raise ValueError("OU[gaussian] requires positive --root-variance.")


def _validate_mutually_exclusive_options(args):
    if (
        getattr(args, "alpha", None) is not None
        and getattr(args, "alpha_bounds", None) is not None
    ):
        raise ValueError("--alpha-bounds cannot be combined with fixed --alpha.")
    if (
        getattr(args, "evolution_parameter", None) is not None
        and getattr(args, "evolution_parameter_bounds", None) is not None
    ):
        raise ValueError(
            "--evolution-parameter-bounds cannot be combined with fixed "
            "--evolution-parameter."
        )
    if (
        getattr(args, "eb_rate", None) is not None
        and getattr(args, "eb_rate_bounds", None) is not None
    ):
        raise ValueError("--eb-rate-bounds cannot be combined with fixed --eb-rate.")
    generic_transform_options = {
        name
        for name in ("evolution_parameter", "evolution_parameter_bounds")
        if getattr(args, name, None) is not None
    }
    eb_transform_options = {
        name
        for name in ("eb_rate", "eb_rate_bounds")
        if getattr(args, name, None) is not None
    }
    if generic_transform_options and eb_transform_options:
        raise ValueError(
            "--evolution-parameter/--evolution-parameter-bounds cannot be combined "
            "with --eb-rate/--eb-rate-bounds."
        )


def _validate_regime_options(args, trait_type):
    if getattr(args, "regime_parameters", None) not in (None, "") and getattr(
        args, "regime_map", None
    ) in (None, ""):
        raise ValueError("--regime-parameters requires --regime-map.")
    if getattr(args, "regime_model", None) not in (None, "") and getattr(
        args, "regime_map", None
    ) in (None, ""):
        raise ValueError("--regime-model requires --regime-map.")
    if trait_type == "continuous" and getattr(args, "regime_model", None) not in (
        None,
        "",
    ):
        raise ValueError("--regime-model is available only for discrete comparison.")
    if trait_type == "discrete" and getattr(args, "regime_parameters", None) not in (
        None,
        "",
    ):
        raise ValueError(
            "--regime-parameters is available only for continuous comparison."
        )


def _validate_regime_parameter_contract(context, candidates):
    if not _option_supplied(context.args, "regime_parameters"):
        return
    contracts = {
        candidate.model: REGIME_PARAMETER_COLUMNS[candidate.model]
        for candidate in candidates
        if candidate.model in REGIME_PARAMETER_COLUMNS
        and _candidate_inapplicability(context, candidate) is None
    }
    distinct = {columns for columns in contracts.values()}
    if len(distinct) > 1:
        details = "; ".join(
            f"{model}={','.join(columns)}"
            for model, columns in sorted(contracts.items())
        )
        raise ValueError(
            "--regime-parameters cannot be shared by selected models with different "
            f"required columns ({details}); select one continuous regime-model "
            "contract."
        )


def _validate_candidate_model_options(context, candidates):
    from nwkit.asr_input import _validate_model_arguments

    for candidate in candidates:
        if _candidate_inapplicability(context, candidate) is not None:
            continue
        candidate_args = _candidate_args(context.args, candidate)
        _validate_model_arguments(
            candidate_args,
            context.trait_type,
            candidate.model,
            candidate.root_prior,
        )


def _validate_comparison_options(context, candidates, *, automatic=False):
    args = context.args
    _validate_option_consumers(context, candidates, automatic=automatic)
    _validate_root_options(context, candidates)
    _validate_mutually_exclusive_options(args)
    _validate_regime_options(args, context.trait_type)
    _validate_regime_parameter_contract(context, candidates)
    alpha = getattr(args, "alpha", None)
    if alpha is not None and float(alpha) <= 0.0:
        raise ValueError("--alpha must be strictly positive for OU models.")
    _validate_transform_values(context, candidates)
    _validate_candidate_model_options(context, candidates)


def _candidate_inapplicability(context, candidate):
    definition = model_definition(candidate.model)
    if candidate.root_prior not in definition.root_priors:
        return f"root prior {candidate.root_prior} is unsupported; choose " + ", ".join(
            definition.root_priors
        )
    multiple = len(context.trait_columns) > 1
    if context.trait_type == "discrete":
        if multiple and candidate.model != "MK-MIXTURE":
            return (
                "multiple character columns are currently compared only by MK-MIXTURE"
            )
        if not multiple and candidate.model == "MK-MIXTURE":
            return "MK-MIXTURE requires at least two character columns"
    else:
        if multiple and candidate.model not in MULTIVARIATE_MODELS:
            return "multiple continuous traits require MV-BM or MV-OU"
        if not multiple and candidate.model in MULTIVARIATE_MODELS:
            return "multivariate models require at least two trait columns"
    if candidate.model in REGIME_MODELS and getattr(
        context.args, "regime_map", None
    ) in (None, ""):
        return "the model requires --regime-map"
    if candidate.model == "CUSTOM" and getattr(context.args, "rate_matrix", None) in (
        None,
        "",
    ):
        return "CUSTOM requires --rate-matrix"
    if candidate.model == "THRESHOLD" and getattr(context.args, "states", None) in (
        None,
        "",
    ):
        return "THRESHOLD requires ordered --states"
    return None


def _cached(context, key, factory):
    if key not in context.cache:
        try:
            context.cache[key] = factory()
        except Exception as exc:
            context.cache[key] = _CachedFailure(exc)
    cached = context.cache[key]
    if isinstance(cached, _CachedFailure):
        raise cached.error
    return cached


def _option_or_default(args, name, default):
    value = getattr(args, name, None)
    return default if value is None else value


def _single_discrete_data(context):
    def build():
        from nwkit.asr import DEFAULT_AMBIGUOUS_SEPARATOR, _read_tip_states

        return _read_tip_states(
            context.args.trait,
            context.trait_columns[0],
            context.leaf_names,
            states_arg=getattr(context.args, "states", None),
            missing_values_arg=getattr(context.args, "missing_values", None),
            ambiguous_separator=_option_or_default(
                context.args, "ambiguous_separator", DEFAULT_AMBIGUOUS_SEPARATOR
            ),
            unmatched=getattr(context.args, "unmatched", "warn"),
            trait_df=context.trait_df,
        )

    return _cached(context, "single_discrete_data", build)


def _discrete_transition_graph(context, states):
    def build():
        from nwkit.discrete_asr_models import read_transition_graph

        return read_transition_graph(
            getattr(context.args, "transition_graph", None),
            states,
            state_source=(
                "--states"
                if getattr(context.args, "states", None) not in (None, "")
                else "--trait"
            ),
        )[0]

    return _cached(context, "discrete_transition_graph", build)


def _regime_assignment(context):
    def build():
        from nwkit.asr_regimes import read_regime_map

        return read_regime_map(context.args.regime_map, context.tree)

    return _cached(context, "regime_assignment", build)


def _custom_discrete_data(context):
    def build():
        from nwkit.asr import (
            DEFAULT_AMBIGUOUS_SEPARATOR,
            _parse_state_argument,
            _read_tip_states,
        )
        from nwkit.discrete_asr_models import read_rate_matrix

        matrix_states, matrix = read_rate_matrix(context.args.rate_matrix)
        explicit = getattr(context.args, "states", None)
        if explicit not in (None, "") and _parse_state_argument(explicit) != list(
            matrix_states
        ):
            raise ValueError(
                "--states must exactly match the state order in --rate-matrix."
            )
        states, observed, likelihoods = _read_tip_states(
            context.args.trait,
            context.trait_columns[0],
            context.leaf_names,
            states_arg=matrix_states,
            missing_values_arg=getattr(context.args, "missing_values", None),
            ambiguous_separator=_option_or_default(
                context.args, "ambiguous_separator", DEFAULT_AMBIGUOUS_SEPARATOR
            ),
            unmatched=getattr(context.args, "unmatched", "warn"),
            trait_df=context.trait_df,
            state_source="--rate-matrix",
        )
        return states, observed, likelihoods, matrix

    return _cached(context, "custom_discrete_data", build)


def _fit_single_discrete(context, candidate):
    from nwkit.asr import (
        _parse_rate_bounds,
        compute_covarion_marginals,
        compute_hrm_marginals,
        compute_mk_marginals,
    )

    candidate_args = _candidate_args(context.args, candidate)
    if candidate.model == "CUSTOM":
        states, observed, likelihoods, fixed_rate_matrix = _custom_discrete_data(
            context
        )
    else:
        states, observed, likelihoods = _single_discrete_data(context)
        fixed_rate_matrix = None
    rate_bounds = _parse_rate_bounds(getattr(candidate_args, "rate_bounds", None))
    graph = (
        None
        if candidate.model == "CUSTOM"
        else _discrete_transition_graph(context, states)
    )
    common = {
        "tree": context.tree,
        "states": states,
        "observed_state_by_leaf": observed,
        "likelihood_by_leaf": likelihoods,
        "rate": getattr(candidate_args, "rate", None),
        "root_prior_mode": candidate.root_prior,
        "rate_bounds": rate_bounds,
        "transition_graph": graph,
    }
    if candidate.model == "HRM":
        return compute_hrm_marginals(
            **common,
            hidden_categories=_option_or_default(
                candidate_args, "hidden_categories", 2
            ),
        )[1]
    if candidate.model == "COVARION":
        return compute_covarion_marginals(
            **common,
            hidden_categories=_option_or_default(
                candidate_args, "hidden_categories", 2
            ),
        )[1]
    return compute_mk_marginals(
        **common,
        model=candidate.model,
        fixed_rate_matrix=fixed_rate_matrix,
        regime_assignment=(
            _regime_assignment(context) if candidate.model == "MK-REGIME" else None
        ),
        regime_model=getattr(candidate_args, "regime_model", None) or "ER",
    )[1]


def _multi_discrete_data(context):
    def build():
        from nwkit.asr import (
            DEFAULT_AMBIGUOUS_SEPARATOR,
            _parse_state_argument,
            _read_tip_states,
        )

        explicit = getattr(context.args, "states", None)
        if explicit in (None, ""):
            states = []
            seen = set()
            for character in context.trait_columns:
                character_states, _observed, _likelihoods = _read_tip_states(
                    context.args.trait,
                    character,
                    context.leaf_names,
                    missing_values_arg=getattr(context.args, "missing_values", None),
                    ambiguous_separator=_option_or_default(
                        context.args,
                        "ambiguous_separator",
                        DEFAULT_AMBIGUOUS_SEPARATOR,
                    ),
                    unmatched=getattr(context.args, "unmatched", "warn"),
                    trait_df=context.trait_df,
                )
                for state in character_states:
                    if state not in seen:
                        states.append(state)
                        seen.add(state)
        else:
            states = _parse_state_argument(explicit)
        observed_by_character = {}
        likelihoods_by_character = {}
        for character in context.trait_columns:
            _states, observed, likelihoods = _read_tip_states(
                context.args.trait,
                character,
                context.leaf_names,
                states_arg=states,
                missing_values_arg=getattr(context.args, "missing_values", None),
                ambiguous_separator=_option_or_default(
                    context.args,
                    "ambiguous_separator",
                    DEFAULT_AMBIGUOUS_SEPARATOR,
                ),
                unmatched=getattr(context.args, "unmatched", "warn"),
                trait_df=context.trait_df,
                state_source="--states"
                if explicit not in (None, "")
                else "joint characters",
            )
            observed_by_character[character] = observed
            likelihoods_by_character[character] = likelihoods
        return states, observed_by_character, likelihoods_by_character

    return _cached(context, "multi_discrete_data", build)


def _fit_mk_mixture(context, candidate):
    from nwkit.asr import _parse_rate_bounds, compute_mk_mixture_marginals

    candidate_args = _candidate_args(context.args, candidate)
    states, observed, likelihoods = _multi_discrete_data(context)
    graph = _discrete_transition_graph(context, states)
    return compute_mk_mixture_marginals(
        context.tree,
        states,
        observed,
        likelihoods,
        base_model=getattr(candidate_args, "mixture_model", None) or "ER",
        mixture=getattr(candidate_args, "rate_mixture", None) or "gamma",
        categories=_option_or_default(candidate_args, "rate_categories", 4),
        rate=getattr(candidate_args, "rate", None),
        root_prior_mode=candidate.root_prior,
        rate_bounds=_parse_rate_bounds(getattr(candidate_args, "rate_bounds", None)),
        transition_graph=graph,
    )[1]


def _continuous_data(context):
    def build():
        if len(context.trait_columns) == 1:
            observed, errors = continuous_tip_values(
                context.trait_df,
                context.trait_columns[0],
                context.leaf_names,
                None if context.error_columns is None else context.error_columns[0],
            )
        else:
            observed = continuous_tip_vectors(
                context.trait_df, context.trait_columns, context.leaf_names
            )
            errors = continuous_tip_vector_errors(
                context.trait_df,
                context.error_columns,
                context.trait_columns,
                context.leaf_names,
            )
        return observed, errors

    return _cached(context, "continuous_data", build)


def _fit_continuous(context, candidate):
    from nwkit.asr import _fit_continuous_model

    observed, errors = _continuous_data(context)
    settings = AsrSettings(
        trait_type="continuous",
        model=candidate.model,
        root_prior=candidate.root_prior,
        output="summary",
        tree_annotation="summary",
        ci_level=0.95,
    )
    candidate_args = _candidate_args(context.args, candidate)
    regime_assignment = (
        _regime_assignment(context) if candidate.model in REGIME_MODELS else None
    )
    return _fit_continuous_model(
        context.tree,
        observed,
        errors,
        context.trait_columns,
        candidate_args,
        settings,
        regime_assignment,
    )[1]


def _fit_candidate(context, candidate):
    if context.trait_type == "continuous":
        return _fit_continuous(context, candidate)
    if candidate.model == "MK-MIXTURE":
        return _fit_mk_mixture(context, candidate)
    return _fit_single_discrete(context, candidate)


def _fit_value(fit, name, default=None):
    return (
        fit.get(name, default) if isinstance(fit, dict) else getattr(fit, name, default)
    )


def _base_row(context, candidate, elapsed, *, status, message=""):
    return {
        "model_id": candidate.model_id,
        "model": candidate.model,
        "trait_type": context.trait_type,
        "trait_type_requested": getattr(context.args, "trait_type", "auto"),
        "trait_columns": ",".join(context.trait_columns),
        "num_traits": len(context.trait_columns),
        "comparison_group": "",
        "likelihood_kind": "",
        "root_prior": candidate.root_prior,
        "root_prior_source": candidate.root_prior_source,
        "equivalent_to": "",
        "status": status,
        "rankable": "no",
        "message": " ".join(str(message).split()),
        "fit_status": "",
        "optimizer_success": "",
        "optimizer_message": "",
        "sample_size": np.nan,
        "num_parameters": np.nan,
        "estimated_parameters": "",
        "fixed_parameters": "",
        "log_likelihood": np.nan,
        "shared_preparation_seconds": float(
            context.cache.get("shared_preparation_seconds", 0.0)
        ),
        "elapsed_seconds": float(elapsed),
    }


def _comparison_group(context, candidate, likelihood_kind):
    if context.trait_type == "continuous":
        data_kind = (
            f"{len(context.trait_columns)}d"
            if len(context.trait_columns) > 1
            else "scalar"
        )
    else:
        data_kind = (
            f"{len(context.trait_columns)}char"
            if len(context.trait_columns) > 1
            else "scalar"
        )
    return f"{context.trait_type}:{data_kind}:{likelihood_kind}:{candidate.root_prior}"


def _classify_fit(candidate, fit, summary):
    fit_status = str(summary.get("fit_status", "ok"))
    likelihood = summary.get("log_likelihood")
    if likelihood is None or not math.isfinite(float(likelihood)):
        return "no_likelihood", "no", "The fit has no finite marginal likelihood."
    if candidate.model == "HRM":
        return (
            "nonregular",
            "no",
            "Hidden-class label switching makes regular-model IC ranking invalid.",
        )
    if "singular" in fit_status:
        return "nonregular", "no", f"Non-regular fit status: {fit_status}."
    if has_nonregular_variance_boundary(fit_status):
        return (
            "nonregular",
            "no",
            f"A fitted variance component reached zero: {fit_status}.",
        )
    if candidate.model == "COVARION" and "boundary" in fit_status:
        return "nonregular", "no", f"Non-regular covarion boundary: {fit_status}."
    optimizer_success = _fit_value(fit, "optimizer_success", None)
    if optimizer_success is False:
        return "nonconverged", "no", "The optimizer did not report convergence."
    if "boundary" in fit_status:
        return "boundary", "yes", f"Boundary fit retained: {fit_status}."
    return "ok", "yes", ""


def _parameter_names(fit, value_name, base_name):
    values = _fit_value(fit, value_name, None)
    if isinstance(values, dict):
        return [f"{base_name}[{key}]" for key in values]
    return [base_name]


def _route_parameters(fit, names, flag_name, estimated, fixed):
    target = estimated if bool(_fit_value(fit, flag_name, False)) else fixed
    target.extend(names)


def _discrete_parameter_contract(candidate, fit):
    estimated: list[str] = []
    fixed: list[str] = []
    model = candidate.model
    if model == "CUSTOM":
        fixed.append("Q")
    elif model == "COVARION":
        target = estimated if _fit_value(fit, "base_rate_estimated", True) else fixed
        target.append("base_rate")
        estimated.extend(("hidden_rate_spread", "switching_rate"))
    elif model == "MK-MIXTURE":
        target = estimated if _fit_value(fit, "base_rate_estimated", True) else fixed
        target.append("base_rates")
        estimated.append(
            "gamma_shape"
            if _fit_value(fit, "rate_mixture") == "gamma"
            else "free_rate_weights"
        )
    else:
        target = estimated if _fit_value(fit, "rate_estimated", False) else fixed
        target.append("transition_rates")
    return estimated, fixed


def _continuous_parameter_contract(candidate, fit):
    estimated: list[str] = []
    fixed: list[str] = []
    model = candidate.model
    if model in MULTIVARIATE_MODELS:
        estimated.append("trait_covariance")
        if model == "MV-OU":
            _route_parameters(fit, ["alpha"], "alpha_estimated", estimated, fixed)
            theta_names = [
                f"theta[{name}]" for name in _fit_value(fit, "trait_names", ())
            ]
            _route_parameters(fit, theta_names, "theta_estimated", estimated, fixed)
    elif model in {"BMS", "BMS-DRIFT"}:
        names = _parameter_names(fit, "sigma2_by_regime", "sigma2")
        _route_parameters(fit, names, "sigma2_estimated", estimated, fixed)
        if model == "BMS-DRIFT":
            names = _parameter_names(fit, "drift_by_regime", "drift")
            _route_parameters(fit, names, "drift_estimated", estimated, fixed)
    elif model in {"OUM", "OUMA", "OUMV", "OUMVA"}:
        alpha_names = (
            _parameter_names(fit, "alpha_by_regime", "alpha")
            if model in {"OUMA", "OUMVA"}
            else ["alpha"]
        )
        sigma_names = (
            _parameter_names(fit, "sigma2_by_regime", "sigma2")
            if model in {"OUMV", "OUMVA"}
            else ["sigma2"]
        )
        _route_parameters(fit, alpha_names, "alpha_estimated", estimated, fixed)
        _route_parameters(
            fit,
            _parameter_names(fit, "theta_by_regime", "theta"),
            "theta_estimated",
            estimated,
            fixed,
        )
        _route_parameters(fit, sigma_names, "sigma2_estimated", estimated, fixed)
    else:
        _route_parameters(fit, ["sigma2"], "sigma2_estimated", estimated, fixed)
        if model in {"LAMBDA", "KAPPA", "DELTA", "EB", "ACDC"}:
            name = _fit_value(fit, "evolution_parameter_name", "transform")
            _route_parameters(
                fit,
                [str(name)],
                "evolution_parameter_estimated",
                estimated,
                fixed,
            )
        elif model == "BM-DRIFT":
            _route_parameters(fit, ["drift"], "drift_estimated", estimated, fixed)
        elif model == "OU":
            _route_parameters(fit, ["alpha"], "alpha_estimated", estimated, fixed)
            _route_parameters(fit, ["theta"], "theta_estimated", estimated, fixed)
    if candidate.root_prior in {"fixed", "gaussian"}:
        fixed.append("root_mean")
    if candidate.root_prior == "gaussian":
        fixed.append("root_variance")
    return estimated, fixed


def _parameter_contract(context, candidate, fit):
    estimated, fixed = (
        _discrete_parameter_contract(candidate, fit)
        if context.trait_type == "discrete"
        else _continuous_parameter_contract(candidate, fit)
    )
    if context.error_columns is not None:
        fixed.append("measurement_error")
    return ",".join(estimated), ",".join(fixed)


def _successful_row(context, candidate, fit, elapsed):
    summary = summarize_fit(
        candidate.model,
        fit,
        trait_type=context.trait_type,
        root_prior=candidate.root_prior,
    )
    likelihood = summary.get("log_likelihood")
    if likelihood is not None and math.isfinite(float(likelihood)):
        if int(summary["num_parameters"]) < 0:
            raise ValueError("the fit returned a negative parameter count")
        if int(summary["sample_size"]) <= 0:
            raise ValueError("the fit returned a non-positive sample size")
    status, rankable, message = _classify_fit(candidate, fit, summary)
    estimated_parameters, fixed_parameters = _parameter_contract(
        context, candidate, fit
    )
    row = _base_row(context, candidate, elapsed, status=status, message=message)
    optimizer_success = _fit_value(fit, "optimizer_success", None)
    row.update(
        {
            "comparison_group": (
                ""
                if status == "no_likelihood"
                else _comparison_group(context, candidate, summary["likelihood_kind"])
            ),
            "likelihood_kind": summary["likelihood_kind"],
            "rankable": rankable,
            "fit_status": str(summary.get("fit_status", "ok")),
            "optimizer_success": (
                ""
                if optimizer_success is None
                else "yes"
                if optimizer_success
                else "no"
            ),
            "optimizer_message": " ".join(
                str(_fit_value(fit, "optimizer_message", "")).split()
            ),
            "sample_size": summary["sample_size"],
            "num_parameters": summary["num_parameters"],
            "estimated_parameters": estimated_parameters,
            "fixed_parameters": fixed_parameters,
            "log_likelihood": summary["log_likelihood"],
        }
    )
    return row


def _automatic_diagnostic_skip(candidate, automatic):
    return automatic and candidate.model == "HRM"


def _prepare_shared_inputs(context, candidates, *, automatic):
    if "shared_preparation_seconds" in context.cache:
        return
    started = context.cache.get("_shared_preparation_started")
    if started is None:
        started = time.perf_counter()
    applicable = [
        candidate
        for candidate in candidates
        if _candidate_inapplicability(context, candidate) is None
        and candidate.model != "THRESHOLD"
        and not _automatic_diagnostic_skip(candidate, automatic)
    ]
    try:
        if context.trait_type == "continuous" and applicable:
            _continuous_data(context)
            if any(candidate.model in REGIME_MODELS for candidate in applicable):
                _regime_assignment(context)
        elif context.trait_type == "discrete" and applicable:
            if len(context.trait_columns) > 1:
                if any(candidate.model == "MK-MIXTURE" for candidate in applicable):
                    states, _observed, _likelihoods = _multi_discrete_data(context)
                    _discrete_transition_graph(context, states)
            else:
                if any(candidate.model != "CUSTOM" for candidate in applicable):
                    states, _observed, _likelihoods = _single_discrete_data(context)
                    if any(
                        candidate.model in DISCRETE_RATE_MODELS
                        for candidate in applicable
                    ):
                        _discrete_transition_graph(context, states)
                if any(candidate.model == "CUSTOM" for candidate in applicable):
                    _custom_discrete_data(context)
                if any(candidate.model == "MK-REGIME" for candidate in applicable):
                    _regime_assignment(context)
    except Exception:
        # _cached stores and re-raises the shared input failure.  Automatic
        # comparison still records it independently for each affected model.
        pass
    context.cache["shared_preparation_seconds"] = time.perf_counter() - started


def _transform_parameter(args, model):
    value = getattr(args, "evolution_parameter", None)
    if model in {"EB", "ACDC"} and value is None:
        value = getattr(args, "eb_rate", None)
    return value


def _transform_bounds(args, model):
    value = getattr(args, "evolution_parameter_bounds", None)
    if model in {"EB", "ACDC"} and value in (None, ""):
        value = getattr(args, "eb_rate_bounds", None)
    return value


def _one_regime(context):
    try:
        assignment = _regime_assignment(context)
    except Exception:
        return None
    return assignment if len(assignment.regimes) == 1 else None


def _continuous_equivalence_contract(context, candidate):
    if context.trait_type != "continuous" or len(context.trait_columns) != 1:
        return None
    model = candidate.model
    root = candidate.root_prior
    args = context.args
    bm_key = ("continuous", "BM", root)
    if model == "BM":
        return bm_key, "under the same Brownian process contract"
    if model in {"LAMBDA", "KAPPA", "DELTA"}:
        parameter = _transform_parameter(args, model)
        if parameter is not None and float(parameter) == 1.0:
            return bm_key, f"with the neutral fixed {model.lower()} parameter"
        return None
    if model in {"EB", "ACDC"}:
        parameter = _transform_parameter(args, model)
        if parameter is not None:
            if float(parameter) == 0.0:
                return bm_key, "with a fixed zero exponential-rate parameter"
            return (
                ("continuous", "exponential-rate", root, "fixed", float(parameter)),
                "under the same fixed exponential-rate parameter",
            )
        bounds = _transform_bounds(args, model)
        if bounds not in (None, ""):
            return (
                ("continuous", "exponential-rate", root, "bounded", str(bounds)),
                "under the same explicitly bounded nonpositive parameter space",
            )
        return None
    if model == "BM-DRIFT":
        drift = getattr(args, "drift", None)
        if drift is not None and float(drift) == 0.0:
            return bm_key, "with fixed zero drift"
        return (
            ("continuous", "BM-DRIFT", root, None if drift is None else float(drift)),
            "under the same Brownian-drift process contract",
        )
    if model in {"BMS", "BMS-DRIFT"}:
        if getattr(args, "regime_parameters", None) not in (None, ""):
            return None
        if _one_regime(context) is None:
            return None
        if model == "BMS":
            return bm_key, "because the regime map contains one regime"
        drift = getattr(args, "drift", None)
        if drift is not None and float(drift) == 0.0:
            return bm_key, "because one regime with fixed zero drift reduces to BM"
        return (
            ("continuous", "BM-DRIFT", root, None if drift is None else float(drift)),
            "because the regime map contains one regime",
        )
    if model == "OU" and root == "stationary":
        return (
            ("continuous", "OU-stationary", root),
            "under the same stationary OU process contract",
        )
    if model in {"OUM", "OUMA", "OUMV", "OUMVA"}:
        if getattr(args, "regime_parameters", None) not in (None, ""):
            return None
        if _one_regime(context) is None:
            return None
        return (
            ("continuous", "OU-stationary", root),
            "because the regime map contains one regime",
        )
    return None


def _discrete_equivalence_contract(context, candidate, states, graph):
    if context.trait_type != "discrete" or len(context.trait_columns) != 1:
        return None
    model = candidate.model
    is_regime = model == "MK-REGIME"
    if is_regime:
        if _one_regime(context) is None:
            return None
        model = getattr(context.args, "regime_model", None) or "ER"
    if model not in {"ER", "SYM", "ARD", "F81", "GTR"}:
        return None
    rate = getattr(context.args, "rate", None)
    if candidate.model == "ER" and rate is not None:
        # For ordinary ER --rate is fixed, whereas it is only an optimizer
        # starting value for every other fitted structured model.
        return None
    # ARD/F81 share direct q_ij bounds for two states. GTR deliberately remains
    # distinct because it bounds exchangeability and frequency-ratio coordinates.
    family = model_equivalence_family(model, states, graph)
    source = "a one-regime MK model" if is_regime else "this transition model"
    return (
        (
            "discrete",
            family,
            candidate.root_prior,
            graph.shape,
            graph.tobytes(),
        ),
        f"for {source} and the same transition/root contract",
    )


def _equivalent_candidate_representatives(context, candidates):
    states = graph = None
    if context.trait_type == "discrete" and len(context.trait_columns) == 1:
        try:
            states, _observed, _likelihoods = _single_discrete_data(context)
            graph = _discrete_transition_graph(context, states)
        except Exception:
            return {}
    representatives: dict[str, tuple[str, str]] = {}
    first_by_key: dict[tuple[Any, ...], str] = {}
    for candidate in candidates:
        contract = (
            _continuous_equivalence_contract(context, candidate)
            if context.trait_type == "continuous"
            else _discrete_equivalence_contract(
                context,
                candidate,
                states,
                graph,  # type: ignore[arg-type]
            )
        )
        if contract is None:
            continue
        key, reason = contract
        if key in first_by_key:
            representatives[candidate.model_id] = (first_by_key[key], reason)
        else:
            first_by_key[key] = candidate.model_id
    return representatives


def _equivalent_parameter_contract(context, candidate):
    args = context.args
    estimated: list[str] = []
    fixed: list[str] = []

    def route(name, value):
        (estimated if value is None else fixed).append(name)

    if context.trait_type == "discrete":
        estimated.append("transition_rates")
    else:
        model = candidate.model
        assignment = _one_regime(context) if model in REGIME_MODELS else None
        regime = None if assignment is None else assignment.regimes[0]
        sigma_name = f"sigma2[{regime}]" if model in {"BMS", "BMS-DRIFT"} else "sigma2"
        route(sigma_name, getattr(args, "sigma2", None))
        if model in TRANSFORMED_MODELS:
            route("evolution_parameter", _transform_parameter(args, model))
        elif model in {"BM-DRIFT", "BMS-DRIFT"}:
            drift_name = f"drift[{regime}]" if regime is not None else "drift"
            route(drift_name, getattr(args, "drift", None))
        elif model in {"OU", "OUM", "OUMA", "OUMV", "OUMVA"}:
            alpha_name = f"alpha[{regime}]" if model in {"OUMA", "OUMVA"} else "alpha"
            theta_name = f"theta[{regime}]" if regime is not None else "theta"
            sigma_name = f"sigma2[{regime}]" if model in {"OUMV", "OUMVA"} else "sigma2"
            estimated.clear()
            fixed.clear()
            route(alpha_name, getattr(args, "alpha", None))
            route(theta_name, getattr(args, "theta", None))
            route(sigma_name, getattr(args, "sigma2", None))
    if context.error_columns is not None:
        fixed.append("measurement_error")
    return ",".join(estimated), ",".join(fixed)


def _equivalent_row(context, candidate, representative, reason):
    if representative.get("comparison_group", "") == "":
        row = _base_row(
            context,
            candidate,
            0.0,
            status="failed",
            message=(
                f"Equivalent representative {representative['model_id']} did not "
                "produce a reusable fit."
            ),
        )
        row["equivalent_to"] = representative["model_id"]
        return row
    row = dict(representative)
    estimated_parameters, fixed_parameters = _equivalent_parameter_contract(
        context, candidate
    )
    row.update(
        {
            "model_id": candidate.model_id,
            "model": candidate.model,
            "root_prior": candidate.root_prior,
            "root_prior_source": candidate.root_prior_source,
            "equivalent_to": representative["model_id"],
            "status": "equivalent",
            "rankable": "no",
            "message": (
                f"Statistically equivalent to {representative['model_id']} {reason}; "
                "excluded from duplicate IC weighting."
            ),
            "fit_status": "equivalent",
            "optimizer_success": "",
            "optimizer_message": "",
            "estimated_parameters": estimated_parameters,
            "fixed_parameters": fixed_parameters,
            "elapsed_seconds": 0.0,
        }
    )
    return row


def evaluate_comparison_candidates(context, candidates, *, automatic):
    _prepare_shared_inputs(context, candidates, automatic=automatic)
    equivalent_to = _equivalent_candidate_representatives(context, candidates)
    rows = []
    rows_by_model_id = {}
    for candidate in candidates:
        started = time.perf_counter()
        reason = _candidate_inapplicability(context, candidate)
        if reason is not None:
            if not automatic:
                raise ValueError(
                    f"ASR comparison model '{candidate.model_id}' is not applicable: {reason}."
                )
            row = _base_row(
                context,
                candidate,
                time.perf_counter() - started,
                status="not_applicable",
                message=reason,
            )
            rows.append(row)
            rows_by_model_id[candidate.model_id] = row
            continue
        if candidate.model == "THRESHOLD":
            row = _base_row(
                context,
                candidate,
                time.perf_counter() - started,
                status="no_likelihood",
                message=(
                    "THRESHOLD uses posterior MCMC diagnostics and does not "
                    "provide a marginal likelihood for IC comparison."
                ),
            )
            rows.append(row)
            rows_by_model_id[candidate.model_id] = row
            continue
        if _automatic_diagnostic_skip(candidate, automatic):
            row = _base_row(
                context,
                candidate,
                time.perf_counter() - started,
                status="not_fitted",
                message=(
                    "Automatic comparison skips HRM because hidden-class label "
                    "switching invalidates regular-model IC ranking; request HRM "
                    "explicitly with --models to fit it diagnostically."
                ),
            )
            rows.append(row)
            rows_by_model_id[candidate.model_id] = row
            continue
        equivalence = equivalent_to.get(candidate.model_id)
        if equivalence is not None:
            representative_id, reason = equivalence
            row = _equivalent_row(
                context,
                candidate,
                rows_by_model_id[representative_id],
                reason,
            )
            rows.append(row)
            rows_by_model_id[candidate.model_id] = row
            continue
        sys.stderr.write(f"ASR comparison: fitting {candidate.model_id}.\n")
        try:
            fit = _fit_candidate(context, candidate)
            row = _successful_row(
                context, candidate, fit, time.perf_counter() - started
            )
        except Exception as exc:
            if not automatic:
                raise ValueError(
                    f"ASR comparison model '{candidate.model_id}' failed: {exc}"
                ) from exc
            row = _base_row(
                context,
                candidate,
                time.perf_counter() - started,
                status="failed",
                message=exc,
            )
            sys.stderr.write(
                f"Warning: ASR comparison model '{candidate.model_id}' failed: {exc}\n"
            )
        rows.append(row)
        rows_by_model_id[candidate.model_id] = row
    return rows


def comparison_table(context, candidates, *, automatic, criterion):
    rows = evaluate_comparison_candidates(context, candidates, automatic=automatic)
    table = grouped_model_comparison_table(rows, criterion=criterion)
    for column in ASR_COMPARE_COLUMNS:
        if column not in table:
            table[column] = np.nan
    table = table.loc[:, ASR_COMPARE_COLUMNS]
    for column in INTEGER_COLUMNS:
        table[column] = pd.array(table[column], dtype="Int64")
    return table


def _write_table(table, path):
    if path == "-":
        print(table.to_csv(sep="\t", index=False, na_rep=""), end="")
    else:
        table.to_csv(path, sep="\t", index=False, na_rep="")


def _validate_output_paths(args):
    figure = getattr(args, "figure_out", None)
    if figure == "-":
        raise ValueError("'--figure-out' must be a PDF file path, not STDOUT.")
    if figure not in (None, "") and os.path.splitext(str(figure))[1].lower() != ".pdf":
        raise ValueError("'--figure-out' must use the .pdf extension.")
    outputs = [("--outfile", args.outfile), ("--figure-out", figure)]
    validate_distinct_output_paths(outputs)
    inputs = [
        ("--infile", args.infile),
        ("--trait", args.trait),
        ("--regime-map", getattr(args, "regime_map", None)),
        ("--regime-parameters", getattr(args, "regime_parameters", None)),
        ("--rate-matrix", getattr(args, "rate_matrix", None)),
    ]
    transition_graph = getattr(args, "transition_graph", None)
    if transition_graph not in (None, "", "complete", "ordered"):
        inputs.append(("--transition-graph", transition_graph))
    validate_outputs_do_not_replace_inputs(
        inputs, outputs, label="ASR comparison output"
    )


def _write_outputs(table, args, criterion, *, include_figure):
    figure = getattr(args, "figure_out", None) if include_figure else None
    outputs = []
    if args.outfile != "-":
        outputs.append(args.outfile)
    if figure not in (None, ""):
        outputs.append(figure)
    if not outputs:
        _write_table(table, "-")
        return
    write_stdout = (lambda: _write_table(table, "-")) if args.outfile == "-" else None
    with output_transaction(outputs, after_install=write_stdout) as staged:
        if args.outfile != "-":
            _write_table(table, staged[args.outfile])
        if figure not in (None, ""):
            draw_comparison_figure(table, staged[figure], criterion)


def _has_completed_fit(table):
    return bool(
        table["status"].isin({"ok", "boundary", "nonregular", "nonconverged"}).any()
    )


def _preflight_comparison_figure(context, candidates, criterion):
    if getattr(context.args, "figure_out", None) in (None, ""):
        return
    predictable_text = [
        "ASR model comparison",
        _criterion_label(criterion),
        "Comparison set Not assigned Rank Model Fit status Notes",
        context.trait_type,
        *context.trait_columns,
        *(candidate.model_id for candidate in candidates),
    ]
    _font_family_for_text("\n".join(predictable_text))


def asr_compare_main(args):
    from nwkit.asr import _validate_tree_for_asr
    from nwkit.asr_input import _validate_mode_arguments

    _validate_output_paths(args)
    shared_preparation_started = time.perf_counter()
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    _validate_tree_for_asr(tree)
    trait_columns = _parse_columns(args.state_column, "--state-column")
    error_columns = _parse_error_columns(
        getattr(args, "standard_error_column", None), trait_columns
    )
    state_input = trait_columns if len(trait_columns) > 1 else trait_columns[0]
    trait_df = read_asr_table(
        args.trait,
        state_input,
        list(tree.leaf_names()),
        missing_values=getattr(args, "missing_values", None),
        unmatched=getattr(args, "unmatched", "warn"),
        standard_error_column=error_columns,
    )
    requested_type = getattr(args, "trait_type", "auto")
    trait_type = resolve_trait_type(requested_type, trait_df, state_input)
    _validate_mode_arguments(args, trait_type)
    candidates, automatic = resolve_comparison_candidates(trait_type, args)
    context = ComparisonContext(
        tree,
        trait_df,
        trait_type,
        trait_columns,
        error_columns,
        args,
    )
    context.cache["_shared_preparation_started"] = shared_preparation_started
    _validate_comparison_options(context, candidates, automatic=automatic)
    criterion = getattr(args, "criterion", "aic")
    _preflight_comparison_figure(context, candidates, criterion)
    sys.stderr.write(
        f"ASR comparison trait type: {trait_type} ({requested_type}); "
        f"{len(candidates)} candidate(s).\n"
    )
    if requested_type == "auto" and trait_type == "continuous":
        sys.stderr.write(
            "ASR comparison auto-detection treats numeric columns as continuous; "
            "use --trait-type discrete for numeric category codes.\n"
        )
    table = comparison_table(
        context,
        candidates,
        automatic=automatic,
        criterion=criterion,
    )
    completed = _has_completed_fit(table)
    _write_outputs(table, args, criterion, include_figure=completed)
    if not completed:
        raise ValueError(
            "Every selected ASR model failed or was not comparable; the summary TSV "
            "was written, and the PDF was left unchanged."
        )
    if not table["criterion_rank"].notna().any():
        sys.stderr.write(
            "Warning: no comparison group contains two finite models for the "
            f"selected {_criterion_label(criterion)} criterion.\n"
        )
    return table
