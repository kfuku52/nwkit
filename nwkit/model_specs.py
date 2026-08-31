"""Lightweight model metadata shared by CLI help and numerical implementations."""

from dataclasses import dataclass


@dataclass(frozen=True)
class EvolutionModelSpec:
    name: str
    parameter_name: str | None
    contrast_supported: bool = True
    branch_lengths_used: bool = True


EVOLUTION_MODEL_SPECS = {
    spec.name: spec
    for spec in [
        EvolutionModelSpec("brownian", None),
        EvolutionModelSpec("lambda", "lambda"),
        EvolutionModelSpec("ou", "alpha"),
        EvolutionModelSpec("kappa", "kappa"),
        EvolutionModelSpec("delta", "delta"),
        EvolutionModelSpec("eb", "rate_change"),
        EvolutionModelSpec("acdc", "rate_change"),
        EvolutionModelSpec("independent", None, branch_lengths_used=False),
        EvolutionModelSpec(
            "custom", None, contrast_supported=False, branch_lengths_used=False
        ),
    ]
}
EVOLUTION_MODELS = tuple(EVOLUTION_MODEL_SPECS)
CONTRAST_EVOLUTION_MODELS = tuple(
    name for name, spec in EVOLUTION_MODEL_SPECS.items() if spec.contrast_supported
)
PARAMETERIZED_EVOLUTION_MODELS = tuple(
    name for name, spec in EVOLUTION_MODEL_SPECS.items() if spec.parameter_name
)


def evolution_model_spec(model: str) -> EvolutionModelSpec:
    try:
        return EVOLUTION_MODEL_SPECS[str(model)]
    except KeyError as exc:
        raise ValueError("Unsupported evolutionary model: {}.".format(model)) from exc


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
