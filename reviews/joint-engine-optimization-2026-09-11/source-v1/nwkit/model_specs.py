"""Lightweight model metadata shared by CLI help and numerical implementations."""

from dataclasses import dataclass


@dataclass(frozen=True)
class EvolutionModelSpec:
    name: str
    parameter_name: str | None
    contrast_supported: bool = True
    branch_lengths_used: bool = True
    asr_name: str | None = None
    asr_default_root_prior: str | None = None
    asr_root_priors: tuple[str, ...] = ()


EVOLUTION_MODEL_SPECS = {
    spec.name: spec
    for spec in [
        EvolutionModelSpec(
            "brownian",
            None,
            asr_name="BM",
            asr_default_root_prior="flat",
            asr_root_priors=("flat",),
        ),
        EvolutionModelSpec(
            "lambda",
            "lambda",
            asr_name="LAMBDA",
            asr_default_root_prior="flat",
            asr_root_priors=("flat",),
        ),
        EvolutionModelSpec(
            "ou",
            "alpha",
            asr_name="OU",
            asr_default_root_prior="stationary",
            asr_root_priors=("stationary", "fixed", "gaussian"),
        ),
        EvolutionModelSpec(
            "kappa",
            "kappa",
            asr_name="KAPPA",
            asr_default_root_prior="flat",
            asr_root_priors=("flat",),
        ),
        EvolutionModelSpec(
            "delta",
            "delta",
            asr_name="DELTA",
            asr_default_root_prior="flat",
            asr_root_priors=("flat",),
        ),
        EvolutionModelSpec(
            "eb",
            "rate_change",
            asr_name="EB",
            asr_default_root_prior="flat",
            asr_root_priors=("flat",),
        ),
        EvolutionModelSpec(
            "acdc",
            "rate_change",
            asr_name="ACDC",
            asr_default_root_prior="flat",
            asr_root_priors=("flat",),
        ),
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
