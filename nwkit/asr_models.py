"""Central ASR model definitions and mode-specific defaults."""

from dataclasses import dataclass

from nwkit.model_specs import evolution_model_spec


@dataclass(frozen=True)
class AsrModelDefinition:
    name: str
    trait_type: str
    default_root_prior: str
    root_priors: tuple[str, ...]
    process_model: str | None = None
    parameter_name: str | None = None


def _evolution_definition(model):
    spec = evolution_model_spec(model)
    if spec.asr_name is None or spec.asr_default_root_prior is None:
        raise ValueError(f"Evolution model '{model}' is not available for ASR.")
    return AsrModelDefinition(
        spec.asr_name,
        "continuous",
        spec.asr_default_root_prior,
        spec.asr_root_priors,
        process_model=spec.name,
        parameter_name=spec.parameter_name,
    )


_MODEL_DEFINITIONS = (
    AsrModelDefinition("ER", "discrete", "equal", ("equal", "empirical", "stationary")),
    AsrModelDefinition(
        "SYM", "discrete", "equal", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition(
        "ARD", "discrete", "equal", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition(
        "F81", "discrete", "stationary", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition(
        "GTR", "discrete", "stationary", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition(
        "MK-REGIME", "discrete", "equal", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition(
        "HRM", "discrete", "equal", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition(
        "COVARION", "discrete", "equal", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition(
        "MK-MIXTURE", "discrete", "equal", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition("THRESHOLD", "discrete", "gaussian", ("gaussian",)),
    AsrModelDefinition(
        "CUSTOM", "discrete", "equal", ("equal", "empirical", "stationary")
    ),
    _evolution_definition("brownian"),
    AsrModelDefinition("BMS", "continuous", "flat", ("flat",)),
    AsrModelDefinition("BMS-DRIFT", "continuous", "flat", ("flat",)),
    _evolution_definition("lambda"),
    _evolution_definition("kappa"),
    _evolution_definition("delta"),
    _evolution_definition("eb"),
    _evolution_definition("acdc"),
    AsrModelDefinition("BM-DRIFT", "continuous", "flat", ("flat",)),
    AsrModelDefinition("MV-BM", "continuous", "flat", ("flat",)),
    AsrModelDefinition("MV-OU", "continuous", "stationary", ("stationary",)),
    _evolution_definition("ou"),
    AsrModelDefinition("OUM", "continuous", "stationary", ("stationary",)),
    AsrModelDefinition("OUMA", "continuous", "stationary", ("stationary",)),
    AsrModelDefinition("OUMV", "continuous", "stationary", ("stationary",)),
    AsrModelDefinition("OUMVA", "continuous", "stationary", ("stationary",)),
)

MODEL_REGISTRY = {definition.name: definition for definition in _MODEL_DEFINITIONS}
_DEFAULT_MODELS = {"discrete": "ER", "continuous": "BM"}


def model_names(trait_type=None):
    return tuple(
        definition.name
        for definition in _MODEL_DEFINITIONS
        if trait_type is None or definition.trait_type == trait_type
    )


def default_model(trait_type):
    try:
        return _DEFAULT_MODELS[trait_type]
    except KeyError as exc:
        raise ValueError(f"Unsupported resolved trait type: {trait_type}.") from exc


def model_definition(name):
    try:
        return MODEL_REGISTRY[name]
    except KeyError as exc:
        raise ValueError(f"Unsupported '--model': {name}") from exc
