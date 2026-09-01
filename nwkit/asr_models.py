"""Central ASR model definitions and mode-specific defaults."""

from dataclasses import dataclass


@dataclass(frozen=True)
class AsrModelDefinition:
    name: str
    trait_type: str
    default_root_prior: str
    root_priors: tuple[str, ...]


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
        "CUSTOM", "discrete", "equal", ("equal", "empirical", "stationary")
    ),
    AsrModelDefinition("BM", "continuous", "flat", ("flat",)),
    AsrModelDefinition("BMS", "continuous", "flat", ("flat",)),
    AsrModelDefinition("EB", "continuous", "flat", ("flat",)),
    AsrModelDefinition("BM-DRIFT", "continuous", "flat", ("flat",)),
    AsrModelDefinition("MV-BM", "continuous", "flat", ("flat",)),
    AsrModelDefinition("OU", "continuous", "stationary", ("stationary",)),
    AsrModelDefinition("OUM", "continuous", "stationary", ("stationary",)),
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
