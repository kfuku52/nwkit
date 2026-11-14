"""Content-checked, atomic restart records for WGD count bootstraps."""

import hashlib
import json
import math
import os
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from nwkit.output_transaction import output_lock_path, write_text_output
from nwkit.shift_native_provenance import native_implementation_sha256
from nwkit.util import acquire_exclusive_lock


def _digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def checkpoint_inputs(args):
    inputs: tuple[str, ...] = ("infile", "counts", "detection_tsv", "rate_groups_tsv")
    if getattr(args, "fixed_events_tsv", None):
        inputs += ("fixed_events_tsv",)
    if any(getattr(args, key) == "-" for key in inputs):
        raise ValueError("WGD bootstrap checkpoint requires file inputs, not STDIN.")
    return {
        key: _digest(path) if (path := getattr(args, key)) else None for key in inputs
    }


def checkpoint_identity(args, result, inputs):
    import numpy as np

    from nwkit.wgd_count import _fit_json

    if checkpoint_inputs(args) != inputs:
        raise ValueError("WGD bootstrap inputs changed during the observed fit.")
    settings = (
        "format",
        "quoted_node_names",
        "input_rooted",
        "missing_values",
        "rate_model",
        "candidate_branches",
        "event_fractions",
        "multiplicity",
        "family_gamma_shape",
        "family_rate_categories",
        "bootstrap",
        "seed",
        "max_states",
        "state_tolerance",
        "max_iterations",
        "ascertainment",
    )
    identity = {
        "inputs": inputs,
        "settings": {key: getattr(args, key) for key in settings},
        "background": _fit_json(result.background),
        "searched_nodes_in_order": [candidate.node for candidate in result.candidates],
        "implementation_sha256": native_implementation_sha256(),
        "numeric_configuration": {
            name: np.__config__.get_info(name)
            for name in ("blas_opt_info", "lapack_opt_info")
        }
        if hasattr(np.__config__, "get_info")
        else np.show_config(mode="dicts"),
        "numeric_environment": {
            name: os.environ.get(name)
            for name in (
                "OPENBLAS_CORETYPE",
                "NPY_DISABLE_CPU_FEATURES",
                "MKL_CBWR",
                "BLIS_NUM_THREADS",
            )
        },
    }
    # Match the JSON representation on disk even when parser values are tuples.
    return json.loads(json.dumps(identity, allow_nan=False, sort_keys=True))


@dataclass
class BootstrapCheckpoint:
    path: Path
    identity: dict
    statistics: list[float | None]

    def save(self, values):
        payload = {
            "schema_version": 1,
            "identity": self.identity,
            "statistics": values,
            "statistics_sha256": _statistics_digest(values),
        }
        write_text_output(
            self.path,
            lambda stream: json.dump(payload, stream, allow_nan=False, sort_keys=True),
        )


def _statistics_digest(values):
    return hashlib.sha256(
        json.dumps(values, allow_nan=False, separators=(",", ":")).encode()
    ).hexdigest()


def _compatible_draw_count(saved_identity, identity, draws):
    if not isinstance(saved_identity, dict):
        return None
    previous = saved_identity.get("settings")
    current = identity.get("settings")
    if not isinstance(previous, dict) or not isinstance(current, dict):
        return None
    saved_draws = previous.get("bootstrap")
    if type(saved_draws) is not int or not 1 <= saved_draws <= draws:
        return None
    updated = {**saved_identity, "settings": {**previous, "bootstrap": draws}}
    return saved_draws if updated == identity else None


@contextmanager
def open_checkpoint(path, identity, draws):
    if path is None:
        yield None
        return
    path = Path(path).resolve()
    lock = output_lock_path(str(path)) + ".bootstrap"
    with acquire_exclusive_lock(lock, lock_label="WGD bootstrap checkpoint"):
        statistics: list[float | None] = [None] * draws
        if path.exists():
            with path.open(encoding="utf-8") as stream:
                payload = json.load(stream)
            if not isinstance(payload, dict) or payload.get("schema_version") != 1:
                raise ValueError("WGD bootstrap checkpoint has an unsupported schema.")
            saved_draws = _compatible_draw_count(
                payload.get("identity"), identity, draws
            )
            if saved_draws is None:
                raise ValueError(
                    "WGD bootstrap checkpoint does not match the inputs, search, or runtime."
                )
            saved = payload.get("statistics")
            if (
                not isinstance(saved, list)
                or len(saved) != saved_draws
                or any(
                    x is not None
                    and (type(x) not in (float, int) or not math.isfinite(x) or x < 0)
                    for x in saved
                )
            ):
                raise ValueError(
                    "WGD bootstrap checkpoint contains invalid statistics."
                )
            statistics = saved
            if payload.get("statistics_sha256") != _statistics_digest(statistics):
                raise ValueError(
                    "WGD bootstrap checkpoint statistics checksum does not match."
                )
            # Simulation consumes the same seeded prefix regardless of how many
            # draws were requested. Existing holes still need to be refitted.
            statistics = [*statistics, *([None] * (draws - saved_draws))]
        checkpoint = BootstrapCheckpoint(path, identity, statistics)
        if not path.exists():
            checkpoint.save(statistics)
        yield checkpoint
