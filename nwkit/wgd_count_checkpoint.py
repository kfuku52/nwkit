"""Content-checked, atomic restart records for WGD count bootstraps."""

import hashlib
import json
import math
import os
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from nwkit.output_transaction import output_lock_path, write_text_output
from nwkit.shift_native_provenance import native_implementation_sha256
from nwkit.util import acquire_exclusive_lock
from nwkit.wgd_count_fit import CountFit, CountScan, ScanCandidate
from nwkit.wgd_count_model import MultiplicationEvent


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
        raise ValueError("WGD checkpoint requires file inputs, not STDIN.")
    return {
        key: _digest(path) if (path := getattr(args, key)) else None for key in inputs
    }


_OBSERVED_SETTINGS = (
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
    "max_states",
    "state_tolerance",
    "max_iterations",
    "ascertainment",
)


def observed_identity(args, inputs):
    identity = {
        "inputs": inputs,
        "settings": {key: getattr(args, key) for key in _OBSERVED_SETTINGS},
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

    return json.loads(json.dumps(identity, allow_nan=False, sort_keys=True))


def checkpoint_identity(args, result, inputs):
    from nwkit.wgd_count import _fit_json

    if checkpoint_inputs(args) != inputs:
        raise ValueError("WGD bootstrap inputs changed during the observed fit.")
    identity = observed_identity(args, inputs)
    identity["settings"].update(bootstrap=args.bootstrap, seed=args.seed)
    identity["background"] = _fit_json(result.background)
    identity["searched_nodes_in_order"] = [
        candidate.node for candidate in result.candidates
    ]
    return identity


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


def _fit_payload(fit):
    return {**asdict(fit), "rates": fit.rates.tolist()}


def _restore_fit(values):
    fit = CountFit(
        **{
            **values,
            "rates": np.asarray(values["rates"], dtype=float),
            "event": MultiplicationEvent(**values["event"])
            if values["event"]
            else None,
        }
    )
    if not fit.converged or not np.all(np.isfinite(fit.rates)):
        raise ValueError("WGD fit checkpoint contains an invalid fit.")
    return fit


def _encode_fit_record(record):
    return json.dumps(record, allow_nan=False, sort_keys=True, separators=(",", ":"))


class ObservedCheckpoint:
    """One atomic record for the background, candidates and unfinished branches."""

    def __init__(self, path, identity):
        self.path = path
        self.identity = identity
        self.background = None
        self.candidates = None
        self.complete = False
        self.components: dict[int, dict[int, CountFit]] = {}
        if path.exists():
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                data = payload["data"]
                if (
                    data["schema_version"] != 1
                    or data["identity"] != identity
                    or payload["sha256"] != self._checksum(data)
                ):
                    raise ValueError(
                        "WGD fit checkpoint does not match the inputs, search, or runtime."
                    )
                self.background = _restore_fit(data["background"])
                self.complete = data["complete"]
                if type(self.complete) is not bool or (
                    self.complete and not data["candidates"]
                ):
                    raise ValueError(
                        "WGD fit checkpoint has an invalid completion state."
                    )
                if data["candidates"] is not None:
                    self.candidates = tuple(
                        ScanCandidate(
                            **{
                                **values,
                                "event_fit": _restore_fit(values["event_fit"]),
                                "burst_fit": _restore_fit(values["burst_fit"]),
                            }
                        )
                        for values in data["candidates"]
                    )
                self.components = {
                    int(node): {
                        int(index): _restore_fit(fit) for index, fit in fits.items()
                    }
                    for node, fits in data.get("components", {}).items()
                }
                if self.complete and self.components:
                    raise ValueError(
                        "A completed WGD fit checkpoint has unfinished components."
                    )
            except (AttributeError, KeyError, TypeError) as exc:
                raise ValueError("Malformed WGD fit checkpoint.") from exc

    @staticmethod
    def _checksum(data):
        return hashlib.sha256(_encode_fit_record(data).encode()).hexdigest()

    def save(self, background, candidates=None, *, complete=False, components=None):
        if components is None:
            components = self.components
        rows = (
            None
            if candidates is None
            else [
                {
                    "node": row.node,
                    "event_fit": _fit_payload(row.event_fit),
                    "burst_fit": _fit_payload(row.burst_fit),
                    "improvement": row.improvement,
                    "burst_aic_difference": row.burst_aic_difference,
                }
                for row in candidates
            ]
        )
        data = {
            "schema_version": 1,
            "identity": self.identity,
            "background": _fit_payload(background),
            "candidates": rows,
            "complete": complete,
        }
        if components:
            data["components"] = {
                str(node): {
                    str(index): _fit_payload(fit) for index, fit in fits.items()
                }
                for node, fits in components.items()
            }
        payload = {"data": data, "sha256": self._checksum(data)}
        write_text_output(
            self.path, lambda stream: stream.write(_encode_fit_record(payload))
        )
        self.background = background
        self.candidates = candidates
        self.complete = complete
        self.components = components

    def save_component(self, node, index, fit):
        self.save(
            self.background,
            self.candidates,
            components={
                **self.components,
                node: {**self.components.get(node, {}), index: fit},
            },
        )

    def save_candidate(self, candidate):
        self.save(
            self.background,
            (*self.candidates, candidate) if self.candidates else (candidate,),
            components={
                node: fits
                for node, fits in self.components.items()
                if node != candidate.node
            },
        )

    def scan(self, fractions, multiplicity):
        if not self.complete or self.background is None or self.candidates is None:
            return None
        return CountScan(
            self.background,
            self.candidates,
            fractions=fractions,
            multiplicity=multiplicity,
        )


@contextmanager
def open_observed_checkpoint(path, args, inputs):
    path = Path(path).resolve()
    lock = output_lock_path(str(path)) + ".observed"
    with acquire_exclusive_lock(lock, lock_label="WGD fit checkpoint"):
        yield ObservedCheckpoint(path, observed_identity(args, inputs))
