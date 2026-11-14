"""Observation mechanisms for censored regression and its bootstrap."""

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from nwkit.model_matrix import ReplicatedObservation, parse_key_values


def expand_observations(values, template, *, default=float("nan")):
    """Expand tip auxiliaries in the same biological-replicate order as y."""
    if values is not None and len(values) != len(template):
        raise ValueError("Observation mechanism must specify every tree tip.")
    expanded: list[float] = []
    for index, original in enumerate(template):
        count = (
            len(original.values) if isinstance(original, ReplicatedObservation) else 1
        )
        value = default if values is None else values[index]
        if isinstance(value, ReplicatedObservation):
            if len(value.values) != count:
                raise ValueError("Observation mechanism replicates are misaligned.")
            expanded.extend(value.values)
        else:
            expanded.extend([value] * count)
    return np.asarray(expanded, dtype=float)


def pack_observations(values, template):
    """Preserve biological replication when returning simulated observations."""
    result = []
    position = 0
    for original in template:
        replicated = isinstance(original, ReplicatedObservation)
        count = len(original.values) if replicated else 1
        selected = tuple(float(value) for value in values[position : position + count])
        result.append(ReplicatedObservation(selected) if replicated else selected[0])
        position += count
    if position != len(values):
        raise ValueError("Generated observations are misaligned with the template.")
    return result


@dataclass(frozen=True)
class CensoringModel:
    """Known, non-informative coarsening rule on the entire sample space.

    Detection limits apply to every observation, including originally exact
    values. Missing limits mean no limit on that side. Interval bins contain
    finite internal cutpoints; the two exterior intervals extend to infinity.
    Observed censor bounds are separate data, never a substitute for this rule.
    """

    kind: str
    detection_lower: Sequence | None = None
    detection_upper: Sequence | None = None
    cutpoints: tuple[float, ...] = ()

    def limits(self, template):
        if self.kind not in {"uncensored", "detection-limits", "interval-bins"}:
            raise ValueError("Unknown censored-Gaussian observation model.")
        lower = expand_observations(self.detection_lower, template)
        upper = expand_observations(self.detection_upper, template)
        if np.isinf(lower).any() or np.isinf(upper).any():
            raise ValueError("Detection limits must be finite or missing.")
        if self.kind != "detection-limits" and (
            self.detection_lower is not None or self.detection_upper is not None
        ):
            raise ValueError("Detection columns require detection-limits mode.")
        if self.kind == "detection-limits":
            if self.detection_lower is None and self.detection_upper is None:
                raise ValueError("Detection-limits mode requires limit columns.")
            if np.any(lower >= upper):
                raise ValueError("Lower detection limits must be below upper limits.")
        cuts = np.asarray(self.cutpoints, dtype=float)
        if self.kind == "interval-bins":
            if (
                cuts.ndim != 1
                or not len(cuts)
                or not np.isfinite(cuts).all()
                or np.any(np.diff(cuts) <= 0)
            ):
                raise ValueError(
                    "Interval cutpoints must be finite and strictly increasing."
                )
        elif len(cuts):
            raise ValueError("Cutpoints require interval-bins mode.")
        return lower, upper

    def observe(self, latent_values, template):
        """Return fresh response values and fresh observed lower/upper bounds."""
        lower_limit, upper_limit = self.limits(template)
        latent_values = np.asarray(latent_values, dtype=float)
        if (
            latent_values.shape != lower_limit.shape
            or not np.isfinite(latent_values).all()
        ):
            raise ValueError("Latent responses must be finite and aligned.")
        values = latent_values.copy()
        lower = np.full(len(values), np.nan)
        upper = np.full(len(values), np.nan)
        if self.kind == "detection-limits":
            left = latent_values <= lower_limit
            right = latent_values >= upper_limit
            upper[left] = lower_limit[left]
            lower[right] = upper_limit[right]
            values[left | right] = np.nan
        elif self.kind == "interval-bins":
            cuts = np.asarray(self.cutpoints)
            indices = np.searchsorted(cuts, latent_values, side="left")
            lower = np.r_[np.nan, cuts][indices]
            upper = np.r_[cuts, np.nan][indices]
            values[:] = np.nan
        return tuple(
            pack_observations(array, template) for array in (values, lower, upper)
        )

    def validate_observed(self, values, lower, upper):
        """Check that the supplied observed data could follow this mechanism."""
        low_limit, high_limit = self.limits(values)
        observed = expand_observations(values, values)
        low = expand_observations(lower, values)
        high = expand_observations(upper, values)
        exact = np.isnan(low) & np.isnan(high)
        if self.kind == "uncensored":
            valid = exact & np.isfinite(observed)
        elif self.kind == "detection-limits":
            valid_exact = exact & np.isfinite(observed)
            valid_exact &= np.isnan(low_limit) | (observed > low_limit)
            valid_exact &= np.isnan(high_limit) | (observed < high_limit)
            valid_left = np.isnan(low) & (high == low_limit) & np.isnan(observed)
            valid_right = (low == high_limit) & np.isnan(high) & np.isnan(observed)
            valid = valid_exact | valid_left | valid_right
        else:
            edges = np.r_[-np.inf, self.cutpoints, np.inf]
            low_numeric = np.where(np.isnan(low), -np.inf, low)
            high_numeric = np.where(np.isnan(high), np.inf, high)
            valid = np.zeros(len(observed), dtype=bool)
            for first, second in zip(edges[:-1], edges[1:], strict=True):
                valid |= (low_numeric == first) & (high_numeric == second)
            valid &= np.isnan(observed)
        if not valid.all():
            raise ValueError(
                "Observed censoring data are inconsistent with the observation model."
            )


def read_censoring_models(args, response_specs, leaf_names, read_auxiliary):
    """Build CLI mechanisms with a caller-provided, replicate-aware reader."""
    names = parse_key_values(
        getattr(args, "response_observation_model", None),
        "--response-observation-model",
    )
    lower_columns = parse_key_values(
        getattr(args, "response_detection_lower", None), "--response-detection-lower"
    )
    upper_columns = parse_key_values(
        getattr(args, "response_detection_upper", None), "--response-detection-upper"
    )
    bins = parse_key_values(
        getattr(args, "response_observation_bins", None), "--response-observation-bins"
    )
    selected = set(names) | set(lower_columns) | set(upper_columns) | set(bins)
    for response in selected:
        if (
            response not in response_specs
            or response_specs[response].family != "censored-gaussian"
        ):
            raise ValueError(
                "Observation models apply only to censored-Gaussian responses."
            )
        if response not in names:
            raise ValueError(
                "Detection limits and bins require --response-observation-model."
            )
    auxiliary = read_auxiliary(
        sorted(set(lower_columns.values()) | set(upper_columns.values()))
    )
    result = {}
    for response, kind in names.items():
        lower = (
            None
            if response not in lower_columns
            else [auxiliary[lower_columns[response]][name] for name in leaf_names]
        )
        upper = (
            None
            if response not in upper_columns
            else [auxiliary[upper_columns[response]][name] for name in leaf_names]
        )
        try:
            cutpoints = (
                tuple(float(value) for value in bins[response].split("|"))
                if response in bins
                else ()
            )
        except ValueError as exc:
            raise ValueError(
                "Observation bins must be separated numeric cutpoints."
            ) from exc
        result[response] = CensoringModel(kind, lower, upper, cutpoints)
    return result
