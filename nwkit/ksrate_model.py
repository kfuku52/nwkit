"""Family-aware focal-lineage Ks correction on a rooted species tree."""

from dataclasses import dataclass

import numpy as np
from scipy.stats import binom


@dataclass(frozen=True)
class KsObservation:
    species_a: str
    species_b: str
    family_id: str
    ks: float


@dataclass(frozen=True)
class KsTrio:
    focal: str
    sister: str
    outgroup: str
    species_event_id: str
    branch_id: int


def corrected_ks(focal_sister, focal_outgroup, sister_outgroup):
    """Twice the focal path from the focal/sister ancestor to the focal tip."""
    values = np.asarray([focal_sister, focal_outgroup, sister_outgroup], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Pairwise Ks distances must be finite and nonnegative.")
    first, second, third = (float(value) for value in values)
    correction = (
        (first - third) + second if first >= third else second - (third - first)
    )
    if not np.isfinite(correction):
        raise ValueError("Corrected Ks must be representable as a finite float.")
    return correction


def _midpoint(low, high):
    low, high = float(low), float(high)
    return (low + high) / 2 if low < 0 < high else low + (high - low) / 2


def _stable_median(values):
    values = np.sort(np.asarray(values, dtype=float))
    if not len(values) or np.any(np.isnan(values)):
        return float("nan")
    return _midpoint(values[(len(values) - 1) // 2], values[len(values) // 2])


def _weighted_median(values, weights):
    order = np.argsort(values, kind="stable")
    values = np.asarray(values)[order]
    weights = np.asarray(weights, dtype=np.int64)[order]
    total = int(weights.sum())
    if total == 0:
        return float("nan")
    cumulative = np.cumsum(weights)
    low = int(np.searchsorted(cumulative, (total - 1) // 2, side="right"))
    high = int(np.searchsorted(cumulative, total // 2, side="right"))
    return _midpoint(values[low], values[high])


class PairKs:
    """One independent family observation per unordered species pair."""

    def __init__(self, observations, species):
        species = set(species)
        self.observations = tuple(observations)
        seen = set()
        grouped: dict[tuple[str, str], list[KsObservation]] = {}
        for record in self.observations:
            if record.species_a not in species or record.species_b not in species:
                raise ValueError(
                    "Ks observations contain species absent from the tree."
                )
            if record.species_a == record.species_b:
                raise ValueError(
                    "Ks correction needs between-species ortholog distances."
                )
            if not record.family_id:
                raise ValueError("Ks observations need nonempty family IDs.")
            if not np.isfinite(record.ks) or record.ks < 0:
                raise ValueError("Ks observations must be finite and nonnegative.")
            pair = tuple(sorted((record.species_a, record.species_b)))
            key = (*pair, record.family_id)
            if key in seen:
                raise ValueError(
                    "Ks input repeats a family for one species pair; select one representative first."
                )
            seen.add(key)
            grouped.setdefault(pair, []).append(record)
        if not self.observations:
            raise ValueError("Ks correction needs at least one observation.")
        self.families = tuple(
            sorted({record.family_id for record in self.observations})
        )
        family_index = {family: index for index, family in enumerate(self.families)}
        self.pairs = {
            pair: (
                np.array([record.ks for record in rows]),
                np.array(
                    [family_index[record.family_id] for record in rows], dtype=int
                ),
            )
            for pair, rows in grouped.items()
        }

    def estimates(self, family_weights=None):
        if family_weights is None:
            family_weights = np.ones(len(self.families), dtype=int)
        weights = np.asarray(family_weights)
        if (
            weights.shape != (len(self.families),)
            or np.any(~np.isfinite(weights))
            or np.any(weights < 0)
            or np.any(weights != np.floor(weights))
        ):
            raise ValueError(
                "Family weights must be finite nonnegative integers with one value per family."
            )
        if sum(int(weight) for weight in weights) > np.iinfo(np.int64).max:
            raise ValueError("Family weights exceed the integer accumulation limit.")
        return {
            pair: _weighted_median(values, weights[indices])
            for pair, (values, indices) in self.pairs.items()
        }

    def count(self, first, second):
        values = self.pairs.get(tuple(sorted((first, second))))
        return 0 if values is None else len(values[0])


def trio_values(trios, estimates):
    values = np.full(len(trios), np.nan)
    for index, trio in enumerate(trios):
        distances = [
            estimates.get(tuple(sorted(pair)), np.nan)
            for pair in (
                (trio.focal, trio.sister),
                (trio.focal, trio.outgroup),
                (trio.sister, trio.outgroup),
            )
        ]
        if np.all(np.isfinite(distances)):
            values[index] = corrected_ks(*distances)
    return values


def grouped_corrections(trios, values):
    grouped: dict[tuple[str, str], list[int]] = {}
    for index, trio in enumerate(trios):
        grouped.setdefault((trio.focal, trio.species_event_id), []).append(index)
    result = {}
    for key, indices in grouped.items():
        selected = values[indices]
        available = selected[np.isfinite(selected)]
        result[key] = _stable_median(available)
    return result


def bootstrap_corrections(pairs, trios, draws, seed):
    """Resample shared family IDs jointly, preserving across-pair dependence."""
    if not isinstance(draws, int) or draws < 0 or not isinstance(seed, int) or seed < 0:
        raise ValueError("Ks bootstrap count and seed must be nonnegative integers.")
    keys = tuple(grouped_corrections(trios, np.zeros(len(trios))))
    results = {key: np.full(draws, np.nan) for key in keys}
    rng = np.random.default_rng(seed)
    num_families = len(pairs.families)
    for draw in range(draws):
        selected = rng.integers(0, num_families, size=num_families)
        weights = np.bincount(selected, minlength=num_families)
        corrections = grouped_corrections(
            trios, trio_values(trios, pairs.estimates(weights))
        )
        for key in keys:
            results[key][draw] = corrections[key]
    return results


def median_order_statistic_interval(values, error_probability):
    """Noninterpolated finite-sample interval for an iid population median."""
    values = np.asarray(values, dtype=float)
    if (
        values.ndim != 1
        or not len(values)
        or np.any(~np.isfinite(values))
        or not np.isfinite(error_probability)
        or not 0 < error_probability < 1
    ):
        raise ValueError("Median intervals need finite samples and error in (0,1).")
    values = np.sort(values)
    n = len(values)
    ranks = np.arange(1, n // 2 + 1)
    admissible = ranks[2 * binom.cdf(ranks - 1, n, 0.5) <= error_probability]
    if not len(admissible):
        return float("nan"), float("nan")
    k = int(admissible[-1])
    return float(values[k - 1]), float(values[n - k])


def simultaneous_median_corrections(pairs, trios, ci_level):
    """Bonferroni pair-median bounds propagated through trio and node medians.

    Across-pair dependence is unrestricted. Families must be independent and
    identically distributed within each pair; population selection remains fixed.
    """
    if not np.isfinite(ci_level) or not 0 < ci_level < 1:
        raise ValueError("Simultaneous median confidence level must be in (0,1).")
    required = {
        tuple(sorted(pair))
        for trio in trios
        for pair in (
            (trio.focal, trio.sister),
            (trio.focal, trio.outgroup),
            (trio.sister, trio.outgroup),
        )
    }
    available = sorted(required & set(pairs.pairs))
    intervals = {
        pair: median_order_statistic_interval(
            pairs.pairs[pair][0], (1 - ci_level) / len(available)
        )
        for pair in available
    }
    grouped: dict[tuple[str, str], list[tuple[float, float]]] = {}
    for trio in trios:
        bounds = [
            intervals.get(tuple(sorted(pair)), (np.nan, np.nan))
            for pair in (
                (trio.focal, trio.sister),
                (trio.focal, trio.outgroup),
                (trio.sister, trio.outgroup),
            )
        ]
        key = trio.focal, trio.species_event_id
        if not all(
            tuple(sorted(pair)) in pairs.pairs
            for pair in (
                (trio.focal, trio.sister),
                (trio.focal, trio.outgroup),
                (trio.sister, trio.outgroup),
            )
        ):
            continue
        # Never discard an unbounded complete trio to manufacture a narrow CI.
        if np.all(np.isfinite(bounds)):
            low = corrected_ks(bounds[0][0], bounds[1][0], bounds[2][1])
            high = corrected_ks(bounds[0][1], bounds[1][1], bounds[2][0])
        else:
            low = high = float("nan")
        grouped.setdefault(key, []).append((low, high))
    return {
        key: (
            _stable_median(np.array(values)[:, 0]),
            _stable_median(np.array(values)[:, 1]),
        )
        for key, values in grouped.items()
    }, len(available)
