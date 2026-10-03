"""Research-only conditional integration of full hidden locus histories.

No CLI integration/default is changed. Large histories may use a predeclared
genealogy draw, but are never truncated, rejected or retried to fit a work cap.
"""

import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, product
from typing import Any

import numpy as np

from nwkit.mul_locus import (
    LocusCoalescent,
    detected_signature,
    finite_number,
    sample_locus_tree,
    signature_size,
    transition,
)
from nwkit.mul_locus_mc import LocusBank, parameter_record

OVERFLOW = ("overflow",)


class WorkBudget:
    def __init__(self, limit):
        if type(limit) is not int or limit < 1:
            raise ValueError("Conditional integration work cap must be positive.")
        self.limit, self.used = limit, 0

    def charge(self, count=1):
        self.used += count
        if self.used > self.limit:
            raise ValueError(
                "Conditional integration exceeded work cap; no history discarded."
            )


def add_log(distribution, key, value):
    distribution[key] = float(np.logaddexp(distribution.get(key, -math.inf), value))


def merger_levels(initial, budget):
    """Random-pair merger law, independent of the branch waiting times."""
    level = {initial: 0.0}
    yield level
    while len(next(iter(level))) > 1:
        following: dict[tuple, float] = {}
        for forest, mass in level.items():
            pairs = len(forest) * (len(forest) - 1) // 2
            for a, b in combinations(range(len(forest)), 2):
                budget.charge()
                merged = ("node", *sorted((forest[a], forest[b])))
                state = tuple(
                    sorted(
                        (*[g for i, g in enumerate(forest) if i not in (a, b)], merged)
                    )
                )
                add_log(following, state, mass - math.log(pairs))
        level = following
        yield level


def advance_forest(initial, time, daughter, budget, *, ancestral=False):
    if not initial:
        return {(): 0.0}
    start = len(initial)
    result = {}
    for level in merger_levels(initial, budget):
        end = len(next(iter(level)))
        if (ancestral or daughter) and end != 1:
            continue
        budget.charge((start - end + 1) ** 2)
        log_probability = 0.0 if ancestral else transition(start, end, time)
        if math.isfinite(log_probability):
            result.update(
                {forest: mass + log_probability for forest, mass in level.items()}
            )
        if time == 0 and not ancestral:
            break
    return result


def genealogy_distribution(locus, ne, *, max_states=100000):
    """Exact all-pair forest DP, retaining undetected tips in daughter bounds."""
    return _genealogy_distribution(locus, ne, WorkBudget(max_states))


def _genealogy_distribution(locus, ne, budget):
    counts = LocusCoalescent(locus, ne, max_states=budget.limit)
    budget.charge(counts.work)
    forests: dict[Any, dict[tuple, float]] = {}
    tip = 0
    for node in counts.order:
        current: dict[tuple, float] = {(): 0.0}
        if node.kind == "tip":
            current = {(("tip", node.species, tip),): 0.0}
            tip += 1
        for child in node.children:
            joined: dict[tuple, float] = {}
            for (first, p), (second, q) in product(
                current.items(), forests.pop(child).items()
            ):
                budget.charge()
                add_log(joined, tuple(sorted(first + second)), p + q)
            current = joined
        following: dict[tuple, float] = {}
        for initial, mass in current.items():
            for state, log_probability in advance_forest(
                initial,
                counts.times[node],
                node.daughter,
                budget,
                ancestral=node is locus,
            ).items():
                budget.charge()
                add_log(following, state, mass + log_probability)
        normalizer = counts.bound_log_normalizers.get(node, 0.0)
        forests[node] = {state: mass - normalizer for state, mass in following.items()}
    result = {
        state[0] if state else None: math.exp(mass)
        for state, mass in forests[locus].items()
    }
    check_distribution(result)
    return result


def check_distribution(distribution):
    if any(
        not finite_number(p) or p < 0 for p in distribution.values()
    ) or not math.isclose(
        math.fsum(distribution.values()), 1.0, rel_tol=2e-11, abs_tol=2e-13
    ):
        raise ArithmeticError("Conditional distribution failed mass validation.")


def detection_distribution(
    gene, detection, *, max_states=100000, max_observed_tips=None
):
    """Integrate all detection masks AFTER the full genealogy has been formed."""
    if max_observed_tips is not None and (
        type(max_observed_tips) is not int or max_observed_tips < 1
    ):
        raise ValueError("Detection observation-size bound must be positive.")
    return _detection_distribution(
        gene, detection, WorkBudget(max_states), max_observed_tips
    )


def _detection_distribution(gene, detection, budget, max_tips=None):
    if any(not finite_number(p) or not 0 <= p <= 1 for p in detection.values()):
        raise ValueError("Detection probabilities must be in [0,1].")
    values: dict[Any, dict[Any, float]] = {}
    pending = [(gene, False)]
    while pending:
        node, visited = pending.pop()
        budget.charge()
        if node is None:
            values[node] = {None: 1.0}
        elif node[0] == "tip":
            p = detection[node[1]]
            values[node] = {
                key: mass
                for key, mass in ((None, 1 - p), (("tip", node[1]), p))
                if mass > 0
            }
        elif not visited:
            pending.append((node, True))
            pending.extend((child, False) for child in node[1:])
        else:
            left, right = (values.pop(child) for child in node[1:])
            joined: dict[Any, float] = defaultdict(float)
            for (first, p), (second, q) in product(left.items(), right.items()):
                budget.charge()
                key = (
                    (second if first is None else first)
                    if first is None or second is None
                    else ("node", *sorted((first, second)))
                )
                if (
                    first == OVERFLOW
                    or second == OVERFLOW
                    or (max_tips is not None and signature_size(key) > max_tips)
                ):
                    key = OVERFLOW
                joined[key] += p * q
            values[node] = dict(joined)
    result = {key: p for key, p in values[gene].items() if p > 0}
    check_distribution(result)
    return result


def conditional_observations(locus, ne, detection, *, max_states=100000):
    result: dict[Any, float] = defaultdict(float)
    budget = WorkBudget(max_states)
    for gene, mass in _genealogy_distribution(locus, ne, budget).items():
        for signature, p in _detection_distribution(gene, detection, budget).items():
            budget.charge()
            result[signature] += mass * p
    check_distribution(result)
    return dict(result)


@dataclass
class Moments:
    n: int = 0
    mean: float = 0.0
    m2: float = 0.0

    def add(self, value):
        if not finite_number(value) or not 0 <= value <= 1:
            raise ValueError("Conditional probability must be in [0,1].")
        self.n += 1
        delta = value - self.mean
        self.mean += delta / self.n
        self.m2 += delta * (value - self.mean)

    def with_zeros(self, total):
        if type(total) is not int or total < self.n:
            raise ValueError("Moment sample count cannot lose observations.")
        if not total:
            return Moments()
        return Moments(
            total,
            self.mean * self.n / total,
            self.m2 + self.mean**2 * self.n * (total - self.n) / total,
        )


def bounded_interval(moments, alpha):
    """Two-sided Maurer-Pontil (2009), Theorem 4, with log(4/alpha)."""
    if not finite_number(alpha) or not 0 < alpha < 1:
        raise ValueError("Bounded-mean error budget must be in (0,1).")
    if type(moments.n) is not int or moments.n < 2:
        raise ValueError("Bounded-mean interval requires at least two IID draws.")
    if (
        not finite_number(moments.mean)
        or not 0 <= moments.mean <= 1
        or not finite_number(moments.m2)
        or moments.m2 < 0
    ):
        raise ArithmeticError("Invalid bounded-variable moments.")
    logarithm = math.log(4) - math.log(alpha)
    radius = math.sqrt(
        2 * moments.m2 / (moments.n - 1) * logarithm / moments.n
    ) + 7 * logarithm / (3 * (moments.n - 1))
    return max(0.0, moments.mean - radius), min(1.0, moments.mean + radius)


def chernoff_interval(moments, alpha):
    """Invert the two-sided unit-interval Chernoff bound (Foong et al., 2022)."""
    bounded_interval(moments, alpha)  # Validate the same IID bounded moments.
    mean = moments.mean
    threshold = (math.log(2) - math.log(alpha)) / moments.n

    def divergence(p):
        if p in (0, 1):
            return 0.0 if p == mean else math.inf
        return math.fsum(
            (
                mean * (math.log(mean) - math.log(p)) if mean else 0.0,
                (1 - mean) * (math.log1p(-mean) - math.log1p(-p)) if mean < 1 else 0.0,
            )
        )

    # Keep the outer bracket: rounding must widen, not shrink, coverage.
    low, inner = 0.0, mean
    for _ in range(1075):
        mid = (low + inner) / 2
        if mid in (low, inner):
            break
        if divergence(mid) > threshold:
            low = mid
        else:
            inner = mid
    inner, high = mean, 1.0
    for _ in range(1075):
        mid = (inner + high) / 2
        if mid in (inner, high):
            break
        if divergence(mid) > threshold:
            high = mid
        else:
            inner = mid
    return low, high


@dataclass
class IntegratedBank:
    candidate: int
    h2: str
    grid: int
    population: object
    parameters: object
    strata: list
    method: str
    interval_method: str = "empirical-bernstein"

    @property
    def samples(self):
        return sum(stratum["samples"] for stratum in self.strata)

    @property
    def attempts(self):
        return self.samples


def integrated_probability(bank, signature, alpha):
    interval = (
        chernoff_interval if bank.interval_method == "chernoff-kl" else bounded_interval
    )
    if bank.interval_method not in ("chernoff-kl", "empirical-bernstein"):
        raise ValueError("Unknown bounded-contribution interval method.")
    numerator = denominator = n_low = n_high = d_low = d_high = 0.0
    for stratum in bank.strata:
        weight = stratum["weight"]
        pattern = (
            stratum["patterns"].get(signature, Moments()).with_zeros(stratum["samples"])
        )
        selected = stratum["selected"]
        lo, hi = interval(pattern, alpha)
        slo, shi = interval(selected, alpha)
        numerator += weight * pattern.mean
        denominator += weight * selected.mean
        n_low += weight * lo
        n_high += weight * hi
        d_low += weight * slo
        d_high += weight * shi
    if denominator <= 0:
        raise ValueError("Conditional integration has no selected mass.")
    return (
        numerator / denominator,
        n_low / d_high,
        min(1.0, n_high / d_low) if d_low else 1.0,
    )


def score_integrated_bank(bank, observations, alpha):
    scores, lower, upper, rows = [], [], [], []
    for signature, multiplicity in Counter(observations).items():
        estimate, lo, hi = integrated_probability(bank, signature, alpha)
        scores.append(multiplicity * math.log(estimate) if estimate else -math.inf)
        lower.append(multiplicity * math.log(lo) if lo else -math.inf)
        upper.append(multiplicity * math.log(hi))
        rows.append(
            {
                "signature": signature,
                "multiplicity": multiplicity,
                "probability_estimate": estimate,
                "probability_lower": lo,
                "probability_upper": hi,
            }
        )
    return {
        "mul.tree": bank.candidate,
        "h2.node": bank.h2,
        "grid": bank.grid,
        "parameters": parameter_record(bank),
        "log_likelihood": math.fsum(scores),
        "mc_lower": math.fsum(lower),
        "mc_upper": math.fsum(upper),
        "patterns": rows,
    }


def accumulate(stratum, distribution, max_tips):
    check_distribution(distribution)
    selected = {
        key: p
        for key, p in distribution.items()
        if key != OVERFLOW and 2 <= signature_size(key) <= max_tips
    }
    mass = math.fsum(selected.values())
    if mass > 1 and mass <= 1 + 2e-13:
        mass = 1.0
    stratum["selected"].add(mass)
    for signature, probability in selected.items():
        stratum["patterns"].setdefault(signature, Moments()).add(min(1.0, probability))


def build_paired_banks(task, model, *, exact_tip_limit=4, max_states=100000):
    """Same IID DL histories/gene draws for histogram, detection and hybrid RB."""
    if type(exact_tip_limit) is not int or exact_tip_limit < 0:
        raise ValueError("Exact hidden-tip limit must be a nonnegative integer.")
    WorkBudget(max_states)
    candidate, h2, grid, population, point = task
    rate, stem = point.duplication + point.loss, model["ancestral_stem"]
    carry = math.exp(-rate * stem)
    birth = point.duplication / rate * -math.expm1(-rate * stem) if rate else 0.0
    if carry == 0 or (point.duplication > 0 and stem > 0 and birth == 0):
        raise ArithmeticError("Conditional integration stratum prior underflowed.")
    active = [("none", carry)] + ([("birth", birth)] if birth else [])
    if model["samples"] < 2 * len(active) or model["samples"] > model["max_attempts"]:
        raise ValueError(
            "Conditional budget requires two draws per active stratum and obeys attempt cap."
        )
    histogram, detection, hybrid = [], [], []
    counts: Counter[tuple] = Counter()
    exact_histories = sampled_histories = 0
    for i, (condition, weight) in enumerate(active):
        samples = model["samples"] // len(active) + (i < model["samples"] % len(active))
        raw = {
            "condition": condition,
            "weight": weight,
            "samples": samples,
            "counts": Counter(),
            "selected": 0,
        }
        rb = [
            {
                "condition": condition,
                "weight": weight,
                "samples": samples,
                "patterns": {},
                "selected": Moments(),
            }
            for _ in range(2)
        ]
        history_rng = np.random.default_rng(
            np.random.SeedSequence([model["seed"], 3, candidate, grid, i])
        )
        for draw in range(samples):
            locus = sample_locus_tree(
                population,
                point,
                stem,
                history_rng,
                max_nodes=model["max_locus_nodes"],
                stem_condition=condition,
            )
            gene_rng = np.random.default_rng(
                np.random.SeedSequence([model["seed"], 4, candidate, grid, i, draw])
            )
            sampler = LocusCoalescent(
                locus, point.ne, max_states=model["max_coalescent_states"]
            )
            gene = sampler.sample(gene_rng)
            detect_rng = np.random.default_rng(
                np.random.SeedSequence([model["seed"], 5, candidate, grid, i, draw])
            )
            observed = detected_signature(gene, model["detection"], detect_rng)
            if 2 <= signature_size(observed) <= model["max_observed_tips"]:
                raw["counts"][observed] += 1
                raw["selected"] += 1
            probabilities = detection_distribution(
                gene,
                model["detection"],
                max_states=max_states,
                max_observed_tips=model["max_observed_tips"],
            )
            accumulate(rb[0], probabilities, model["max_observed_tips"])
            hidden_tips = sum(node.kind == "tip" for node in sampler.order)
            if hidden_tips <= exact_tip_limit:
                integrated = conditional_observations(
                    locus, point.ne, model["detection"], max_states=max_states
                )
                exact_histories += 1
            else:
                integrated = probabilities
                sampled_histories += 1
            accumulate(rb[1], integrated, model["max_observed_tips"])
        histogram.append(raw)
        detection.append(rb[0])
        hybrid.append(rb[1])
        counts.update(raw["counts"])
    banks = {
        "histogram": LocusBank(
            candidate,
            h2,
            grid,
            population,
            point,
            counts,
            model["samples"],
            model["samples"],
            histogram,
        ),
        "detection": IntegratedBank(
            candidate, h2, grid, population, point, detection, "detection-RB"
        ),
        "hybrid": IntegratedBank(
            candidate,
            h2,
            grid,
            population,
            point,
            hybrid,
            "predeclared-small-history-genealogy-detection-RB",
        ),
    }
    return banks, {
        "exact_histories": exact_histories,
        "sampled_histories": sampled_histories,
    }


def build_integrated_bank(task, model):
    method = model["integration"]
    if method not in ("detection-rb", "hybrid-rb"):
        raise ValueError("Unknown conditional integration method.")
    banks, _ = build_paired_banks(
        task,
        model,
        exact_tip_limit=4 if method == "hybrid-rb" else 0,
        max_states=model["max_coalescent_states"],
    )
    bank = banks["hybrid" if method == "hybrid-rb" else "detection"]
    bank.interval_method = "chernoff-kl"
    return bank
