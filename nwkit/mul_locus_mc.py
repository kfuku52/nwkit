"""Explicit finite-grid locus-model comparison and search-wide null bootstrap."""

import math
from collections import Counter
from dataclasses import asdict, dataclass

import numpy as np
from scipy.stats import beta

from nwkit.mul_locus import (
    LocusParameters,
    finite_number,
    sample_observation,
    sample_selected,
    signature_size,
)
from nwkit.mul_msc_model import dated_candidates


@dataclass
class LocusBank:
    candidate: int
    h2: str
    grid: int
    population: object
    parameters: LocusParameters
    counts: Counter
    samples: int
    attempts: int
    strata: list | None = None


def validate_model(model, species):
    required = {
        "schema",
        "copy_role",
        "root_locus_count",
        "species_time_unit",
        "ancestral_stem",
        "detection",
        "max_observed_tips",
        "parameter_grid",
        "samples",
        "seed",
        "confidence",
        "max_attempts",
        "max_locus_nodes",
        "max_coalescent_states",
    }
    if (
        not isinstance(model, dict)
        or set(model) - {"integration"} != required
        or model["schema"] != "nwkit-mul-locus-mc-model-v1"
    ):
        raise ValueError("Locus model requires exactly the documented v1 keys.")
    if model.get("integration", "selected-histogram") not in (
        "selected-histogram",
        "ancestral-stratified",
        "detection-rb",
        "hybrid-rb",
    ):
        raise ValueError("Unknown locus integration method.")
    if (
        model["copy_role"] != "distinct-loci"
        or type(model["root_locus_count"]) is not int
        or model["root_locus_count"] != 1
        or model["species_time_unit"] != "generations"
    ):
        raise ValueError(
            "Locus v1 requires distinct loci, one ancestral origin locus and generations."
        )
    if not finite_number(model["ancestral_stem"]) or model["ancestral_stem"] < 0:
        raise ValueError("Ancestral DL stem must be finite/nonnegative.")
    if (
        not isinstance(model["detection"], dict)
        or set(model["detection"]) != set(species.leaf_names())
        or any(
            not finite_number(p) or not 0 < p <= 1 for p in model["detection"].values()
        )
    ):
        raise ValueError(
            "Known detection probabilities must cover every species in (0,1]."
        )
    for name in (
        "max_observed_tips",
        "samples",
        "max_attempts",
        "max_locus_nodes",
        "max_coalescent_states",
    ):
        if type(model[name]) is not int or model[name] < (
            2 if name == "max_observed_tips" else 1
        ):
            raise ValueError(f"Locus {name} must be a positive integer (tips >=2).")
    if (
        type(model["seed"]) is not int
        or model["seed"] < 0
        or not finite_number(model["confidence"])
        or not 0 < model["confidence"] < 1
    ):
        raise ValueError("Locus seed must be nonnegative; confidence in (0,1).")
    if not isinstance(model["parameter_grid"], list) or not model["parameter_grid"]:
        raise ValueError("Locus parameter grid must be a nonempty JSON array.")
    parameters = []
    for point in model["parameter_grid"]:
        if not isinstance(point, dict) or set(point) != {
            "duplication",
            "loss",
            "ne",
            "hybridization_age",
        }:
            raise ValueError(
                "Locus grid points require exactly duplication/loss/ne/hybridization_age."
            )
        parameters.append(LocusParameters(**point))
    if len(set(parameters)) != len(parameters):
        raise ValueError("Locus grid points must be unique.")
    return parameters


def make_tasks(species, h1, h2, parameters, max_candidates):
    tasks, excluded = [], []
    null_points = set()
    for grid, point in enumerate(parameters):
        candidates, _ = dated_candidates(
            species,
            h1,
            h2,
            point.hybridization_age,
            max_candidates=max_candidates,
            require_evaluated=False,
        )
        null_key = point.duplication, point.loss, point.ne
        if null_key not in null_points:
            tasks.append((0, "NA", grid, candidates[0].tree, point))
            null_points.add(null_key)
        for candidate in candidates[1:]:
            if candidate.status == "evaluated":
                tasks.append((candidate.id, candidate.h2, grid, candidate.tree, point))
            else:
                excluded.append(
                    {
                        "mul.tree": candidate.id,
                        "grid": grid,
                        "h2.node": candidate.h2,
                        "reason": candidate.reason,
                    }
                )
    if not any(task[0] for task in tasks):
        raise ValueError("No temporally compatible allopolyploid grid points.")
    return tasks, excluded


def build_bank(task, model):
    candidate, h2, grid, population, point = task
    if model.get("integration") in ("detection-rb", "hybrid-rb"):
        from nwkit.mul_locus_integral import build_integrated_bank

        return build_integrated_bank(task, model)
    if model.get("integration") == "ancestral-stratified":
        return build_stratified_bank(task, model)
    rng = np.random.default_rng(
        np.random.SeedSequence([model["seed"], 0, candidate, grid])
    )
    values, attempts = sample_selected(
        population, point, model, rng, count=model["samples"]
    )
    return LocusBank(
        candidate, h2, grid, population, point, Counter(values), len(values), attempts
    )


def build_stratified_bank(task, model):
    candidate, h2, grid, population, point = task
    rate, stem = point.duplication + point.loss, model["ancestral_stem"]
    if not math.isfinite(rate):
        raise ValueError("Combined locus rate must be finite.")
    carry = math.exp(-rate * stem)
    birth = point.duplication / rate * -math.expm1(-rate * stem) if rate else 0.0
    if carry == 0:
        raise ArithmeticError("Ancestral no-event stratum prior underflowed.")
    if point.duplication > 0 and stem > 0 and birth == 0:
        raise ArithmeticError("Ancestral birth stratum prior underflowed.")
    active = [("none", carry)] + ([("birth", birth)] if birth else [])
    if model["samples"] < len(active) or model["samples"] > model["max_attempts"]:
        raise ValueError(
            "Stratified budget must cover active strata and obey attempt cap."
        )
    strata = []
    counts: Counter[tuple] = Counter()
    for i, (condition, weight) in enumerate(active):
        samples = model["samples"] // len(active) + (i < model["samples"] % len(active))
        rng = np.random.default_rng(
            np.random.SeedSequence([model["seed"], 0, candidate, grid, i])
        )
        selected: Counter[tuple] = Counter()
        for _ in range(samples):
            signature = sample_observation(
                population, point, model, rng, stem_condition=condition
            )
            if 2 <= signature_size(signature) <= model["max_observed_tips"]:
                selected[signature] += 1
        strata.append(
            {
                "condition": condition,
                "weight": weight,
                "samples": samples,
                "selected": sum(selected.values()),
                "counts": selected,
            }
        )
        counts.update(selected)
    return LocusBank(
        candidate,
        h2,
        grid,
        population,
        point,
        counts,
        model["samples"],
        model["samples"],
        strata,
    )


def integration_alpha(model, bank_count):
    multiplier = (
        3
        if model.get("integration")
        in ("ancestral-stratified", "detection-rb", "hybrid-rb")
        else 1
    )
    budget = (1 - model["confidence"], bank_count * multiplier)
    categories = category_bound(
        len(model["detection"]), model["max_observed_tips"], confidence_budget=budget
    )
    return budget[0] / (categories * budget[1])


def pattern_probability(bank, signature, alpha):
    if bank.strata is None:
        hits = bank.counts[signature]
        low, high = probability_interval(hits, bank.samples, alpha)
        return hits / bank.samples, low, high
    numerator, denominator = 0.0, 0.0
    n_low, n_high, d_low, d_high = 0.0, 0.0, 0.0, 0.0
    for stratum in bank.strata:
        weight, samples = stratum["weight"], stratum["samples"]
        hits, selected = stratum["counts"][signature], stratum["selected"]
        lower, upper = probability_interval(hits, samples, alpha)
        s_lower, s_upper = probability_interval(selected, samples, alpha)
        numerator += weight * hits / samples
        denominator += weight * selected / samples
        n_low += weight * lower
        n_high += weight * upper
        d_low += weight * s_lower
        d_high += weight * s_upper
    if denominator <= 0:
        raise ValueError("Stratified MC has no selected probability mass.")
    return (
        numerator / denominator,
        n_low / d_high,
        min(1.0, n_high / d_low) if d_low else 1.0,
    )


def category_bound(species_count, max_tips, *, confidence_budget=None):
    total, topologies = 0, 1
    for n in range(2, max_tips + 1):
        if n > 2:
            topologies *= 2 * n - 3
        total += topologies * species_count**n
        if confidence_budget is not None:
            alpha = confidence_budget[0] / (total * confidence_budget[1])
            if alpha <= 0 or 1 - alpha / 2 == 1:
                raise ValueError(
                    "Simultaneous MC confidence precision exceeds floating-point range."
                )
    return total


def probability_interval(hits, samples, alpha):
    lower = 0.0 if hits == 0 else float(beta.ppf(alpha / 2, hits, samples - hits + 1))
    upper = (
        1.0
        if hits == samples
        else float(beta.ppf(1 - alpha / 2, hits + 1, samples - hits))
    )
    if not math.isfinite(lower + upper) or (hits < samples and upper == 1):
        raise ArithmeticError("MC confidence bound is numerically unrepresentable.")
    return lower, upper


def score_bank(bank, observations, alpha):
    scores, lower, upper, rows = [], [], [], []
    for signature, multiplicity in Counter(observations).items():
        hits = bank.counts[signature]
        estimate, lo, hi = pattern_probability(bank, signature, alpha)
        score = math.log(estimate) if estimate else -math.inf
        low = math.log(lo) if lo else -math.inf
        high = math.log(hi)
        scores.append(multiplicity * score)
        lower.append(multiplicity * low)
        upper.append(multiplicity * high)
        rows.append(
            {
                "signature": signature,
                "multiplicity": multiplicity,
                "hits": hits,
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


def parameter_record(bank):
    values = asdict(bank.parameters)
    if bank.candidate == 0:
        values["hybridization_age"] = None
    return values


def search_banks(banks, observations, alpha, *, scorer=None):
    observations = tuple(observations)
    scorer = score_bank if scorer is None else scorer
    rows = [scorer(bank, observations, alpha) for bank in banks]
    null = sorted(
        (r for r in rows if r["mul.tree"] == 0),
        key=lambda r: (-r["log_likelihood"], r["grid"]),
    )
    alt = sorted(
        (r for r in rows if r["mul.tree"] != 0),
        key=lambda r: (-r["log_likelihood"], r["mul.tree"], r["grid"]),
    )
    if (
        not null
        or not alt
        or not all(math.isfinite(r[0]["log_likelihood"]) for r in (null, alt))
    ):
        raise ValueError(
            "Locus MC has no finite null/alternative score; increase samples, no pseudocount fallback."
        )
    point = 2 * (alt[0]["log_likelihood"] - null[0]["log_likelihood"])
    low = 2 * (max(r["mc_lower"] for r in alt) - max(r["mc_upper"] for r in null))
    high = 2 * (max(r["mc_upper"] for r in alt) - max(r["mc_lower"] for r in null))
    return {
        "rows": rows,
        "null": null[0],
        "alternative": alt[0],
        "contrast": point,
        "contrast_lower": low,
        "contrast_upper": high,
    }


def null_replicates(
    banks,
    null,
    observations,
    model,
    alpha,
    replicates,
    sampler,
    *,
    grid,
    scorer=None,
):
    rows = []
    for replicate in range(replicates):
        seed = (
            [model["seed"], 2, null.grid, replicate]
            if grid
            else [model["seed"], 1, replicate]
        )
        rng = np.random.default_rng(np.random.SeedSequence(seed))
        try:
            data, attempts = sampler(
                null.population, null.parameters, model, rng, count=len(observations)
            )
            data = tuple(data)
            if len(data) != len(observations):
                raise ValueError(
                    f"Null sampler expected {len(observations)} families, got {len(data)}."
                )
            options = {} if scorer is None else {"scorer": scorer}
            fit = search_banks(banks, data, alpha, **options)
        except (ValueError, ArithmeticError) as error:
            source = f" generating grid {null.grid}" if grid else ""
            raise type(error)(
                f"Null bootstrap{source} replicate {replicate + 1}: {error}"
            ) from error
        rows.append(
            {
                "replicate": replicate + 1,
                "null_grid": fit["null"]["grid"],
                "alternative_grid": fit["alternative"]["grid"],
                "alternative_candidate": fit["alternative"]["mul.tree"],
                "contrast": fit["contrast"],
                "contrast_lower": fit["contrast_lower"],
                "contrast_upper": fit["contrast_upper"],
                "attempts": attempts,
            }
        )
        if grid:
            rows[-1]["generating_null_grid"] = null.grid
    return rows


def calibration_probabilities(observed, rows):
    return {
        key: (1 + sum(r[first] >= observed[second] for r in rows)) / (len(rows) + 1)
        for key, first, second in (
            ("p_value", "contrast", "contrast"),
            ("mc_p_lower", "contrast_lower", "contrast_upper"),
            ("mc_p_upper", "contrast_upper", "contrast_lower"),
        )
    }


def calibrate(
    banks,
    observations,
    model,
    alpha,
    replicates,
    *,
    sampler=None,
    null_calibration="plug-in",
    scorer=None,
):
    if null_calibration not in ("plug-in", "grid-supremum"):
        raise ValueError("Unknown locus null calibration method.")
    if (
        type(replicates) is not int
        or replicates < 0
        or (null_calibration == "grid-supremum" and not replicates)
    ):
        raise ValueError(
            "Locus bootstrap replicates must be a nonnegative integer; grid-supremum requires a positive count."
        )
    sampler = sample_selected if sampler is None else sampler
    banks = tuple(banks)
    observations = tuple(observations)
    options = {} if scorer is None else {"scorer": scorer}
    observed = search_banks(banks, observations, alpha, **options)
    if not replicates:
        return observed, None
    nulls = sorted((b for b in banks if b.candidate == 0), key=lambda b: b.grid)
    if null_calibration == "plug-in":
        null = next(b for b in nulls if b.grid == observed["null"]["grid"])
        rows = null_replicates(
            banks,
            null,
            observations,
            model,
            alpha,
            replicates,
            sampler,
            grid=False,
            scorer=scorer,
        )
        return observed, {
            "method": "plug-in-finite-grid-parametric-bootstrap",
            "replicates": replicates,
            **calibration_probabilities(observed, rows),
            "null_generating_parameters": parameter_record(null),
            "rows": rows,
            "limitations": "Conditional on supplied model/grid and numerical fitted null; not uniform composite-null coverage or a gene-tree bootstrap. MC bounds cover score-bank error, not uncertainty in the generating null grid point, bootstrap sampling error, or biological uncertainty.",
        }
    return observed, grid_calibration(
        banks,
        nulls,
        observations,
        observed,
        model,
        alpha,
        replicates,
        sampler,
        scorer=scorer,
    )


def grid_calibration(
    banks,
    nulls,
    observations,
    observed,
    model,
    alpha,
    replicates,
    sampler,
    *,
    scorer=None,
):
    rows, summaries = [], []
    for null in nulls:
        draws = null_replicates(
            banks,
            null,
            observations,
            model,
            alpha,
            replicates,
            sampler,
            grid=True,
            scorer=scorer,
        )
        summaries.append(
            {
                "generating_null_grid": null.grid,
                "null_generating_parameters": parameter_record(null),
                "replicates": replicates,
                **calibration_probabilities(observed, draws),
            }
        )
        rows.extend(draws)
    return {
        "method": "finite-grid-supremum-Monte-Carlo-test",
        "replicates": replicates,
        "num_null_grid_points": len(nulls),
        "total_replicates": len(rows),
        **{
            key: max(s[key] for s in summaries)
            for key in ("p_value", "mc_p_lower", "mc_p_upper")
        },
        "least_favorable_grid": max(summaries, key=lambda s: s["p_value"])[
            "generating_null_grid"
        ],
        "least_favorable_mc_upper_grid": max(summaries, key=lambda s: s["mc_p_upper"])[
            "generating_null_grid"
        ],
        "null_grid_calibrations": summaries,
        "rows": rows,
        "limitations": "Supremum over every supplied finite null grid point, not a continuous nuisance space. Rank-test validity requires the true null on that grid, independently frozen banks, identical observation/selection laws and complete IID replicates. No model-misspecification or empirical gene-tree-error guarantee. MC-overlap bounds cover score-bank error, not bootstrap sampling uncertainty or biological confidence.",
    }
