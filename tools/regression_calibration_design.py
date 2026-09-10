"""Independent data generation for fixed-model regression calibration.

This is a validation tool, not a production estimator. In particular, the RSC
event-balanced objective is composite; its working covariance is not a claim
that it defines an ordinary Gaussian likelihood.
"""

import hashlib
import json
import re
from dataclasses import asdict, dataclass, replace

import numpy as np
from scipy.special import expit


@dataclass(frozen=True)
class Case:
    name: str
    engine: str = "rsc"
    size: int = 20
    copies: int = 1
    concentrated: bool = False
    predictors: int = 1
    correlation: float = 0.0
    beta: float = 0.0
    event_variance: float = 0.0
    lineage_variance: float = 0.0
    sampling_variance: float = 0.0
    biological_replicates: int = 1
    technical_replicates: int = 1
    estimated_se: bool = False
    predictor_variance: float = 0.0
    missing: float = 0.0
    missing_mechanism: str = "mcar"
    working_covariance_dgm: bool = False
    family: str = "binomial"
    phylogenetic_variance: float = 0.5
    baseline: float = 0.5
    # Generator NB size r: Var(Y|eta) = mu + mu^2/r. NWKIT reports alpha=1/r.
    dispersion: float = 2.0
    zero_inflation: float = 0.0
    count_error: float = 0.0
    tree_shape: str = "balanced"
    gene_evolution_model: str = "brownian"

    def __post_init__(self):
        if not self.name or self.engine not in {"rsc", "rsc-tips", "glmm"}:
            raise ValueError("A case needs a name and a recognized engine")
        for field in (
            "size",
            "copies",
            "predictors",
            "biological_replicates",
            "technical_replicates",
        ):
            value = getattr(self, field)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{field} must be a positive integer")
        if self.size < 2:
            raise ValueError(
                "At least two tips/events are required in a generation case"
            )
        for field in (
            "event_variance",
            "lineage_variance",
            "sampling_variance",
            "predictor_variance",
            "phylogenetic_variance",
        ):
            value = getattr(self, field)
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{field} must be finite and nonnegative")
        if not np.isfinite(self.beta) or not np.isfinite(self.correlation):
            raise ValueError("Effect and correlation must be finite")
        if (
            self.predictors > 1
            and not -1 / (self.predictors - 1) <= self.correlation <= 1
        ):
            raise ValueError("Predictor correlation is not positive semidefinite")
        for field in ("missing", "zero_inflation", "count_error"):
            if not 0 <= getattr(self, field) < 1:
                raise ValueError(f"{field} must be in [0,1)")
        if self.family not in {"binomial", "poisson", "negative-binomial"}:
            raise ValueError("Unsupported response family")
        if (
            not np.isfinite(self.baseline)
            or self.baseline <= 0
            or (self.family == "binomial" and self.baseline >= 1)
        ):
            raise ValueError("Invalid response baseline")
        if not np.isfinite(self.dispersion) or self.dispersion <= 0:
            raise ValueError("NB size must be positive")
        if self.tree_shape not in {"balanced", "comb"}:
            raise ValueError("Unsupported tree shape")
        if self.missing_mechanism not in {"mcar", "mar", "mnar", "clade"}:
            raise ValueError("Unsupported missingness mechanism")
        if self.gene_evolution_model not in {"brownian", "lambda"}:
            raise ValueError(
                "Raw-tip validation currently supports Brownian or fitted lambda"
            )


def seed_for(seed, case_name, replicate, stream):
    """Stable streams independent of methods, ordering, workers and sharding."""
    payload = json.dumps([int(seed), case_name, int(replicate), stream]).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def cases():
    result = [Case(f"rsc-e{n}", size=n) for n in (2, 3, 5, 10, 20, 50)]
    for n in (5, 20):
        base = Case(f"rsc-e{n}-copies5", size=n, copies=5, event_variance=0.5)
        result.extend(
            [
                base,
                replace(base, name=base.name + "-boundary", event_variance=0.0),
                replace(base, name=base.name + "-concentrated", concentrated=True),
                replace(base, name=base.name + "-working", working_covariance_dgm=True),
                replace(base, name=base.name + "-lineage", lineage_variance=0.5),
            ]
        )
    base = Case("rsc-e20-replicates", sampling_variance=1.0, biological_replicates=5)
    result.extend(
        [
            base,
            replace(base, name=base.name + "-estimated-se", estimated_se=True),
            replace(base, name=base.name + "-predictor-error", predictor_variance=0.5),
            Case("rsc-e20-joint-rho09", predictors=2, correlation=0.9),
            Case("rsc-e20-joint-rho099", predictors=2, correlation=0.99),
            Case("rsc-e2-rank-ineligible", size=2, predictors=2),
        ]
    )
    for mechanism in ("mcar", "mar", "mnar", "clade"):
        result.append(
            Case(
                f"rsc-e20-missing-{mechanism}", missing=0.3, missing_mechanism=mechanism
            )
        )
    for n in (8, 30, 60):
        for baseline in (0.5, 0.05):
            result.append(
                Case(
                    f"binomial-n{n}-p{baseline}",
                    engine="glmm",
                    size=n,
                    baseline=baseline,
                )
            )
    for family in ("binomial", "poisson", "negative-binomial"):
        base = Case(
            f"{family}-n30",
            engine="glmm",
            size=30,
            family=family,
            baseline=0.5 if family == "binomial" else 2.0,
        )
        if family != "binomial":
            result.append(base)
        result.extend(
            [
                replace(base, name=base.name + "-boundary", phylogenetic_variance=0.0),
                replace(
                    base,
                    name=base.name + "-strong-phylogeny",
                    phylogenetic_variance=2.0,
                ),
                replace(
                    base,
                    name=base.name + "-missing-mar",
                    missing=0.3,
                    missing_mechanism="mar",
                ),
                replace(
                    base, name=base.name + "-count-error", beta=0.5, count_error=0.3
                ),
            ]
        )
    result.extend(
        [
            Case("binomial-n30-comb", engine="glmm", size=30, tree_shape="comb"),
            Case(
                "poisson-n30-zero-inflation",
                engine="glmm",
                size=30,
                family="poisson",
                baseline=2.0,
                zero_inflation=0.3,
            ),
            Case(
                "negative-binomial-n30-poisson-limit",
                engine="glmm",
                size=30,
                family="negative-binomial",
                baseline=2.0,
                dispersion=1e6,
            ),
            Case("rsc-e20-alternative", beta=0.5),
            Case("binomial-n30-alternative", engine="glmm", size=30, beta=0.5),
        ]
    )
    raw = Case("rsc-tips-n8", engine="rsc-tips", size=8)
    result.extend(
        [
            raw,
            replace(raw, name="rsc-tips-n8-copies3", copies=3, event_variance=0.5),
            replace(
                raw,
                name="rsc-tips-n8-replicates",
                biological_replicates=5,
                sampling_variance=1.0,
            ),
            replace(
                raw,
                name="rsc-tips-n8-technical",
                biological_replicates=5,
                technical_replicates=2,
                sampling_variance=1.0,
            ),
            replace(
                raw,
                name="rsc-tips-n8-predictor-replicates",
                biological_replicates=5,
                sampling_variance=1.0,
                predictor_variance=0.5,
            ),
            replace(
                raw,
                name="rsc-tips-n8-partial-replicates",
                biological_replicates=5,
                sampling_variance=1.0,
                missing=0.3,
            ),
            replace(raw, name="rsc-tips-n8-lambda-auto", gene_evolution_model="lambda"),
        ]
    )
    assert len({case.name for case in result}) == len(result)
    return result


def tree_covariance(n, shape="balanced"):
    """Brownian covariance computed from shared edges, without NWKIT helpers."""
    covariance = np.zeros((n, n))

    def visit(indices):
        if len(indices) == 1:
            covariance[indices[0], indices[0]] += 1.0
            return f"S{indices[0]}:1"
        cut = 1 if shape == "comb" else len(indices) // 2
        children = []
        for child in (indices[:cut], indices[cut:]):
            if len(child) > 1:
                covariance[np.ix_(child, child)] += 1.0
            children.append(visit(child))
        return "(" + ",".join(children) + "):1"

    newick = visit(list(range(n))).rsplit(":1", 1)[0] + ";"
    covariance /= np.mean(np.diag(covariance))
    return newick, covariance


def observed_mask(rng, case, x, y, groups):
    if case.missing == 0:
        return np.ones(len(y), dtype=bool)
    if case.missing_mechanism == "mcar":
        probability = np.full(len(y), case.missing)
    elif case.missing_mechanism in {"mar", "mnar"}:
        score = x[:, 0] if case.missing_mechanism == "mar" else y
        score = (score - np.mean(score)) / max(np.std(score), 1e-12)
        probability = expit(np.log(case.missing / (1 - case.missing)) + score)
    elif case.missing_mechanism == "clade":
        # Deliberate contiguous group loss, not an independent Bernoulli sample.
        return groups >= int(np.ceil(case.size * case.missing))
    else:
        raise ValueError(case.missing_mechanism)
    return rng.random(len(y)) >= probability


def generate_rsc(case, rng):
    counts = np.full(case.size, case.copies)
    if case.concentrated:
        counts[:] = 1
        counts[: max(2, case.size // 5)] = case.copies
    events = np.repeat(np.arange(case.size), counts)
    lineages = np.concatenate([np.arange(c) for c in counts])
    corr = np.full((case.predictors, case.predictors), case.correlation)
    np.fill_diagonal(corr, 1.0)
    event_x = rng.multivariate_normal(np.zeros(case.predictors), corr, case.size)
    x = event_x[events]
    beta = np.zeros(case.predictors)
    beta[0] = case.beta
    # Known conditional predictor uncertainty, shared by paralogs in an event.
    true_x = event_x + rng.normal(size=event_x.shape) * np.sqrt(case.predictor_variance)
    diagonal = np.ones(len(events))
    if case.working_covariance_dgm:
        diagonal *= counts[events]
    covariance = np.diag(diagonal)
    covariance += case.event_variance * (events[:, None] == events[None, :])
    if case.lineage_variance:
        for column in range(case.predictors):
            normalized = x[:, column] / np.sqrt(np.mean(event_x[:, column] ** 2))
            covariance += (
                case.lineage_variance
                * np.outer(normalized, normalized)
                * (lineages[:, None] == lineages[None, :])
            )
    sampling = np.full(len(events), case.sampling_variance / case.biological_replicates)
    covariance += np.diag(
        sampling * (counts[events] if case.working_covariance_dgm else 1)
    )
    y = true_x[events] @ beta + rng.multivariate_normal(
        np.zeros(len(events)), covariance
    )
    true_covariance = covariance + np.dot(beta, beta) * case.predictor_variance * (
        events[:, None] == events[None, :]
    )
    estimated_sampling = sampling.copy()
    if case.estimated_se:
        df = case.biological_replicates - 1
        if df < 1:
            raise ValueError("Estimated SE requires at least two biological replicates")
        estimated_sampling *= rng.chisquare(df, len(events)) / df
    keep = observed_mask(rng, case, x, y, events)
    return {
        "x": x[keep],
        "y": y[keep],
        "events": events[keep],
        "lineages": lineages[keep],
        "sampling": estimated_sampling[keep],
        "true_covariance": true_covariance[np.ix_(keep, keep)],
        "generated_rows": len(y),
        "generated_events": case.size,
        "removed_rows": int((~keep).sum()),
        "beta": beta,
    }


def generate_glmm(case, rng):
    tree, covariance = tree_covariance(case.size, case.tree_shape)
    copy_latent = rng.multivariate_normal(np.zeros(case.size), covariance)
    true_counts = rng.poisson(np.exp(np.log(2.0) + 0.5 * copy_latent))
    counts = true_counts.copy()
    if case.count_error:
        counts = rng.binomial(counts, 1 - case.count_error) + rng.poisson(
            case.count_error, case.size
        )
    x = np.log1p(counts).reshape(-1, 1)
    latent = rng.multivariate_normal(
        np.zeros(case.size), case.phylogenetic_variance * covariance
    )
    intercept = (
        np.log(case.baseline / (1 - case.baseline))
        if case.family == "binomial"
        else np.log(case.baseline)
    )
    eta = intercept + case.beta * np.log1p(true_counts) + latent
    if case.family == "binomial":
        y = rng.binomial(1, expit(eta))
    elif case.family == "poisson":
        y = rng.poisson(np.exp(eta))
    elif case.family == "negative-binomial":
        y = rng.negative_binomial(
            case.dispersion, case.dispersion / (case.dispersion + np.exp(eta))
        )
    else:
        raise ValueError(case.family)
    if case.zero_inflation:
        y[rng.random(case.size) < case.zero_inflation] = 0
    keep = observed_mask(rng, case, x, y, np.arange(case.size))
    return {
        "x": x[keep],
        "y": y[keep],
        "covariance": covariance[np.ix_(keep, keep)],
        "tree": tree,
        "retained_tips": np.flatnonzero(keep),
        "counts": counts[keep],
        "true_counts": true_counts[keep],
        "generated_rows": case.size,
        "removed_rows": int((~keep).sum()),
        "beta": np.array([case.beta]),
    }


def generate(case, seed):
    rng = np.random.default_rng(seed)
    if case.engine == "rsc-tips":
        return generate_raw_rsc(case, rng)
    return generate_rsc(case, rng) if case.engine == "rsc" else generate_glmm(case, rng)


def generate_raw_rsc(case, rng):
    """Tip-level Brownian generation followed by production reconciliation/PIC.

    Technical replicates are exact duplicates of biological observations in
    this invariance design; they must not increase the biological sample size.
    Missingness affects biological observations and is not silently regenerated.
    """
    tree, covariance = tree_covariance(case.size, case.tree_shape)
    species_tree = re.sub(r"S(\d+)", r"Genus_s\1", tree)
    x = rng.multivariate_normal(np.zeros(case.size), covariance)
    shared = rng.multivariate_normal(
        np.zeros(case.size), covariance * case.event_variance
    )
    expressions, traits, copy_trees, values = [], [], [], []
    for copy in range(case.copies):
        copy_trees.append(re.sub(r"S(\d+)", rf"Genus_s\1_g{copy}", tree).rstrip(";"))
        gene_noise = rng.multivariate_normal(np.zeros(case.size), covariance)
        means = case.beta * x + shared + gene_noise
        values.extend(means)
        for tip, value in enumerate(means):
            for biological in range(case.biological_replicates):
                observed = value + rng.normal() * np.sqrt(case.sampling_variance)
                if rng.random() < case.missing:
                    observed = None
                for technical in range(case.technical_replicates):
                    expressions.append(
                        {
                            "leaf_name": f"Genus_s{tip}_g{copy}",
                            "y": observed,
                            "biological": str(biological),
                            "technical": str(technical),
                        }
                    )
    while len(copy_trees) > 1:
        first, second, *remaining = copy_trees
        copy_trees = [f"({first}:1,{second}:1)", *remaining]
    gene_tree = copy_trees[0] + ";"
    for tip, value in enumerate(x):
        for biological in range(
            case.biological_replicates if case.predictor_variance else 1
        ):
            traits.append(
                {
                    "leaf_name": f"Genus_s{tip}",
                    "x": value + rng.normal() * np.sqrt(case.predictor_variance),
                    "biological": str(biological),
                }
            )
    return {
        "species_tree": species_tree,
        "gene_tree": gene_tree,
        "expression_records": expressions,
        "trait_records": traits,
        "x": x[:, None],
        "y": np.asarray(values),
        "beta": np.array([case.beta]),
        "generated_rows": len(expressions),
        "generated_species": case.size,
        "removed_rows": sum(row["y"] is None for row in expressions),
    }


def case_dict(case):
    return asdict(case)
