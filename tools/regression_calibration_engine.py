"""Exercise existing NWKIT estimators; keep reference tests out of production."""

import io
import signal
import tempfile
import time
import warnings
from contextlib import ExitStack, contextmanager, redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm

from nwkit import phylogenetic_glmm as glmm
from nwkit import regress
from nwkit.cli import main as nwkit_main
from nwkit.gaussian import draw_from_factor


def sidecar(ids, trait, matrix, tree_id="OG"):
    return pd.DataFrame(
        [
            {
                "tree_id": tree_id,
                "trait": trait,
                "contrast_id_1": first,
                "contrast_id_2": second,
                "sampling_covariance": matrix[i, j],
            }
            for i, first in enumerate(ids)
            for j, second in enumerate(ids)
            if j >= i
        ]
    )


def rsc_tables(case, data):
    events = np.asarray(data["events"])
    rows = []
    predictor_rows = []
    names = [f"x{j}" for j in range(case.predictors)]
    for i, (event, lineage, value) in enumerate(
        zip(events, data["lineages"], data["y"], strict=True)
    ):
        rows.append(
            {
                "tree_id": "OG",
                "gene_clade_id": f"g{i}",
                "lineage_clade_id": f"l{lineage}",
                "event_type": "speciation",
                "eligible": "yes",
                "coverage_status": "complete",
                "species_event_id": f"e{event}",
                "species_event_taxa": f"t{event}",
                "species_numerator_event_id": f"a{event}",
                "species_denominator_event_id": f"b{event}",
                "trait": "y",
                "evolution_model": "brownian",
                "evolution_parameter_name": "",
                "evolution_parameter": "",
                "branch_length_mode": "original",
                "raw_contrast": value,
                "contrast_variance": 1.0,
            }
        )
    for event in np.unique(events):
        index = np.flatnonzero(events == event)[0]
        for j, name in enumerate(names):
            predictor_rows.append(
                {
                    "tree_id": "species",
                    "branch_clade_id": f"e{event}",
                    "descendant_taxa": f"t{event}",
                    "numerator_clade_id": f"a{event}",
                    "denominator_clade_id": f"b{event}",
                    "trait": name,
                    "evolution_model": "brownian",
                    "evolution_parameter_name": "",
                    "evolution_parameter": "",
                    "branch_length_mode": "original",
                    "raw_contrast": data["x"][index, j],
                    "contrast_variance": 1.0,
                }
            )
    response = pd.DataFrame(rows)
    predictor = pd.DataFrame(predictor_rows)
    sampling = sidecar(response.gene_clade_id, "y", np.diag(data["sampling"]))
    predictor_sampling = None
    if case.predictor_variance:
        # NWKIT consumes contrast-scale conditional predictor uncertainty here.
        ids = [f"e{event}" for event in np.unique(events)]
        predictor_sampling = pd.concat(
            [
                sidecar(
                    ids, name, np.eye(len(ids)) * case.predictor_variance, "species"
                )
                for name in names
            ],
            ignore_index=True,
        )
    return response, predictor, names, sampling, predictor_sampling


@contextmanager
def record_refits(engine, enabled, diagnostics):
    """Observe unmodified calls, including discarded attempts, within one worker.

    The first call is the original fit; subsequent calls are bootstrap attempts.
    These hooks are limited to the direct fixed-model engines exercised here.
    """
    if not enabled:
        yield
        return
    attempts = []

    def observer(original):
        def call(*args, **kwargs):
            row = {"accepted": False}
            attempts.append(row)
            try:
                fit = original(*args, **kwargs)
                if engine == "rsc":
                    row["accepted"] = bool(fit["optimizer_converged"])
                    row["variance_components"] = fit["component_variances"]
                else:
                    row["accepted"] = bool(
                        fit.optimizer_converged
                        and np.isfinite(fit.log_likelihood)
                        and np.isfinite(fit.coefficients).all()
                    )
                    row["boundary"] = bool(fit.boundary_warning)
                    row["separation"] = bool(fit.separation_warning)
                    row["covariance_status"] = fit.coefficient_covariance_status
                return fit
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                row["error"] = f"{type(exc).__name__}: {exc}"
                raise

        return call

    with ExitStack() as stack:
        if engine == "rsc":
            for name in ("_profile_covariance_fit", "fit_conditional_eiv_gaussian"):
                stack.enter_context(
                    patch.object(regress, name, observer(getattr(regress, name)))
                )
        else:
            stack.enter_context(
                patch.object(
                    glmm,
                    "_call_phylogenetic_glmm",
                    observer(glmm._call_phylogenetic_glmm),
                )
            )
        try:
            yield
        finally:
            refits = attempts[1:]
            diagnostics.update(
                bootstrap_attempts=len(refits),
                bootstrap_successes=sum(row["accepted"] for row in refits),
                bootstrap_attempt_records=refits,
            )


def fit_rsc(case, data, method, replicates, seed):
    response, predictor, names, sampling, predictor_sampling = rsc_tables(case, data)
    fitted = regress.fit_reconciled_pgls(
        response,
        predictor,
        ["y"],
        names,
        model="hierarchical",
        event_weighting="event",
        inference=method,
        bootstrap_replicates=replicates,
        seed=seed,
        response_sampling_covariance=sampling,
        predictor_sampling_covariance=predictor_sampling,
        event_random_effect="auto",
        lineage_random_slope="auto",
        reml=True,
    )
    return fitted.loc[fitted.term == "x0"].iloc[0].to_dict()


def call_glmm(case, data, method="wald", replicates=2, seed=1, design=None, y=None):
    if design is None:
        design = np.column_stack([np.ones(len(data["y"])), data["x"]])
    return glmm.fit_phylogenetic_glmm(
        data["y"] if y is None else y,
        design,
        lambda _: data["covariance"],
        family=case.family,
        levels=["0", "1"] if case.family == "binomial" else None,
        reference="0" if case.family == "binomial" else None,
        coefficient_penalty="none",
        inference=method,
        bootstrap_replicates=replicates,
        seed=seed,
    )


def fit_glmm(case, data, method, replicates, seed):
    fit = call_glmm(case, data, method, replicates, seed)
    coefficient = float(fit.coefficients.reshape(-1)[1])
    se, statistic, p_value, lower, upper, status = glmm.summarize_glmm_coefficient(
        fit, 1, coefficient, norm.ppf(0.975)
    )
    return {
        "coefficient": coefficient,
        "standard_error": se,
        "statistic": statistic,
        "p_value": p_value,
        "confidence_interval_lower": lower,
        "confidence_interval_upper": upper,
        "inference_status": status,
        "optimizer_converged": fit.optimizer_converged,
        "optimizer_message": fit.optimizer_message,
        "boundary_warning": fit.boundary_warning,
        "separation_warning": fit.separation_warning,
        "coefficient_covariance_status": fit.coefficient_covariance_status,
        "component_variances": fit.component_variances,
        "log_likelihood": fit.log_likelihood,
        "response_dispersion": fit.dispersion,
    }


def oracle(data):
    x, y, covariance = data["x"], data["y"], data["true_covariance"]
    inverse_x = np.linalg.solve(covariance, x)
    information = x.T @ inverse_x
    beta_covariance = np.linalg.inv(information)
    coefficient = float((beta_covariance @ inverse_x.T @ y)[0])
    se = np.sqrt(beta_covariance[0, 0])
    return {
        "coefficient": coefficient,
        "standard_error": se,
        "p_value": float(2 * norm.sf(abs(coefficient / se))),
        "confidence_interval_lower": coefficient - norm.ppf(0.975) * se,
        "confidence_interval_upper": coefficient + norm.ppf(0.975) * se,
        "inference_status": "ok",
    }


def null_tail(statistic, simulated_statistics):
    """Fixed B attempts; failed refits yield bounds, never silently replaced."""
    successes = [value for value in simulated_statistics if value is not None]
    exceedances = sum(value >= statistic for value in successes)
    failures = len(simulated_statistics) - len(successes)
    denominator = len(simulated_statistics) + 1
    lower = (1 + exceedances) / denominator
    upper = (1 + exceedances + failures) / denominator
    return {
        "p_value": lower if not failures else None,
        "p_value_failure_lower": lower,
        "p_value_failure_upper": upper,
        "bootstrap_attempts": len(simulated_statistics),
        "bootstrap_successes": len(successes),
        "inference_status": "ok" if not failures else "bootstrap-refit-failed",
        "statistic": statistic,
    }


def checked_statistic(null_objective, full_objective):
    statistic = 2 * (null_objective - full_objective)
    if not np.isfinite(statistic) or statistic < -1e-6:
        raise ValueError(f"Invalid nested objective difference: {statistic}")
    return max(0.0, float(statistic))


def rsc_null_bootstrap(case, data, replicates, seed):
    """Reference composite-objective test; covariance terms fixed across models.

    This does not claim a chi-squared reference distribution. Its unconditional
    calibration must be measured like that of the production procedures.
    """
    x, y = data["x"], data["y"]
    events, event_inverse, counts = np.unique(
        data["events"], return_inverse=True, return_counts=True
    )
    lineages, lineage_inverse, lineage_counts = np.unique(
        data["lineages"], return_inverse=True, return_counts=True
    )
    fixed, components, factors, *_ = regress._build_covariance_components(
        x,
        np.ones(len(y)),
        np.diag(data["sampling"]),
        event_inverse,
        counts,
        lineage_inverse,
        lineage_counts,
        events,
        lineages,
        [f"x{j}" for j in range(x.shape[1])],
        n_events=len(events),
        num_parameters=x.shape[1],
        event_weighting="event",
        model="hierarchical",
        event_random_effect="auto",
        lineage_random_slope="auto",
    )

    def fit(response, design):
        value = regress._profile_covariance_fit(
            response,
            design,
            fixed,
            components,
            reml=False,
            component_factors=factors,
            likelihood_observations=len(events),
            likelihood_groups=event_inverse,
        )
        if not value["optimizer_converged"]:
            raise ValueError("Variance optimizer did not converge")
        return value

    null = fit(y, x[:, 1:])
    full = fit(y, x)
    statistic = checked_statistic(null["objective"], full["objective"])
    rng = np.random.default_rng(seed)
    values, errors = [], []
    for _ in range(replicates):
        simulated = x[:, 1:] @ null["beta"] + draw_from_factor(
            null["cholesky"], rng.normal(size=len(y)), rng=rng
        )
        try:
            null_refit, full_refit = fit(simulated, x[:, 1:]), fit(simulated, x)
            values.append(
                checked_statistic(null_refit["objective"], full_refit["objective"])
            )
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            values.append(None)
            errors.append(str(exc))
    return {
        **null_tail(statistic, values),
        "coefficient": float(full["beta"][0]),
        "bootstrap_errors": errors,
        "reference_objective": "event-balanced-composite-ML",
        "interval_method": "not-computed",
    }


def glmm_null_bootstrap(case, data, replicates, seed):
    x = np.column_stack([np.ones(len(data["y"])), data["x"]])

    def fit(y, design):
        value = call_glmm(case, data, design=design, y=y)
        if not value.optimizer_converged:
            raise ValueError("GLMM optimizer did not converge")
        return value

    null, full = fit(data["y"], x[:, :1]), fit(data["y"], x)
    statistic = checked_statistic(-null.log_likelihood, -full.log_likelihood)
    covariance = data["covariance"] / np.mean(np.diag(data["covariance"]))
    rng = np.random.default_rng(seed)
    values, errors = [], []
    for _ in range(replicates):
        latent = rng.multivariate_normal(
            np.zeros(len(x)), null.component_variances["phylogenetic"] * covariance
        )
        eta = null.coefficients.reshape(-1)[0] + latent
        if case.family == "binomial":
            y = rng.binomial(1, expit(eta))
        elif case.family == "poisson":
            y = rng.poisson(np.exp(eta))
        else:
            # NWKIT's NB2 dispersion is alpha, while NumPy takes size=1/alpha.
            dispersion = 1.0 / null.dispersion
            y = rng.negative_binomial(
                dispersion, dispersion / (dispersion + np.exp(eta))
            )
        try:
            null_refit, full_refit = fit(y, x[:, :1]), fit(y, x)
            values.append(
                checked_statistic(
                    -null_refit.log_likelihood, -full_refit.log_likelihood
                )
            )
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            values.append(None)
            errors.append(str(exc))
    return {
        **null_tail(statistic, values),
        "coefficient": float(full.coefficients.reshape(-1)[1]),
        "bootstrap_errors": errors,
        "reference_objective": "laplace-ML",
        "interval_method": "not-computed",
    }


def is_applicable(case, method):
    if case.engine == "rsc-tips":
        return method in {"wald", "parametric-bootstrap"}
    if case.engine == "rsc":
        return method in {"wald", "parametric-bootstrap", "oracle"} or (
            method == "null-bootstrap" and case.predictor_variance == 0
        )
    return method in {
        "wald",
        "parametric-bootstrap",
        "null-bootstrap",
        "profile-likelihood",
        "likelihood-ratio",
    }


@contextmanager
def time_limit(seconds):
    def timeout(signum, frame):
        raise TimeoutError(f"Fit exceeded {seconds}s")

    previous = signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def evaluate(case, data, method, replicates, seed, timeout):
    started = time.monotonic()
    result = {"method": method, "status": "completed", "target_beta": case.beta}
    if not is_applicable(case, method):
        return {**result, "status": "not_applicable"}
    n = len(data["y"])
    if case.engine == "rsc-tips":
        observed_by_leaf = {}
        for row in data["expression_records"]:
            observed_by_leaf.setdefault(row["leaf_name"], False)
            observed_by_leaf[row["leaf_name"]] |= row["y"] is not None
        eligible = all(observed_by_leaf.values())
        reason = "missing_all_response_replicates"
    elif case.engine == "rsc":
        eligible = len(np.unique(data["events"])) > case.predictors
        reason = "too_few_events"
    else:
        eligible = (
            n >= 4 and len(np.unique(data["y"])) > 1 and len(np.unique(data["x"])) > 1
        )
        reason = "too_few_species" if n < 4 else "invariant_response_or_predictor"
    if not eligible:
        return {**result, "status": "ineligible", "reason": reason, "seconds": 0.0}
    try:
        with time_limit(timeout), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with record_refits(
                case.engine,
                method == "parametric-bootstrap" and case.engine != "rsc-tips",
                result,
            ):
                if case.engine == "rsc-tips":
                    fitted = fit_raw_rsc(case, data, method, replicates, seed)
                elif method == "oracle":
                    fitted = oracle(data)
                elif method == "null-bootstrap":
                    operation = (
                        rsc_null_bootstrap
                        if case.engine == "rsc"
                        else glmm_null_bootstrap
                    )
                    fitted = operation(case, data, replicates, seed)
                else:
                    operation = fit_rsc if case.engine == "rsc" else fit_glmm
                    fitted = operation(case, data, method, replicates, seed)
                result.update(fitted)
            result["warnings"] = sorted({str(w.message) for w in caught})
    except (ValueError, RuntimeError, np.linalg.LinAlgError, TimeoutError) as exc:
        result.update(status="fit_failed", error=f"{type(exc).__name__}: {exc}")
    result["seconds"] = time.monotonic() - started
    return result


def fit_raw_rsc(case, data, method, replicates, seed):
    with tempfile.TemporaryDirectory(prefix="nwkit-calibration-tips-") as temporary:
        root = Path(temporary)
        (root / "species.nwk").write_text(data["species_tree"])
        (root / "gene.nwk").write_text(data["gene_tree"])
        pd.DataFrame(data["expression_records"]).to_csv(
            root / "expression.tsv", sep="\t", index=False
        )
        pd.DataFrame(data["trait_records"]).to_csv(
            root / "traits.tsv", sep="\t", index=False
        )
        argv = [
            "regress",
            "--gene-tree",
            str(root / "gene.nwk"),
            "--reconciliation-tree",
            str(root / "gene.nwk"),
            "--species-tree",
            str(root / "species.nwk"),
            "--expression",
            str(root / "expression.tsv"),
            "--species-traits",
            str(root / "traits.tsv"),
            "--responses",
            "y",
            "--predictors",
            "x",
            "--tree-id",
            "OG",
            "--out-prefix",
            str(root / "result"),
            "--event-source",
            "lca",
            "--species-parser",
            "legacy",
            "--reconciled-model",
            "hierarchical",
            "--event-weighting",
            "event",
            "--speciation-coverage",
            "complete",
            "--gene-evolution-model",
            case.gene_evolution_model,
            "--species-evolution-model",
            "brownian",
            "--reml",
            "yes",
            "--inference",
            method,
            "--bootstrap-replicates",
            str(replicates),
            "--seed",
            str(seed),
        ]
        if case.gene_evolution_model != "brownian":
            argv.extend(["--gene-evolution-parameter", "auto"])
        if case.biological_replicates > 1:
            argv.extend(
                [
                    "--response-biological-id",
                    "biological",
                    "--response-within-variance",
                    "pooled",
                ]
            )
        if case.technical_replicates > 1:
            argv.extend(
                [
                    "--response-technical-id",
                    "technical",
                    "--response-technical-aggregation",
                    "mean",
                ]
            )
        if case.predictor_variance:
            argv.extend(
                [
                    "--predictor-biological-id",
                    "biological",
                    "--predictor-within-variance",
                    "pooled",
                ]
            )
        log = io.StringIO()
        with redirect_stdout(log), redirect_stderr(log):
            nwkit_main(argv)
        frame = pd.read_csv(
            root / "result.regression.tsv", sep="\t", keep_default_na=False
        )
        selected = frame.loc[frame.term == "x"]
        if len(selected) != 1:
            raise ValueError(
                f"Raw RSC produced {len(selected)} target coefficient rows: {log.getvalue()}"
            )
        row = selected.iloc[0].to_dict()
        row["raw_cli_log"] = log.getvalue()
        row["replicate_summaries"] = {
            path.name: pd.read_csv(path, sep="\t", keep_default_na=False).to_dict(
                "records"
            )
            for path in root.glob("*tip-summary.tsv")
        }
        row["bootstrap_attempt_accounting"] = (
            "not-instrumented-in-raw-pipeline"
            if method == "parametric-bootstrap"
            else "not-applicable"
        )
        return row
