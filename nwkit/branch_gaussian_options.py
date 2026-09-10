"""ASR options and validation for fixed branch-specific Gaussian models."""

import math

MODEL = "BRANCH-GAUSSIAN"
BRANCH_OPTIONS = (
    "branch_models",
    "branch_regimes",
    "regime_models",
    "branch_models_out",
    "process_out",
    "prior_samples",
    "prior_root_value",
)


def register_branch_gaussian_options(parser):
    for name, help_text in (
        (
            "branch-models",
            "BRANCH-GAUSSIAN: complete non-root branch_id/model/parameter TSV.",
        ),
        ("branch-regimes", "BRANCH-GAUSSIAN: complete non-root branch_id/regime TSV."),
        (
            "regime-models",
            "BRANCH-GAUSSIAN: regime/model/parameter TSV with --branch-regimes.",
        ),
        (
            "branch-models-out",
            "BRANCH-GAUSSIAN: normalized direct model TSV for reuse.",
        ),
        (
            "process-out",
            "BRANCH-GAUSSIAN: JSON containing the tree, fixed parameters and run settings.",
        ),
    ):
        parser.add_argument("--" + name, "--" + name.replace("-", "_"), help=help_text)
    parser.add_argument(
        "--prior-samples",
        "--prior_samples",
        type=int,
        help="BRANCH-GAUSSIAN with --output prior-samples: latent prior draws (default 1000, 1–10000).",
    )
    parser.add_argument(
        "--prior-root-value",
        "--prior_root_value",
        type=float,
        help="BRANCH-GAUSSIAN prior-samples with a flat root: required fixed starting value.",
    )


def branch_root(args):
    from nwkit.gaussian_tree import GaussianRootPrior

    mode = getattr(args, "root_prior", None)
    if mode not in {"fixed", "flat", "gaussian", "stationary"}:
        raise ValueError(
            "BRANCH-GAUSSIAN requires an explicit --root-prior fixed, flat, gaussian, or stationary."
        )
    mean, variance = (
        getattr(args, "root_mean", None),
        getattr(args, "root_variance", None),
    )
    if mode == "flat":
        if mean is not None or variance is not None:
            raise ValueError("A flat root cannot take --root-mean or --root-variance.")
        return GaussianRootPrior("flat", variance=None)
    if mean is None:
        raise ValueError("BRANCH-GAUSSIAN proper roots require --root-mean.")
    if mode == "fixed" and variance is not None:
        raise ValueError("A fixed root has zero variance; omit --root-variance.")
    return GaussianRootPrior(mode, mean, 0.0 if mode == "fixed" else variance)


def _validate_prior_options(args, root):
    output = getattr(args, "output", None) or "summary"
    count, root_value = (
        getattr(args, "prior_samples", None),
        getattr(args, "prior_root_value", None),
    )
    if output != "prior-samples":
        if count is not None or root_value is not None:
            raise ValueError(
                "--prior-samples and --prior-root-value require --output prior-samples."
            )
        return
    if count is not None and (
        isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 10000
    ):
        raise ValueError("--prior-samples must be between 1 and 10000.")
    if root.mode == "flat":
        if root_value is None or not math.isfinite(root_value):
            raise ValueError(
                "Flat-root prior sampling requires a finite --prior-root-value."
            )
    elif root_value is not None:
        raise ValueError("--prior-root-value is only valid with a flat root.")
    if (
        getattr(args, "trait", None) is not None
        or getattr(args, "standard_error_column", None) is not None
    ):
        raise ValueError(
            "Prior samples are not conditioned on data; omit --trait and --standard-error-column."
        )
    if getattr(args, "figure_simulation_mode", None) == "conditional":
        raise ValueError(
            "Prior-only output cannot have a conditional simulation figure."
        )
    if getattr(args, "figure_tip_heatmap", "no") == "yes":
        raise ValueError("Prior-only output has no observed traits for a tip heatmap.")


def validate_branch_options(args, model):
    if model != MODEL:
        if any(getattr(args, name, None) is not None for name in BRANCH_OPTIONS):
            raise ValueError(
                "Branch-model and prior-sample options require --model BRANCH-GAUSSIAN."
            )
        return
    direct, regimes, definitions = (
        getattr(args, name, None)
        for name in ("branch_models", "branch_regimes", "regime_models")
    )
    if not (
        (direct and not regimes and not definitions)
        or (not direct and regimes and definitions)
    ):
        raise ValueError(
            "Use --branch-models alone, or --branch-regimes with --regime-models."
        )
    if any(value == "" for value in (direct, regimes, definitions)):
        raise ValueError("Branch assignment paths must not be empty.")
    root = branch_root(args)
    _validate_prior_options(args, root)
    forbidden = [
        name
        for name in (
            "sigma2",
            "compare_models",
            "tree_ensemble",
            "tree_ensemble_out",
            "bootstrap_out",
            "bootstrap_intervals_out",
            "measurement_covariance",
        )
        if getattr(args, name, None) not in (None, "")
    ]
    if forbidden:
        raise ValueError(
            "Options not defined for fixed branch assignments: "
            + ", ".join("--" + name.replace("_", "-") for name in forbidden)
        )
    if getattr(args, "seed", None) is not None and args.seed < 0:
        raise ValueError("--seed must be nonnegative.")
    output = getattr(args, "output", None)
    if output in {"likelihood", "prior-samples"}:
        unused: tuple[str, ...] = (
            "posterior_samples_out",
            "posterior_predictive_out",
            "cross_validation_out",
            "replicate_observations",
        )
        if output == "likelihood":
            unused += ("figure_out", "tree_out")
        elif getattr(args, "tree_out", None):
            raise ValueError(
                "--tree-out annotates posterior ASR results, not prior-only samples."
            )
        if any(getattr(args, name, None) not in (None, "") for name in unused):
            raise ValueError(
                f"--output {output} cannot use posterior diagnostic, tree or incompatible figure outputs."
            )
        if getattr(args, "target", "all") != "all":
            raise ValueError(f"--output {output} requires --target all.")
