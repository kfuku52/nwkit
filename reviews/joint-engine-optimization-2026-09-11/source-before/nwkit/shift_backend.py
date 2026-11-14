"""File-based adapter to the separately installed GPL kfl1ou program.

This adapter only invokes public R APIs; no inference code is vendored.
"""

import shutil
import subprocess
from pathlib import Path

from nwkit.shift_backend_probe import R_PBIC_PROBE, collect_pbic_validation
from nwkit.shift_convergence import R_CONVERGENCE

# Labels and user settings are passed as files/arguments, never interpolated R code.
R_SCRIPT = r"""
__PBIC_PREFLIGHT__
args <- commandArgs(trailingOnly=TRUE)
if (!requireNamespace("kfl1ou", quietly=TRUE) ||
    utils::packageVersion("kfl1ou") < "3.0.9") {
    stop("Install kfl1ou >= 3.0.9 in this R environment.")
}
if (identical(args[2], "pBIC")) nwkit_pbic_preflight()
options(digits=17)
tr <- ape::read.tree("tree.nwk")
y <- utils::read.delim("trait.tsv", check.names=FALSE, stringsAsFactors=FALSE)
Y <- matrix(y$value, ncol=1, dimnames=list(y$leaf_name, "trait"))
dat <- kfl1ou::adjust_data(tr, Y, normalize=FALSE, quietly=TRUE,
                          repair.tree=FALSE, drop.all.missing=FALSE,
                          drop.invariant=FALSE)
input.error <- if ("observation_variance" %in% names(y))
    matrix(y$observation_variance, ncol=1, dimnames=list(y$leaf_name, "trait")) else NULL
# All-zero known variances are exactly the no-error model. Canonicalize them
# so adding an explicit zero-SE column does not select a different optimizer.
if (!is.null(input.error) && all(input.error == 0)) input.error <- NULL
fit <- kfl1ou::estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=as.integer(args[1]), criterion=args[2],
    root.model=args[3], search.strategy=args[4],
    exhaustive.max.configurations=as.integer(args[5]),
    ensemble.seed=as.integer(args[6]), nCores=1, quietly=TRUE,
    input_error=input.error, measurement_error=FALSE)
tr <- fit$tree
# Derive clades from topology, independent of ape's edge ordering.
clades <- vector("list", length(tr$tip.label) + tr$Nnode)
for (i in seq_along(tr$tip.label)) clades[[i]] <- tr$tip.label[i]
post <- ape::reorder.phylo(tr, "postorder")
for (i in seq_len(nrow(post$edge))) {
    parent <- post$edge[i,1]; child <- post$edge[i,2]
    clades[[parent]] <- c(clades[[parent]], clades[[child]])
}
keys <- vapply(tr$edge[,2], function(i) paste(sort(clades[[i]]), collapse="/"), "")
write_tsv <- function(x, path) utils::write.table(
    x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
__CONVERGENCE_SETUP__
write_tsv(data.frame(clade=keys[fit$shift.configuration],
    mean_effect=as.numeric(fit$shift.means),
    optimum_effect=as.numeric(fit$shift.values)), "shifts.tsv")
write_tsv(data.frame(token=tr$tip.label,
    observed=as.numeric(fit$Y[tr$tip.label,1]),
    predicted=as.numeric(fit$mu[tr$tip.label,1]),
    residual=as.numeric(fit$residuals[tr$tip.label,1]),
    optimum=as.numeric(fit$optima[tr$tip.label,1])), "tips.tsv")
write_tsv(data.frame(backend_version=as.character(utils::packageVersion("kfl1ou")),
    alpha=fit$alpha, sigma2=fit$sigma2, intercept=fit$intercept,
    log_likelihood=sum(fit$logLik), score=fit$score), "model.tsv")
profile <- search.fit$profile
write_tsv(data.frame(score=profile$scores,
    clades=vapply(profile$configurations, function(edges)
        paste(keys[edges], collapse=";"), "")), "candidates.tsv")
diagnostic <- search.fit$search.diagnostics
scalar <- function(x) if (is.null(x) || length(x) != 1L) NA else x
write_tsv(data.frame(
    strategy=scalar(diagnostic$strategy),
    configuration_space_size=scalar(diagnostic$configuration.space.size),
    evaluated_configurations=scalar(diagnostic$evaluated.configurations),
    coverage=scalar(diagnostic$coverage),
    globally_optimal=scalar(diagnostic$globally.optimal),
    ensemble_attempted=scalar(diagnostic$ensemble$attempted),
    ensemble_successful=scalar(diagnostic$ensemble$successful),
    ensemble_failed=scalar(diagnostic$ensemble$failed),
    alpha_lower=scalar(search.fit$l1ou.options$alpha.lower.bound),
    alpha_upper=scalar(search.fit$l1ou.options$alpha.upper.bound)), "search.tsv")
writeLines(capture.output(dput(diagnostic)), "diagnostics.txt")
saveRDS(fit, "fit.rds")
if (length(args) >= 7L && as.integer(args[7]) > 0L) {
    boot <- if (convergence.enabled) convergent_bootstrap() else kfl1ou::l1ou_bootstrap_support(
        fit, nItrs=as.integer(args[7]), type="parametric",
        seed=as.integer(args[8]), multicore=FALSE, nCores=1, quietly=TRUE)
    write_tsv(data.frame(attempted=boot$attempted,
        successful=boot$successful, failed=boot$failed), "bootstrap.tsv")
    write_tsv(data.frame(success_index=seq_along(boot$all.shifts),
        clades=vapply(boot$all.shifts, function(edges)
        paste(keys[edges], collapse=";"), "")), "bootstrap-configurations.tsv")
    utils::write.table(data.frame(message=as.character(names(boot$failure.messages)),
        count=as.integer(boot$failure.messages)), "bootstrap-failures.tsv",
        sep="\t", quote=TRUE, qmethod="double", row.names=FALSE)
}
""".replace("__CONVERGENCE_SETUP__", R_CONVERGENCE).replace(
    "__PBIC_PREFLIGHT__", R_PBIC_PROBE
)


def run_backend(directory, args):
    executable = shutil.which(args.rscript)
    if executable is None:
        raise ValueError(
            "Rscript was not found; install R and kfl1ou >= 3.0.9, or use --rscript."
        )
    executable = str(Path(executable).resolve())
    script = Path(directory) / "runner.R"
    script.write_text(R_SCRIPT, encoding="utf-8")
    result = subprocess.run(
        [
            executable,
            "--vanilla",
            str(script),
            str(args.max_shifts),
            args.criterion,
            args.root_model,
            args.search_strategy,
            str(args.exhaustive_max_configurations),
            str(args.seed),
            str(args.bootstrap),
            str(args.bootstrap_seed),
            "TRUE" if args.convergence else "FALSE",
        ],
        cwd=directory,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if result.returncode:
        raise RuntimeError(
            f"kfl1ou failed (exit {result.returncode}):\n{result.stderr[-8000:]}"
        )
    return {
        "rscript": executable,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "pbic_validation": collect_pbic_validation(directory)
        if args.criterion == "pBIC"
        else {"status": "not_applicable", "criterion": args.criterion},
    }
