"""Public kfl1ou fixed-configuration fits for a small enumerated reference."""

from nwkit.shift_backend_probe import R_PBIC_PROBE

R_SCRIPT = r"""
__PBIC_PREFLIGHT__
args <- commandArgs(trailingOnly=TRUE)
if (!requireNamespace("kfl1ou", quietly=TRUE) || utils::packageVersion("kfl1ou") < "3.0.9")
    stop("Install kfl1ou >= 3.0.9.")
options(digits=17)
baseline <- readRDS("two-stage.rds")
initial <- attr(baseline, "nwkit.unconstrained.fit")
tr <- initial$tree
Y <- initial$Y
opt <- initial$l1ou.options
if (identical(opt$criterion, "pBIC")) nwkit_pbic_preflight()
opt$use.saved.scores <- FALSE
opt$fixed.alpha <- FALSE
opt$compute.hessian <- FALSE
clades <- vector("list", length(tr$tip.label) + tr$Nnode)
for (i in seq_along(tr$tip.label)) clades[[i]] <- tr$tip.label[i]
post <- ape::reorder.phylo(tr, "postorder")
for (i in seq_len(nrow(post$edge))) {
    parent <- post$edge[i,1]; child <- post$edge[i,2]
    clades[[parent]] <- c(clades[[parent]], clades[[child]])
}
keys <- vapply(tr$edge[,2], function(i) paste(sort(clades[[i]]), collapse="/"), "")
map <- read.delim("branches.tsv", stringsAsFactors=FALSE)
edge <- match(map$clade, keys)
if (anyNA(edge)) stop("Branch mapping failed.")
edge <- setNames(edge, as.character(map$branch_id))
edge <- c("0"=0L, edge)
candidates <- read.delim("candidates.tsv", stringsAsFactors=FALSE, colClasses="character")
parse_ids <- function(value) if (!nzchar(value)) integer() else as.integer(strsplit(value, ";", fixed=TRUE)[[1]])
write_tsv <- function(x, path) write.table(x, path, sep="\t", quote=TRUE, qmethod="double", row.names=FALSE, na="NA")
anchor <- tryCatch({
    refit <- kfl1ou::fit_OU(tr, Y, unname(baseline$shift.configuration),
        cr.regimes=baseline$convergent.regimes, l1ou.options=opt)
    data.frame(score=refit$score, log_likelihood=sum(refit$logLik), error="")
}, error=function(e) data.frame(score=NA_real_, log_likelihood=NA_real_, error=conditionMessage(e)))
write_tsv(anchor, "baseline-refit.tsv")
best <- NULL; best.id <- NA_integer_
records <- vector("list", nrow(candidates))
for (i in seq_len(nrow(candidates))) {
    row <- candidates[i,]
    selected <- unname(edge[as.character(parse_ids(row$shifts))])
    groups <- lapply(strsplit(row$groups, "|", fixed=TRUE)[[1]], function(g)
        unname(edge[as.character(parse_ids(g))]))
    warnings <- character()
    record <- tryCatch(withCallingHandlers({
        fit <- kfl1ou::fit_OU(tr, Y, shift.configuration=selected,
            cr.regimes=groups, l1ou.options=opt)
        if (!identical(fit$tree$edge, tr$edge) || !identical(fit$tree$tip.label, tr$tip.label))
            stop("Candidate changed tree coordinates.")
        if (length(fit$alpha) != 1L || !is.finite(fit$alpha) || fit$alpha <= 0)
            stop("Candidate has unidentifiable optima.")
        if (!all(is.finite(c(fit$score, fit$logLik, fit$sigma2))) || fit$sigma2 <= 0)
            stop("Candidate has nonfinite score or invalid variance.")
        free.score <- free.ll <- NA_real_
        if (length(groups) == length(selected) + 1L) {
            free <- kfl1ou::fit_OU(tr, Y, shift.configuration=selected, l1ou.options=opt)
            free.score <- free$score; free.ll <- sum(free$logLik)
        }
        if (is.null(best) || fit$score < best$score) { best <- fit; best.id <- as.integer(row$candidate_id) }
        data.frame(candidate_id=row$candidate_id, status="completed", score=fit$score,
            log_likelihood=sum(fit$logLik), alpha=fit$alpha, sigma2=fit$sigma2,
            free_score=free.score, free_log_likelihood=free.ll, error="", warnings="")
    }, warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error=function(e) data.frame(candidate_id=row$candidate_id, status="failed", score=NA_real_,
        log_likelihood=NA_real_, alpha=NA_real_, sigma2=NA_real_, free_score=NA_real_, free_log_likelihood=NA_real_, error=conditionMessage(e), warnings=""))
    record$warnings <- paste(unique(warnings), collapse=" | ")
    records[[i]] <- record
    cat("candidate", i, "of", nrow(candidates), record$status, "\n")
}
write_tsv(do.call(rbind, records), "candidate-results.tsv")
if (is.null(best)) stop("No joint candidate returned an identifiable finite fit.")
saveRDS(best, "joint.rds")
write_tsv(data.frame(candidate_id=best.id, score=best$score, log_likelihood=sum(best$logLik),
    alpha=best$alpha, sigma2=best$sigma2, intercept=best$intercept,
    backend_version=as.character(utils::packageVersion("kfl1ou"))), "joint.tsv")
write_tsv(data.frame(token=tr$tip.label, observed=as.numeric(best$Y[,1]),
    predicted=as.numeric(best$mu[,1]), optimum=as.numeric(best$optima[,1])), "joint-tips.tsv")
write_tsv(data.frame(branch_id=as.integer(names(edge)[match(unname(best$shift.configuration), edge)]),
    mean_effect=as.numeric(best$shift.means), optimum_effect=as.numeric(best$shift.values)), "joint-effects.tsv")
""".replace("__PBIC_PREFLIGHT__", R_PBIC_PROBE)
