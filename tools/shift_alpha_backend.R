# Public-API research driver. BIC() uses kfl1ou's public logLik/nobs methods.
# No inference implementation from the external GPL backend is copied here.
args <- commandArgs(trailingOnly=TRUE)
if (!requireNamespace("kfl1ou", quietly=TRUE)) stop("Install the corrected kfl1ou backend.")
if (!exists("nwkit_pbic_preflight", mode="function"))
    stop("Generate this driver through validate_shift_alpha.py to include the pBIC capability contract.")
options(digits=17)
write_tsv <- function(x, path) utils::write.table(x, path, sep="\t", quote=TRUE,
    qmethod="double", row.names=FALSE, na="NA")
parse_ids <- function(x) if (!nzchar(x)) integer() else as.integer(strsplit(x, ";", fixed=TRUE)[[1]])
compact <- function(x) paste(x, collapse=";")

if (identical(args[[1]], "probe")) {
    nwkit_pbic_preflight(dirname(args[[2]]))
    write_tsv(data.frame(version=as.character(utils::packageVersion("kfl1ou")),
        library=normalizePath(find.package("kfl1ou")),
        contract="ou-optimum-information-v1"), args[[2]])
    quit(status=0)
}

nwkit_pbic_preflight()
tr <- ape::read.tree("tree.nwk")
y <- read.delim("traits.tsv", stringsAsFactors=FALSE, check.names=FALSE)
Y <- matrix(y$value, ncol=1, dimnames=list(y$leaf_name, "trait"))
dat <- kfl1ou::adjust_data(tr, Y, normalize=FALSE, quietly=TRUE,
    repair.tree=FALSE, drop.all.missing=FALSE, drop.invariant=FALSE)
tr <- dat$tree; Y <- dat$Y
input.error <- matrix(y$se^2, ncol=1, dimnames=list(y$leaf_name, "trait"))
if (all(input.error == 0)) input.error <- NULL
settings <- read.delim("settings.tsv", stringsAsFactors=FALSE)
candidates <- read.delim("candidates.tsv", colClasses="character")
map <- read.delim("branches.tsv", stringsAsFactors=FALSE)
clades <- vector("list", length(tr$tip.label) + tr$Nnode)
for (i in seq_along(tr$tip.label)) clades[[i]] <- tr$tip.label[i]
post <- ape::reorder.phylo(tr, "postorder")
for (i in seq_len(nrow(post$edge))) {
    parent <- post$edge[i,1]; child <- post$edge[i,2]
    clades[[parent]] <- c(clades[[parent]], clades[[child]])
}
keys <- vapply(tr$edge[,2], function(i) paste(sort(clades[[i]]), collapse="/"), "")
edge <- match(map$clade, keys)
if (anyNA(edge)) stop("Cannot map branch coordinates.")
edge <- c("0"=0L, setNames(edge, as.character(map$branch_id)))
to_nwkit <- function(x) as.integer(names(edge)[match(as.integer(x), edge)])
dir.create("models")
records <- list(); ledger <- list()
model_row <- function(fit, setting, method, criterion, warnings="") {
    label <- paste(setting$floor_id, method, criterion, sep="-")
    selected <- to_nwkit(fit$shift.configuration)
    groups <- lapply(fit$convergent.regimes, to_nwkit)
    if (is.null(fit$convergent.regimes)) groups <- as.list(c(0L, selected))
    saveRDS(fit, file.path("models", paste0(label, ".rds")))
    write_tsv(data.frame(leaf_name=tr$tip.label, predicted=as.numeric(fit$mu[,1]),
        optimum=as.numeric(fit$optima[,1])), file.path("models", paste0(label, "-tips.tsv")))
    write_tsv(data.frame(branch_id=selected, mean_effect=as.numeric(fit$shift.means)),
        file.path("models", paste0(label, "-effects.tsv")))
    data.frame(floor_id=setting$floor_id, method=method, criterion=criterion,
        status="completed", score=if(criterion=="BIC") stats::BIC(fit) else fit$score,
        reported_score=fit$score, alpha=fit$alpha, sigma2=fit$sigma2,
        intercept=fit$intercept, log_likelihood=sum(fit$logLik),
        shift_branch_ids=compact(selected), groups=paste(vapply(groups, compact, ""),collapse="|"),
        error="", warnings=warnings)
}
failed_row <- function(setting, method, criterion, error, warnings="") data.frame(
    floor_id=setting$floor_id, method=method, criterion=criterion, status="failed",
    score=NA_real_, reported_score=NA_real_, alpha=NA_real_, sigma2=NA_real_, intercept=NA_real_,
    log_likelihood=NA_real_, shift_branch_ids="", groups="", error=error, warnings=warnings)

for (j in seq_len(nrow(settings))) {
    setting <- settings[j,]
    for (criterion in c("BIC", "pBIC")) {
        warnings <- character()
        row <- tryCatch(withCallingHandlers({
            initial <- kfl1ou::estimate_shift_configuration(tr, Y,
                max.nShifts=2L, criterion=criterion, root.model=setting$root_model,
                search.strategy="exhaustive", exhaustive.max.configurations=5000L,
                alpha.lower=setting$lower, alpha.upper=setting$upper,
                alpha.starting.value=setting$starting, input_error=input.error,
                measurement_error=FALSE, optimizer.starts=1L, compute.hessian=FALSE,
                nCores=1L, quietly=TRUE)
            fit <- kfl1ou::estimate_convergent_regimes(initial, criterion=criterion,
                method="backward", fixed.alpha=FALSE, nCores=1L)
            model_row(fit, setting, "two_stage", criterion)
        }, warning=function(w){warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
        error=function(e) failed_row(setting, "two_stage", criterion, conditionMessage(e)))
        row$warnings <- paste(unique(warnings), collapse=" | ")
        records[[length(records)+1L]] <- row
    }
    best <- list(BIC=NULL, pBIC=NULL)
    for (i in seq_len(nrow(candidates))) {
        candidate <- candidates[i,]
        selected <- unname(edge[as.character(parse_ids(candidate$shifts))])
        groups <- lapply(strsplit(candidate$groups, "|", fixed=TRUE)[[1]], function(g)
            unname(edge[as.character(parse_ids(g))]))
        warnings <- character()
        row <- tryCatch(withCallingHandlers({
            fit <- kfl1ou::fit_OU(tr, Y, selected, cr.regimes=groups,
                criterion="pBIC", root.model=setting$root_model,
                alpha.lower=setting$lower, alpha.upper=setting$upper,
                alpha.starting.value=setting$starting, input_error=input.error,
                measurement_error=FALSE, optimizer.starts=1L, compute.hessian=FALSE,
                search.max.nShifts=2L)
            scores <- c(BIC=stats::BIC(fit), pBIC=fit$score)
            if (any(!is.finite(c(scores, fit$alpha, fit$sigma2, fit$logLik))) || fit$alpha<=0 || fit$sigma2<=0)
                stop("Nonfinite or unidentifiable candidate fit.")
            if (!identical(fit$tree$edge, tr$edge) || !identical(fit$tree$tip.label, tr$tip.label))
                stop("Candidate changed tree coordinates.")
            for (criterion in names(scores)) {
                if (is.null(best[[criterion]]) || scores[[criterion]] < best[[criterion]]$score)
                    best[[criterion]] <- list(fit=fit, score=scores[[criterion]], candidate_id=as.integer(candidate$candidate_id))
            }
            data.frame(floor_id=setting$floor_id, candidate_id=candidate$candidate_id,
                status="completed", BIC=scores[["BIC"]], pBIC=scores[["pBIC"]],
                alpha=fit$alpha, sigma2=fit$sigma2, log_likelihood=sum(fit$logLik),
                singleton_free_score=if(length(groups)==length(selected)+1L) fit$unconstrained.score else NA_real_,
                error="", warnings="")
        }, warning=function(w){warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
        error=function(e) data.frame(floor_id=setting$floor_id, candidate_id=candidate$candidate_id,
            status="failed", BIC=NA_real_, pBIC=NA_real_, alpha=NA_real_, sigma2=NA_real_,
            log_likelihood=NA_real_, singleton_free_score=NA_real_, error=conditionMessage(e),warnings=""))
        row$warnings <- paste(unique(warnings), collapse=" | ")
        ledger[[length(ledger)+1L]] <- row
    }
    for (criterion in names(best)) {
        row <- if (is.null(best[[criterion]])) failed_row(setting, "joint", criterion, "No successful candidate") else
            model_row(best[[criterion]]$fit, setting, "joint", criterion)
        records[[length(records)+1L]] <- row
    }
    write_tsv(do.call(rbind, ledger), "candidates-results.tsv")
    write_tsv(do.call(rbind, records), "models.tsv")
    cat(setting$floor_id, "completed\n")
}
