#!/usr/bin/env Rscript
# Input is one complete post-burnin sample file per independent chain.
args <- commandArgs(trailingOnly = TRUE)
workers <- 1L
thresholds <- NULL
if (length(args) >= 2L && identical(args[[1L]], "--jobs")) {
  workers <- suppressWarnings(as.integer(args[[2L]]))
  if (is.na(workers) || workers < 1L) stop("Expected a positive worker count")
  args <- args[-c(1L, 2L)]
}
if (length(args) >= 3L && identical(args[[1L]], "--stop-at")) {
  thresholds <- suppressWarnings(as.numeric(args[2:3]))
  if (any(!is.finite(thresholds) | thresholds <= 0)) stop("Expected positive stopping thresholds")
  args <- args[-(1:3)]
}
if (length(args) < 3L) stop("Expected output TSV and at least two chains")
if (!requireNamespace("posterior", quietly = TRUE)) stop("R package posterior is required")
tables <- lapply(args[-1L], function(p) {
  read.table(p, header = TRUE, check.names = FALSE, comment.char = "", quote = "", colClasses = "numeric")
})
if (!all(vapply(tables, function(x) identical(names(x), names(tables[[1L]])), logical(1)))) stop("Different sample headers")
parameters <- names(tables[[1L]])
if (length(parameters) < 2L || parameters[[1L]] != "Gen" ||
    anyDuplicated(parameters) || any(!nzchar(parameters))) {
  stop("Invalid sample headers: expected unique Gen and parameter columns")
}
if (length(unique(vapply(tables, nrow, integer(1)))) != 1L) stop("Unequal chain lengths")
if (nrow(tables[[1L]]) < 8L) stop("Insufficient samples")
columns <- seq.int(2L, length(parameters))
workers <- if (.Platform$OS.type == "windows") 1L else min(workers, length(columns))

by_parameter <- function(summarize) {
  if (workers == 1L) return(lapply(columns, summarize))
  # Read-only chain tables are shared by forked workers. Preserve column order.
  rows <- parallel::mclapply(columns, summarize, mc.cores = workers, mc.set.seed = FALSE)
  failed <- which(!vapply(rows, is.data.frame, logical(1)))
  if (length(failed)) {
    failed_parameter <- parameters[columns[failed[[1L]]]]
    message <- rows[[failed[[1L]]]]
    if (is.null(message)) message <- "worker returned no result"
    stop(sprintf("Diagnostic worker failed for %s: %s", failed_parameter, as.character(message)))
  }
  rows
}

chain_matrix <- function(j) do.call(cbind, lapply(tables, function(t) t[[j]]))

mcse <- function(x) data.frame(
  mcse_mean = posterior::mcse_mean(x),
  mcse_q025 = posterior::mcse_quantile(x, probs = 0.025),
  mcse_q975 = posterior::mcse_quantile(x, probs = 0.975)
)
empty_mcse <- data.frame(mcse_mean = NA_real_, mcse_q025 = NA_real_, mcse_q975 = NA_real_)

result <- do.call(rbind, by_parameter(function(j) {
  x <- chain_matrix(j)
  if (any(!is.finite(x))) stop(paste("Nonfinite or nonnumeric samples:", parameters[j]))
  measures <- data.frame(parameter = parameters[j], rhat = posterior::rhat(x),
    ess_bulk = posterior::ess_bulk(x), ess_tail = posterior::ess_tail(x))
  cbind(measures, if (is.null(thresholds)) mcse(x) else empty_mcse)
}))

if (!is.null(thresholds)) {
  # Allow for decimal rounding in the TSV used by the stopping decision.
  slack <- 128 * .Machine$double.eps
  may_stop <- with(result, rhat <= thresholds[1L] * (1 + slack) &
    ess_bulk >= thresholds[2L] * (1 - slack) & ess_tail >= thresholds[2L] * (1 - slack))
  if (isTRUE(all(may_stop))) {
    result[names(empty_mcse)] <- do.call(rbind, by_parameter(function(j) mcse(chain_matrix(j))))
  }
}
write.table(result, args[1L], sep = "\t", row.names = FALSE, quote = FALSE)

cat("Posterior diagnostic workers:", workers, "\n")
print(sessionInfo())
