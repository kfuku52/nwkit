#!/usr/bin/env Rscript
# Input is one complete post-burnin sample file per independent chain.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Expected output TSV and at least two chains")
if (!requireNamespace("posterior", quietly = TRUE)) stop("R package posterior is required")
tables <- lapply(args[-1L], function(p) read.table(p, header = TRUE, check.names = FALSE, comment.char = "", quote = ""))
if (!all(vapply(tables, function(x) identical(names(x), names(tables[[1L]])), logical(1)))) stop("Different sample headers")
parameters <- names(tables[[1L]])
if (length(parameters) < 2L || parameters[[1L]] != "Gen" ||
    anyDuplicated(parameters) || any(!nzchar(parameters))) {
  stop("Invalid sample headers: expected unique Gen and parameter columns")
}
if (length(unique(vapply(tables, nrow, integer(1)))) != 1L) stop("Unequal chain lengths")
if (nrow(tables[[1L]]) < 8L) stop("Insufficient samples")
rows <- lapply(seq.int(2L, ncol(tables[[1L]])), function(j) {
  x <- do.call(cbind, lapply(tables, function(t) t[[j]]))
  if (!is.numeric(x) || any(!is.finite(x))) stop(paste("Nonfinite or nonnumeric samples:", parameters[j]))
  data.frame(parameter = names(tables[[1L]])[j], rhat = posterior::rhat(x),
    ess_bulk = posterior::ess_bulk(x), ess_tail = posterior::ess_tail(x),
    mcse_mean = posterior::mcse_mean(x),
    mcse_q025 = posterior::mcse_quantile(x, probs = 0.025),
    mcse_q975 = posterior::mcse_quantile(x, probs = 0.975))
})
write.table(do.call(rbind, rows), args[1L], sep = "\t", row.names = FALSE, quote = FALSE)

print(sessionInfo())
