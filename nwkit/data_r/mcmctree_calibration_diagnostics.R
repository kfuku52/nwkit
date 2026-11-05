#!/usr/bin/env Rscript
# Diagnostics for separate calibration experiments, using posterior's tested
# rank-normalized/folded split R-hat and bulk/tail ESS implementation.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: script summary.tsv chain1.tsv chain2.tsv [chainN.tsv]")
if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Package posterior is required. Rebuild the GeneGalleon container with current dependencies.")
}
chains <- lapply(args[-1L], function(path) {
  read.table(path, header = TRUE, check.names = FALSE, comment.char = "", quote = "")
})
if (any(vapply(chains, nrow, integer(1)) < 8L)) stop("At least 8 retained draws per chain are required.")
if (length(unique(vapply(chains, nrow, integer(1)))) != 1L) stop("Chain lengths differ.")
if (!all(vapply(chains, function(x) identical(names(x), names(chains[[1L]])), logical(1)))) {
  stop("Chain parameter columns differ.")
}
parameters <- setdiff(names(chains[[1L]]), "Gen")
if (!length(parameters) || anyDuplicated(parameters)) stop("Invalid parameter columns.")
rows <- lapply(parameters, function(parameter) {
  x <- vapply(chains, function(chain) as.numeric(chain[[parameter]]), numeric(nrow(chains[[1L]])))
  if (any(!is.finite(x))) stop(paste("Nonfinite samples:", parameter))
  limits <- quantile(as.vector(x), c(0.025, 0.5, 0.975), names = FALSE)
  data.frame(parameter = parameter, mean = mean(x), median = limits[2L],
             eti_low = limits[1L], eti_high = limits[3L], sd = sd(as.vector(x)),
             rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
             ess_tail = posterior::ess_tail(x), mcse_mean = posterior::mcse_mean(x))
})
write.table(do.call(rbind, rows), args[1L], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
