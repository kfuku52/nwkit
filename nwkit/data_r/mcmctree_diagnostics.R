#!/usr/bin/env Rscript
# Input is one complete post-burnin sample file per independent chain.
args <- commandArgs(trailingOnly = TRUE)
workers <- 1L
thresholds <- NULL
candidates <- integer()
confirm <- 0L
final <- FALSE
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
if (length(args) >= 4L && identical(args[[1L]], "--try-burnin")) {
  candidates <- suppressWarnings(as.integer(strsplit(args[[2L]], ",", fixed = TRUE)[[1L]]))
  confirm <- suppressWarnings(as.integer(args[[3L]]))
  final <- args[[4L]] == "1"
  if (is.null(thresholds) || !length(candidates) ||
      any(is.na(candidates) | candidates < 1L) || anyDuplicated(candidates) ||
      is.na(confirm) || !(confirm %in% c(0L, candidates))) {
    stop("Invalid burn-in candidate counts")
  }
  args <- args[-(1:4)]
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
draws <- nrow(tables[[1L]])
if (draws < 8L || any(candidates > draws - 8L)) stop("Insufficient samples")
columns <- seq.int(2L, length(parameters))
workers <- if (.Platform$OS.type == "windows") 1L else min(workers, length(columns))

by_parameter <- function(summarize, indices = columns) {
  cores <- min(workers, length(indices))
  if (cores <= 1L) return(lapply(indices, summarize))
  # Read-only chain tables are shared by forked workers. Preserve column order.
  rows <- parallel::mclapply(indices, summarize, mc.cores = cores, mc.set.seed = FALSE)
  failed <- which(!vapply(rows, is.data.frame, logical(1)))
  if (length(failed)) {
    failed_parameter <- parameters[indices[failed[[1L]]]]
    message <- rows[[failed[[1L]]]]
    if (is.null(message)) message <- "worker returned no result"
    stop(sprintf("Diagnostic worker failed for %s: %s", failed_parameter, as.character(message)))
  }
  rows
}

window_rows <- function(discard) {
  if (discard == 0L) NULL else seq.int(discard + 1L, draws)
}

chain_matrix <- function(j, rows) {
  if (is.null(rows)) return(do.call(cbind, lapply(tables, function(t) t[[j]])))
  do.call(cbind, lapply(tables, function(t) t[[j]][rows]))
}

mcse <- function(x) data.frame(
  mcse_mean = posterior::mcse_mean(x),
  mcse_q025 = posterior::mcse_quantile(x, probs = 0.025),
  mcse_q975 = posterior::mcse_quantile(x, probs = 0.975)
)
empty_mcse <- data.frame(mcse_mean = NA_real_, mcse_q025 = NA_real_, mcse_q975 = NA_real_)

summarize_window <- function(discard, indices = columns) {
  rows <- window_rows(discard)
  do.call(rbind, by_parameter(function(j) {
    x <- chain_matrix(j, rows)
    if (discard == 0L && any(!is.finite(x))) {
      stop(paste("Nonfinite or nonnumeric samples:", parameters[j]))
    }
    measures <- data.frame(parameter = parameters[j], rhat = posterior::rhat(x),
      ess_bulk = posterior::ess_bulk(x), ess_tail = posterior::ess_tail(x))
    cbind(measures, if (is.null(thresholds)) mcse(x) else empty_mcse)
  }, indices))
}

valid <- function(result) {
  with(result, all(is.finite(rhat) & rhat > 0 &
    is.finite(ess_bulk) & ess_bulk > 0 & is.finite(ess_tail) & ess_tail > 0))
}

passes <- function(result) {
  valid(result) && with(result, all(rhat < thresholds[1L] &
    ess_bulk >= thresholds[2L] & ess_tail >= thresholds[2L]))
}

complete_mcse <- function(result, discard) {
  rows <- window_rows(discard)
  result[names(empty_mcse)] <- do.call(rbind, by_parameter(function(j) mcse(chain_matrix(j, rows))))
  result
}

result <- summarize_window(0L)

if (length(candidates)) {
  selected <- 0L
  full_passes <- passes(result)
  if (!full_passes && valid(result)) {
    # A single failed parameter rejects a window. Test the limiting parameter
    # of the full trace first and reuse its row when the window passes.
    distance <- with(result, pmax(rhat / thresholds[1L], thresholds[2L] / pmin(ess_bulk, ess_tail)))
    limiting <- columns[[which.max(distance)]]
    for (discard in candidates) {
      first <- summarize_window(discard, limiting)
      if (!passes(first)) next
      candidate <- rbind(first, summarize_window(discard, columns[columns != limiting]))
      if (passes(candidate)) {
        candidate <- candidate[match(result$parameter, candidate$parameter), , drop = FALSE]
        selected <- discard
        break
      }
    }
  }
  if (full_passes || (final && selected == 0L)) result <- complete_mcse(result, 0L)
  result$discarded_draws <- 0L
  if (selected > 0L) {
    if (selected == confirm) candidate <- complete_mcse(candidate, selected)
    candidate$discarded_draws <- selected
    result <- rbind(result, candidate)
  }
  # Keep the strict crossing decision identical in Python and R.
  diagnostic_columns <- c("rhat", "ess_bulk", "ess_tail")
  result[diagnostic_columns] <- lapply(result[diagnostic_columns], function(x) sprintf("%.17g", x))
} else if (!is.null(thresholds)) {
  # Allow for decimal rounding in the TSV used by the stopping decision.
  slack <- 128 * .Machine$double.eps
  may_stop <- with(result, rhat <= thresholds[1L] * (1 + slack) &
    ess_bulk >= thresholds[2L] * (1 - slack) & ess_tail >= thresholds[2L] * (1 - slack))
  if (isTRUE(all(may_stop))) result <- complete_mcse(result, 0L)
}
write.table(result, args[1L], sep = "\t", row.names = FALSE, quote = FALSE)

cat("Posterior diagnostic workers:", workers, "\n")
print(sessionInfo())
