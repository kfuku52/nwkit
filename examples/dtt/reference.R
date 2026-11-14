# Run from repository root. Reference: geiger 2.0.11; no NWKIT/R RNG equivalence assumed.
suppressPackageStartupMessages(library(geiger))
tree <- read.tree("examples/dtt/tree.nwk")
x <- as.matrix(read.delim("examples/dtt/traits.tsv", row.names=1)[, c("size", "shape")])
options(digits=17)
reference <- dtt(tree, x, index="avg.sq", nsim=0, plot=FALSE)
print(reference)
print(ratematrix(tree, x))
# Keep the first (pre-split) root point, combine simultaneous post-split events,
# then extend the zero tail to the present. These are NWKIT's time rows.
keep <- c(TRUE, !duplicated(reference$times[-1], fromLast=TRUE))
times <- c(reference$times[keep], 1)
values <- c(reference$dtt[keep], 0)
print(data.frame(relative_time=times, relative_disparity=values))
# If the example NWKIT tables exist, independently verify DTT and MDI integration.
if (file.exists("/tmp/dtt.tsv") && file.exists("/tmp/dtt-summary.tsv")) {
  nw <- read.delim("/tmp/dtt.tsv")
  summary <- read.delim("/tmp/dtt-summary.tsv")
  stopifnot(isTRUE(all.equal(nw$relative_time, unname(times), tolerance=1e-10)))
  stopifnot(isTRUE(all.equal(nw$relative_disparity, unname(values), tolerance=1e-10)))
  delta <- nw$relative_disparity - nw$bm_median
  mdi <- sum(diff(nw$relative_time) * (head(delta,-1) + tail(delta,-1)) / 2)
  stopifnot(abs(mdi - summary$mdi) < 1e-10)
  cat("NWKIT observed DTT and MDI quadrature agree.\n")
}
