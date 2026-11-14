# Run from the repository root. Independent reference used by Python tests.
# Verified with Rphylopars 0.3.10 and ape; neither is a NWKIT dependency.
library(Rphylopars)
a <- read.delim("examples/individual_asr/individuals.tsv",
                colClasses=c("character", "character", "character", "numeric"))
a <- reshape(a, direction="wide", timevar="trait",
             idvar=c("leaf_name", "individual_id"))
a <- a[, c("leaf_name", "value.x", "value.y")]
names(a) <- c("species", "x", "y")
tree <- ape::read.tree("examples/individual_asr/tree.nwk")
for (reml in c(TRUE, FALSE)) {
  for (full in c(TRUE, FALSE)) {
    fit <- phylopars(a, tree, model="BM", pheno_error=TRUE,
                    pheno_correlated=full, REML=reml,
                    repeat_optim_limit=3, repeat_optim_tol=1e-9)
    # Diagonal-W ML is close to a singular evolutionary covariance here.
    # Refine Rphylopars' stopping tolerance with warm-started optimizations.
    if (!reml && !full) {
      for (restart in 1:5) {
        fit <- phylopars(a, tree, model="BM", pheno_error=TRUE,
                        pheno_correlated=full, REML=reml, skip_EM=TRUE,
                        phylocov_start=fit$pars$phylocov,
                        phenocov_start=fit$pars$phenocov,
                        repeat_optim_limit=10, repeat_optim_tol=1e-14)
      }
    }
    cat("\nREML:", reml, "full W:", full, "\n")
    print(fit$pars, digits=13)
    print(fit$logLik, digits=13)
    print(fit$anc_recon, digits=13)
  }
}
