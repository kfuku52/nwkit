"""Behavioral contract for the external, single-trait pBIC implementation.

The R fragment uses public fits and an independent dense Gaussian calculation.
It must execute inside the same R process and namespace as the requested fit.
Package versions alone cannot distinguish the original and corrected 3.0.9.
"""

import csv
import hashlib
import math
from pathlib import Path

PBIC_CONTRACT = "ou-optimum-information-v1"
PBIC_TOLERANCE = 2e-6

R_PBIC_PROBE = r"""
nwkit_pbic_preflight <- function(directory=".") {
    seed.existed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    saved.seed <- if (seed.existed) get(".Random.seed", envir=.GlobalEnv) else NULL
    saved.kind <- RNGkind()
    on.exit({
        do.call(RNGkind, as.list(saved.kind))
        if (seed.existed) assign(".Random.seed", saved.seed, envir=.GlobalEnv) else
            if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))
                rm(".Random.seed", envir=.GlobalEnv)
    }, add=TRUE)
    fail <- function(message) stop(paste0(
        "kfl1ou pBIC capability check failed: ", message,
        ". Install a corrected backend; original 3.0.9 is incompatible. ",
        "No alternative criterion has been substituted."), call.=FALSE)
    if (!requireNamespace("kfl1ou", quietly=TRUE)) fail("package unavailable")
    probe <- function() {
        tr <- ape::read.tree(text=paste0(
            "(((a:1,b:1):1,(c:1,d:1):1):1,",
            "((e:1,f:1):1,(g:1,h:1):1):1);"))
        tr <- ape::reorder.phylo(tr, "postorder")
        names <- letters[1:8]
        y <- matrix(c(1.55,1.31,-.21,-.24,1.70,1.85,.72,.25), ncol=1,
                    dimnames=list(names, "trait"))
        shared.y <- matrix(c(.65,-.90,-.09,.75,-.89,.45,-.73,-.26), ncol=1,
                           dimnames=list(names, "trait"))
        edges <- match(match(c("a", "e"), tr$tip.label), tr$edge[,2])
        rows <- list()
        check <- function(label, actual, expected) {
            if (length(actual) != 1L || length(expected) != 1L ||
                !all(is.finite(c(actual, expected))) ||
                abs(actual - expected) > 2e-6) fail(label)
            rows[[length(rows)+1L]] <<- data.frame(
                check=label, actual=actual, expected=expected,
                absolute_error=abs(actual-expected))
        }
        # These two disconnected terminal shifts have one unit of exposure.
        # X maps absolute regime optima to tip means; no backend design or
        # information-matrix helper participates in this reference calculation.
        dense <- function(fit, shared, fixed, root) {
            labels <- fit$tree$tip.label
            alpha <- as.numeric(fit$alpha)
            variance <- as.numeric(fit$sigma2)
            if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 ||
                length(variance) != 1L || !is.finite(variance) || variance <= 0)
                fail("invalid fitted covariance parameters")
            D <- ape::cophenetic.phylo(fit$tree)[labels, labels]
            V <- variance / (2*alpha) * exp(-alpha*D)
            if (root == "OUfixedRoot")
                V <- V * (-expm1(-alpha*(6-D)))
            A <- cbind(as.numeric(labels == "a"), as.numeric(labels == "e"))
            exposure <- -expm1(-alpha) * A
            X <- if (shared) cbind(1-rowSums(exposure), rowSums(exposure)) else
                cbind(1-rowSums(exposure), exposure)
            L <- chol(V)
            design <- forwardsolve(t(L), X)
            observed <- as.numeric(fit$Y[labels,1])
            response <- forwardsolve(t(L), observed)
            qr.fit <- qr(design)
            if (qr.fit$rank != ncol(X)) fail("rank deficient reference fixture")
            residual <- qr.resid(qr.fit, response)
            ll <- -.5*(8*log(2*pi)+2*sum(log(diag(L)))+sum(residual^2))
            information <- 2*sum(log(abs(diag(qr.R(qr.fit))))) +
                ncol(X)*log(stats::var(observed))
            # Two locations, one process variance and optionally estimated alpha.
            score <- -2*ll + 4*log(13) + (2-as.integer(fixed))*log(8) + information
            c(log_likelihood=ll, score=score)
        }
        for (root in c("OUfixedRoot", "OUrandomRoot")) {
            public.fit <- function(shared=FALSE, fixed=FALSE, singleton=FALSE) {
                groups <- if (shared) list(0L, edges) else
                    if (singleton) as.list(c(0L, edges)) else NULL
                suppressWarnings(kfl1ou::fit_OU(
                    tr, if (shared) shared.y else y, edges,
                    cr.regimes=groups, criterion="pBIC", root.model=root,
                    alpha.lower=if (fixed) .4 else .03,
                    alpha.upper=if (fixed) .4 else 3,
                    compute.hessian=FALSE, optimizer.starts=1L))
            }
            for (fixed in c(TRUE, FALSE)) {
                free <- public.fit(fixed=fixed)
                singleton <- public.fit(fixed=fixed, singleton=TRUE)
                tag <- paste(root, if (fixed) "fixed" else "estimated", sep="/")
                for (representation in c("free", "singleton")) {
                    fitted <- if (representation == "free") free else singleton
                    reference <- dense(fitted, FALSE, fixed, root)
                    check(paste(tag, representation, "likelihood", sep="/"),
                          sum(fitted$logLik), reference[["log_likelihood"]])
                    check(paste(tag, representation, "pBIC", sep="/"),
                          fitted$score, reference[["score"]])
                }
                # Fixed-alpha fits share the exact same model and covariance.
                # Estimated-alpha fits are each checked at their own fitted alpha.
                if (fixed) check(paste(tag, "coordinate_equivalence", sep="/"),
                                 free$score, singleton$score)
            }
            merged <- public.fit(shared=TRUE)
            reference <- dense(merged, TRUE, FALSE, root)
            check(paste(root, "shared/likelihood", sep="/"),
                  sum(merged$logLik), reference[["log_likelihood"]])
            check(paste(root, "shared/pBIC", sep="/"),
                  merged$score, reference[["score"]])
        }
        library <- normalizePath(find.package("kfl1ou"), mustWork=TRUE)
        files <- sort(list.files(library, recursive=TRUE, full.names=TRUE,
                                 all.files=TRUE, no..=TRUE))
        files <- files[!file.info(files)$isdir]
        hashes <- tools::md5sum(files)
        if (anyNA(hashes)) fail("cannot fingerprint installed package")
        write_tsv <- function(value, name) utils::write.table(value,
            file.path(directory, name), sep="\t", quote=TRUE, qmethod="double",
            row.names=FALSE, na="NA")
        write_tsv(do.call(rbind, rows), "pbic-checks.tsv")
        write_tsv(data.frame(contract="ou-optimum-information-v1",
            version=as.character(utils::packageVersion("kfl1ou")), library=library),
            "pbic-identity.tsv")
        write_tsv(data.frame(path=substring(files, nchar(library)+2L),
                            md5=unname(hashes)), "pbic-files.tsv")
        invisible(TRUE)
    }
    tryCatch(probe(), error=function(e) {
        message <- conditionMessage(e)
        if (startsWith(message, "kfl1ou pBIC capability check failed:"))
            stop(message, call.=FALSE)
        fail(message)
    })
}
"""


def collect_pbic_validation(directory):
    """Read attestation, rejecting changed library files or invalid probe output."""
    directory = Path(directory)

    def rows(name):
        with (directory / name).open(encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    identities = rows("pbic-identity.tsv")
    if len(identities) != 1 or identities[0]["contract"] != PBIC_CONTRACT:
        raise ValueError("Missing or incompatible pBIC capability attestation")
    identity = identities[0]
    checks = rows("pbic-checks.tsv")
    expected_checks = {
        f"{root}/{phase}/{representation}/{quantity}"
        for root in ("OUfixedRoot", "OUrandomRoot")
        for phase in ("fixed", "estimated")
        for representation in ("free", "singleton")
        for quantity in ("likelihood", "pBIC")
    } | {
        f"{root}/{suffix}"
        for root in ("OUfixedRoot", "OUrandomRoot")
        for suffix in (
            "fixed/coordinate_equivalence",
            "shared/likelihood",
            "shared/pBIC",
        )
    }
    if (
        len(checks) != len(expected_checks)
        or {row["check"] for row in checks} != expected_checks
    ):
        raise ValueError("Incomplete pBIC capability checks")
    for row in checks:
        actual, expected, error = (
            float(row[key]) for key in ("actual", "expected", "absolute_error")
        )
        if (
            not all(math.isfinite(value) for value in (actual, expected, error))
            or not math.isclose(error, abs(actual - expected), abs_tol=1e-12)
            or not 0 <= error <= PBIC_TOLERANCE
        ):
            raise ValueError("Invalid pBIC capability check result")
        row.update(actual=actual, expected=expected, absolute_error=error)
    library = Path(identity["library"]).resolve(strict=True)
    files = rows("pbic-files.tsv")
    if not files or len({row["path"] for row in files}) != len(files):
        raise ValueError("Incomplete pBIC package fingerprint")
    hashes = {}
    for row in files:
        path = (library / row["path"]).resolve(strict=True)
        if not path.is_relative_to(library):
            raise ValueError("pBIC package fingerprint escapes its library")
        data = path.read_bytes()
        if hashlib.md5(data, usedforsecurity=False).hexdigest() != row["md5"]:
            raise ValueError("kfl1ou package changed during inference")
        hashes[row["path"]] = hashlib.sha256(data).hexdigest()
    return {
        **identity,
        "status": "passed",
        "tolerance": PBIC_TOLERANCE,
        "checks": checks,
        "installed_files_sha256": hashes,
        "scope": "single-trait optimum-coordinate pBIC; not selection calibration",
    }
