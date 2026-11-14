"""Public R convergence orchestration and validated shared-optimum groups."""

from typing import Any

from nwkit.shift_bootstrap import _table
from nwkit.shift_math import close_in_units

R_CONVERGENCE = r"""
convergence.enabled <- length(args) >= 9L && args[9] == "TRUE"
converge <- function(model) {
    if (!is.finite(model$alpha) || model$alpha <= 0)
        stop("Convergence requires identifiable OU optima (alpha > 0).")
    result <- kfl1ou::estimate_convergent_regimes(
        model, criterion=args[2], method="backward", fixed.alpha=FALSE, nCores=1)
    if (!is.finite(result$alpha) || result$alpha <= 0)
        stop("Convergence refit has unidentifiable OU optima (alpha <= 0).")
    result
}
group_keys <- function(model) {
    if (!identical(model$tree$edge, tr$edge) ||
        !identical(model$tree$tip.label, tr$tip.label))
        stop("Convergence changed tree edge coordinates.")
    vapply(model$convergent.regimes, function(group)
        paste(vapply(group, function(edge)
            if (edge == 0) "@root" else keys[edge], ""), collapse=";"), "")
}
write_model <- function(model, path) write_tsv(data.frame(
    backend_version=as.character(utils::packageVersion("kfl1ou")),
    alpha=model$alpha, sigma2=model$sigma2, intercept=model$intercept,
    log_likelihood=sum(model$logLik), score=model$score), path)
search.fit <- fit
if (convergence.enabled) {
    write_model(search.fit, "unconstrained-model.tsv")
    fit <- converge(search.fit)
    attr(fit, "nwkit.unconstrained.fit") <- search.fit
    write_tsv(data.frame(clades=group_keys(fit)), "convergence.tsv")
}
convergent_bootstrap <- function() {
    simulations <- stats::simulate(fit, nsim=as.integer(args[7]),
        seed=as.integer(args[8]), preserve.missing=TRUE, engine="tree")
    records <- lapply(simulations, function(Ystar) tryCatch({
        initial <- kfl1ou::estimate_shift_configuration(
            tr, Ystar, max.nShifts=as.integer(args[1]), criterion=args[2],
            root.model=args[3], search.strategy=args[4],
            exhaustive.max.configurations=as.integer(args[5]),
            ensemble.seed=as.integer(args[6]), nCores=1, quietly=TRUE,
            input_error=input.error, measurement_error=FALSE)
        refit <- converge(initial)
        list(ok=TRUE, shifts=unname(refit$shift.configuration),
             groups=paste(group_keys(refit), collapse="|"), error="")
    }, error=function(e) list(ok=FALSE, error=conditionMessage(e))))
    good <- vapply(records, function(record) record$ok, logical(1))
    if (!any(good)) stop("all bootstrap replicates failed.")
    errors <- vapply(records[!good], function(record) record$error, "")
    write_tsv(data.frame(success_index=seq_len(sum(good)),
        groups=vapply(records[good], function(record) record$groups, "")),
        "bootstrap-convergence.tsv")
    list(attempted=length(records), successful=sum(good), failed=sum(!good),
         all.shifts=lapply(records[good], function(record) record$shifts),
         failure.messages=table(errors[nzchar(errors)]))
}
"""


def decode_groups(values, mapping, selected):
    groups = []
    for value in values:
        keys = value.split(";")
        if any(key != "@root" and key not in mapping for key in keys):
            raise ValueError("Invalid convergence clade.")
        group = [0 if key == "@root" else mapping[key] for key in keys]
        if any(key != "@root" and mapping[key] == 0 for key in keys):
            raise ValueError(
                "Convergence background requires its explicit root marker."
            )
        groups.append(sorted(group))
    flat = [branch for group in groups for branch in group]
    if len(flat) != len(set(flat)) or set(flat) != {0, *selected}:
        raise ValueError(
            "Convergence groups must partition background and selected shifts."
        )
    return sorted(groups)


def collect_convergence(directory, mapping, selected, enabled):
    if not enabled:
        return None, {branch: f"shift_{branch}" for branch in selected}
    groups = decode_groups(
        [row["clades"] for row in _table(directory, "convergence.tsv", ["clades"])],
        mapping,
        selected,
    )
    aliases = {
        branch: "baseline" if 0 in group else f"shift_{min(group)}"
        for group in groups
        for branch in group
    }
    return {
        "method": "backward",
        "alpha_refitted": True,
        "globally_optimal": False,
        "groups": [
            {"regime": aliases[group[0]], "branch_ids": group} for group in groups
        ],
        "merges": 1 + len(selected) - len(groups),
    }, aliases


def shared_regime_rows(rows, scales):
    unique: dict[str, Any] = {}
    for row in rows:
        previous = unique.get(row["regime"])
        if previous is not None:
            if not close_in_units(
                previous["optimum"],
                row["optimum"],
                operands=(scales[previous["branch_id"]], scales[row["branch_id"]]),
            ):
                raise ValueError("Constrained regime optima disagree.")
            if row["branch_id"] >= previous["branch_id"]:
                continue
        unique[row["regime"]] = row
    return sorted(unique.values(), key=lambda row: row["branch_id"])


def convergence_bootstrap(directory, configurations, mapping, tree, ids):
    from collections import Counter

    rows = _table(directory, "bootstrap-convergence.tsv", ["success_index", "groups"])
    if len(rows) != len(configurations):
        raise ValueError("Convergence bootstrap counts disagree.")
    partitions: Counter[tuple[tuple[str, ...], ...]] = Counter()
    groups_by_replicate = []
    for index, (row, selected) in enumerate(zip(rows, configurations, strict=True), 1):
        if row["success_index"] != str(index):
            raise ValueError("Invalid convergence bootstrap success index.")
        groups = decode_groups(row["groups"].split("|"), mapping, selected)
        groups_by_replicate.append(groups)
        aliases = {branch: min(group) for group in groups for branch in group}
        regimes: dict[Any, int] = {}
        tips: dict[int, list[str]] = {}
        for node in tree.traverse("preorder"):
            regimes[node] = (
                aliases[ids[node]] if ids[node] in aliases else regimes[node.up]
            )
            if node.is_leaf:
                tips.setdefault(regimes[node], []).append(node.name)
        partition = tuple(sorted(tuple(sorted(names)) for names in tips.values()))
        partitions[partition] += 1
    return {
        "selection": "shift_and_convergence",
        "successful_convergence_groups": groups_by_replicate,
        "shared_optimum_partition_frequencies": [
            {
                "groups": [list(group) for group in partition],
                "count": count,
                "frequency": count / len(rows),
            }
            for partition, count in sorted(
                partitions.items(), key=lambda item: (-item[1], item[0])
            )
        ],
    }
