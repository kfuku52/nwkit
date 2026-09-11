"""Resolve automatic search caps without changing explicit integer requests."""

NATIVE_SEARCH_DEFAULTS = {
    "candidate_pool": 128,
    "refit_budget": 256,
    "screening_budget": 100000,
    "beam_width": 2,
    "lasso_iterations": 150,
    "search_memory_mb": 512,
}


def native_shift_limit(data, args, *, budgeted=False):
    requested = args.max_shifts
    tree_limit = len(data.tree.leaf_names) - 2
    if requested != "auto":
        if (
            isinstance(requested, bool)
            or not isinstance(requested, int)
            or not 0 <= requested <= tree_limit
        ):
            raise ValueError(
                "Maximum shifts must be a nonnegative integer smaller than tips minus one, or auto."
            )
        return requested, {
            "requested": requested,
            "resolved": requested,
            "constraints": {"explicit": requested, "tree": tree_limit},
            "budget_limited": False,
        }
    constraints = {"tree": tree_limit}
    if budgeted:
        if args.refit_budget < 1:
            raise ValueError(
                "Automatic maximum shifts requires a positive refit budget."
            )
        constraints["refit_budget"] = args.refit_budget - 1
        if args.search_strategy != "native-path":
            if args.candidate_pool < 1:
                raise ValueError(
                    "Automatic maximum shifts requires a positive candidate pool."
                )
            constraints["candidate_pool"] = args.candidate_pool
    resolved = min(constraints.values())
    return resolved, {
        "requested": "auto",
        "resolved": resolved,
        "constraints": constraints,
        "budget_limited": resolved < tree_limit,
    }
