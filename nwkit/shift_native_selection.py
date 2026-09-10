"""Replayable native search with bootstrap or information-criterion selection."""

from nwkit.shift_native_bootstrap import (
    calibrate_native_search,
    gate_native_aic,
    native_selection_support,
)
from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_path import sparse_native_search
from nwkit.shift_native_search import (
    NativeLayoutEvaluator,
    enumerate_native_layouts,
)


class NativeSearchRunner:
    def __init__(self, data, args, fit_arguments):
        if args.search_strategy not in {"auto", "exhaustive", "lasso", "native-path"}:
            raise ValueError(
                "Native strategy must be auto, exhaustive, lasso, or native-path."
            )
        if args.search_strategy == "native-path" and (
            args.criterion not in {"AIC", "AICc"} or args.convergence
        ):
            raise ValueError(
                "Native path search requires AIC or AICc without convergence."
            )
        self.use_path = args.search_strategy == "native-path"
        self.criterion = args.criterion
        self.fit_arguments = fit_arguments
        self.layouts = None
        self.metadata: dict = {}
        self.options = NativeSearchOptions(
            max_shifts=args.max_shifts,
            convergence=args.convergence,
            candidate_pool=args.candidate_pool,
            refit_budget=args.refit_budget,
            screening_budget=args.screening_budget,
            beam_width=args.beam_width,
            lasso_iterations=args.lasso_iterations,
            memory_limit=args.search_memory_mb * 1024**2,
        )
        if args.search_strategy not in {"lasso", "native-path"}:
            try:
                self.layouts, self.metadata = enumerate_native_layouts(
                    data,
                    args.max_shifts,
                    convergence=args.convergence,
                    limit=args.exhaustive_max_configurations,
                )
            except ValueError as exc:
                if args.search_strategy != "auto" or not any(
                    message in str(exc)
                    for message in ("traversal budget", "candidate limit")
                ):
                    raise
        if self.layouts is None:
            self.options.validate(data, uses_candidate_pool=not self.use_path)

    def __call__(self, data):
        if self.layouts is None:
            search = sparse_native_search if self.use_path else heuristic_native_search
            return search(
                data,
                options=self.options,
                fit_arguments=self.fit_arguments,
                criterion=self.criterion,
            )
        evaluator = NativeLayoutEvaluator(data, self.fit_arguments, self.criterion)
        for layout in self.layouts:
            evaluator.evaluate(layout)
        return evaluator.finish(
            {
                "strategy": "exhaustive",
                **self.metadata,
                "continuous_global_optimum_certified": False,
            }
        )


def select_native(data, args, arguments):
    if args.bootstrap < 0 or args.seed < 0 or args.bootstrap_seed < 0:
        raise ValueError("Native bootstrap count and seeds must be nonnegative.")
    gate = getattr(args, "global_null_gate", False)
    if gate and (args.criterion != "AIC" or getattr(args, "regime_map", None)):
        raise ValueError(
            "--global-null-gate requires native AIC search without --regime-map."
        )
    if (args.criterion is None or gate) and (
        args.calibration_replicates < 1
        or not 0 < args.calibration_level < 1
        or 1 / (args.calibration_replicates + 1) > args.calibration_level
    ):
        raise ValueError("Calibration draws cannot resolve the requested test level.")
    runner = NativeSearchRunner(data, args, arguments)
    search = runner(data)

    def choose(searched, sample, seed):
        if gate:
            return gate_native_aic(
                sample,
                searched,
                runner,
                replicates=args.calibration_replicates,
                seed=seed,
                level=args.calibration_level,
            )
        if args.criterion is not None:
            if searched.best_information is None:
                raise ValueError(
                    "No candidate has finite native information criterion."
                )
            return searched.best_information, None
        return calibrate_native_search(
            sample,
            searched,
            runner,
            replicates=args.calibration_replicates,
            seed=seed,
            level=args.calibration_level,
        )

    result, calibration = choose(search, data, args.seed)
    support = None
    if args.bootstrap:

        def select(sample, seed):
            searched = runner(sample)
            return choose(searched, sample, seed)[0]

        support = native_selection_support(
            data, result, select, replicates=args.bootstrap, seed=args.bootstrap_seed
        )
    return result, {
        "inference_role": "research_joint_shift_selection",
        "adoption_status": "not_validated_for_production",
        "selection_calibration": calibration,
        "information_criterion": result.get("information_criterion"),
        "selection_support": support,
        "search": search.metadata,
        "candidates": search.records,
    }
