"""Budgeted joint location/regime search with explicit incomplete coverage."""

import itertools
import math
from dataclasses import dataclass

from nwkit.shift_joint_screen import JointQuickProfile, joint_group_lasso_screen
from nwkit.shift_native_limits import NATIVE_SEARCH_DEFAULTS
from nwkit.shift_native_model import ShiftLayout
from nwkit.shift_native_quick import NativeQuickProfile
from nwkit.shift_native_screen import group_lasso_screen
from nwkit.shift_native_search import (
    NativeLayoutEvaluator,
    layout_complexity,
    observable_layout,
)


@dataclass(frozen=True)
class NativeSearchOptions:
    max_shifts: int = 2
    convergence: bool = False
    candidate_pool: int = NATIVE_SEARCH_DEFAULTS["candidate_pool"]
    refit_budget: int = NATIVE_SEARCH_DEFAULTS["refit_budget"]
    screening_budget: int = NATIVE_SEARCH_DEFAULTS["screening_budget"]
    beam_width: int = NATIVE_SEARCH_DEFAULTS["beam_width"]
    lasso_iterations: int = NATIVE_SEARCH_DEFAULTS["lasso_iterations"]
    memory_limit: int = NATIVE_SEARCH_DEFAULTS["search_memory_mb"] * 1024**2

    def validate(self, data, *, uses_candidate_pool=True):
        if not 0 <= self.max_shifts < len(data.tree.leaf_names) - 1:
            raise ValueError("Maximum shifts must be smaller than tips minus one.")
        if self.refit_budget <= self.max_shifts or (
            uses_candidate_pool and self.candidate_pool < self.max_shifts
        ):
            raise ValueError(
                "Search needs at least max_shifts+1 refits and a candidate pool covering max_shifts."
            )
        if (
            min(
                self.screening_budget,
                self.beam_width,
                self.lasso_iterations,
                self.memory_limit,
            )
            < 1
        ):
            raise ValueError("Native search budgets must be positive.")


def _build(data, shifts, groups):
    try:
        layout = ShiftLayout.build(data.tree, shifts, groups)
    except ValueError:
        return None
    return layout if observable_layout(data, layout) else None


def added_layouts(data, layout, pool, convergence):
    for branch in pool:
        if branch in layout.shifts:
            continue
        groups = [list(g) for g in layout.groups]
        proposed = [groups + [[branch]]]
        if convergence:
            for i in range(len(groups)):
                modified = [g.copy() for g in groups]
                modified[i].append(branch)
                proposed.append(modified)
        for group in proposed:
            candidate = _build(data, (*layout.shifts, branch), group)
            if candidate is not None:
                yield candidate


def neighboring_layouts(data, layout, pool, convergence):
    for branch in layout.shifts:
        shifts = [b for b in layout.shifts if b != branch]
        groups = [[b for b in g if b != branch] for g in layout.groups]
        groups = [g for g in groups if g]
        candidate = _build(data, shifts, groups)
        if candidate is not None:
            yield candidate
        for replacement in pool:
            if replacement in layout.shifts:
                continue
            moved = [
                [replacement if b == branch else b for b in g] for g in layout.groups
            ]
            candidate = _build(data, (*shifts, replacement), moved)
            if candidate is not None:
                yield candidate
    if not convergence:
        return
    for first, second in itertools.combinations(range(len(layout.groups)), 2):
        groups = [
            list(g) for i, g in enumerate(layout.groups) if i not in (first, second)
        ]
        groups.append([*layout.groups[first], *layout.groups[second]])
        candidate = _build(data, layout.shifts, groups)
        if candidate is not None:
            yield candidate
    for group in layout.groups:
        if len(group) < 2:
            continue
        for branch in group:
            groups = [[b for b in g if b != branch] for g in layout.groups]
            groups = [g for g in groups if g] + [[branch]]
            candidate = _build(data, layout.shifts, groups)
            if candidate is not None:
                yield candidate


class _HeuristicSearch:
    def __init__(self, data, options, fit_arguments, criterion=None):
        self.data, self.options = data, options
        self.evaluator = NativeLayoutEvaluator(data, fit_arguments, criterion)
        self.null = ShiftLayout.build(data.tree)
        self.evaluator.evaluate(self.null)
        null_fit = self.evaluator.best[0]
        screen = group_lasso_screen
        profile: type[NativeQuickProfile] | type[JointQuickProfile] = NativeQuickProfile
        if "joint_fit" in null_fit:
            screen, profile = joint_group_lasso_screen, JointQuickProfile
        self.pool, self.screen_metadata = screen(
            data,
            null_fit,
            pool_size=options.candidate_pool,
            iterations=options.lasso_iterations,
            memory_limit=options.memory_limit,
        )
        self.profile = profile(
            data,
            null_fit,
            self.pool,
            fit_arguments.get("alpha_height"),
            **(
                {"memory_limit": options.memory_limit}
                if "joint_fit" in null_fit
                else {}
            ),
        )
        self.quick_scores = {}

    def rank(self, candidates):
        accepted = set()
        for layout in candidates:
            if layout in self.evaluator.scores:
                continue
            if layout not in self.quick_scores:
                if self.profile.evaluations >= self.options.screening_budget:
                    break
                self.quick_scores[layout] = self.profile.score(layout)
            if math.isfinite(self.quick_scores[layout]):
                accepted.add(layout)
        return sorted(
            accepted,
            key=lambda layout: (
                -self.quick_scores[layout],
                layout.shifts,
                layout.groups,
            ),
        )

    def refit(self, ranked, limit):
        if not ranked:
            return []
        by_complexity: dict[int, ShiftLayout] = {}
        for layout in ranked:
            by_complexity.setdefault(layout_complexity(layout), layout)
        # Keep both the least-constrained and simplest promising layouts. This
        # supplies several nested complexity families without enumerating all.
        representatives = []
        costs = sorted(by_complexity)
        while costs:
            representatives.append(by_complexity[costs.pop()])
            if costs:
                representatives.append(by_complexity[costs.pop(0)])
        candidates = list(dict.fromkeys([*representatives, *ranked]))
        budget = min(limit, self.options.refit_budget - len(self.evaluator.records))
        chosen = candidates[:budget]
        for layout in chosen:
            self.evaluator.evaluate(layout)
        return sorted(
            chosen,
            key=lambda layout: (
                -self.evaluator.scores[layout],
                layout.shifts,
                layout.groups,
            ),
        )

    def forward(self):
        frontier = [self.null]
        for size in range(1, self.options.max_shifts + 1):
            proposed = itertools.chain.from_iterable(
                added_layouts(self.data, layout, self.pool, self.options.convergence)
                for layout in frontier
            )
            ranked = self.rank(proposed)
            remaining = self.options.max_shifts - size + 1
            budget = max(
                1,
                (self.options.refit_budget - len(self.evaluator.records))
                // (remaining + 1),
            )
            chosen = self.refit(ranked, min(2 * self.options.beam_width, budget))
            if not chosen:
                break
            unconstrained = next(
                (
                    candidate
                    for candidate in chosen
                    if len(candidate.groups) == len(candidate.shifts) + 1
                ),
                None,
            )
            convergent = next(
                (
                    candidate
                    for candidate in chosen
                    if len(candidate.groups) < len(candidate.shifts) + 1
                ),
                None,
            )
            frontier = list(
                dict.fromkeys(
                    [
                        candidate
                        for candidate in (unconstrained, convergent, *chosen)
                        if candidate is not None
                    ]
                )
            )[: self.options.beam_width]

    def refine(self):
        seeds = [
            result["layout"]
            for _, result in sorted(self.evaluator.best.items(), reverse=True)
        ]
        for layout in seeds:
            if len(self.evaluator.records) >= self.options.refit_budget:
                break
            ranked = self.rank(
                neighboring_layouts(
                    self.data, layout, self.pool, self.options.convergence
                )
            )
            self.refit(ranked, self.options.beam_width)


def heuristic_native_search(data, *, options=None, fit_arguments=None, criterion=None):
    options = NativeSearchOptions() if options is None else options
    options.validate(data)
    search = _HeuristicSearch(data, options, fit_arguments or {}, criterion)
    search.forward()
    search.refine()
    metadata = {
        "strategy": "group_lasso_beam",
        "complete_discrete_enumeration": False,
        "continuous_global_optimum_certified": False,
        "screening": search.screen_metadata,
        "quick_evaluations": search.profile.evaluations,
        "refitted_candidates": len(search.evaluator.records),
        "refit_budget": options.refit_budget,
        "screening_budget": options.screening_budget,
        "budget_exhausted": len(search.evaluator.records) >= options.refit_budget
        or search.profile.evaluations >= options.screening_budget,
        "refinement": "one_pass_drop_move_merge_split",
        "candidate_ranking": "fixed_joint_null_covariance_profile; final candidates refitted without penalties"
        if "joint_fit" in search.evaluator.best[0]
        else "fixed_null_covariance_profiles; final candidates refitted without penalties",
    }
    return search.evaluator.finish(metadata)
