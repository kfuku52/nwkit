"""Prespecified paired OU selection and alpha-bound sensitivity design."""

import math

import numpy as np
from shift_simulation_cases import SCENARIOS, simulate

ROOTS = ("OUfixedRoot", "OUrandomRoot")
FLOORS = (("small", 1e-7), ("raised", 0.1))


def protocol(replicates=50, extension_replicates=10, seed=20260911):
    if replicates < 1 or extension_replicates < 0 or seed < 0:
        raise ValueError("Invalid replication counts or seed")
    cells = []
    specifications = [("primary", 8, 2.1, 0.0, SCENARIOS, replicates)]
    specifications += [
        ("weak_pull", 8, 0.2, 0.0, ("null", "convergent"), extension_replicates),
        ("known_error", 8, 2.1, 0.2, ("null", "convergent"), extension_replicates),
        ("sixteen_tips", 16, 2.1, 0.0, ("null", "convergent"), extension_replicates),
    ]
    for family, tips, alpha_height, se, scenarios, count in specifications:
        if not count:
            continue
        for scenario in scenarios:
            for root in ROOTS:
                cells.append(
                    {
                        "cell_id": len(cells),
                        "family": family,
                        "tips": tips,
                        "scenario": scenario,
                        "root_model": root,
                        "alpha_height": alpha_height,
                        "sigma2_height": 0.25,
                        "standard_error": se,
                        "effect": 2.0,
                        "replicates": count,
                    }
                )
    return {
        "schema_version": 1,
        "master_seed": seed,
        "cells": cells,
        "floors": [{"floor_id": key, "alpha_height": value} for key, value in FLOORS],
        "alpha_upper_height": 10.0,
        "alpha_start_height": 1.0,
        "criteria": ["BIC", "pBIC"],
        "methods": ["two_stage", "joint"],
        "max_shifts": 2,
        "workers": 4,
        "boundary_relative_tolerance": 1e-4,
        "comparison_score_tolerance": 1e-4,
        "interpretation": "Independent paired simulation; rates are conditional on these finite balanced-tree cells, not nominal error guarantees.",
    }


def cases(specification):
    for cell in specification["cells"]:
        for replicate in range(cell["replicates"]):
            seed = int(
                np.random.SeedSequence(
                    [specification["master_seed"], cell["cell_id"], replicate]
                ).generate_state(1)[0]
            )
            yield {
                **cell,
                "replicate": replicate,
                "seed": seed,
                "case_id": f"c{cell['cell_id']:02d}-r{replicate:03d}",
            }


def generate_case(case):
    height = math.log2(case["tips"])
    return simulate(
        tips=case["tips"],
        scenario=case["scenario"],
        root_model=case["root_model"],
        se=case["standard_error"],
        effect=case["effect"],
        seed=case["seed"],
        alpha=case["alpha_height"] / height,
        sigma2=case["sigma2_height"] / height,
    )
