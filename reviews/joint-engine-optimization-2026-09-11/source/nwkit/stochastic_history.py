"""Retain conditional CTMC bridge histories without changing count-only seeds."""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

from nwkit.util import assign_branch_ids

_MAX_HISTORY_SEGMENTS = 1_000_000


@dataclass(frozen=True)
class MapDraw:
    simulation: int
    node_states: np.ndarray
    # Each branch holds (local start, local end, projected state) segments.
    branches: dict[int, tuple]
    counts: dict


@dataclass(frozen=True)
class MapSample:
    states: tuple
    draws: tuple
    branch_ids: dict
    depths: dict
    lengths: dict


def _segments(path, length, projection, rng):
    if len(path) == 1:
        return ((0.0, length, int(projection[path[0]])),)
    times = np.sort(rng.random(len(path) - 1)) * length
    if np.any(np.diff(times) <= 0) or times[0] <= 0 or times[-1] >= length:
        raise ValueError(
            "Stochastic-map event times cannot be represented distinctly inside the branch; reduce rate/time scale."
        )
    result = []
    start, state = 0.0, int(projection[path[0]])
    for time, next_state in zip(times, path[1:], strict=True):
        projected = int(projection[next_state])
        if projected != state:
            result.append((start, float(time), state))
            start, state = float(time), projected
    result.append((start, length, state))
    return tuple(result)


def _draw_history(spec, simulation, seed_sequence):
    from nwkit.asr import (
        _sample_bridge_transition_counts,
        _sample_node_states_from_spec,
    )

    rng = np.random.default_rng(seed_sequence)
    # Times are independent conditional on N. Never consume the legacy stream.
    time_seed = np.random.SeedSequence(
        seed_sequence.entropy, spawn_key=(*seed_sequence.spawn_key, 918273)
    )
    times_rng = np.random.default_rng(time_seed)
    node_states = _sample_node_states_from_spec(spec, rng)
    projection = spec.get("state_projection")
    if projection is None:
        projection = np.arange(spec["num_states"])
    branches = {}
    counts: dict[tuple[int, int, int], int] = {}
    for index, parent in enumerate(spec["parent_indices"]):
        if parent < 0:
            continue
        path = [int(node_states[parent])]
        _sample_bridge_transition_counts(
            int(node_states[parent]),
            int(node_states[index]),
            spec["branch_lengths"][index],
            spec["rate_matrices"][index],
            rng,
            uniformization_context=spec["uniformization_contexts"][index],
            history=path,
        )
        if path[-1] != node_states[index]:
            raise ValueError("Conditional stochastic-map bridge missed its endpoint.")
        branch = spec["branch_ids"][index]
        segments = _segments(path, spec["branch_lengths"][index], projection, times_rng)
        branches[branch] = segments
        for first, second in zip(segments, segments[1:], strict=False):
            key = (branch, first[2], second[2])
            counts[key] = counts.get(key, 0) + 1
    return MapDraw(simulation, np.asarray(projection)[node_states], branches, counts)


def _draw_chunk(payload):
    spec, indexed_seeds = payload
    return [_draw_history(spec, index + 1, seed) for index, seed in indexed_seeds]


def _history_work(fit, branches, count):
    from nwkit.asr import _uniformization_parameters, _validate_stochastic_map_work

    if count * (len(branches) + 1) > 100_000:
        raise ValueError(
            "Map histories exceed 100,000 draw/node records; reduce --n-sim."
        )
    size = fit["rate_matrix"].shape[0]
    # The legacy specification builds Q, P and log(P), plus inside likelihoods.
    dense_bytes = (len(branches) + 1) * (3 * size * size + size) * 8
    if dense_bytes > 256 * 1024**2:
        raise ValueError(
            "Map simulation arrays exceed 256 MiB before histories; reduce tree or expanded-state size."
        )
    _validate_stochastic_map_work(fit, branches, count)
    matrices = fit.get("rate_matrix_by_node", {})
    terms = {}
    per_draw = 0
    for node in branches:
        matrix = matrices.get(node, fit["rate_matrix"])
        key = (id(matrix), float(node.dist))
        if key not in terms:
            terms[key] = _uniformization_parameters(matrix, node.dist)[2]
        per_draw += 1 + terms[key]
    bound = per_draw * count
    if bound > _MAX_HISTORY_SEGMENTS:
        raise ValueError(
            "Stochastic-map histories exceed the 1,000,000-segment storage budget; reduce --n-sim or the rate/time scale."
        )


def sample_histories(tree, states, fit, count, *, seed=None, threads=1):
    from nwkit.asr import (
        _build_stochastic_map_spec,
        _get_process_pool_context,
        _simulation_seed_sequence,
        _stochastic_uniformization_contexts,
        _validated_simulation_threads,
    )

    count, threads = _validated_simulation_threads(count, threads)
    identifiers = assign_branch_ids(tree)
    nodes = list(tree.traverse("preorder"))
    branches = [node for node in nodes if not node.is_root]
    depths = {tree: 0.0}
    for node in nodes[1:]:
        depths[node] = depths[node.up] + float(node.dist)
        if float(node.dist) > 0 and not np.isclose(
            depths[node] - depths[node.up], float(node.dist), rtol=1e-8, atol=0
        ):
            raise ValueError(
                "Positive map branch length loses precision on the root-time axis; reduce the root-depth/branch-length scale disparity."
            )
        if not np.isfinite(depths[node]):
            raise ValueError(
                "Stochastic-map time depths overflow; rescale branch lengths."
            )
    _history_work(fit, branches, count)
    contexts = _stochastic_uniformization_contexts(fit, branches)
    spec = _build_stochastic_map_spec(tree, fit, identifiers, contexts)
    indexed = list(enumerate(_simulation_seed_sequence(seed, count)))
    if threads == 1 or count == 1:
        draws = _draw_chunk((spec, indexed))
    else:
        workers = min(threads, count)
        chunks = [indexed[index::workers] for index in range(workers)]
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=_get_process_pool_context()
        ) as executor:
            draws = [
                draw
                for chunk in executor.map(
                    _draw_chunk, ((spec, chunk) for chunk in chunks)
                )
                for draw in chunk
            ]
    return MapSample(
        tuple(states),
        tuple(sorted(draws, key=lambda draw: draw.simulation)),
        identifiers,
        depths,
        {identifiers[node]: float(node.dist) for node in branches},
    )
