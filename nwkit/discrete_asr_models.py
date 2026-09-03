"""Transition structures and fixed generators for discrete ASR."""

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

_ROUNDOFF_FACTOR = 128.0
FREQUENCY_RATIO_BOUNDS = (1e-8, 1e8)


@dataclass(frozen=True, slots=True)
class RateDesign:
    """A fitted direct-rate partition over directed CTMC edges."""

    states: tuple[str, ...]
    class_names: tuple[str, ...]
    class_by_edge: np.ndarray
    source: str

    @property
    def graph(self):
        return self.class_by_edge >= 0


def _roundoff_close(residual, scale):
    residual = np.asarray(residual, dtype=float)
    scale = np.asarray(scale, dtype=float)
    tolerance = _ROUNDOFF_FACTOR * np.finfo(float).eps * scale
    return np.abs(residual) <= tolerance


def _validated_graph(graph, num_states):
    graph = np.asarray(graph, dtype=bool)
    if graph.shape != (num_states, num_states):
        raise ValueError("Transition graph dimensions do not match the state space.")
    if np.any(np.diag(graph)):
        raise ValueError("A transition graph cannot contain self edges.")
    return graph


def complete_transition_graph(num_states):
    graph = np.ones((num_states, num_states), dtype=bool)
    np.fill_diagonal(graph, False)
    return graph


def ordered_transition_graph(num_states):
    graph = np.zeros((num_states, num_states), dtype=bool)
    for index in range(max(0, num_states - 1)):
        graph[index, index + 1] = True
        graph[index + 1, index] = True
    return graph


def _canonical_rate_partition(class_by_edge):
    groups = []
    for class_index in sorted(set(int(value) for value in class_by_edge.ravel())):
        if class_index < 0:
            continue
        edges = tuple(
            sorted(
                (int(source), int(target))
                for source, target in zip(
                    *np.nonzero(class_by_edge == class_index), strict=True
                )
            )
        )
        groups.append(edges)
    return tuple(sorted(groups))


def _direct_rate_partition(model, states, graph, rate_design):
    size = len(states)
    classes = np.full((size, size), -1, dtype=int)
    if model == "MK-DESIGN":
        if rate_design is None:
            return None
        return _canonical_rate_partition(rate_design.class_by_edge)
    if model == "ER":
        classes[graph] = 0
    elif model == "SYM":
        if not np.array_equal(graph, graph.T):
            return None
        class_index = 0
        for source in range(size):
            for target in range(source + 1, size):
                if graph[source, target]:
                    classes[source, target] = classes[target, source] = class_index
                    class_index += 1
    elif model == "ARD":
        for class_index, (source, target) in enumerate(
            zip(*np.nonzero(graph), strict=True)
        ):
            classes[source, target] = class_index
    elif model == "F81":
        if not np.array_equal(graph, complete_transition_graph(size)):
            return None
        for target in range(size):
            classes[:, target] = target
        np.fill_diagonal(classes, -1)
    else:
        return None
    return _canonical_rate_partition(classes)


def model_equivalence_family(model, states, graph=None, rate_design=None):
    """Return the exact structured-rate family used for duplicate-fit detection."""

    if model == "MK-DESIGN":
        if rate_design is None:
            return model
        graph = rate_design.graph
    else:
        graph = (
            complete_transition_graph(len(states))
            if graph is None
            else _validated_graph(graph, len(states))
        )
    partition = _direct_rate_partition(model, states, graph, rate_design)
    if partition is not None:
        return ("direct-rate-partition", partition)
    return model


def rate_design_from_edges(states, edges, *, source):
    """Build and validate a rate design from ``(from, to, class)`` triples."""

    states = tuple(str(state) for state in states)
    if not states or len(states) != len(set(states)):
        raise ValueError("A rate design requires a non-empty unique state space.")
    state_to_index = {state: index for index, state in enumerate(states)}
    class_names: list[str] = []
    class_to_index = {}
    class_by_edge = np.full((len(states), len(states)), -1, dtype=int)
    seen_edges = set()
    for raw_source, raw_target, raw_class in edges:
        from_state = str(raw_source)
        to_state = str(raw_target)
        class_name = str(raw_class).strip()
        if from_state not in state_to_index or to_state not in state_to_index:
            unknown = from_state if from_state not in state_to_index else to_state
            raise ValueError(
                "State in the rate design is absent from the model state space: "
                f"{unknown}"
            )
        if from_state == to_state:
            raise ValueError("A rate design cannot contain self edges.")
        if class_name == "":
            raise ValueError("A rate design cannot contain an empty rate_class.")
        edge = (from_state, to_state)
        if edge in seen_edges:
            raise ValueError(
                f"A rate design contains a duplicated edge: {from_state} -> {to_state}"
            )
        seen_edges.add(edge)
        if class_name not in class_to_index:
            class_to_index[class_name] = len(class_names)
            class_names.append(class_name)
        class_by_edge[state_to_index[from_state], state_to_index[to_state]] = (
            class_to_index[class_name]
        )
    if not seen_edges:
        raise ValueError("A rate design must contain at least one directed edge.")
    return RateDesign(
        states,
        tuple(class_names),
        class_by_edge,
        str(source),
    )


def read_rate_design(path, states):
    """Read a fitted edge/rate-class design TSV."""

    source = Path(str(path))
    try:
        table = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError) as exc:
        raise ValueError(f"Failed to read '--rate-design': {source}") from exc
    required = ("from_state", "to_state", "rate_class")
    if tuple(table.columns) != required:
        raise ValueError(
            "'--rate-design' must contain exactly these TSV columns in order: "
            "from_state, to_state, rate_class."
        )
    return rate_design_from_edges(
        states,
        table.itertuples(index=False, name=None),
        source=str(source),
    )


def read_transition_graph(specification, states, *, state_source="--states"):
    """Return an allowed directed-edge mask and a reproducible source label."""

    if specification in (None, "", "complete"):
        return complete_transition_graph(len(states)), "complete"
    if specification == "ordered":
        return ordered_transition_graph(len(states)), "ordered"

    path = Path(str(specification))
    try:
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError) as exc:
        raise ValueError(f"Failed to read '--transition-graph': {path}") from exc
    required = {"from_state", "to_state"}
    if not required.issubset(table.columns):
        raise ValueError(
            "'--transition-graph' must be a TSV edge list with from_state and "
            "to_state columns."
        )
    state_to_index = {state: index for index, state in enumerate(states)}
    graph = np.zeros((len(states), len(states)), dtype=bool)
    seen = set()
    for row in table.itertuples(index=False):
        source = str(row.from_state)
        target = str(row.to_state)
        if source not in state_to_index or target not in state_to_index:
            unknown = source if source not in state_to_index else target
            raise ValueError(
                "State in '--transition-graph' is absent from the state space "
                f"defined by '{state_source}': {unknown}"
            )
        if source == target:
            raise ValueError("'--transition-graph' cannot contain self edges.")
        edge = (source, target)
        if edge in seen:
            raise ValueError(
                f"'--transition-graph' contains a duplicated edge: {source} -> {target}"
            )
        seen.add(edge)
        graph[state_to_index[source], state_to_index[target]] = True
    if not np.any(graph):
        raise ValueError("'--transition-graph' must contain at least one edge.")
    return graph, str(path)


def parameter_labels(model, states, graph=None, rate_design=None):
    if model == "MK-DESIGN":
        if rate_design is None:
            raise ValueError("--model MK-DESIGN requires --rate-design.")
        if tuple(states) != rate_design.states:
            raise ValueError("Rate-design states do not match the model state space.")
        if graph is not None:
            raise ValueError("MK-DESIGN cannot use a separate transition graph.")
        return [("rate_class", name) for name in rate_design.class_names]
    graph = (
        complete_transition_graph(len(states))
        if graph is None
        else _validated_graph(graph, len(states))
    )
    if model == "ER":
        return [("all", "all")] if np.any(graph) else []
    if model == "SYM":
        if not np.array_equal(graph, graph.T):
            raise ValueError("SYM requires a symmetric '--transition-graph'.")
        return [
            (states[i], states[j])
            for i in range(len(states))
            for j in range(i + 1, len(states))
            if graph[i, j]
        ]
    if model == "ARD":
        return [
            (states[i], states[j])
            for i in range(len(states))
            for j in range(len(states))
            if graph[i, j]
        ]
    if model == "F81":
        _require_complete_graph(model, graph)
        return [("target", state) for state in states] if len(states) > 1 else []
    if model == "GTR":
        _require_complete_graph(model, graph)
        exchangeabilities = [
            ("exchangeability", f"{states[i]}<->{states[j]}")
            for i in range(len(states))
            for j in range(i + 1, len(states))
        ]
        frequency_ratios = [("frequency_ratio", state) for state in states[1:]]
        return exchangeabilities + frequency_ratios
    if model == "CUSTOM":
        return []
    raise ValueError(f"Unsupported '--model': {model}")


def _require_complete_graph(model, graph):
    expected = complete_transition_graph(graph.shape[0])
    if not np.array_equal(graph, expected):
        raise ValueError(f"{model} requires a complete '--transition-graph'.")


def parameter_kinds(model, states, graph=None, rate_design=None):
    """Return optimizer-bound classes parallel to :func:`parameter_labels`."""

    labels = parameter_labels(model, states, graph, rate_design)
    if model == "GTR":
        return [
            "frequency_ratio" if source == "frequency_ratio" else "rate"
            for source, _ in labels
        ]
    return ["rate"] * len(labels)


def initial_parameters(model, states, initial_rate, graph=None, rate_design=None):
    """Return a homogeneous, neutral-frequency optimizer starting point."""

    kinds = parameter_kinds(model, states, graph, rate_design)
    return np.asarray(
        [1.0 if kind == "frequency_ratio" else initial_rate for kind in kinds],
        dtype=float,
    )


def build_rate_matrix(model, states, rates, graph=None, rate_design=None):
    if model == "MK-DESIGN":
        labels = parameter_labels(model, states, graph, rate_design)
        rates = np.asarray(rates, dtype=float)
        if len(rates) != len(labels):
            raise ValueError(
                f"Unexpected number of rate parameters for model '{model}'."
            )
        if np.any(~np.isfinite(rates)) or np.any(rates < 0.0):
            raise ValueError("Mk rates must be non-negative finite numbers.")
        assert rate_design is not None
        matrix = np.zeros((len(states), len(states)), dtype=float)
        for class_index, rate in enumerate(rates):
            matrix[rate_design.class_by_edge == class_index] = rate
        np.fill_diagonal(matrix, -matrix.sum(axis=1))
        return matrix
    graph = (
        complete_transition_graph(len(states))
        if graph is None
        else _validated_graph(graph, len(states))
    )
    labels = parameter_labels(model, states, graph)
    rates = np.asarray(rates, dtype=float)
    if len(rates) != len(labels):
        raise ValueError(f"Unexpected number of rate parameters for model '{model}'.")
    if np.any(~np.isfinite(rates)) or np.any(rates < 0.0):
        raise ValueError("Mk rates must be non-negative finite numbers.")
    matrix = np.zeros((len(states), len(states)), dtype=float)
    if model == "ER":
        if labels:
            matrix[graph] = rates[0]
    elif model == "SYM":
        state_to_index = {state: index for index, state in enumerate(states)}
        for rate, (source, target) in zip(rates, labels, strict=True):
            i, j = state_to_index[source], state_to_index[target]
            matrix[i, j] = matrix[j, i] = rate
    elif model == "ARD":
        state_to_index = {state: index for index, state in enumerate(states)}
        for rate, (source, target) in zip(rates, labels, strict=True):
            matrix[state_to_index[source], state_to_index[target]] = rate
    elif model == "F81":
        for target_index, target_rate in enumerate(rates):
            matrix[:, target_index] = target_rate
        np.fill_diagonal(matrix, 0.0)
    elif model == "GTR":
        num_exchangeabilities = len(states) * (len(states) - 1) // 2
        exchangeabilities = rates[:num_exchangeabilities]
        frequency_weights = np.concatenate(
            (np.ones(1, dtype=float), rates[num_exchangeabilities:])
        )
        if not np.any(frequency_weights > 0.0):
            raise ValueError("GTR stationary-frequency weights cannot all be zero.")
        frequencies = frequency_weights / frequency_weights.sum()
        parameter_index = 0
        for i in range(len(states)):
            for j in range(i + 1, len(states)):
                exchangeability = exchangeabilities[parameter_index]
                matrix[i, j] = exchangeability * frequencies[j]
                matrix[j, i] = exchangeability * frequencies[i]
                parameter_index += 1
    elif model != "CUSTOM":
        raise ValueError(f"Unsupported '--model': {model}")
    np.fill_diagonal(matrix, -matrix.sum(axis=1))
    return matrix


def read_rate_matrix(path):
    """Read a labelled TSV generator, deriving zero diagonals when requested."""

    try:
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError) as exc:
        raise ValueError(f"Failed to read '--rate-matrix': {path}") from exc
    if len(table.columns) < 2 or table.columns[0] != "state":
        raise ValueError(
            "'--rate-matrix' must have 'state' as its first column followed by "
            "one column per state."
        )
    states = [str(column) for column in table.columns[1:]]
    row_states = [str(value) for value in table["state"]]
    if len(states) == 0 or len(states) != len(set(states)):
        raise ValueError("'--rate-matrix' state columns must be non-empty and unique.")
    if row_states != states:
        raise ValueError(
            "'--rate-matrix' row states must exactly match its state columns in order."
        )
    try:
        matrix = table[states].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("'--rate-matrix' entries must be numeric.") from exc
    if matrix.shape != (len(states), len(states)) or np.any(~np.isfinite(matrix)):
        raise ValueError("'--rate-matrix' must contain a finite square matrix.")
    off_diagonal = matrix.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    if np.any(off_diagonal < 0.0):
        raise ValueError("'--rate-matrix' off-diagonal rates must be non-negative.")
    expected_diagonal = -off_diagonal.sum(axis=1)
    diagonal = np.diag(matrix)
    if np.all(diagonal == 0.0):
        matrix = off_diagonal
        np.fill_diagonal(matrix, expected_diagonal)
    else:
        diagonal_scale = np.abs(diagonal) + np.abs(expected_diagonal)
        valid_diagonal = _roundoff_close(diagonal - expected_diagonal, diagonal_scale)
        if np.all(valid_diagonal):
            return states, validate_rate_matrix(matrix)
        raise ValueError(
            "'--rate-matrix' diagonals must be zero (to derive them) or the "
            "negative row sums of off-diagonal rates."
        )
    return states, validate_rate_matrix(matrix)


def validate_rate_matrix(rate_matrix, *, num_states=None):
    """Validate and copy a finite continuous-time Markov generator."""

    try:
        matrix = np.asarray(rate_matrix, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("The rate matrix must be numeric.") from exc
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] == 0:
        raise ValueError("A non-empty square rate matrix is required.")
    if num_states is not None and matrix.shape != (num_states, num_states):
        raise ValueError("'--rate-matrix' dimensions do not match the state space.")
    if np.any(~np.isfinite(matrix)):
        raise ValueError("The rate matrix must be finite.")
    off_diagonal = matrix.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    if np.any(off_diagonal < 0.0) or np.any(np.diag(matrix) > 0.0):
        raise ValueError(
            "A rate matrix requires non-negative off-diagonals and non-positive diagonals."
        )
    row_scale = np.sum(np.abs(matrix), axis=1)
    if not np.all(_roundoff_close(matrix.sum(axis=1), row_scale)):
        raise ValueError("Rate-matrix rows must sum to zero.")
    return matrix.copy()


def stationary_distribution(rate_matrix):
    """Return the unique stationary distribution of a finite CTMC generator."""

    matrix = validate_rate_matrix(rate_matrix)
    system = matrix.T.copy()
    system[-1, :] = 1.0
    target = np.zeros(matrix.shape[0], dtype=float)
    target[-1] = 1.0
    try:
        distribution = np.linalg.solve(system, target)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "--root-prior stationary requires a unique stationary distribution."
        ) from exc
    probability_tolerance = _ROUNDOFF_FACTOR * np.finfo(float).eps
    residual_scale = float(np.linalg.norm(distribution, ord=1)) * float(
        np.linalg.norm(matrix, ord=np.inf)
    )
    residual_tolerance = _ROUNDOFF_FACTOR * np.finfo(float).eps * residual_scale
    residual = float(np.linalg.norm(distribution @ matrix, ord=np.inf))
    if (
        np.any(distribution < -probability_tolerance)
        or not math.isclose(
            float(distribution.sum()),
            1.0,
            rel_tol=0.0,
            abs_tol=probability_tolerance,
        )
        or residual > residual_tolerance
    ):
        raise ValueError(
            "--root-prior stationary requires a unique valid stationary distribution."
        )
    distribution = np.maximum(distribution, 0.0)
    return distribution / distribution.sum()
