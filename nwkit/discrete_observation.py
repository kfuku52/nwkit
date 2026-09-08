"""Explicit observation likelihoods and known discrete misclassification."""

import csv

import numpy as np

UNSUPPORTED_MODELS = frozenset(
    {"THRESHOLD", "MK-MIXTURE", "PAGEL-INDEPENDENT", "PAGEL-DEPENDENT"}
)


def validate_discrete_observation_options(args, model):
    paths = [
        getattr(args, name, None)
        for name in ("tip_likelihoods", "misclassification_matrix")
    ]
    if all(path not in (None, "") for path in paths):
        raise ValueError(
            "--tip-likelihoods and --misclassification-matrix are mutually exclusive."
        )
    if model in UNSUPPORTED_MODELS and any(path not in (None, "") for path in paths):
        raise ValueError(
            f"Explicit observation likelihoods are not supported for {model}."
        )


def _matrix_table(path, key, states):
    with open(path, newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, None)
        expected = [key, *states]
        if header != expected:
            raise ValueError(
                f"Observation table columns must exactly match {expected} in state order."
            )
        result = {}
        for row in reader:
            if len(row) != len(expected) or not row[0] or row[0] in result:
                raise ValueError(
                    "Observation tables require unique nonempty keys and complete rows."
                )
            try:
                vector = np.array(row[1:], dtype=float)
            except ValueError as exc:
                raise ValueError("Observation likelihoods must be numeric.") from exc
            if (
                not np.isfinite(vector).all()
                or np.any(vector < 0)
                or np.any(vector > 1)
                or not np.any(vector > 0)
            ):
                raise ValueError(
                    "Observation likelihoods must be finite, in [0, 1], and not all zero."
                )
            result[row[0]] = vector
    return result


def apply_discrete_observation_model(states, likelihoods, args):
    """Replace exact coding for supplied tips, or marginalize observed labels."""
    path = getattr(args, "tip_likelihoods", None)
    matrix_path = getattr(args, "misclassification_matrix", None)
    if path not in (None, "") and matrix_path not in (None, ""):
        raise ValueError(
            "--tip-likelihoods and --misclassification-matrix are mutually exclusive."
        )
    result = {
        name: np.asarray(value, dtype=float).copy()
        for name, value in likelihoods.items()
    }
    if path not in (None, ""):
        overrides = _matrix_table(path, "leaf_name", states)
        unknown = set(overrides) - set(result)
        if unknown:
            raise ValueError(
                "Unknown tips in --tip-likelihoods: " + ", ".join(sorted(unknown))
            )
        if not overrides:
            raise ValueError("--tip-likelihoods requires at least one tip row.")
        result.update(overrides)
    if matrix_path not in (None, ""):
        rows = _matrix_table(matrix_path, "state", states)
        if set(rows) != set(states):
            raise ValueError(
                "Misclassification rows must contain every state exactly once."
            )
        matrix = np.array([rows[state] for state in states])
        if not np.allclose(matrix.sum(axis=1), 1, atol=1e-12, rtol=1e-12):
            raise ValueError(
                "Misclassification rows (true state to observed state) must sum to one."
            )
        result = {name: matrix @ vector for name, vector in result.items()}
        if any(not np.any(vector > 0) for vector in result.values()):
            raise ValueError(
                "A coded observation is impossible under the misclassification matrix."
            )
    return result
