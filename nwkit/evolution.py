"""Evolutionary covariance models shared by PGLS and contrast workflows."""

import math
from dataclasses import dataclass
from io import StringIO
from typing import Any, Callable

import numpy as np
import pandas as pd

from nwkit.clade_index import LcaIndex
from nwkit.util import read_input_text


@dataclass(frozen=True)
class EvolutionModelSpec:
    name: str
    parameter_name: str | None
    contrast_supported: bool = True
    branch_lengths_used: bool = True


EVOLUTION_MODEL_SPECS = {
    spec.name: spec
    for spec in [
        EvolutionModelSpec("brownian", None),
        EvolutionModelSpec("lambda", "lambda"),
        EvolutionModelSpec("ou", "alpha"),
        EvolutionModelSpec("kappa", "kappa"),
        EvolutionModelSpec("delta", "delta"),
        EvolutionModelSpec("eb", "rate_change"),
        EvolutionModelSpec("acdc", "rate_change"),
        EvolutionModelSpec("independent", None, branch_lengths_used=False),
        EvolutionModelSpec(
            "custom", None, contrast_supported=False, branch_lengths_used=False
        ),
    ]
}
EVOLUTION_MODELS = tuple(EVOLUTION_MODEL_SPECS)
CONTRAST_EVOLUTION_MODELS = tuple(
    name for name, spec in EVOLUTION_MODEL_SPECS.items() if spec.contrast_supported
)
PARAMETERIZED_EVOLUTION_MODELS = tuple(
    name for name, spec in EVOLUTION_MODEL_SPECS.items() if spec.parameter_name
)


def evolution_model_spec(model: str) -> EvolutionModelSpec:
    try:
        return EVOLUTION_MODEL_SPECS[str(model)]
    except KeyError as exc:
        raise ValueError("Unsupported evolutionary model: {}.".format(model)) from exc


def _base_edge_lengths(tree, branch_length: str) -> dict[Any, float]:
    if branch_length not in {"original", "unit"}:
        raise ValueError("Unsupported branch-length mode: {}.".format(branch_length))
    lengths: dict[Any, float] = {}
    for node in tree.traverse():
        if node.is_root:
            lengths[node] = 0.0
            continue
        if branch_length == "unit":
            lengths[node] = 1.0
            continue
        try:
            value = float(node.dist)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("Tree branch lengths must be numeric.") from exc
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(
                "Evolutionary models require positive finite non-root branch lengths; "
                "use unit branch lengths to ignore input lengths."
            )
        lengths[node] = value
    return lengths


def _depths_from_edges(tree, edge_lengths: dict[Any, float]) -> dict[Any, float]:
    depths: dict[Any, float] = {}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            depths[node] = 0.0
        else:
            depths[node] = depths[node.up] + edge_lengths[node]
    return depths


def tree_depths(tree, branch_length: str = "original") -> dict[Any, float]:
    """Return root-to-node depths under original or unit branch lengths."""
    return _depths_from_edges(tree, _base_edge_lengths(tree, branch_length))


def _tree_height(tree, depths: dict[Any, float]) -> float:
    height = max(depths[leaf] for leaf in tree.leaves())
    if not math.isfinite(height) or height <= 0.0:
        raise ValueError("The evolutionary-model tree must have positive depth.")
    return float(height)


def _require_ultrametric(tree, depths: dict[Any, float], model: str) -> float:
    tip_depths = np.asarray([depths[leaf] for leaf in tree.leaves()], dtype=float)
    height = float(np.max(tip_depths))
    tolerance = max(
        1e-8 * max(1.0, height),
        np.finfo(float).eps * max(1.0, height) * max(100, len(tip_depths)),
    )
    if float(np.max(tip_depths) - np.min(tip_depths)) > tolerance:
        raise ValueError(
            "Evolutionary model '{}' requires an ultrametric tree.".format(model)
        )
    return height


def validate_evolution_parameter(model: str, parameter: float | None) -> float | None:
    spec = evolution_model_spec(model)
    if spec.parameter_name is None:
        if parameter is not None:
            raise ValueError(
                "Evolutionary model '{}' does not take a parameter.".format(model)
            )
        return None
    if parameter is None:
        raise ValueError(
            "Evolutionary model '{}' requires parameter '{}'.".format(
                model, spec.parameter_name
            )
        )
    if isinstance(parameter, (bool, np.bool_)):
        raise ValueError("Evolutionary-model parameters must be numeric, not boolean.")
    try:
        value = float(parameter)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Evolutionary-model parameters must be numeric.") from exc
    if not math.isfinite(value):
        raise ValueError("Evolutionary-model parameters must be finite.")
    if model == "lambda" and not 0.0 <= value <= 1.0:
        raise ValueError("Pagel's lambda must be between zero and one.")
    if model == "ou" and value <= 0.0:
        raise ValueError("OU alpha must be positive.")
    if model == "kappa" and value < 0.0:
        raise ValueError("Pagel's kappa must be non-negative.")
    if model == "delta" and value <= 0.0:
        raise ValueError("Pagel's delta must be positive.")
    if model == "eb" and value > 0.0:
        raise ValueError(
            "Early-burst rate change must be non-positive; use ACDC for acceleration."
        )
    return value


def _exponential_time_transform(depth: float, rate_change: float) -> float:
    if abs(rate_change * depth) < 1e-8:
        return depth * (1.0 + 0.5 * rate_change * depth)
    try:
        value = math.expm1(rate_change * depth) / rate_change
    except OverflowError as exc:
        raise ValueError(
            "The rate-change parameter overflows the tree time scale."
        ) from exc
    if not math.isfinite(value):
        raise ValueError(
            "The rate-change parameter produces non-finite branch lengths."
        )
    return value


def _ou_ultrametric_depth(depth: float, height: float, alpha: float) -> float:
    first = math.exp(-2.0 * alpha * (height - depth))
    root = math.exp(-2.0 * alpha * height)
    return (first - root) / (2.0 * alpha)


def _transformed_node_depths(
    tree,
    model: str,
    parameter: float,
    base_depths: dict[Any, float],
) -> dict[Any, float]:
    height = _tree_height(tree, base_depths)
    if model == "delta":
        _require_ultrametric(tree, base_depths, model)
        return {
            node: height * (depth / height) ** parameter
            for node, depth in base_depths.items()
        }
    if model in {"eb", "acdc"}:
        return {
            node: _exponential_time_transform(depth, parameter)
            for node, depth in base_depths.items()
        }
    if model == "ou":
        height = _require_ultrametric(tree, base_depths, model)
        return {
            node: _ou_ultrametric_depth(depth, height, parameter)
            for node, depth in base_depths.items()
        }
    raise ValueError("Model '{}' is not a node-depth transform.".format(model))


def transformed_edge_variances(
    tree,
    *,
    model: str = "brownian",
    parameter: float | None = None,
    branch_length: str = "original",
) -> dict[Any, float]:
    """Return Brownian edge variances equivalent to a supported model."""
    if branch_length not in {"original", "unit"}:
        raise ValueError("Unsupported branch-length mode: {}.".format(branch_length))
    spec = evolution_model_spec(model)
    if not spec.contrast_supported:
        raise ValueError(
            "Evolutionary model '{}' cannot be represented by local contrasts.".format(
                model
            )
        )
    parameter = validate_evolution_parameter(model, parameter)
    if model == "independent":
        return {
            node: 0.0 if node.is_root or not node.is_leaf else 1.0
            for node in tree.traverse()
        }
    base_edges = _base_edge_lengths(tree, branch_length)
    if model == "brownian":
        return base_edges
    if model == "kappa":
        assert parameter is not None
        try:
            edges = {
                node: 0.0 if node.is_root else base_edges[node] ** parameter
                for node in tree.traverse()
            }
        except OverflowError as exc:
            raise ValueError(
                "Pagel's kappa produces non-finite transformed branches."
            ) from exc
        if not all(math.isfinite(value) for value in edges.values()):
            raise ValueError("Pagel's kappa produces non-finite transformed branches.")
        return edges
    base_depths = _depths_from_edges(tree, base_edges)
    if model == "lambda":
        assert parameter is not None
        transformed = {}
        for node in tree.traverse():
            if node.is_root:
                transformed[node] = 0.0
            elif node.is_leaf:
                transformed[node] = (
                    base_edges[node] + (1.0 - parameter) * base_depths[node.up]
                )
            else:
                transformed[node] = parameter * base_edges[node]
        return transformed
    assert parameter is not None
    transformed_depths = _transformed_node_depths(
        tree,
        model,
        parameter,
        base_depths,
    )
    edges = {
        node: (
            0.0
            if node.is_root
            else transformed_depths[node] - transformed_depths[node.up]
        )
        for node in tree.traverse()
    }
    scale = max(1.0, max(abs(value) for value in edges.values()))
    tolerance = np.finfo(float).eps * scale * max(100, len(edges))
    for node, value in edges.items():
        if value < -tolerance or not math.isfinite(value):
            raise ValueError(
                "Evolutionary model '{}' produced invalid transformed branches.".format(
                    model
                )
            )
        if value < 0.0:
            edges[node] = 0.0
    return edges


def _direct_ou_covariance(
    tree,
    leaves,
    base_depths: dict[Any, float],
    alpha: float,
) -> np.ndarray:
    lca = LcaIndex(tree)
    size = len(leaves)
    covariance = np.zeros((size, size), dtype=float)
    for first, first_leaf in enumerate(leaves):
        for second in range(first, size):
            second_leaf = leaves[second]
            shared = base_depths[lca.common_ancestor(first_leaf, second_leaf)]
            ancestral_variance = -math.expm1(-2.0 * alpha * shared) / (2.0 * alpha)
            independent_distance = (
                base_depths[first_leaf] + base_depths[second_leaf] - 2.0 * shared
            )
            value = ancestral_variance * math.exp(-alpha * independent_distance)
            covariance[first, second] = value
            covariance[second, first] = value
    return covariance


def validate_custom_covariance(matrix, leaf_names) -> np.ndarray:
    leaf_names = [str(name) for name in leaf_names]
    if isinstance(matrix, pd.DataFrame):
        matrix = matrix.copy()
        matrix.index = matrix.index.map(str)
        matrix.columns = matrix.columns.map(str)
        if matrix.index.duplicated().any() or matrix.columns.duplicated().any():
            raise ValueError("Custom covariance contains duplicated species names.")
        missing_rows = sorted(set(leaf_names) - set(matrix.index))
        missing_columns = sorted(set(leaf_names) - set(matrix.columns))
        extra_rows = sorted(set(matrix.index) - set(leaf_names))
        extra_columns = sorted(set(matrix.columns) - set(leaf_names))
        if missing_rows or missing_columns or extra_rows or extra_columns:
            raise ValueError(
                "Custom covariance species must exactly match tree tips "
                "(missing rows={}; missing columns={}; extra rows={}; extra columns={}).".format(
                    ",".join(missing_rows),
                    ",".join(missing_columns),
                    ",".join(extra_rows),
                    ",".join(extra_columns),
                )
            )
        matrix = matrix.loc[leaf_names, leaf_names]
    try:
        values = (
            matrix.to_numpy(dtype=float)
            if isinstance(matrix, pd.DataFrame)
            else np.asarray(matrix, dtype=float)
        )
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Custom covariance must contain numeric values.") from exc
    if values.shape != (len(leaf_names), len(leaf_names)):
        raise ValueError("Custom covariance has the wrong dimensions.")
    if not np.isfinite(values).all():
        raise ValueError("Custom covariance must contain finite values.")
    symmetric = (values + values.T) / 2.0
    scale = max(1.0, float(np.max(np.abs(symmetric))))
    tolerance = np.finfo(float).eps * scale * max(100, len(values))
    if float(np.max(np.abs(values - values.T))) > tolerance:
        raise ValueError("Custom covariance must be symmetric.")
    try:
        np.linalg.cholesky(symmetric)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Custom covariance must be positive definite.") from exc
    return symmetric


def read_custom_covariance(path: str, leaf_names) -> np.ndarray:
    text = read_input_text(path)
    if text.strip() == "":
        raise ValueError("'--evolution-covariance' is empty.")
    header = text.splitlines()[0].split("\t")
    if not header or header[0] != "leaf_name":
        raise ValueError(
            "'--evolution-covariance' must be a wide TSV whose first column is 'leaf_name'."
        )
    duplicated = sorted({name for name in header if header.count(name) > 1})
    if duplicated:
        raise ValueError(
            "'--evolution-covariance' contains duplicated column(s): {}.".format(
                ", ".join(duplicated)
            )
        )
    try:
        dataframe = pd.read_csv(
            StringIO(text),
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_filter=False,
        )
    except Exception as exc:
        raise ValueError("'--evolution-covariance' is not a valid TSV.") from exc
    if dataframe["leaf_name"].duplicated().any():
        raise ValueError("'--evolution-covariance' contains duplicated leaf_name rows.")
    dataframe = dataframe.set_index("leaf_name")
    try:
        numeric = dataframe.apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError("'--evolution-covariance' values must be numeric.") from exc
    return validate_custom_covariance(numeric, leaf_names)


def build_evolutionary_covariance(
    tree,
    leaf_names,
    *,
    model: str = "brownian",
    parameter: float | None = None,
    branch_length: str = "original",
    custom_covariance=None,
) -> np.ndarray:
    """Build a positive-definite tip covariance in the requested tip order."""
    if branch_length not in {"original", "unit"}:
        raise ValueError("Unsupported branch-length mode: {}.".format(branch_length))
    spec = evolution_model_spec(model)
    leaf_names = [str(name) for name in leaf_names]
    leaf_by_name = {str(leaf.name): leaf for leaf in tree.leaves()}
    missing = sorted(set(leaf_names) - set(leaf_by_name))
    if missing:
        raise ValueError(
            "Evolutionary covariance requested absent tree tips: {}.".format(
                ", ".join(missing)
            )
        )
    if model == "custom":
        if custom_covariance is None:
            raise ValueError("Custom evolution model requires a covariance matrix.")
        return validate_custom_covariance(custom_covariance, leaf_names)
    if custom_covariance is not None:
        raise ValueError("A custom covariance is only valid with model 'custom'.")
    parameter = validate_evolution_parameter(model, parameter)
    leaves = [leaf_by_name[name] for name in leaf_names]
    if model == "ou":
        assert parameter is not None
        base_depths = tree_depths(tree, branch_length)
        covariance = _direct_ou_covariance(tree, leaves, base_depths, parameter)
    else:
        edge_variances = transformed_edge_variances(
            tree,
            model=model,
            parameter=parameter,
            branch_length=branch_length,
        )
        depths = _depths_from_edges(tree, edge_variances)
        lca = LcaIndex(tree)
        size = len(leaves)
        covariance = np.zeros((size, size), dtype=float)
        for first, first_leaf in enumerate(leaves):
            for second in range(first, size):
                value = depths[lca.common_ancestor(first_leaf, leaves[second])]
                covariance[first, second] = value
                covariance[second, first] = value
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Evolutionary model '{}' produced a non-positive-definite covariance matrix.".format(
                spec.name
            )
        ) from exc
    return covariance


def optimization_parameterization(
    tree,
    model: str,
    *,
    branch_length: str = "original",
) -> tuple[tuple[float, float], Callable[[float], float]]:
    """Return bounded optimizer coordinates and their decoder."""
    spec = evolution_model_spec(model)
    if spec.parameter_name is None:
        raise ValueError(
            "Evolutionary model '{}' has no shape parameter.".format(model)
        )
    depths = tree_depths(tree, branch_length)
    height = _tree_height(tree, depths)
    if model == "lambda":
        return (0.0, 1.0), float
    if model == "kappa":
        return (0.0, 3.0), float
    if model == "eb":
        return (-50.0 / height, 0.0), float
    if model == "acdc":
        return (-50.0 / height, 50.0 / height), float
    if model == "ou":
        bounds = (math.log(1e-6 / height), math.log(1e3 / height))
        return bounds, lambda value: math.exp(float(value))
    if model == "delta":
        return (math.log(1e-4), math.log(1e4)), lambda value: math.exp(float(value))
    raise ValueError("No optimizer parameterization for model '{}'.".format(model))


def encoded_evolution_parameter(model: str, parameter: float) -> float:
    if model in {"ou", "delta"}:
        return math.log(parameter)
    return parameter


def parameter_near_boundary(
    tree,
    model: str,
    parameter: float,
    *,
    branch_length: str = "original",
) -> bool:
    bounds, _ = optimization_parameterization(
        tree,
        model,
        branch_length=branch_length,
    )
    encoded = encoded_evolution_parameter(model, parameter)
    tolerance = max(1e-5, (bounds[1] - bounds[0]) * 1e-5)
    return encoded <= bounds[0] + tolerance or encoded >= bounds[1] - tolerance
