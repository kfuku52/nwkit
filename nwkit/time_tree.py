"""MCMCtree-specific time-tree parsing, summaries, and annotations."""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass
from io import StringIO

import numpy as np
import pandas as pd

from nwkit.rooting_state import require_rooted
from nwkit.util import read_input_text

AGE_PROPERTY_NAMES = (
    "age",
    "age_mean",
    "age_median",
    "age_ci_low",
    "age_ci_high",
    "age_ci_kind",
    "age_ci_level",
    "mcmctree_node_id",
)

CALIBRATION_PROPERTY_NAMES = (
    "calibration_raw",
    "calibration_type",
    "calibration_lower",
    "calibration_upper",
    "calibration_lower_tail",
    "calibration_upper_tail",
    "calibration_offset",
    "calibration_scale",
)

_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
_POINT_PATTERN = re.compile(r"^@\s*({})$".format(_NUMBER))
_BOUND_PATTERN = re.compile(
    r"^B\(\s*({0})\s*,\s*({0})\s*,\s*({0})\s*,\s*({0})\s*\)$".format(_NUMBER),
    flags=re.IGNORECASE,
)
_LOWER_PATTERN = re.compile(
    r"^L\(\s*({0})\s*,\s*({0})\s*,\s*({0})\s*,\s*({0})\s*\)$".format(_NUMBER),
    flags=re.IGNORECASE,
)
_UPPER_PATTERN = re.compile(
    r"^U\(\s*({0})\s*,\s*({0})\s*\)$".format(_NUMBER),
    flags=re.IGNORECASE,
)
_LEGACY_BOUNDS_PATTERN = re.compile(
    r"^>\s*({0})\s*<\s*({0})$".format(_NUMBER),
)
_LEGACY_LOWER_PATTERN = re.compile(r"^>\s*({})$".format(_NUMBER))
_LEGACY_UPPER_PATTERN = re.compile(r"^<\s*({})$".format(_NUMBER))


def _validated_calibration(record):
    numeric_values = [
        value for key, value in record.items() if key not in {"raw", "type"}
    ]
    if not all(math.isfinite(float(value)) for value in numeric_values):
        raise ValueError("MCMCtree calibration values must be finite.")
    for key in ("lower", "upper"):
        if key in record and float(record[key]) < 0.0:
            raise ValueError("MCMCtree calibration ages must be non-negative.")
    if (
        "lower" in record
        and "upper" in record
        and float(record["lower"]) > float(record["upper"])
    ):
        raise ValueError(
            "MCMCtree calibration lower age must not exceed its upper age."
        )
    for key in ("lower_tail", "upper_tail"):
        if key in record and not 0.0 < float(record[key]) < 1.0:
            raise ValueError(
                "MCMCtree calibration tail probabilities must be between 0 and 1."
            )
    if "offset" in record and float(record["offset"]) < 0.0:
        raise ValueError("MCMCtree calibration offsets must be non-negative.")
    if "scale" in record and float(record["scale"]) <= 0.0:
        raise ValueError("MCMCtree calibration scales must be positive.")
    return record


@dataclass(frozen=True)
class MCMCtreePosterior:
    """Validated MCMCtree node-age samples tied to one fixed topology."""

    node_ids: tuple
    values: np.ndarray
    generations: np.ndarray | None

    @property
    def sample_count(self):
        return int(self.values.shape[0])


def parse_mcmctree_calibration(value):
    """Return a normalized calibration record, or ``None`` for a normal label."""

    if value in (None, ""):
        return None
    text = str(value).strip().strip("'").strip('"')
    match = _POINT_PATTERN.fullmatch(text)
    if match:
        age = float(match.group(1))
        return _validated_calibration(
            {
                "raw": text,
                "type": "point",
                "lower": age,
                "upper": age,
            }
        )
    match = _BOUND_PATTERN.fullmatch(text)
    if match:
        return _validated_calibration(
            {
                "raw": text,
                "type": "bounded",
                "lower": float(match.group(1)),
                "upper": float(match.group(2)),
                "lower_tail": float(match.group(3)),
                "upper_tail": float(match.group(4)),
            }
        )
    match = _LOWER_PATTERN.fullmatch(text)
    if match:
        return _validated_calibration(
            {
                "raw": text,
                "type": "lower",
                "lower": float(match.group(1)),
                "offset": float(match.group(2)),
                "scale": float(match.group(3)),
                "lower_tail": float(match.group(4)),
            }
        )
    match = _UPPER_PATTERN.fullmatch(text)
    if match:
        return _validated_calibration(
            {
                "raw": text,
                "type": "upper",
                "upper": float(match.group(1)),
                "upper_tail": float(match.group(2)),
            }
        )
    match = _LEGACY_BOUNDS_PATTERN.fullmatch(text)
    if match:
        return _validated_calibration(
            {
                "raw": text,
                "type": "bounded",
                "lower": float(match.group(1)),
                "upper": float(match.group(2)),
            }
        )
    match = _LEGACY_LOWER_PATTERN.fullmatch(text)
    if match:
        return _validated_calibration(
            {
                "raw": text,
                "type": "lower",
                "lower": float(match.group(1)),
            }
        )
    match = _LEGACY_UPPER_PATTERN.fullmatch(text)
    if match:
        return _validated_calibration(
            {
                "raw": text,
                "type": "upper",
                "upper": float(match.group(1)),
            }
        )
    return None


def annotate_mcmctree_calibrations(tree):
    """Attach machine-readable properties to calibration-labelled nodes."""

    annotated = 0
    for node in tree.traverse():
        calibration = parse_mcmctree_calibration(node.name)
        if calibration is None:
            continue
        annotated += 1
        for key, value in calibration.items():
            node.props["calibration_{}".format(key)] = value
    return annotated


def paml_node_mapping(tree):
    """Map PAML node numbers to ETE nodes for an MCMCtree species topology.

    PAML numbers tips in their input order from 1, followed by internal nodes
    in root-first preorder.  This mapping agrees with the numbered topology
    printed in MCMCtree's species-tree output block.
    """

    leaves = list(tree.leaves())
    internal = [node for node in tree.traverse(strategy="preorder") if not node.is_leaf]
    mapping = {index: node for index, node in enumerate(leaves, start=1)}
    mapping.update(
        {len(leaves) + index: node for index, node in enumerate(internal, start=1)}
    )
    return mapping


def _read_mcmctree_dataframe(source):
    is_file = source != "-" and os.path.isfile(source)
    if is_file:
        if os.path.getsize(source) == 0:
            raise ValueError("MCMCtree posterior file is empty.")
        table_source = source
    else:
        text = read_input_text(source)
        if text.strip() == "":
            raise ValueError("MCMCtree posterior file is empty.")
        table_source = StringIO(text)
    try:
        dataframe = pd.read_csv(table_source, sep=r"\s+", comment="#")
    except Exception as exc:
        raise ValueError("Failed to parse the MCMCtree posterior table.") from exc
    if dataframe.empty:
        raise ValueError("MCMCtree posterior table contains no samples.")
    if dataframe.columns.duplicated().any():
        duplicates = dataframe.columns[dataframe.columns.duplicated()].tolist()
        raise ValueError(
            "MCMCtree posterior table contains duplicate columns: {}".format(
                ", ".join(str(value) for value in duplicates),
            )
        )
    return dataframe


def read_mcmctree_posterior(source, tree, burnin=0, thin=1):
    """Read and validate MCMCtree ``mcmc.txt`` node-age samples."""

    if any(
        len(node.get_children()) != 2 for node in tree.traverse() if not node.is_leaf
    ):
        raise ValueError("MCMCtree posterior topology must be rooted and binary.")
    require_rooted(tree, "MCMCtree posterior topology must be rooted and binary.")
    burnin = int(burnin)
    thin = int(thin)
    if burnin < 0:
        raise ValueError("'--posterior-burnin' must be zero or greater.")
    if thin < 1:
        raise ValueError("'--posterior-thin' must be at least 1.")
    dataframe = _read_mcmctree_dataframe(source)
    mapping = paml_node_mapping(tree)
    internal_ids = tuple(
        node_id for node_id, node in mapping.items() if not node.is_leaf
    )
    expected_columns = ["t_n{}".format(node_id) for node_id in internal_ids]
    missing = [column for column in expected_columns if column not in dataframe.columns]
    if missing:
        raise ValueError(
            "MCMCtree posterior is missing node-age column(s): {}".format(
                ", ".join(missing)
            )
        )
    retained = dataframe.iloc[burnin::thin]
    if retained.empty:
        raise ValueError(
            "No MCMCtree posterior samples remain after burn-in and thinning."
        )
    try:
        values = retained.loc[:, expected_columns].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "MCMCtree node-age columns must contain numeric values."
        ) from exc
    if not np.isfinite(values).all():
        raise ValueError("MCMCtree node-age samples must be finite.")
    if (values < 0.0).any():
        raise ValueError("MCMCtree node-age samples must be non-negative.")
    id_to_column = {node_id: index for index, node_id in enumerate(internal_ids)}
    node_to_id = {node: node_id for node_id, node in mapping.items()}
    for child in tree.traverse():
        if child.is_root:
            continue
        parent = child.up
        parent_ages = values[:, id_to_column[node_to_id[parent]]]
        child_ages = (
            np.zeros(values.shape[0], dtype=float)
            if child.is_leaf
            else values[:, id_to_column[node_to_id[child]]]
        )
        invalid = np.flatnonzero(parent_ages + 1e-12 < child_ages)
        if invalid.size:
            raise ValueError(
                "MCMCtree sample {} makes child node {} older than its parent {}.".format(
                    burnin + (int(invalid[0]) * thin) + 1,
                    node_to_id[child],
                    node_to_id[parent],
                )
            )
    generations = None
    if "Gen" in retained.columns:
        try:
            generations = retained["Gen"].to_numpy(dtype=float)
        except (TypeError, ValueError):
            generations = None
    return MCMCtreePosterior(
        node_ids=internal_ids,
        values=values,
        generations=generations,
    )


def _hpd_interval(values, level):
    ordered = np.sort(np.asarray(values, dtype=float))
    count = len(ordered)
    alpha = 1.0 - float(level)
    initial_low = int(count * alpha / 2.0 + 1e-12)
    initial_high = int(count * (1.0 - alpha / 2.0) + 1e-12)
    index_width = max(0, min(count - 1, initial_high - initial_low))
    if index_width >= count - 1:
        return float(ordered[0]), float(ordered[-1])
    widths = ordered[index_width:] - ordered[: count - index_width]
    start = int(np.argmin(widths))
    return float(ordered[start]), float(ordered[start + index_width])


def _credible_interval(values, level, kind):
    kind = str(kind).strip().lower()
    if kind == "hpd":
        return _hpd_interval(values, level)
    if kind == "equal-tail":
        tail = (1.0 - float(level)) / 2.0
        ordered = np.sort(np.asarray(values, dtype=float))
        low_index = min(len(ordered) - 1, int(len(ordered) * tail + 1e-12))
        high_index = min(
            len(ordered) - 1,
            int(len(ordered) * (1.0 - tail) + 1e-12),
        )
        return float(ordered[low_index]), float(ordered[high_index])
    raise ValueError("'--posterior-ci' must be hpd or equal-tail.")


def summarize_mcmctree_posterior(tree, posterior, point="mean", ci="hpd", level=0.95):
    """Attach posterior node ages and convert them to dated branch lengths."""

    point = str(point).strip().lower()
    ci = str(ci).strip().lower()
    level = float(level)
    if point not in {"mean", "median"}:
        raise ValueError("'--posterior-point' must be mean or median.")
    if ci not in {"hpd", "equal-tail"}:
        raise ValueError("'--posterior-ci' must be hpd or equal-tail.")
    if not 0.0 < level < 1.0:
        raise ValueError("'--posterior-ci-level' must be between 0 and 1.")
    mapping = paml_node_mapping(tree)
    column_by_id = {node_id: index for index, node_id in enumerate(posterior.node_ids)}
    age_by_node = {}
    for node_id, node in mapping.items():
        node.props["mcmctree_node_id"] = int(node_id)
        if node.is_leaf:
            age_by_node[node] = 0.0
            node.props["age"] = 0.0
            continue
        values = posterior.values[:, column_by_id[node_id]]
        mean = float(np.mean(values))
        median = float(np.median(values))
        low, high = _credible_interval(values, level=level, kind=ci)
        selected = mean if point == "mean" else median
        age_by_node[node] = selected
        node.props.update(
            {
                "age": selected,
                "age_mean": mean,
                "age_median": median,
                "age_ci_low": low,
                "age_ci_high": high,
                "age_ci_kind": "HPD" if ci == "hpd" else "equal-tail",
                "age_ci_level": level,
            }
        )
    for node in tree.traverse():
        if node.is_root:
            node.dist = None
            continue
        distance = age_by_node[node.up] - age_by_node[node]
        if distance < -1e-10:
            raise ValueError(
                "Posterior point estimates produce a negative dated branch length."
            )
        node.dist = max(0.0, float(distance))
    return tree


def posterior_sample_tree(tree, posterior, sample_index):
    """Return one branch-length tree reconstructed from a posterior row."""

    sample_index = int(sample_index)
    copied = tree.copy()
    mapping = paml_node_mapping(copied)
    column_by_id = {node_id: index for index, node_id in enumerate(posterior.node_ids)}
    age_by_node = {}
    for node_id, node in mapping.items():
        age_by_node[node] = (
            0.0
            if node.is_leaf
            else float(posterior.values[sample_index, column_by_id[node_id]])
        )
    for node in copied.traverse():
        if node.is_root:
            node.dist = None
        else:
            node.dist = max(0.0, age_by_node[node.up] - age_by_node[node])
    return copied


def infer_node_ages_from_branch_lengths(tree, tolerance=1e-5):
    """Annotate an ultrametric dated tree with ages measured back from tips."""

    age_by_node = {}
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            age_by_node[node] = 0.0
            node.props.setdefault("age", 0.0)
            continue
        child_ages = [
            age_by_node[child] + (0.0 if child.dist is None else float(child.dist))
            for child in node.get_children()
        ]
        if max(child_ages) - min(child_ages) > max(
            tolerance, tolerance * max(child_ages)
        ):
            raise ValueError(
                "Dated-tree branch lengths are not ultrametric at an internal node."
            )
        age = float(sum(child_ages) / len(child_ages))
        age_by_node[node] = age
        node.props.setdefault("age", age)
        if "age_ci_low" in node.props and "age_mean" not in node.props:
            node.props["age_mean"] = age
    return tree


def prepare_time_tree_annotations(tree):
    """Recognize constraints and infer ages for PAML/FigTree dated trees."""

    calibration_count = annotate_mcmctree_calibrations(tree)
    ci_count = sum("age_ci_low" in node.props for node in tree.traverse())
    if ci_count:
        infer_node_ages_from_branch_lengths(tree)
    return {
        "calibration_count": calibration_count,
        "credible_interval_count": ci_count,
    }
