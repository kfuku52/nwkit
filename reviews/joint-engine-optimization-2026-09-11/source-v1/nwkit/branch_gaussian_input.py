"""Strict TSV assignments for fixed scalar branch-specific Gaussian models."""

import csv
import math
import sys
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from nwkit.branch_gaussian import (
    BranchGaussianModel,
    BrownianBranch,
    GaussianJump,
    OUBranch,
)
from nwkit.util import assign_branch_ids

_PARAMETER_COLUMNS = frozenset(
    {"sigma2", "alpha", "theta", "jump_mean", "jump_variance"}
)


@dataclass(frozen=True)
class BranchGaussianAssignment:
    models_by_branch_id: Mapping[int, BranchGaussianModel]
    regime_by_branch_id: Mapping[int, str] | None = None

    def model_rows(self):
        """Return normalized direct assignments suitable for a model TSV."""
        rows = []
        for identifier, model in sorted(self.models_by_branch_id.items()):
            diffusion = model.diffusion
            rows.append(
                {
                    "branch_id": identifier,
                    "model": "JUMP"
                    if diffusion is None
                    else ("OU" if isinstance(diffusion, OUBranch) else "BM"),
                    "sigma2": None if diffusion is None else diffusion.variance_rate,
                    "alpha": diffusion.alpha
                    if isinstance(diffusion, OUBranch)
                    else None,
                    "theta": diffusion.optimum
                    if isinstance(diffusion, OUBranch)
                    else None,
                    "jump_mean": None if model.jump is None else model.jump.mean,
                    "jump_variance": None
                    if model.jump is None
                    else model.jump.variance,
                }
            )
        return rows


def _read_rows(path, *, required, allowed):
    context = (
        nullcontext(sys.stdin)
        if str(path) == "-"
        else Path(path).open(encoding="utf-8-sig", newline="")
    )
    with context as stream:
        reader = csv.reader(stream, delimiter="\t", strict=True)
        try:
            header = next(reader, None)
            if not header:
                raise ValueError(f"Empty model/assignment TSV: {path}.")
            # Also accept a BOM on stdin, where no utf-8-sig decoder is applied.
            header[0] = header[0].removeprefix("\ufeff")
            if len(set(header)) != len(header):
                raise ValueError(f"Duplicate TSV column names in {path}.")
            missing, extra = set(required) - set(header), set(header) - set(allowed)
            if missing or extra:
                raise ValueError(
                    f"Invalid TSV columns in {path}: missing={sorted(missing)}, extra={sorted(extra)}."
                )
            rows = []
            for values in reader:
                if not values:
                    raise ValueError(
                        f"Empty TSV row in {path}, line {reader.line_num}."
                    )
                if len(values) != len(header):
                    raise ValueError(
                        f"TSV row width differs from its header in {path}, line {reader.line_num}."
                    )
                row = {
                    key: value.strip()
                    for key, value in zip(header, values, strict=True)
                }
                rows.append((reader.line_num, row))
            return rows
        except csv.Error as exc:
            raise ValueError(
                f"Malformed TSV in {path}, line {reader.line_num}."
            ) from exc


def _branch_id(text, context):
    try:
        value = int(text)
    except ValueError as exc:
        raise ValueError(f"Invalid branch_id '{text}' in {context}.") from exc
    if value <= 0 or text not in {str(value), f"+{value}"}:
        raise ValueError(
            f"branch_id must be a positive integer (root 0 is excluded) in {context}."
        )
    return value


def _number(row, key, context):
    text = row.get(key, "")
    if not text:
        raise ValueError(f"Missing {key} in {context}.")
    try:
        value = float(text)
    except ValueError as exc:
        raise ValueError(f"Invalid numeric {key} in {context}: {text}.") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {key} in {context}.")
    return value


def _model(row, context):
    name = row["model"]
    if name not in {"BM", "OU", "JUMP"}:
        raise ValueError(f"model must be BM, OU, or JUMP in {context}.")
    unused = (
        {"alpha", "theta"}
        if name == "BM"
        else ({"sigma2", "alpha", "theta"} if name == "JUMP" else set())
    )
    if any(row.get(key, "") for key in unused):
        raise ValueError(
            f"Parameters not used by {name} must be blank in {context}: {sorted(unused)}."
        )
    try:
        diffusion: BrownianBranch | OUBranch | None = None
        if name == "BM":
            diffusion = BrownianBranch(_number(row, "sigma2", context))
        elif name == "OU":
            diffusion = OUBranch(
                _number(row, "alpha", context),
                _number(row, "sigma2", context),
                _number(row, "theta", context),
            )
        jump = None
        if name == "JUMP" or row.get("jump_mean", "") or row.get("jump_variance", ""):
            jump = GaussianJump(
                _number(row, "jump_mean", context),
                _number(row, "jump_variance", context),
            )
        return BranchGaussianModel(diffusion, jump)
    except ValueError as exc:
        raise ValueError(f"Invalid model in {context}: {exc}") from exc


def _read_models(path, key):
    result = {}
    rows = _read_rows(
        path, required={key, "model"}, allowed={key, "model", *_PARAMETER_COLUMNS}
    )
    for line, row in rows:
        context = f"{path}, line {line}"
        identifier = _branch_id(row[key], context) if key == "branch_id" else row[key]
        if not identifier:
            raise ValueError(f"Empty {key} in {context}.")
        if identifier in result:
            raise ValueError(f"Duplicate {key} '{identifier}' in {context}.")
        result[identifier] = _model(row, context)
    return result


def _read_regimes(path):
    result = {}
    rows = _read_rows(
        path, required={"branch_id", "regime"}, allowed={"branch_id", "regime"}
    )
    for line, row in rows:
        context = f"{path}, line {line}"
        identifier = _branch_id(row["branch_id"], context)
        if identifier in result:
            raise ValueError(f"Duplicate branch_id '{identifier}' in {context}.")
        if not row["regime"]:
            raise ValueError(f"Empty regime in {context}.")
        result[identifier] = row["regime"]
    return result


def load_branch_gaussian_assignment(
    tree, *, branch_models=None, branch_regimes=None, regime_models=None
):
    """Load direct models or two regime tables; cover every non-root branch."""
    if not tree.is_root:
        raise ValueError("tree must be a root node.")
    paths = (branch_models, branch_regimes, regime_models)
    if sum(str(path) == "-" for path in paths) > 1:
        raise ValueError("STDIN can supply only one model/assignment table.")
    if any(path == "" for path in paths):
        raise ValueError("Model/assignment TSV paths must not be empty.")
    if branch_models is not None:
        if branch_regimes is not None or regime_models is not None:
            raise ValueError(
                "Use --branch-models alone, or --branch-regimes with --regime-models."
            )
        models, regimes = _read_models(branch_models, "branch_id"), None
    else:
        if branch_regimes is None or regime_models is None:
            raise ValueError("Both --branch-regimes and --regime-models are required.")
        regimes = _read_regimes(branch_regimes)
        definitions = _read_models(regime_models, "regime")
        if set(regimes.values()) != set(definitions):
            raise ValueError(
                "Regime definitions must exactly match the used regimes; missing or unused definitions found."
            )
        models = {
            identifier: definitions[regime] for identifier, regime in regimes.items()
        }
    expected = {
        identifier
        for node, identifier in assign_branch_ids(tree).items()
        if not node.is_root
    }
    if set(models) != expected:
        raise ValueError(
            f"Model assignments must cover every non-root branch ID: missing={sorted(expected - set(models))}, extra={sorted(set(models) - expected)}."
        )
    return BranchGaussianAssignment(models, regimes)
