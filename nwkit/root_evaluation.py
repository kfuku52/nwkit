import csv
import hashlib
import math
from dataclasses import dataclass, field
from io import StringIO
from typing import Any

from nwkit.clade_mapping import canonical_split

SPLIT_ID_PREFIX = "split-sha256:"


@dataclass(frozen=True, slots=True)
class RootingCandidate:
    split: tuple[frozenset[str], frozenset[str]]
    position_kind: str
    position_fraction_from_side_a: float | None
    edge_length: float | None
    score: Any = None
    metrics: dict[str, Any] = field(default_factory=dict, compare=False)
    equivalent_splits: tuple[tuple[frozenset[str], frozenset[str]], ...] = field(
        default_factory=tuple, compare=False
    )


@dataclass(frozen=True, slots=True)
class RootingEvaluation:
    method: str
    selection_basis: str
    score_name: str
    score_unit: str
    candidates: tuple[RootingCandidate, ...]
    evaluated_edges: int | None = None
    tie_rule: str = "not_applicable"
    source: str | None = None
    status: str = "ok"
    message: str = ""

    @classmethod
    def failed(cls, method, message, source=None):
        return cls(
            method=method,
            selection_basis="",
            score_name="",
            score_unit="",
            candidates=tuple(),
            source=source,
            status="failed",
            message=str(message),
        )


def root_split_id(split):
    split = canonical_split(split[0], split[1])
    digest = hashlib.sha256(b"nwkit-split-v1\0")
    for side in split:
        names = tuple(sorted(str(name) for name in side))
        digest.update(len(names).to_bytes(8, byteorder="big"))
        for name in names:
            encoded = name.encode("utf-8")
            digest.update(len(encoded).to_bytes(8, byteorder="big"))
            digest.update(encoded)
    return SPLIT_ID_PREFIX + digest.hexdigest()


def format_taxa(taxa):
    output = StringIO()
    csv.writer(output, lineterminator="").writerow(sorted(str(name) for name in taxa))
    return output.getvalue()


def _finite_root_edge_values(children):
    values = []
    for child in children:
        if child.dist is None:
            return None
        value = float(child.dist)
        if not math.isfinite(value) or value < 0.0:
            return None
        values.append(value)
    return values


def candidate_from_rooted_tree(
    tree,
    *,
    position_kind,
    score=None,
    metrics=None,
):
    children = tree.get_children()
    if len(children) != 2:
        raise ValueError(
            "A rooting candidate must be represented by exactly two root children."
        )
    all_taxa = frozenset(str(name) for name in tree.leaf_names())
    child_taxa = [
        frozenset(str(name) for name in child.leaf_names()) for child in children
    ]
    split = canonical_split(child_taxa[0], child_taxa[1])
    values = _finite_root_edge_values(children)
    edge_length = None
    fraction = None
    if values is not None:
        scale = max(values, default=0.0)
        if scale == 0.0:
            fraction = 0.5
            edge_length = 0.0
        else:
            scaled = [value / scale for value in values]
            scaled_total = math.fsum(scaled)
            side_a_index = 0 if child_taxa[0] == split[0] else 1
            fraction = scaled[side_a_index] / scaled_total
            try:
                candidate_length = math.fsum(values)
            except OverflowError:
                candidate_length = math.inf
            if math.isfinite(candidate_length):
                edge_length = candidate_length
    if position_kind == "edge_unspecified":
        fraction = None
    elif position_kind == "exact" and fraction is None:
        position_kind = "edge_unspecified"
    if (split[0] | split[1]) != all_taxa:
        raise ValueError("Rooting candidate split does not cover every tree tip.")
    return RootingCandidate(
        split=split,
        position_kind=position_kind,
        position_fraction_from_side_a=fraction,
        edge_length=edge_length,
        score=score,
        metrics=dict(metrics or {}),
    )
