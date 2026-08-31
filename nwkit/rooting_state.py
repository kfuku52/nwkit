"""Tree rooting declarations, independent of labels and branching degree.

Only the root owns these properties.  ``nwkit_rooted`` is the portable NHX
representation; private properties retain the source of the interpretation.
"""

import re
import sys
from dataclasses import dataclass
from functools import wraps

ROOTED_PROP = "nwkit_rooted"
ROOTING_SOURCE_PROP = "_nwkit_rooting_source"
ROOTING_DECLARED_PROP = "_nwkit_rooting_declared"
ROOTING_PROPERTIES = frozenset(
    (ROOTED_PROP, ROOTING_SOURCE_PROP, ROOTING_DECLARED_PROP)
)
_ROOTING_TOKEN = re.compile(r"^\s*\[\s*&([RU])\s*\]\s*", re.IGNORECASE)
_STATE_VALUES = {"yes": True, "no": False, "unknown": None}

# Primary inputs share one mode; every auxiliary tree has its own mode.
ROOTING_INPUT_OPTIONS = {
    "infile": "input-rooted",
    "tree": "input-rooted",
    "gene_tree": "input-rooted",
    "gene_tree_ensemble": "input-rooted",
    **{
        name: name.replace("_", "-") + "-rooted"
        for name in (
            "infile2",
            "reference",
            "species_tree",
            "reconciliation_tree",
            "root_source",
            "name_source",
            "support_source",
            "length_source",
            "property_source",
            "densitree_trees",
        )
    },
}


def rooting_option_for_input(option_name):
    role = str(option_name).removeprefix("--").replace("-", "_")
    return "--" + ROOTING_INPUT_OPTIONS.get(role, "input-rooted")


@dataclass(frozen=True)
class RootingInfo:
    rooted: bool | None
    source: str
    declared: str = ""

    @property
    def state(self):
        return (
            "unknown"
            if self.rooted is None
            else ("rooted" if self.rooted else "unrooted")
        )


def _state_text(rooted):
    return "unknown" if rooted is None else ("yes" if rooted else "no")


def _parse_state(value):
    text = str(value).strip().lower()
    if text not in _STATE_VALUES:
        raise ValueError("Root NHX 'nwkit_rooted' must be yes, no, or unknown.")
    return _STATE_VALUES[text]


def topology_rooting(tree):
    """Legacy inference, without claiming an unmarked polytomy is unrooted."""
    return True if len(tree.get_children()) <= 2 else None


def get_rooting_info(tree):
    if ROOTED_PROP in tree.props:
        rooted = _parse_state(tree.props[ROOTED_PROP])
        return RootingInfo(
            rooted,
            str(tree.props.get(ROOTING_SOURCE_PROP, "nhx")),
            str(tree.props.get(ROOTING_DECLARED_PROP, "")),
        )
    rooted = topology_rooting(tree)
    return RootingInfo(rooted, "topology" if rooted else "unknown")


def is_rooted(tree):
    return get_rooting_info(tree).rooted is True


def set_rooting_info(tree, rooted, source="operation", declared=""):
    if rooted is not None and not isinstance(rooted, bool):
        raise ValueError("Rooting state must be True, False, or None.")
    tree.props[ROOTED_PROP] = _state_text(rooted)
    tree.props[ROOTING_SOURCE_PROP] = source
    tree.props[ROOTING_DECLARED_PROP] = declared


def copy_rooting_info(source, target):
    info = get_rooting_info(source)
    set_rooting_info(target, info.rooted, info.source, info.declared)


def inherit_subtree_rooting(source, target):
    """A standalone subtree copy keeps its containing tree's interpretation."""
    if source.up is not None:
        copy_rooting_info(source.root, target)
    return target


def rooting_operation(function):
    """Mark a successfully rooted tree, also for the Python rooting API."""

    @wraps(function)
    def wrapped(*args, **kwargs):
        result = function(*args, **kwargs)
        tree = result[0] if isinstance(result, tuple) else result
        for node in tree.traverse():
            for prop in ROOTING_PROPERTIES:
                node.props.pop(prop, None)
        set_rooting_info(tree, True)
        return result

    return wrapped


def extract_rooting_token(text):
    """Strip only leading rooting tokens, never tokens inside node labels."""
    declared = None
    while match := _ROOTING_TOKEN.match(text):
        value = match.group(1).upper() == "R"
        if declared is not None and value != declared:
            raise ValueError("Conflicting [&R] and [&U] rooting declarations.")
        declared = value
        text = text[match.end() :]
    return text, declared


def apply_input_rooting(tree, rooted="auto", declared=None, warn=True):
    mode = str(rooted).lower()
    if mode not in {"auto", "yes", "no"}:
        raise ValueError("Input rootedness must be auto, yes, or no.")
    info = get_rooting_info(tree)
    has_nhx = ROOTED_PROP in tree.props
    for node in tree.traverse():
        if node is not tree and ROOTED_PROP in node.props:
            raise ValueError("NHX 'nwkit_rooted' is allowed only on the root node.")
    if declared is not None:
        if has_nhx and info.rooted is not declared:
            raise ValueError("Conflicting rooting token and root NHX declaration.")
        info = RootingInfo(declared, "marker", _state_text(declared))
    elif has_nhx:
        info = RootingInfo(info.rooted, "nhx", _state_text(info.rooted))
    if mode != "auto":
        forced = mode == "yes"
        if (
            warn
            and info.source in {"marker", "nhx"}
            and info.rooted is not None
            and forced != info.rooted
        ):
            sys.stderr.write(
                "Warning: forcing input to {} overrides its explicit {} declaration.\n".format(
                    "rooted" if forced else "unrooted", info.state
                )
            )
        info = RootingInfo(forced, "override", info.declared)
    set_rooting_info(tree, info.rooted, info.source, info.declared)
    return info


def rooting_needs_serialization(tree):
    info = get_rooting_info(tree)
    return info.source in {
        "marker",
        "nhx",
        "override",
    } or info.rooted is not topology_rooting(tree)


def require_rooted(tree, message, option="--input-rooted"):
    if not is_rooted(tree):
        raise ValueError(
            "{} Rooting state is {}. Supply [&R] or use '{} yes' to treat the "
            "current top-level node as the root; this does not reroot the tree.".format(
                message, get_rooting_info(tree).state, option
            )
        )
