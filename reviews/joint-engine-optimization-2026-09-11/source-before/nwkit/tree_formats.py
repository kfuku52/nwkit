"""Loss-aware tree container and annotation syntax shared by NWKIT readers.

The lexer keeps quoted names and comments distinct from branch lengths. ETE
remains the tree parser; this module translates container/annotation syntax.
"""

import re
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from typing import Any

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
TREE_ASSIGNMENT = re.compile(
    r"""^\s*(?:UTREE|TREE)\s+(?:\*\s+)?(?:'(?:[^']|'')*'|"(?:[^"]|"")*"|[^\s=]+)\s*=\s*""",
    re.IGNORECASE,
)
MCMCTREE_MARKER = re.compile(r"^\s*Species tree for FigTree\b[^\n]*$", re.MULTILINE)
CI_KEY = re.compile(rf"(?:height_)?({NUMBER})%(?:_?(HPD))?", re.IGNORECASE)
AGE_FIELDS = frozenset({"age", "age_mean", "age_median", "age_ci_low", "age_ci_high"})
CI_FIELDS = frozenset({"age_ci_low", "age_ci_high", "age_ci_kind", "age_ci_level"})
HEIGHT_FIELDS = {
    "height": "age",
    "height_mean": "age_mean",
    "height_median": "age_median",
}


@dataclass(frozen=True)
class Token:
    kind: str
    text: str


def tokens(text):
    """Tokenize punctuation only outside quoted names and nested comments."""
    result = []
    i = 0
    while i < len(text):
        char = text[i]
        if char.isspace():
            i += 1
            continue
        start = i
        if char in "'\"":
            i = _quote_end(text, i)
            kind = "name"
        elif char == "[":
            i = _comment_end(text, i)
            kind = "comment"
        elif char in "(),:;":
            i += 1
            kind = char
        else:
            i += 1
            while (
                i < len(text) and not text[i].isspace() and text[i] not in "(),'\"[]:;"
            ):
                i += 1
            kind = "name"
        result.append(Token(kind, text[start:i]))
    return result


def _quote_end(text, start):
    quote = text[start]
    i = start + 1
    while i < len(text):
        if text[i] == quote:
            if i + 1 < len(text) and text[i + 1] == quote:
                i += 2
                continue
            return i + 1
        i += 1
    raise ValueError("Unterminated quoted tree label.")


def _comment_end(text, start):
    depth = 1
    for i in range(start + 1, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i + 1
    raise ValueError("Unterminated tree comment.")


def statement_spans(text):
    """Find semicolon-terminated records without splitting names/comments."""
    start = i = 0
    while i < len(text):
        if text[i] in "'\"":
            i = _quote_end(text, i)
        elif text[i] == "[":
            i = _comment_end(text, i)
        elif text[i] == ";":
            yield start, i + 1
            i += 1
            start = i
        else:
            i += 1
    if text[start:].strip():
        raise ValueError("Input ended before a terminal semicolon.")


def _nexus_trees(text):
    # The leading report prose of a recovered public summary is not NEXUS.
    marker = re.search(r"(?im)^\s*#NEXUS\b", text)
    if marker:
        text = text[marker.end() :]
    result = []
    for start, end in statement_spans(text):
        statement = text[start:end].strip()
        while statement.startswith("["):
            statement = statement[_comment_end(statement, 0) :].lstrip()
        if re.match(r"TRANSLATE\b", statement, re.IGNORECASE):
            raise ValueError(
                "NEXUS TRANSLATE tables are not supported; provide direct tip labels."
            )
        match = TREE_ASSIGNMENT.match(statement)
        if match:
            result.append(statement[match.end() :])
    if not result:
        raise ValueError(
            "NEXUS input did not contain a complete TREE or UTREE statement."
        )
    return result


def _mcmctree_trees(text):
    marker = MCMCTREE_MARKER.search(text)
    if marker is None:
        raise ValueError(
            "MCMCtree output did not contain its FigTree species-tree block."
        )
    if MCMCTREE_MARKER.search(text, marker.end()):
        raise ValueError(
            "Multiple MCMCtree FigTree blocks are ambiguous; provide one block."
        )
    result: list[str] = []
    pending: list[str] = []
    for line in text[marker.end() :].splitlines():
        if not pending and not line.strip():
            continue
        if not pending and not line.lstrip().startswith("("):
            break
        pending.append(line)
        statement = "\n".join(pending)
        try:
            spans = list(statement_spans(statement))
        except ValueError:
            continue
        if len(spans) != 1:
            raise ValueError("Expected one tree per MCMCtree FigTree statement.")
        result.append(statement.strip())
        pending = []
    if pending or not result:
        raise ValueError(
            "MCMCtree output did not contain a complete FigTree species tree."
        )
    return result


def read_container(text, source="auto"):
    """Return the detected format and every tree in its original input order."""
    if source not in {"auto", "newick", "nhx", "figtree", "mcmctree-output"}:
        raise ValueError(f"Unknown input format: {source}")
    text = str(text).strip()
    if source == "auto":
        if MCMCTREE_MARKER.search(text):
            source = "mcmctree-output"
        elif re.search(r"(?im)^\s*(?:#NEXUS\b|(?:UTREE|TREE)\s+[^=]+?=)", text):
            source = "figtree"
        else:
            source = "nhx" if "[&&NHX:" in text else "newick"
    if source == "mcmctree-output":
        return source, _mcmctree_trees(text)
    if source == "figtree":
        return source, _nexus_trees(text)
    trees = [text[start:end].strip() for start, end in statement_spans(text)]
    if not trees:
        raise ValueError("No input tree was found.")
    return source, trees


def finite_decimal(value, *, positive=False):
    try:
        number = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"Invalid numeric time value: {value}") from exc
    if not number.is_finite() or (positive and number <= 0):
        raise ValueError(
            "Time factors must be finite and positive."
            if positive
            else "Time values must be finite."
        )
    return number


def decimal_text(number):
    if number == 0:
        return "0"
    return format(number.normalize(), "f")


def scaled_time(value, factor):
    number = finite_decimal(value)
    if number < 0:
        raise ValueError("Time values must be non-negative.")
    if factor == 1:
        return str(value)
    with localcontext() as context:
        context.prec = max(
            28, len(number.as_tuple().digits) + len(factor.as_tuple().digits)
        )
        return decimal_text(number * factor)


def _annotation_fields(body):
    """Split NEXUS attribute fields, respecting quoted values and intervals."""
    result = []
    start = i = depth = 0
    while i < len(body):
        if body[i] in "'\"":
            i = _quote_end(body, i)
            continue
        if body[i] == "{":
            depth += 1
        elif body[i] == "}":
            depth -= 1
            if depth < 0:
                raise ValueError("Malformed tree annotation interval.")
        elif body[i] == "," and depth == 0:
            result.append(body[start:i].strip())
            start = i + 1
        i += 1
    if depth:
        raise ValueError("Malformed tree annotation interval.")
    result.append(body[start:].strip())
    return result


def _add_attribute(attributes, key, value):
    if key in attributes:
        raise ValueError(f"Duplicate tree annotation: {key}")
    attributes[key] = value


def annotation_attributes(comment):
    """Return normalized attributes and any uninterpreted comment fields."""
    attributes: dict[str, str] = {}
    if comment.startswith("[&&NHX:"):
        for field in comment[7:-1].split(":"):
            if "=" not in field:
                raise ValueError("Malformed NHX property.")
            key, value = field.split("=", 1)
            _add_attribute(attributes, key, value)
        return attributes, []
    if not comment.startswith("[&") or comment.upper() in {"[&R]", "[&U]"}:
        return attributes, [comment]
    other = []
    for field in _annotation_fields(comment[2:-1]):
        key, sep, value = field.partition("=")
        key, value = key.strip(), value.strip()
        ci = CI_KEY.fullmatch(key)
        if ci and sep:
            interval = re.fullmatch(r"\{\s*([^,{}]+)\s*,\s*([^,{}]+)\s*\}", value)
            if interval is None:
                raise ValueError("Malformed node-age credible interval.")
            fields = {
                "age_ci_low": interval[1].strip(),
                "age_ci_high": interval[2].strip(),
                "age_ci_kind": "HPD" if ci[2] else "equal-tail",
                "age_ci_level": decimal_text(finite_decimal(ci[1]) / 100),
            }
            for name, number in fields.items():
                _add_attribute(attributes, name, number)
        elif sep and (key in AGE_FIELDS or key in HEIGHT_FIELDS or key in CI_FIELDS):
            _add_attribute(attributes, HEIGHT_FIELDS.get(key, key), value)
        else:
            other.append(field)
    return attributes, (["[&" + ",".join(other) + "]"] if other else [])


def validate_age_attributes(attributes):
    for key in AGE_FIELDS & attributes.keys():
        if finite_decimal(attributes[key]) < 0:
            raise ValueError(f"Node-age annotation {key} must be non-negative.")
    present = CI_FIELDS & attributes.keys()
    if not present:
        return
    if present != CI_FIELDS:
        raise ValueError(
            "Node-age intervals require low, high, kind and level attributes."
        )
    if finite_decimal(attributes["age_ci_low"]) > finite_decimal(
        attributes["age_ci_high"]
    ):
        raise ValueError("Node-age interval lower bound exceeds its upper bound.")
    if not 0 < finite_decimal(attributes["age_ci_level"]) < 1:
        raise ValueError("Node-age interval level must be between zero and one.")
    if attributes["age_ci_kind"].lower() not in {"hpd", "equal-tail"}:
        raise ValueError("Node-age interval kind must be HPD or equal-tail.")


def _render_attributes(attributes, output):
    if not attributes:
        return ""
    if output == "newick":
        raise ValueError(
            "Plain Newick cannot retain these properties: "
            + ", ".join(sorted(attributes))
            + ". Use --to nhx/figtree to retain them, --properties drop to discard "
            "all properties, or --age-ci drop for age intervals only."
        )
    if output == "nhx":
        return (
            "[&&NHX:"
            + ":".join(f"{key}={value}" for key, value in attributes.items())
            + "]"
        )
    attributes = dict(attributes)
    comments = []
    if CI_FIELDS <= attributes.keys():
        level = decimal_text(finite_decimal(attributes.pop("age_ci_level")) * 100)
        kind = attributes.pop("age_ci_kind")
        key = level + ("%HPD" if kind.lower() == "hpd" else "%")
        low, high = attributes.pop("age_ci_low"), attributes.pop("age_ci_high")
        comments.append(f"[&{key}={{{low},{high}}}]")
    ages = [
        f"{key}={attributes.pop(key)}" for key in sorted(AGE_FIELDS & attributes.keys())
    ]
    if ages:
        comments.append("[&" + ",".join(ages) + "]")
    # NHX is itself a Newick comment and remains valid within NEXUS. Preserve
    # non-time attributes verbatim rather than guessing FigTree equivalents.
    comments.append(_render_attributes(attributes, "nhx"))
    return "".join(comments)


def _node_suffix(
    parts,
    output,
    factor,
    age_ci,
    parser,
    node_label="",
    properties="keep",
    internal=False,
):
    attributes: dict[str, str] = {}
    comments: list[str] = []
    fields: list[Token] = []
    for part in parts:
        if part.kind == "comment":
            values, other = annotation_attributes(part.text)
            for key, value in values.items():
                _add_attribute(attributes, key, value)
            comments.extend(other)
        else:
            fields.append(part)
    kinds = [field.kind for field in fields]
    if kinds not in ([], ["name"], [":", "name"], ["name", ":", "name"]):
        raise ValueError("Malformed Newick node label or branch length.")
    validate_age_attributes(attributes)
    if internal and node_label and node_label in attributes:
        value = attributes[node_label]
        label = Token("name", "'" + value.replace("'", "''") + "'")
        if fields and fields[0].kind == "name":
            fields[0] = label
        else:
            fields.insert(0, label)
    if age_ci == "drop":
        attributes = {
            key: value for key, value in attributes.items() if key not in CI_FIELDS
        }
    for key in AGE_FIELDS & attributes.keys():
        attributes[key] = scaled_time(attributes[key], factor)
    if properties == "drop":
        attributes = {}
    rendered = []
    for index, field in enumerate(fields):
        if index and fields[index - 1].kind == ":":
            # ETE verifies unscaled lengths, including read_tree's diagnostic
            # allow_non_finite mode. Scaling must always validate its operands.
            rendered.append(
                scaled_time(field.text, factor) if factor != 1 else field.text
            )
        else:
            rendered.append(field.text)
    if parser:
        comments = [
            comment for comment in comments if comment.upper() in {"[&R]", "[&U]"}
        ]
    # ETE's NHX grammar requires annotations after the branch length.
    return (
        "".join(rendered) + _render_attributes(attributes, output) + "".join(comments)
    )


def transform_annotations(
    text,
    *,
    output="nhx",
    factor=Decimal(1),
    age_ci="keep",
    parser=False,
    node_label="",
    properties="keep",
):
    """Convert annotations without reserializing node names or support values."""
    factor = finite_decimal(factor, positive=True)
    if age_ci not in {"keep", "drop"}:
        raise ValueError("--age-ci must be keep or drop.")
    if properties not in {"keep", "drop"}:
        raise ValueError("--properties must be keep or drop.")
    result: list[str] = []
    suffix: list[Token] = []
    internal = False
    for token in tokens(text):
        if token.kind in {"(", ")", ",", ";"}:
            result.append(
                _node_suffix(
                    suffix,
                    output,
                    factor,
                    age_ci,
                    parser,
                    node_label,
                    properties,
                    internal,
                )
            )
            result.append(token.text)
            suffix = []
            internal = token.kind == ")"
        else:
            suffix.append(token)
    result.append(
        _node_suffix(
            suffix, output, factor, age_ci, parser, node_label, properties, internal
        )
    )
    return "".join(result)


def _dated_tree(statement):
    from ete4 import Tree

    from nwkit.rooting_state import extract_rooting_token

    text = transform_annotations(statement, parser=True)
    text, _ = extract_rooting_token(text)
    tree = Tree(text, parser=1)
    names = list(tree.leaf_names())
    if (
        len(names) < 2
        or any(not name for name in names)
        or len(names) != len(set(names))
    ):
        raise ValueError("Dated trees require at least two distinct, named tips.")
    distances = [node.dist for node in tree.traverse() if not node.is_root]
    if all(value is None for value in distances):
        return None
    if any(value is None or finite_decimal(value) < 0 for value in distances):
        raise ValueError(
            "Dated trees require complete non-negative finite branch lengths."
        )
    if tree.dist is not None and finite_decimal(tree.dist) < 0:
        raise ValueError("Dated-tree root length must be non-negative and finite.")
    return tree


def select_mcmctree_tree(statements):
    """Recognize PAML's topology/plain/annotated views of one dated tree."""
    import hashlib

    candidates = []
    for statement in statements:
        tree = _dated_tree(statement)
        if tree is None:
            continue
        signatures: dict[Any, str] = {}
        for node in tree.traverse("postorder"):
            identity = (
                node.name if node.is_leaf else "",
                node.dist,
                sorted(signatures[child] for child in node.children),
            )
            signatures[node] = hashlib.sha256(repr(identity).encode()).hexdigest()
        ci_count = sum("age_ci_low" in node.props for node in tree.traverse())
        candidates.append((statement, signatures[tree], ci_count))
    if len(candidates) == 1:
        return candidates[0][0]
    # A plain and a CI-annotated rendering of exactly the same branch lengths
    # are the documented paired PAML views, not two posterior tree samples.
    annotated = [entry for entry in candidates if entry[2]]
    if (
        len(candidates) == 2
        and len(annotated) == 1
        and len({entry[1] for entry in candidates}) == 1
    ):
        return annotated[0][0]
    if not candidates:
        raise ValueError("No valid dated tree was found in the MCMCtree FigTree block.")
    raise ValueError(
        "Multiple dated trees are ambiguous; select one with --tree-index."
    )
