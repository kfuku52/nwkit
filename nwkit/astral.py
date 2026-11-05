"""ASTRAL extended-support annotation dialect."""

from nwkit.tree_formats import tokens


def parse_annotation(annotation):
    parsed = {}
    for field in annotation.split(";"):
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        key, value = key.strip(), value.strip()
        if key in parsed:
            raise ValueError(f"Duplicate ASTRAL support field: {key}")
        parsed[key] = value
    return parsed


def extract_support_labels(text, label_key):
    output, missing, internal = [], 0, False
    for token in tokens(text):
        raw = token.text
        content = raw[1:-1] if raw[:1] in {"'", '"'} else raw
        if (
            internal
            and content.startswith("[")
            and content.endswith("]")
            and not content.startswith("[&")
        ):
            fields = parse_annotation(content[1:-1])
            if fields:
                if label_key in fields:
                    raw = fields[label_key]
                else:
                    missing += 1
        output.append(raw)
        if token.kind == ")":
            internal = True
        elif token.kind in {"(", ",", ":", ";"}:
            internal = False
    return "".join(output), missing


def has_label_annotations(text):
    internal = False
    for token in tokens(text):
        if internal and token.kind == "name" and token.text[:1] in {"'", '"'}:
            content = token.text[1:-1]
            if content.startswith("[") and content.endswith("]"):
                fields = parse_annotation(content[1:-1])
                if any(
                    key in fields for key in ("q1", "q2", "q3", "pp1", "pp2", "pp3")
                ):
                    return True
        if token.kind == ")":
            internal = True
        elif token.kind in {"(", ",", ":", ";"}:
            internal = False
    return False
