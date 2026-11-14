"""XML serialization without a Python recursion limit for deep saved trees."""

from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape


def serialize_xml(root, *, declaration=False):
    # Preserve existing formatting/namespace output on ordinary source files.
    try:
        return ET.tostring(root, encoding="unicode", xml_declaration=declaration)
    except RecursionError:
        pass
    out = ["<?xml version='1.0' encoding='utf-8'?>\n"] if declaration else []
    stack: list[tuple[ET.Element, str | None]] = [(root, None)]
    while stack:
        element, closing = stack.pop()
        if closing is not None:
            out.append(closing)
            out.append(escape(element.tail or ""))
        elif not len(element):
            out.append(ET.tostring(element, encoding="unicode"))
        else:
            # A shallow element lets ElementTree handle attribute quoting and
            # namespace declarations without recursing into its children.
            shallow = ET.Element(element.tag, element.attrib)
            serialized = ET.tostring(
                shallow, encoding="unicode", short_empty_elements=False
            )
            opening, closing = serialized.split("></", 1)
            out.extend((opening + ">", escape(element.text or "")))
            stack.append((element, "</" + closing))
            stack.extend((child, None) for child in reversed(element))
    return "".join(out)
