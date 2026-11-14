#!/usr/bin/env python3
"""Choose workflow time units and scale MCMCtree calibrations/control files."""

from __future__ import annotations

import argparse
import re
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path

NUMERIC_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
CALIBRATION_RE = re.compile(r"(?P<kind>\b[BUL])\((?P<body>[^)]*)\)")


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def parse_decimal(text: str) -> Decimal | None:
    stripped = text.strip()
    if not stripped or stripped == "-":
        return None
    try:
        return Decimal(stripped)
    except InvalidOperation:
        return None


def format_decimal(value: Decimal) -> str:
    if value.is_zero():
        return "0"
    normalized = value.normalize()
    text = format(normalized, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    if text == "-0":
        return "0"
    return text


def scale_decimal(value: Decimal, scale: Decimal, direction: str) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = 28
        if direction == "down":
            return value / scale
        if direction == "up":
            return value * scale
    raise ValueError(f"invalid direction: {direction}")


def scale_number_text(text: str, scale: Decimal, direction: str) -> str:
    value = parse_decimal(text)
    if value is None:
        return text
    return format_decimal(scale_decimal(value, scale, direction))


def calibration_age_field_count(kind: str, field_count: int) -> int:
    if kind == "B":
        return min(2, field_count)
    if kind in {"L", "U"}:
        return min(1, field_count)
    return 0


def split_field_spacing(field: str) -> tuple[str, str, str]:
    leading_len = len(field) - len(field.lstrip())
    trailing_len = len(field) - len(field.rstrip())
    leading = field[:leading_len]
    trailing = field[len(field) - trailing_len :] if trailing_len else ""
    core_end = len(field) - trailing_len if trailing_len else len(field)
    return leading, field[leading_len:core_end], trailing


def calibration_ages(text: str) -> list[Decimal]:
    ages: list[Decimal] = []
    for match in CALIBRATION_RE.finditer(text):
        kind = match.group("kind")
        fields = match.group("body").split(",")
        for idx in range(calibration_age_field_count(kind, len(fields))):
            value = parse_decimal(fields[idx])
            if value is not None:
                if not value.is_finite():
                    raise ValueError("Calibration ages must be finite numbers")
                ages.append(value)
    return ages


def choose_scale_factor(text: str, target_max: Decimal) -> Decimal:
    ages = calibration_ages(text)
    if not ages:
        return Decimal(1)
    max_age = max(abs(age) for age in ages)
    factor = Decimal(1)
    while max_age / factor > target_max:
        factor *= Decimal(10)
    return factor


def scale_calibration_labels(text: str, scale: Decimal, direction: str) -> str:
    if scale == 1:
        return text

    def replace(match: re.Match[str]) -> str:
        kind = match.group("kind")
        fields = match.group("body").split(",")
        for idx in range(calibration_age_field_count(kind, len(fields))):
            leading, core, trailing = split_field_spacing(fields[idx])
            fields[idx] = leading + scale_number_text(core, scale, direction) + trailing
        return f"{kind}({','.join(fields)})"

    return CALIBRATION_RE.sub(replace, text)


def scale_rootage_line(line: str, scale: Decimal, direction: str) -> str:
    if not re.match(r"^\s*RootAge\s*=", line):
        return line
    if "#" in line:
        value_part, comment = line.split("#", 1)
        comment = "#" + comment
    else:
        value_part = line
        comment = ""
    if CALIBRATION_RE.search(value_part):
        return scale_calibration_labels(value_part, scale, direction) + comment
    return (
        re.sub(
            NUMERIC_RE,
            lambda match: scale_number_text(match.group(0), scale, direction),
            value_part,
        )
        + comment
    )


def scale_ctl_rootage_text(text: str, scale: Decimal, direction: str) -> str:
    if scale == 1:
        return text
    return "".join(
        scale_rootage_line(line, scale, direction)
        for line in text.splitlines(keepends=True)
    )


def parse_scale(value: str) -> Decimal:
    scale = parse_decimal(value)
    if scale is None or not scale.is_finite() or scale <= 0:
        raise argparse.ArgumentTypeError("scale must be a positive number")
    return scale


def cmd_factor(args: argparse.Namespace) -> int:
    target_max = parse_scale(args.target_max)
    factor = choose_scale_factor(read_text(args.infile), target_max)
    print(format_decimal(factor))
    return 0


def cmd_scale_calibrations(args: argparse.Namespace) -> int:
    scale = parse_scale(args.scale)
    text = scale_calibration_labels(read_text(args.infile), scale, args.direction)
    write_text(args.outfile, text)
    return 0


def cmd_scale_ctl_rootage(args: argparse.Namespace) -> int:
    scale = parse_scale(args.scale)
    text = scale_ctl_rootage_text(read_text(args.infile), scale, args.direction)
    write_text(args.outfile, text)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Scale MCMCTree time values between public and internal units."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    factor = subparsers.add_parser(
        "factor",
        help="Choose a power-of-ten divisor that keeps calibration ages <= target.",
    )
    factor.add_argument("--infile", type=Path, required=True)
    factor.add_argument("--target-max", default="10")
    factor.set_defaults(func=cmd_factor)

    scale_cal = subparsers.add_parser(
        "scale-calibrations",
        help="Scale B/L/U calibration labels in a Newick file.",
    )
    scale_cal.add_argument("--infile", type=Path, required=True)
    scale_cal.add_argument("--outfile", type=Path, required=True)
    scale_cal.add_argument("--scale", required=True)
    scale_cal.add_argument("--direction", choices=["down", "up"], required=True)
    scale_cal.set_defaults(func=cmd_scale_calibrations)

    scale_ctl = subparsers.add_parser(
        "scale-ctl-rootage",
        help="Scale RootAge values in a MCMCTree control file.",
    )
    scale_ctl.add_argument("--infile", type=Path, required=True)
    scale_ctl.add_argument("--outfile", type=Path, required=True)
    scale_ctl.add_argument("--scale", required=True)
    scale_ctl.add_argument("--direction", choices=["down", "up"], required=True)
    scale_ctl.set_defaults(func=cmd_scale_ctl_rootage)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
