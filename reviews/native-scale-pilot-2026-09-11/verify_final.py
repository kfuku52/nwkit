"""Recheck every retained final pair, including optimizer and nested outputs."""

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def scientific(record):
    """Remove only measured cost and the output destination."""
    result = dict(record)
    result.pop("peak_rss_bytes")
    result["configuration"] = dict(result["configuration"])
    result["configuration"].pop("output")
    for key in list(result):
        if key.endswith("_seconds"):
            result.pop(key)
    result["search_calls"] = [
        {key: value for key, value in row.items() if key != "seconds"}
        for row in result["search_calls"]
    ]
    return result


def finite(value):
    if isinstance(value, dict):
        return all(finite(item) for item in value.values())
    if isinstance(value, list):
        return all(finite(item) for item in value)
    return not isinstance(value, float) or math.isfinite(value)


def verify(first, second):
    before, after = [json.loads((ROOT / path).read_text()) for path in (first, second)]
    left, right = scientific(before), scientific(after)
    if not finite(left) or not finite(right) or left != right:
        raise ValueError(f"Scientific output differs: {first} vs {second}")
    if after["status"] != "ok":
        raise ValueError(f"Incomplete run: {second}")
    workflow = None
    if "support" in after:
        calibration = after["calibration"]
        nested = after["nested_calibrations"]

        def draws(value):
            return sum(test["replicates"] for test in value["tests"])

        expected = 1 + draws(calibration) + sum(1 + draws(value) for value in nested)
        if len(nested) != after["support"]["replicates"]:
            raise ValueError("Missing nested support calibration")
        if len(after["search_calls"]) != expected or any(
            row["status"] != "ok" for row in after["search_calls"]
        ):
            raise ValueError("Missing or incomplete workflow search")
        workflow = {"expected_search_calls": expected, "all_completed": True}
    return {"before": first, "after": second, "exact": True, "workflow": workflow}


if __name__ == "__main__":
    pairs = [
        (
            str(path.relative_to(ROOT)),
            str(path.relative_to(ROOT)).replace("-before-", "-after-"),
        )
        for path in sorted((ROOT / "final-comparisons").glob("*-before-*.json"))
    ]
    if len(pairs) != 12:
        raise ValueError("Expected twelve final repeated pairs")
    pairs += [
        (
            "final-followups/workflow512-before.json",
            "final-followups/workflow512-after.json",
        ),
        ("results/workflow.json", "final-followups/workflow32-after.json"),
        (
            "results/fit100-estimated.json",
            "final-followups/fit100-estimated-after.json",
        ),
        ("results/search100.json", "final-followups/search100-after.json"),
    ]
    print(json.dumps({"checks": [verify(*pair) for pair in pairs]}, indent=2))
