"""Enforce NWKIT's ratcheted cyclomatic-complexity budget."""

import json
import subprocess
import sys

MAX_BLOCK_COMPLEXITY = 40
MAX_AVERAGE_COMPLEXITY = 6.821

# These renderer entry points predate the complexity check. Keep their exact
# current ceilings as explicit, non-increasable debt while enforcing the normal
# budget for every other block. Removing or splitting one requires deleting its
# stale entry here, so the ratchet can only move downward.
LEGACY_BLOCK_LIMITS = {
    ("nwkit/draw.py", "_draw_tree"): 363,
    ("nwkit/draw.py", "draw_main"): 50,
}


def main() -> int:
    completed = subprocess.run(
        [sys.executable, "-m", "radon", "cc", "nwkit", "--json"],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    results = json.loads(completed.stdout)
    blocks = [block for file_blocks in results.values() for block in file_blocks]
    if not blocks:
        raise RuntimeError("Radon did not find any code blocks to analyze.")
    observed_legacy = set()
    violations = []
    for path, file_blocks in results.items():
        for block in file_blocks:
            key = (path, block["name"])
            legacy_limit = LEGACY_BLOCK_LIMITS.get(key)
            if legacy_limit is not None:
                observed_legacy.add(key)
                if block["complexity"] > legacy_limit:
                    violations.append(
                        "{}:{} complexity {} exceeds legacy ceiling {}".format(
                            path,
                            block["name"],
                            block["complexity"],
                            legacy_limit,
                        )
                    )
            elif block["complexity"] > MAX_BLOCK_COMPLEXITY:
                violations.append(
                    "{}:{} complexity {} exceeds {}".format(
                        path,
                        block["name"],
                        block["complexity"],
                        MAX_BLOCK_COMPLEXITY,
                    )
                )
    stale_limits = set(LEGACY_BLOCK_LIMITS) - observed_legacy
    if stale_limits:
        stale = ", ".join("{}:{}".format(*key) for key in sorted(stale_limits))
        violations.append("remove stale legacy complexity limits: " + stale)
    if violations:
        raise RuntimeError(
            "Cyclomatic-complexity budget failed:\n- " + "\n- ".join(violations)
        )

    grade_f = [block for block in blocks if block["rank"] == "F"]
    maximum = max(block["complexity"] for block in blocks)
    average = sum(block["complexity"] for block in blocks) / len(blocks)
    if average > MAX_AVERAGE_COMPLEXITY:
        raise RuntimeError(
            f"Average complexity {average:.2f} exceeds {MAX_AVERAGE_COMPLEXITY:.2f}."
        )
    print(
        "Radon baseline: {:.2f} average, {} grade-F blocks, {} maximum.".format(
            average,
            len(grade_f),
            maximum,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
