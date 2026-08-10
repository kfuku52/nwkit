"""Enforce NWKIT's ratcheted cyclomatic-complexity budget."""

import json
import subprocess
import sys

MAX_GRADE_F_BLOCKS = 0
MAX_BLOCK_COMPLEXITY = 40
MAX_AVERAGE_COMPLEXITY = 6.5


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
    grade_f = [block for block in blocks if block["rank"] == "F"]
    maximum = max(block["complexity"] for block in blocks)
    average = sum(block["complexity"] for block in blocks) / len(blocks)
    if len(grade_f) > MAX_GRADE_F_BLOCKS:
        raise RuntimeError(
            f"Grade-F block count {len(grade_f)} exceeds {MAX_GRADE_F_BLOCKS}."
        )
    if maximum > MAX_BLOCK_COMPLEXITY:
        raise RuntimeError(
            f"Maximum complexity {maximum} exceeds {MAX_BLOCK_COMPLEXITY}."
        )
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
