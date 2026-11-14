"""Hard timeout kills the worker and permits a clean subsequent calculation."""

import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from calibration_timeout import close_worker, run_with_timeout  # noqa: E402


def test_worker_reuse_and_timeout_recovery():
    try:
        assert run_with_timeout(sum, ([1, 2],), 15) == 3
        assert run_with_timeout(sum, ([3, 4],), 15) == 7
        with pytest.raises(TimeoutError, match="Fit exceeded"):
            run_with_timeout(time.sleep, (10,), 0.1)
        assert run_with_timeout(sum, ([5, 6],), 15) == 11
    finally:
        close_worker()
