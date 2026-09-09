"""Protocol failures and loaded-state reuse, including real IQ-TREE evaluation."""

import sys

import numpy as np
import pytest

from nwkit.radte_iqtree_session import IQTreeSession

READY = """print('IQTREE_SESSION_V1 READY 3', flush=True)
print('IQTREE_SESSION_V1 EDGE 0.1 A', flush=True)
print('IQTREE_SESSION_V1 EDGE 0.2 B', flush=True)
print('IQTREE_SESSION_V1 EDGE 0.3 C', flush=True)
print('IQTREE_SESSION_V1 END', flush=True)
"""


def command(source):
    return [sys.executable, "-u", "-c", source]


def test_protocol_state_reuse_and_close():
    q = IQTreeSession(
        command(
            READY
            + """import sys
for line in sys.stdin:
    if line.strip() == 'QUIT': break
    print('IQTREE_SESSION_V1 VALUE -6 1 2 3 -1 -2 -3', flush=True)
"""
        )
    )
    pid = q.process.pid
    for scale in [1, 2, 3]:
        value, gradient, diagonal = q.evaluate(np.array([0.1, 0.2, 0.3]) * scale)
        assert value == 6
        np.testing.assert_equal(gradient, [-1, -2, -3])
        np.testing.assert_equal(diagonal, [1, 2, 3])
        assert q.process.pid == pid and q.process.poll() is None
    q.close()
    assert q.process.poll() == 0
    assert not q.reader.is_alive()
    q.close()
    with pytest.raises(ValueError, match="closed"):
        q.evaluate([0.1, 0.2, 0.3])


@pytest.mark.parametrize(
    "response",
    [
        "VALUE -6 1 2",
        "VALUE nan 1 2 3 4 5 6",
        "ERROR nonfinite likelihood or branch derivative",
    ],
)
def test_protocol_errors_close_worker(response):
    q = IQTreeSession(
        command(
            READY
            + f"""import sys
sys.stdin.readline()
print('IQTREE_SESSION_V1 {response}', flush=True)
sys.stdin.readline()
"""
        )
    )
    with pytest.raises(ValueError):
        q.evaluate([0.1, 0.2, 0.3])
    assert q.process.poll() is not None


def test_dead_worker_is_not_restarted():
    q = IQTreeSession(command(READY + "import sys\nsys.stdin.readline()\n"))
    with pytest.raises(ValueError, match="exited"):
        q.evaluate([0.1, 0.2, 0.3])
    assert q.process.poll() is not None


def test_timeout_closes_worker():
    q = IQTreeSession(
        command(
            READY
            + 'import sys\nfor line in sys.stdin:\n if line.strip() == "QUIT": break\n'
        )
    )
    q.timeout = 0.05
    with pytest.raises(ValueError, match="timed out"):
        q.evaluate([0.1, 0.2, 0.3])
    assert q.process.poll() is not None


def test_missing_feature_is_explicit():
    with pytest.raises(ValueError, match="requires an IQ-TREE build"):
        IQTreeSession(command("print('Invalid --likelihood-session')"))


def test_invalid_request_does_not_change_session():
    q = IQTreeSession(
        command(
            READY
            + 'import sys\nfor line in sys.stdin:\n if line.strip() == "QUIT": break\n'
        )
    )
    try:
        for lengths in [[0, 1, 1], [float("nan"), 1, 1], [1, 2]]:
            with pytest.raises(ValueError, match="one positive length"):
                q.evaluate(lengths)
        assert q.process.poll() is None
    finally:
        q.close()


def test_releasing_likelihood_session_reaps_worker():
    import gc
    import weakref

    q = IQTreeSession(
        command(
            READY
            + 'import sys\nfor line in sys.stdin:\n if line.strip() == "QUIT": break\n'
        )
    )
    worker = q.process
    reference = weakref.ref(q)
    del q
    gc.collect()
    assert reference() is None
    assert worker.poll() == 0


def test_unknown_protocol_version_is_rejected():
    with pytest.raises(ValueError, match="Unsupported IQ-TREE session protocol"):
        IQTreeSession(command("print('IQTREE_SESSION_V2 READY 3', flush=True)"))


def test_score_request_has_no_placeholder_diagonal():
    q = IQTreeSession(
        command(
            READY
            + """import sys
for line in sys.stdin:
    if line.strip() == 'QUIT': break
    if not line.startswith('SCORE '): raise RuntimeError('Expected SCORE')
    print('IQTREE_SESSION_V1 SCORE -6 1 2 3', flush=True)
"""
        )
    )
    try:
        value, gradient, diagonal = q.evaluate(
            [0.1, 0.2, 0.3], second_derivatives=False
        )
        assert value == 6
        np.testing.assert_array_equal(gradient, [-1, -2, -3])
        assert diagonal is None
    finally:
        q.close()


@pytest.mark.parametrize(
    "response", ["VALUE -6 1 2 3", "SCORE -6 1 2", "SCORE -6 1 2 nan"]
)
def test_score_response_failure_closes_worker(response):
    q = IQTreeSession(
        command(
            READY
            + f"""import sys
sys.stdin.readline()
print('IQTREE_SESSION_V1 {response}', flush=True)
sys.stdin.readline()
"""
        )
    )
    with pytest.raises(ValueError):
        q.evaluate([0.1, 0.2, 0.3], second_derivatives=False)
    assert q.process.poll() is not None


def test_statistics_do_not_evaluate_or_close_worker():
    q = IQTreeSession(
        command(
            READY
            + """import sys
for line in sys.stdin:
    if line.strip() == 'QUIT': break
    if line.strip() != 'STATS': raise RuntimeError('Unexpected evaluation')
    print('IQTREE_SESSION_V1 STATS 2 1 2 0 2 0.01 0.03', flush=True)
"""
        )
    )
    try:
        first = q.statistics()
        assert q.statistics() == first
        assert first["requests"] == 2
        assert first["derivative_seconds"] == 0.03
        assert q.process.poll() is None
    finally:
        q.close()


@pytest.mark.parametrize(
    "response", ["STATS 1", "STATS -1 0 0 0 0 0 0", "STATS 1 0 0 0 0 nan 0"]
)
def test_invalid_statistics_close_worker(response):
    q = IQTreeSession(
        command(
            READY
            + f"""import sys
sys.stdin.readline()
print('IQTREE_SESSION_V1 {response}', flush=True)
sys.stdin.readline()
"""
        )
    )
    with pytest.raises(ValueError):
        q.statistics()
    assert q.process.poll() is not None
