"""Protocol failures and loaded-state reuse, including real IQ-TREE evaluation."""

import sys

import numpy as np
import pytest

from nwkit.iqtree_worker import IQTreeWorker

READY = """print('NWKIT_IQTREE_V1 READY 3', flush=True)
print('NWKIT_IQTREE_V1 EDGE 0.1 A', flush=True)
print('NWKIT_IQTREE_V1 EDGE 0.2 B', flush=True)
print('NWKIT_IQTREE_V1 EDGE 0.3 C', flush=True)
print('NWKIT_IQTREE_V1 END', flush=True)
"""


def command(source):
    return [sys.executable, "-u", "-c", source]


def test_protocol_state_reuse_and_close():
    q = IQTreeWorker(
        command(
            READY
            + """import sys
for line in sys.stdin:
    if line.strip() == 'QUIT': break
    print('NWKIT_IQTREE_V1 VALUE -6 1 2 3 -1 -2 -3', flush=True)
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
    q = IQTreeWorker(
        command(
            READY
            + f"""import sys
sys.stdin.readline()
print('NWKIT_IQTREE_V1 {response}', flush=True)
sys.stdin.readline()
"""
        )
    )
    with pytest.raises(ValueError):
        q.evaluate([0.1, 0.2, 0.3])
    assert q.process.poll() is not None


def test_dead_worker_is_not_restarted():
    q = IQTreeWorker(command(READY + "import sys\nsys.stdin.readline()\n"))
    with pytest.raises(ValueError, match="exited"):
        q.evaluate([0.1, 0.2, 0.3])
    assert q.process.poll() is not None


def test_timeout_closes_worker():
    q = IQTreeWorker(
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
    with pytest.raises(ValueError, match="exited before completing"):
        IQTreeWorker(command("print('Not a library worker')"))


def test_invalid_request_does_not_change_session():
    q = IQTreeWorker(
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

    q = IQTreeWorker(
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
    with pytest.raises(ValueError, match="Unsupported IQ-TREE library worker protocol"):
        IQTreeWorker(command("print('NWKIT_IQTREE_V2 READY 3', flush=True)"))
