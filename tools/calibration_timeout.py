"""Reusable spawned worker with a hard timeout on platforms without SIGALRM."""

import atexit
import multiprocessing
import threading

_pool = None
_lock = threading.RLock()


def close_worker():
    global _pool
    with _lock:
        if _pool is not None:
            _pool.terminate()
            _pool.join()
            _pool = None


def run_with_timeout(function, args, seconds):
    global _pool
    with _lock:
        if _pool is None:
            _pool = multiprocessing.get_context("spawn").Pool(1)
        pending = _pool.apply_async(function, args)
        try:
            return pending.get(timeout=seconds)
        except multiprocessing.TimeoutError as exc:
            close_worker()
            raise TimeoutError(f"Fit exceeded {seconds}s") from exc
        except BaseException:
            close_worker()
            raise


atexit.register(close_worker)
