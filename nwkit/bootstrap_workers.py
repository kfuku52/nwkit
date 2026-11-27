"""Ordered, bounded process execution for independent bootstrap refits."""

from collections import deque
from concurrent.futures import Future, ProcessPoolExecutor
from itertools import islice
from multiprocessing import get_context
from typing import Any

_worker_function: Any = None


def validate_bootstrap_workers(workers):
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError("Bootstrap workers must be a positive integer.")


def _init_bootstrap_worker(function):
    global _worker_function
    _worker_function = function


def _run_draw(draw):
    return _worker_function(draw)


def _submit_draws(pool, draws):
    while True:
        try:
            draw = next(draws)
        except StopIteration:
            return
        except Exception as error:
            # Re-raise generation failures at their original draw position.
            failed: Future[Any] = Future()
            failed.set_exception(error)
            yield failed
            return
        yield pool.submit(_run_draw, draw)


def bootstrap_map(function, draws, *, workers, replicates):
    """Generate draws in the coordinator and return results in draw order.

    Send the common fitting context once per process and keep at most two
    draws per worker in flight. The single-worker path preserves local callers
    and uses no extra process.
    """
    validate_bootstrap_workers(workers)
    workers = min(workers, replicates)
    if workers <= 1:
        yield from map(function, draws)
        return
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=get_context("spawn"),
        initializer=_init_bootstrap_worker,
        initargs=(function,),
    ) as pool:
        submitted = _submit_draws(pool, iter(draws))
        pending = deque(islice(submitted, 2 * workers))
        try:
            while pending:
                future = pending.popleft()
                yield future.result()
                pending.extend(islice(submitted, 1))
        finally:
            for future in pending:
                future.cancel()
