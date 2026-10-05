"""Run RDKit work in worker processes with a wall-clock timeout, memory cap and bounded queue.

RDKit calls are synchronous C++; running them on the event loop would let one
expensive query (for example a pathological SMARTS) stall every other request.
"""
import asyncio
import logging
import time
import concurrent.futures as cf
import multiprocessing
import os
import resource

TIMEOUT_S = float(os.getenv("MOLTALK_TIMEOUT_S", "20"))
WORKERS = int(os.getenv("MOLTALK_WORKERS", "2"))
MAX_QUEUED = int(os.getenv("MOLTALK_MAX_QUEUED", "16"))
WORKER_MEMORY_MB = int(os.getenv("MOLTALK_WORKER_MEMORY_MB", "2048"))


def _init_worker():
    # Workers must never write to stdout: in stdio mode it carries the MCP JSON-RPC stream.
    os.dup2(2, 1)
    if WORKER_MEMORY_MB > 0:
        limit = WORKER_MEMORY_MB * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))


class Runner:
    def __init__(self):
        self._pool = None
        self._sem = None
        self._sem_loop = None
        self._waiting = 0

    def _executor(self):
        if self._pool is None:
            os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
            self._pool = cf.ProcessPoolExecutor(WORKERS, mp_context=multiprocessing.get_context("spawn"),
                                                initializer=_init_worker)
        return self._pool

    def reset(self):
        pool, self._pool = self._pool, None
        if pool is not None:
            for process in list((pool._processes or {}).values()):
                process.kill()
            pool.shutdown(wait=False, cancel_futures=True)

    async def run(self, fn, *args):
        loop = asyncio.get_running_loop()
        if self._sem_loop is not loop:
            self._sem, self._sem_loop = asyncio.Semaphore(WORKERS), loop
        if self._waiting >= MAX_QUEUED:
            raise ValueError("Server is busy; too many chemistry requests are queued. Try again shortly.")
        self._waiting += 1
        try:
            async with self._sem:
                started = time.monotonic()
                future = loop.run_in_executor(self._executor(), fn, *args)
                future.add_done_callback(lambda _: log.warning("moltalk %s took %.2f s", fn.__name__, time.monotonic() - started))
                try:
                    return await asyncio.wait_for(future, TIMEOUT_S)
                except (asyncio.TimeoutError, cf.process.BrokenProcessPool) as exc:
                    self.reset()
                    if isinstance(exc, asyncio.TimeoutError):
                        raise ValueError(f"Computation exceeded {TIMEOUT_S:g} s and was stopped; simplify the query.") from None
                    raise ValueError("Computation worker failed (possibly the memory limit); simplify the query.") from None
        finally:
            self._waiting -= 1


log = logging.getLogger("moltalk")
runner = Runner()
