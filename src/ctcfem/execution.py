"""Fault-tolerant parallel execution.

* One BLAS/OpenMP thread is enforced in the CLI and again inside every worker (threadpoolctl), so an
  inherited OMP_NUM_THREADS cannot oversubscribe the machine.
* If a worker dies (out-of-memory killer, native crash) the pool is rebuilt; tasks that were in flight are
  retried in isolation, one per pool, so only the task that kills its worker is recorded as failed.
* Per-task wall-time limit (SIGALRM) and per-worker memory limit (RLIMIT_DATA).
* Submission pauses when available memory falls below a threshold; never more tasks in flight than workers.
* If workers die without any task finishing, the run aborts (SystematicFailure) without marking tasks.
"""
from __future__ import annotations
import os, signal, resource
from collections import deque
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from concurrent.futures.process import BrokenProcessPool
import multiprocessing as mp

THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
_TPL = None


def force_single_thread():
    """Call BEFORE numpy is imported (CLI) and again inside each worker."""
    for v in THREAD_VARS:
        os.environ[v] = "1"


def available_memory_gb() -> float:
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 2 ** 20
    except OSError:
        pass
    return float("inf")


def process_threads() -> int:
    try:
        for line in open("/proc/self/status"):
            if line.startswith("Threads:"):
                return int(line.split()[1])
    except OSError:
        pass
    return -1


def _init_worker(memory_limit_gb, nice):
    global _TPL
    force_single_thread()
    try:
        from threadpoolctl import threadpool_limits
        _TPL = threadpool_limits(limits=1)
    except ImportError:
        pass
    if memory_limit_gb:
        lim = int(memory_limit_gb * 2 ** 30)
        _, hard = resource.getrlimit(resource.RLIMIT_DATA)
        if hard != resource.RLIM_INFINITY:
            lim = min(lim, hard)
        resource.setrlimit(resource.RLIMIT_DATA, (lim, hard))
    if nice:
        try:
            os.nice(int(nice))
        except OSError:
            pass


class TaskTimeout(Exception):
    pass


class SystematicFailure(RuntimeError):
    """Workers die without any task finishing: an environment problem, not a task problem."""


def _alarm(signum, frame):
    raise TaskTimeout("maximum time per task exceeded")


def call_with_limit(fn, arg, timeout_s):
    if timeout_s:
        signal.signal(signal.SIGALRM, _alarm)
        signal.setitimer(signal.ITIMER_REAL, float(timeout_s))
    try:
        return fn(arg)
    finally:
        if timeout_s:
            signal.setitimer(signal.ITIMER_REAL, 0)


def _wrapped(pack):
    fn, arg, timeout_s = pack
    return call_with_limit(fn, arg, timeout_s)


def run_batch(items, fn, workers, *, timeout_s=3600, memory_limit_gb=4.0, nice=0, min_free_gb=4.0,
              tasks_per_child=25, on_death=None, log=print):
    """items: iterable of (key, argument); fn must be importable ('spawn' context).
    Yields (key, result | ('died', reason) | ('error', repr))."""
    ctx = mp.get_context("spawn")
    pending, suspects = deque(items), deque()
    successes, consecutive_deaths = 0, 0

    def pool(n):
        return ProcessPoolExecutor(max_workers=n, mp_context=ctx, initializer=_init_worker,
                                   initargs=(memory_limit_gb, nice), max_tasks_per_child=tasks_per_child)

    while pending:
        ex = pool(workers); in_flight = {}
        try:
            while pending or in_flight:
                while pending and len(in_flight) < workers:
                    if available_memory_gb() < min_free_gb and in_flight:
                        break
                    k, a = pending.popleft()
                    in_flight[ex.submit(_wrapped, (fn, a, timeout_s))] = (k, a)
                done, _ = wait(list(in_flight), timeout=30, return_when=FIRST_COMPLETED)
                for f in done:
                    k, a = in_flight[f]
                    try:
                        res = f.result()
                    except BrokenProcessPool:
                        raise                       # f stays in in_flight: it is not lost
                    except Exception as e:
                        res = ("error", repr(e)[:300])
                    del in_flight[f]; successes += 1
                    yield k, res
        except BrokenProcessPool:
            log(f"a worker died; {len(in_flight)} in-flight tasks will be retried in isolation; rebuilding the pool")
            suspects.extend(in_flight.values()); ex.shutdown(wait=False, cancel_futures=True)
            continue
        ex.shutdown(wait=True)

    while suspects:
        k, a = suspects.popleft(); ex = pool(1)
        try:
            res = ex.submit(_wrapped, (fn, a, timeout_s)).result()
            successes += 1; consecutive_deaths = 0
            yield k, res
        except BrokenProcessPool:
            consecutive_deaths += 1
            if successes == 0 or consecutive_deaths >= 3:
                ex.shutdown(wait=False, cancel_futures=True)
                raise SystematicFailure(f"{consecutive_deaths} workers died in a row ({successes} successes in total): "
                                        "an environment problem (memory, libraries, spawn context), not a task problem. "
                                        "No task was marked; fix the environment and restart.")
            reason = "the worker process died (possible out-of-memory or native-library crash)"
            log(f"task {k}: {reason}; recorded and skipped")
            if on_death:
                on_death(k, a, reason)
            yield k, ("died", reason)
        except Exception as e:
            yield k, ("error", repr(e)[:300])
        finally:
            ex.shutdown(wait=False, cancel_futures=True)


def safe_workers(requested, memory_per_task_gb=1.0, reserve_gb=8.0, log=print):
    cpu = os.cpu_count() or 1
    mem = available_memory_gb()
    by_mem = int(max(1, (mem - reserve_gb) // memory_per_task_gb)) if mem != float("inf") else requested
    w = max(1, min(requested, cpu, by_mem))
    if w < requested:
        log(f"workers reduced from {requested} to {w} (logical CPUs {cpu}; available memory {mem:.0f} GB)")
    return w
