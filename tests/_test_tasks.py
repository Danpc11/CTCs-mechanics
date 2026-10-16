"""Tasks for tests/test_execution.py (must be importable: spawn context)."""
import os, signal, time


def task(arg):
    kind, x = arg
    if kind == "normal":
        return x * 2
    if kind == "die":                       # simulates the out-of-memory killer
        os.kill(os.getpid(), signal.SIGKILL)
    if kind == "hang":
        time.sleep(30)
    if kind == "threads":
        from threadpoolctl import threadpool_info
        import numpy  # noqa: F401
        return dict(env=os.environ.get("OMP_NUM_THREADS"), blas=[i["num_threads"] for i in threadpool_info()])
