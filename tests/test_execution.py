"""The parallel executor neither hangs nor oversubscribes."""
import os, time
import pytest
import _test_tasks as T
from ctcfem import execution as E


def test_dead_worker_and_hanging_task_do_not_hang_the_batch():
    dead = []
    items = [(f"n{i}", ("normal", i)) for i in range(6)] + [("killer", ("die", 0)), ("slow", ("hang", 0))]
    t0 = time.time()
    res = dict(E.run_batch(items, T.task, 2, timeout_s=2, memory_limit_gb=0, min_free_gb=0,
                           on_death=lambda k, a, m: dead.append(k), log=lambda m: None))
    assert time.time() - t0 < 90
    assert all(res[f"n{i}"] == 2 * i for i in range(6))        # no healthy task is lost
    assert res["killer"][0] == "died" and dead == ["killer"]   # only the culprit is marked
    assert res["slow"][0] == "error" and "TaskTimeout" in res["slow"][1]


def test_one_thread_per_worker_even_with_inherited_environment():
    old = os.environ.get("OMP_NUM_THREADS"); os.environ["OMP_NUM_THREADS"] = "8"
    try:
        r = dict(E.run_batch([("h", ("threads", 0))], T.task, 1, memory_limit_gb=0, min_free_gb=0, log=lambda m: None))["h"]
    finally:
        os.environ.pop("OMP_NUM_THREADS") if old is None else os.environ.__setitem__("OMP_NUM_THREADS", old)
    assert r["env"] == "1" and all(n == 1 for n in r["blas"])


def test_systematic_failure_aborts_without_marking_tasks():
    marked = []
    with pytest.raises(E.SystematicFailure):
        list(E.run_batch([(f"k{i}", ("die", 0)) for i in range(4)], T.task, 2, memory_limit_gb=0, min_free_gb=0,
                         on_death=lambda k, a, m: marked.append(k), log=lambda m: None))
    assert marked == []


def test_workers_capped_by_cpu():
    assert 1 <= E.safe_workers(10_000, 1.0, 0.0, log=lambda m: None) <= (os.cpu_count() or 1)
