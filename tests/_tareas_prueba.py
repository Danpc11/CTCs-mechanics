"""Tareas de prueba para tests/test_ejecucion.py (deben ser importables: contexto spawn)."""
import os, signal, time

def tarea(arg):
    tipo, x = arg
    if tipo == "normal":
        return x * 2
    if tipo == "morir":                       # simula el OOM killer: SIGKILL al propio worker
        os.kill(os.getpid(), signal.SIGKILL)
    if tipo == "colgar":
        time.sleep(30)
    if tipo == "hilos":
        from threadpoolctl import threadpool_info
        import numpy  # noqa: F401  (carga BLAS en el worker)
        return dict(env=os.environ.get("OMP_NUM_THREADS"), blas=[i["num_threads"] for i in threadpool_info()])
