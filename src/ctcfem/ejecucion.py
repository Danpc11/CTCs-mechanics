"""Ejecucion paralela que no se desborda.

Traslada a Python las soluciones del repositorio T2D_landscape
(R_scripts_sub_networks/03_bootstrap_nulls_reference.R):

[FIX-threads] Un hilo por worker FORZADO dentro de cada worker (threadpoolctl, el equivalente de
              RhpcBLASctl), no solo con variables de entorno del proceso padre: un OMP_NUM_THREADS
              heredado (p. ej. de .bashrc) hacia que N workers x M hilos sobresuscribieran la maquina.
[FIX-worker]  Si un worker muere (OOM killer, segfault de una libreria nativa), un Pool clasico se
              cuelga para siempre. Aqui el pool se recrea; las tareas que estaban en vuelo se
              reintentan AISLADAS (una por pool) para identificar a la culpable, que se registra como
              error y se salta, sin perder las demas.

Ademas:
* limite de tiempo por tarea (SIGALRM) y de memoria por worker (RLIMIT_DATA): una tarea desbocada
  falla con un error registrado en vez de tumbar la maquina;
* se pausa el envio de tareas si la memoria disponible baja de un umbral;
* nunca hay mas tareas en vuelo que workers (el padre no acumula colas enormes).
"""
from __future__ import annotations
import os, signal, time, resource
from collections import deque
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from concurrent.futures.process import BrokenProcessPool
import multiprocessing as mp

HILOS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
         "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
_TPL = None


def forzar_un_hilo():
    """Se llama ANTES de importar numpy (cli) y otra vez dentro de cada worker."""
    for v in HILOS:
        os.environ[v] = "1"


def memoria_disponible_gb() -> float:
    try:
        for l in open("/proc/meminfo"):
            if l.startswith("MemAvailable:"):
                return int(l.split()[1]) / 2 ** 20
    except OSError:
        pass
    return float("inf")


def hilos_del_proceso() -> int:
    try:
        for l in open("/proc/self/status"):
            if l.startswith("Threads:"):
                return int(l.split()[1])
    except OSError:
        pass
    return -1


def _init_worker(limite_mem_gb, nice):
    global _TPL
    forzar_un_hilo()
    try:
        from threadpoolctl import threadpool_limits
        _TPL = threadpool_limits(limits=1)            # se mantiene la referencia: el limite persiste
    except ImportError:
        pass
    if limite_mem_gb:
        lim = int(limite_mem_gb * 2 ** 30)
        _, duro = resource.getrlimit(resource.RLIMIT_DATA)
        if duro != resource.RLIM_INFINITY:
            lim = min(lim, duro)
        resource.setrlimit(resource.RLIMIT_DATA, (lim, duro))
    if nice:
        try:
            os.nice(int(nice))
        except OSError:
            pass


class TiempoAgotado(Exception):
    pass


class FalloSistematico(RuntimeError):
    """Los workers mueren sin que ninguna tarea termine: problema de entorno, no de una tarea."""


def _alarma(signum, frame):
    raise TiempoAgotado("tiempo maximo por tarea agotado")


def llamar_con_limite(fn, arg, timeout_s):
    """Ejecuta fn(arg) en el worker con un limite de tiempo de pared."""
    if timeout_s:
        signal.signal(signal.SIGALRM, _alarma)
        signal.setitimer(signal.ITIMER_REAL, float(timeout_s))
    try:
        return fn(arg)
    finally:
        if timeout_s:
            signal.setitimer(signal.ITIMER_REAL, 0)


def _envolver(paquete):
    fn, arg, timeout_s = paquete
    return llamar_con_limite(fn, arg, timeout_s)


def ejecutar_lote(items, fn, workers, *, timeout_s=3600, limite_mem_gb=4.0, nice=0,
                  min_libre_gb=4.0, tareas_por_proceso=25, al_morir=None, log=print):
    """items: iterable de (clave, argumento). fn debe ser importable (contexto 'spawn').

    Devuelve un generador de (clave, resultado | ('muerto', motivo) | ('error', repr)).
    al_morir(clave, argumento, motivo): se llama para la tarea que mato a su worker.
    """
    ctx = mp.get_context("spawn")
    pendientes = deque(items)
    sospechosas = deque()
    exitos = 0
    muertes_seguidas = 0

    def pool(n):
        return ProcessPoolExecutor(max_workers=n, mp_context=ctx, initializer=_init_worker,
                                   initargs=(limite_mem_gb, nice), max_tasks_per_child=tareas_por_proceso)

    # ---------------- fase principal: N workers
    while pendientes:
        ex = pool(workers)
        en_vuelo = {}
        try:
            while pendientes or en_vuelo:
                while pendientes and len(en_vuelo) < workers:
                    if memoria_disponible_gb() < min_libre_gb and en_vuelo:
                        break                                   # pausa: no se envia mas hasta liberar memoria
                    k, a = pendientes.popleft()
                    en_vuelo[ex.submit(_envolver, (fn, a, timeout_s))] = (k, a)
                hechos, _ = wait(list(en_vuelo), timeout=30, return_when=FIRST_COMPLETED)
                for f in hechos:
                    k, a = en_vuelo[f]
                    try:
                        res = f.result()
                    except BrokenProcessPool:
                        raise                                   # f sigue en en_vuelo: no se pierde
                    except Exception as e:                      # error de la tarea: se registra, se sigue
                        res = ("error", repr(e)[:300])
                    del en_vuelo[f]
                    exitos += 1
                    yield k, res
        except BrokenProcessPool:
            log(f"[FIX-worker] un worker murio; {len(en_vuelo)} tareas en vuelo pasan a revision aislada; se recrea el pool")
            sospechosas.extend(en_vuelo.values())
            ex.shutdown(wait=False, cancel_futures=True)
            continue
        ex.shutdown(wait=True)

    # ---------------- revision aislada: una tarea por pool, para identificar a la culpable
    while sospechosas:
        k, a = sospechosas.popleft()
        ex = pool(1)
        try:
            res = ex.submit(_envolver, (fn, a, timeout_s)).result()
            exitos += 1; muertes_seguidas = 0
            yield k, res
        except BrokenProcessPool:
            muertes_seguidas += 1
            if exitos == 0 or muertes_seguidas >= 3:
                ex.shutdown(wait=False, cancel_futures=True)
                raise FalloSistematico(
                    f"{muertes_seguidas} workers murieron seguidos sin terminar ninguna tarea ({exitos} exitos en total): "
                    "es un problema del entorno (memoria, librerias, contexto spawn), no de una tarea. "
                    "No se marco ninguna tarea; corrige el entorno y relanza.")
            motivo = "el proceso murio (posible falta de memoria o fallo de libreria nativa)"
            log(f"[FIX-worker] tarea {k}: {motivo}; se registra y se salta")
            if al_morir:
                al_morir(k, a, motivo)
            yield k, ("muerto", motivo)
        except Exception as e:
            yield k, ("error", repr(e)[:300])
        finally:
            ex.shutdown(wait=False, cancel_futures=True)


def workers_seguros(pedidos, mem_por_tarea_gb=1.0, reserva_gb=8.0, log=print):
    """Limita los workers por CPU y por memoria disponible."""
    cpu = os.cpu_count() or 1
    mem = memoria_disponible_gb()
    por_mem = int(max(1, (mem - reserva_gb) // mem_por_tarea_gb)) if mem != float("inf") else pedidos
    w = max(1, min(pedidos, cpu, por_mem))
    if w < pedidos:
        log(f"workers reducidos de {pedidos} a {w} (CPU logicas {cpu}; memoria disponible {mem:.0f} GB)")
    return w
