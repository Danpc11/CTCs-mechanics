"""[FIX-worker] y [FIX-threads]: el ejecutor no se cuelga ni sobresuscribe."""
import os, time
import _tareas_prueba as T
from ctcfem import ejecucion as E


def test_worker_muerto_y_tarea_colgada_no_cuelgan_el_lote():
    muertos = []
    items = [(f"n{i}", ("normal", i)) for i in range(6)] + [("killer", ("morir", 0)), ("lenta", ("colgar", 0))]
    t0 = time.time()
    res = dict(E.ejecutar_lote(items, T.tarea, 2, timeout_s=2, limite_mem_gb=0, min_libre_gb=0,
                               al_morir=lambda k, a, m: muertos.append(k), log=lambda m: None))
    assert time.time() - t0 < 90
    assert all(res[f"n{i}"] == 2 * i for i in range(6))          # ninguna tarea sana se pierde
    assert res["killer"][0] == "muerto" and muertos == ["killer"]  # solo la culpable se marca
    assert res["lenta"][0] == "error" and "TiempoAgotado" in res["lenta"][1]


def test_un_hilo_por_worker_aunque_el_entorno_herede_otro_valor():
    viejo = os.environ.get("OMP_NUM_THREADS")
    os.environ["OMP_NUM_THREADS"] = "8"                          # p. ej. exportado en .bashrc
    try:
        r = dict(E.ejecutar_lote([("h", ("hilos", 0))], T.tarea, 1, limite_mem_gb=0, min_libre_gb=0, log=lambda m: None))["h"]
    finally:
        if viejo is None: os.environ.pop("OMP_NUM_THREADS")
        else: os.environ["OMP_NUM_THREADS"] = viejo
    assert r["env"] == "1" and all(n == 1 for n in r["blas"])


def test_workers_limitados_por_cpu():
    assert 1 <= E.workers_seguros(10_000, 1.0, 0.0, log=lambda m: None) <= (os.cpu_count() or 1)


def test_fallo_sistematico_aborta_sin_marcar_tareas():
    import pytest
    marcadas = []
    items = [(f"k{i}", ("morir", 0)) for i in range(4)]
    with pytest.raises(E.FalloSistematico):
        list(E.ejecutar_lote(items, T.tarea, 2, limite_mem_gb=0, min_libre_gb=0,
                             al_morir=lambda k, a, m: marcadas.append(k), log=lambda m: None))
    assert marcadas == []
