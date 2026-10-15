"""Linea de comandos:  python -m ctcfem.cli <orden> --config configs/x.yaml

ordenes:
  muestrear   escribe muestras.csv y procedencia_parametros.csv
  planificar  cuenta tareas, celulas que requieren FEM y estimacion de coste
  correr      ejecuta (reanuda si se interrumpio)  [--workers N] [--limite N]
  analizar    agrega crudos -> tablas, ventanas, informe 8/10/12 um y figuras
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(v, "1")        # 1 hilo por proceso: 40 procesos = 40 nucleos, sin sobresuscripcion
import argparse
from . import corrida as C
from . import analisis as A


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ctcfem")
    ap.add_argument("orden", choices=["muestrear", "planificar", "correr", "analizar"])
    ap.add_argument("--config", required=True)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--limite", type=int, default=None)
    a = ap.parse_args(argv)
    cfg = C.cargar_config(a.config)
    if a.orden == "muestrear":
        df = C.muestrear(cfg); print(df.groupby("poblacion").size().to_string())
    elif a.orden == "planificar":
        mu = C.muestrear(cfg); ts = C.tareas(cfg, mu)
        Dmin = min(C.aperturas(cfg))
        fem = [t for t in ts if t["celula"]["D"] > Dmin]
        print(f"tareas: {len(ts)}  con FEM (D_celula > {Dmin} um): {len(fem)}  sin contacto: {len(ts)-len(fem)}")
        w = a.workers or cfg["ejecucion"]["workers"]
        for s in (30, 120, 400):
            print(f"  si cada tarea FEM tarda {s:4d} s -> {len(fem)*s/w/3600:6.1f} h con {w} procesos")
    elif a.orden == "correr":
        n, err = C.correr(cfg, workers=a.workers, limite=a.limite); print(f"hechas {n}, errores {err}")
    elif a.orden == "analizar":
        T, R, V = A.analizar(cfg, C.RAIZ); print(f"filas {len(T)}; ver {cfg['salida']}/analisis/")


if __name__ == "__main__":
    main()
