"""Paso 5-A: diseno FEM para el sustituto (surrogate) de la energia.

Por analisis dimensional la energia almacenada se factoriza como
      W = Gs * a^2 * f(c, an/a, C, Pi, r_nuc, r_env)
con   c = a/R  (confinamiento),  Pi = G_cyt * a / Gs  (razon citoplasma/cortex).
Se muestrea f en 6 dimensiones por hipercubo latino y se ajusta un polinomio
de segundo grado en log(f); el Sobol se corre despues sobre el sustituto.
"""
import sys, os, json, time
import numpy as np
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L

# rangos tomados de parametros_compartimento.csv (todos los nucleares y de
# membrana son TRANSFERIDOS; por eso entran aqui como variables)
RANGOS = {
    "c":      (1.20, 2.43),    # envolvente de validez del solver
    "an_a":   (0.40, 0.60),    # transferido (Wang 2023/2025)
    "C":      (2.40, 37.0),    # Ks/Gs = 1+2C ; extremos de PC-3 (2025)
    "Pi":     (0.030, 1.00),   # log ; G_cyt a / Gs ; rango VERIFICADO
                               # por la Fase 1 (0.3x-10x sobre mu_cyt=0.1 Gs/a).
                               # Por debajo de 0.03 el problema de cascara
                               # rigida sobre interior casi vacio se vuelve mal
                               # puesto (ver resultados_fase2.md).
    "r_nuc":  (1.40, 4.00),    # transferido (condrocito, in situ vs aislado)
    "r_env":  (2.00, 8.00),    # transferido (Wang: 2x, maximo explorado 8x)
}
A_REF, GS_REF = 8.74, 8.881
ZT, LTAP = 4.0, 5.0
NDIS = 192


def muestra(n=NDIS, semilla=7):
    s = qmc.LatinHypercube(d=6, seed=semilla).random(n)
    X = {}
    for i, (k, (lo, hi)) in enumerate(RANGOS.items()):
        if k == "Pi":
            X[k] = np.exp(np.log(lo) + s[:, i] * (np.log(hi) - np.log(lo)))
        else:
            X[k] = lo + s[:, i] * (hi - lo)
    return X


def tarea(p):
    c, an_a, C, Pi, r_nuc, r_env = p
    a, Gs = A_REF, GS_REF
    an = an_a * a
    G_cyt = Pi * Gs / a
    Rt = a / c
    t0 = time.time()
    try:
        rec = L.resolver_constriccion(
            a, an, Gs, C, G_cyt, r_nuc, r_env, Rt, zt=ZT, Ltap=LTAP,
            cont_kw=dict(dlam_min=2e-3, max_intentos=55))
    except Exception as exc:
        rec = dict(converge=False, error=repr(exc)[:200])
    rec.update(c=c, an_a=an_a, C=C, Pi=Pi, r_nuc=r_nuc, r_env=r_env,
               a=a, Gs=Gs, G_cyt=G_cyt, segundos=round(time.time() - t0, 2))
    if rec.get("converge"):
        rec["f"] = rec["W_total"] / (Gs * a**2)
        for cc in ("cortex", "citoplasma", "nucleo", "envoltura"):
            rec["f_" + cc] = rec["W_" + cc] / (Gs * a**2)
    return rec


if __name__ == "__main__":
    shard, nsh = 0, 1
    if "--shard" in sys.argv:
        shard = int(sys.argv[sys.argv.index("--shard") + 1])
        nsh = int(sys.argv[sys.argv.index("--de") + 1])
    X = muestra()
    pts = list(zip(*[X[k] for k in RANGOS]))
    mios = pts[shard::nsh]
    t0, res = time.time(), []
    for i, p in enumerate(mios):
        res.append(tarea(p))
        if (i + 1) % 7 == 0:
            print("  shard %d: %3d/%d (%.0f s)"
                  % (shard, i + 1, len(mios), time.time() - t0), flush=True)
    json.dump({"RANGOS": RANGOS, "A_REF": A_REF, "GS_REF": GS_REF,
               "ZT": ZT, "LTAP": LTAP, "registros": res},
              open(os.path.join(HERE, "paso5b_diseno_s%d.json" % shard), "w"),
              indent=1)
    print("shard %d listo: %d/%d convergen, %.0f s"
          % (shard, sum(1 for r in res if r.get("converge")), len(res),
             time.time() - t0), flush=True)
