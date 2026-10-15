"""Fase 6: diseno FEM ampliado para el sustituto.

Igual que paso5_diseno.py pero con an_a extendido a [0.40, 0.90] para cubrir
la razon nucleo/celula del linfocito (muy por encima de 0.60) y asi permitir
un Sobol con muestreo INDEPENDIENTE por poblacion sin extrapolar el sustituto.
"""
import sys, os, json, time
import numpy as np
from scipy.stats import qmc
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L
from paso5_diseno import tarea, A_REF, GS_REF, ZT, LTAP

RANGOS6 = {
    "c":      (1.20, 2.43),
    "an_a":   (0.40, 0.90),    # AMPLIADO (linfocito)
    "C":      (2.40, 37.0),
    "Pi":     (0.050, 1.00),   # log; se evita la esquina Pi<0.05 con an_a alto,
                               # donde el problema se vuelve mal puesto y cada
                               # punto quema las 55 continuaciones sin converger
    "r_nuc":  (1.40, 4.00),
    "r_env":  (2.00, 8.00),
}
NDIS = 160

def muestra(n=NDIS, semilla=11):
    s = qmc.LatinHypercube(d=6, seed=semilla).random(n)
    X = {}
    for i, (k, (lo, hi)) in enumerate(RANGOS6.items()):
        if k == "Pi":
            X[k] = np.exp(np.log(lo) + s[:, i] * (np.log(hi) - np.log(lo)))
        else:
            X[k] = lo + s[:, i] * (hi - lo)
    return X

def tarea_rapida(p):
    """tarea() de paso5_diseno pero con menos continuaciones: un punto que no
    converge en 20 intentos no va a converger, y quemaba 55."""
    c, an_a, C, Pi, r_nuc, r_env = p
    a, Gs = A_REF, GS_REF
    t0 = time.time()
    try:
        rec = L.resolver_constriccion(
            a, an_a * a, Gs, C, Pi * Gs / a, r_nuc, r_env, a / c,
            zt=ZT, Ltap=LTAP, cont_kw=dict(dlam_min=5e-3, max_intentos=20))
    except Exception as exc:
        rec = dict(converge=False, error=repr(exc)[:200])
    rec.update(c=c, an_a=an_a, C=C, Pi=Pi, r_nuc=r_nuc, r_env=r_env,
               a=a, Gs=Gs, G_cyt=Pi * Gs / a, segundos=round(time.time() - t0, 2))
    if rec.get("converge"):
        rec["f"] = rec["W_total"] / (Gs * a**2)
        for cc in ("cortex", "citoplasma", "nucleo", "envoltura"):
            rec["f_" + cc] = rec["W_" + cc] / (Gs * a**2)
    return rec


if __name__ == "__main__":
    X = muestra()
    pts = list(zip(*[X[k] for k in RANGOS6]))
    # anade tambien el diseno ORIGINAL (semilla 7, an_a<=0.60) para reusar
    # cobertura en la zona ya validada
    from paso5_diseno import muestra as m0, RANGOS as R0
    X0 = m0()
    pts += list(zip(*[X0[k] for k in R0]))
    t0 = time.time()
    out = os.path.join(HERE, "diseno_f6.jsonl")
    nok = 0
    with open(out, "w") as fh, Pool(8) as po:
        for i, rec in enumerate(po.imap_unordered(tarea_rapida, pts, chunksize=1)):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            nok += int(bool(rec.get("converge")))
            if (i + 1) % 20 == 0:
                print("  %3d/%d  (%d convergen, %.0f s)"
                      % (i + 1, len(pts), nok, time.time() - t0), flush=True)
    res = [json.loads(l) for l in open(out)]
    json.dump({"RANGOS": RANGOS6, "A_REF": A_REF, "GS_REF": GS_REF,
               "ZT": ZT, "LTAP": LTAP, "registros": res},
              open(os.path.join(HERE, "diseno_f6.json"), "w"), indent=1)
    print("listo: %d/%d convergen, %.0f s" % (nok, len(res), time.time() - t0))
