"""dP_crit de MCF-7 a 2a/D nominal = 1.70 en las cinco geometrias de entrada:
reduccion original frente a corregida (re-ejecucion FEM directa)."""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "01_solver"), os.path.join(HERE, "..", "03_modalidades")]
import reduccion as R, fase2_lib as L
D = os.path.join(HERE, "..", "datos")
recs = [json.loads(l) for f in ["fem_mcf7_referencia_reejecucion.jsonl", "fem_mcf7_geometrias_reejecucion.jsonl"]
        for l in open(os.path.join(D, f))]
df = pd.DataFrame(recs); df = df[df.converge]
a = 8.74; Rw = 1.6 * a; c = 1.70; filas = []
for geo, g in df.groupby("geometria"):
    zt, Lt = float(g.zt.iloc[0]), float(g.Ltap.iloc[0])
    spO = R.spline_original(g.Rt.values, g.W_total.values, a)
    spC = R.spline_por_cuello(g.Rt.values, g.W_total.values, a, Rw, zt, Lt)
    rt = a / c
    pO = L.presion_critica(a, rt, Rw, zt, Lt, spO, spO.derivative())["dP_crit"]
    pC = R.presion_critica(a, rt, Rw, zt, Lt, spC)["dP_crit"]
    filas.append(dict(geometria=geo, zt=zt, Ltap=Lt, g0=round(R.g0(zt, Lt), 3),
                      conf_efectivo=round(a / R.R_cuello(rt, Rw, zt, Lt), 3),
                      angulo=round(L.angulo_de_Ltap(Lt, rt, Rw, zt), 1),
                      original=round(pO, 4), corregido=round(pC, 4),
                      factor=round(pC / pO, 2) if pO > 0 else np.inf))
T = pd.DataFrame(filas).sort_values("corregido"); print(T.to_string(index=False))
T.to_csv(os.path.join(D, "geometrias_170_corregido.csv"), index=False)
v = T[T.corregido > 0].corregido; print("rango corregido:", round(v.max() / v.min(), 1), "x")
