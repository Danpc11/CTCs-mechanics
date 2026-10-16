"""Tabla 2 corregida a partir de re-ejecuciones FEM DIRECTAS por linea
(geometria de referencia), sin escalado entre lineas."""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "01_solver"), os.path.join(HERE, "..", "03_modalidades")]
import reduccion as R
import fase2_lib as L
D = os.path.join(HERE, "..", "datos")
recs = []
for f in ["fem_mcf7_referencia_reejecucion.jsonl", "fem_lineas_referencia_reejecucion.jsonl"]:
    p = os.path.join(D, f)
    if os.path.exists(p):
        recs += [json.loads(l) for l in open(p)]
df = pd.DataFrame(recs)
df = df[df.converge]
zt, Lt = 4.0, 5.0
filas = []
for nm, g in df.groupby("estirpe"):
    a = float(g.a.iloc[0]); Rw = 1.6 * a
    spO = R.spline_original(g.Rt.values, g.W_total.values, a)
    spC = R.spline_por_cuello(g.Rt.values, g.W_total.values, a, Rw, zt, Lt)
    for c in [1.20, 1.42, 1.70, 1.90, 2.10, 2.43]:
        rt = a / c
        pO = L.presion_critica(a, rt, Rw, zt, Lt, spO, spO.derivative())["dP_crit"]
        pC = R.presion_critica(a, rt, Rw, zt, Lt, spC, "quiebre")["dP_crit"]
        pS = R.presion_critica(a, rt, Rw, zt, Lt, spC, "suave")["dP_crit"]
        filas.append(dict(estirpe=nm, conf_nominal=c, original=pO, corregido=pC, cuello_suave=pS))
T = pd.DataFrame(filas)
T.to_csv(os.path.join(D, "tabla2_corregida_directa.csv"), index=False)
for col in ["original", "corregido", "cuello_suave"]:
    print(f"\n== {col} (kPa) ==")
    print(T.pivot(index="conf_nominal", columns="estirpe", values=col).round(4).to_string())
p = T.pivot(index="conf_nominal", columns="estirpe", values="corregido")
s = T.pivot(index="conf_nominal", columns="estirpe", values="cuello_suave")
print("\nselectividad tumoral/leucocito a 1.70 (corregido / suave / publicado):")
for nm in [c for c in p.columns if c != "leucocito"]:
    print(f"  {nm:11s} {p.loc[1.70,nm]/p.loc[1.70,'leucocito']:.2f} / {s.loc[1.70,nm]/s.loc[1.70,'leucocito']:.2f}")
E = {"MCF-7": 2.650, "MDA-MB-231": 1.527, "SK-BR-3": 2.772, "leucocito": 0.620}
print("\nP5 dP_crit/E_app a 1.70 (corregido):", {k: round(p.loc[1.70, k] / v, 3) for k, v in E.items() if k in p})
