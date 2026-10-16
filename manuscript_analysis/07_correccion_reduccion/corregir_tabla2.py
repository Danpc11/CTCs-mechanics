"""
Pruebas de la reduccion corregida y Tabla 2 / beta_c corregidas.

Datos: re-ejecucion FEM de MCF-7, geometria de referencia (zt = 4, Ltap = 5 um),
con el MISMO codigo de produccion (paso1_fem.tarea), en
datos/fem_mcf7_referencia_reejecucion.jsonl  (reejecutar_fem_mcf7.py).

Propagacion a las demas lineas: la colapsabilidad W = Gs a^2 f(c) esta verificada
en el manuscrito (dispersion 1-6 %), de modo que el factor de correccion es comun a
todas las lineas a igual confinamiento; se aplica el cociente publicado
linea/MCF-7 a cada confinamiento.  Es una aproximacion DECLARADA (+-6 %); la cifra
definitiva exige re-ejecutar paso123_analisis.py con reduccion.py.
"""
import os, sys, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "01_solver"), os.path.join(HERE, "..", "03_modalidades")]
import reduccion as R
import fase2_lib as L

DATOS = os.path.join(HERE, "..", "datos")
recs = [json.loads(l) for l in open(os.path.join(DATOS, "fem_mcf7_referencia_reejecucion.jsonl"))]
recs = [r for r in recs if r["converge"]]
a, zt, Lt = 8.74, 4.0, 5.0
Rwide = 1.6 * a
Rt = np.array([r["Rt"] for r in recs]); W = np.array([r["W_total"] for r in recs])
spO = R.spline_original(Rt, W, a)
spC = R.spline_por_cuello(Rt, W, a, Rwide, zt, Lt)

# ---------------- T1: exactitud en el cuello -----------------------------------
rmin = R.R_cuello(Rt, Rwide, zt, Lt)
errC = np.max(np.abs(spC(rmin) - W) / np.maximum(W, 1e-12)) if (W > 0).any() else 0
m = W > 0
errO = np.max(np.abs(spO(rmin[m]) - W[m]) / W[m])
print(f"T1 W en el cuello del propio estado FEM: corregido err={errC:.1e}  original err={errO:.2f}")
assert errC < 1e-9 and errO > 0.5

# ---------------- T2: cero artificial a 2a/D = 1.20 ----------------------------
Wq = lambda x: np.where(x < a, (a - x) ** 2, 0.0)
dWq = lambda x: np.where(x < a, -2 * (a - x), 0.0)
p_orig = L.presion_critica(a, a / 1.20, Rwide, zt, Lt, Wq, dWq, legacy=True)["dP_crit"]
p_corr = R.presion_critica(a, a / 1.20, Rwide, zt, Lt, R.spline_por_cuello([a / 1.2, a / 1.4], [Wq(R.R_cuello(a / 1.2, Rwide, zt, Lt)), Wq(R.R_cuello(a / 1.4, Rwide, zt, Lt))], a, Rwide, zt, Lt))["dP_crit"]
print(f"T2 2a/D=1.20 con W cualquiera: original dP={p_orig}  corregido dP={p_corr:.3g}")
assert p_orig == 0.0 and p_corr > 0

# ---------------- T3: el pipeline original reproduce la Tabla 2 publicada -----
PUB = {  # Tabla 2 del manuscrito, geometria de referencia (kPa)
    1.20: dict(MCF7=0.000, SKBR3=0.000, MDA=0.000, leu=0.000),
    1.42: dict(MCF7=0.007, SKBR3=0.008, MDA=0.004, leu=0.002),
    1.70: dict(MCF7=0.160, SKBR3=0.188, MDA=0.091, leu=0.037),
    1.90: dict(MCF7=0.363, SKBR3=0.423, MDA=0.207, leu=0.083),
    2.10: dict(MCF7=0.712, SKBR3=0.823, MDA=0.406, leu=0.163),
    2.43: dict(MCF7=1.368, SKBR3=1.594, MDA=0.784, leu=0.317)}
filas = []
for c in PUB:
    rt = a / c
    pO = L.presion_critica(a, rt, Rwide, zt, Lt, spO, spO.derivative(), legacy=True)["dP_crit"]
    pC = R.presion_critica(a, rt, Rwide, zt, Lt, spC, "quiebre")
    pC2 = R.presion_critica(a, rt, Rwide, zt, Lt, spC, "quiebre", nz=20001)["dP_crit"]
    pS = R.presion_critica(a, rt, Rwide, zt, Lt, spC, "suave")
    if c >= 1.42:
        assert abs(pO - PUB[c]["MCF7"]) / PUB[c]["MCF7"] < 0.15, (c, pO)
    assert abs(pC2 - pC["dP_crit"]) / pC["dP_crit"] < 0.01          # T4: rejilla
    filas.append(dict(conf_nominal=c, conf_efectivo=round(a / R.R_cuello(rt, Rwide, zt, Lt), 3),
                      MCF7_publicado=PUB[c]["MCF7"], MCF7_original_reproducido=round(pO, 4),
                      MCF7_corregido=round(pC["dP_crit"], 4),
                      MCF7_corregido_cuello_suave=round(pS["dP_crit"], 4),
                      factor_corregido=(pC["dP_crit"] / pO) if pO > 0 else np.inf,
                      factor_suave=(pS["dP_crit"] / pO) if pO > 0 else np.inf))
print("T3 Tabla 2 publicada reproducida con el pipeline original;  T4 rejilla convergida (<1 %)")

T = pd.DataFrame(filas)
# propagacion a las demas lineas con el cociente publicado linea/MCF-7 (1.42 para 1.20)
for col in ["SKBR3", "MDA", "leu"]:
    for var, base in [("corregido", "MCF7_corregido"), ("cuello_suave", "MCF7_corregido_cuello_suave")]:
        T[f"{col}_{var}"] = [round(r[base] * (PUB[c if c > 1.2 else 1.42][col] / PUB[c if c > 1.2 else 1.42]["MCF7"]), 4)
                             for c, (_, r) in zip(PUB, T.iterrows())]
T.to_csv(os.path.join(DATOS, "presion_critica_corregida.csv"), index=False)
pd.set_option("display.width", 220)
print(T.round(4).to_string(index=False))

# ---------------- consecuencias para beta_c y las predicciones -----------------
k = float(T.loc[T.conf_nominal == 1.70, "factor_corregido"].iloc[0])
ks = float(T.loc[T.conf_nominal == 1.70, "factor_suave"].iloc[0])
u = pd.read_csv(os.path.join(DATOS, "umbral_dispositivo.csv"))
u["beta_c_corregido_kPa_a_1mm_s"] = u.beta_c_kPa_a_1mm_s * k
u["beta_c_cuello_suave_kPa_a_1mm_s"] = u.beta_c_kPa_a_1mm_s * ks
u.to_csv(os.path.join(DATOS, "umbral_dispositivo_corregido.csv"), index=False)
print(f"\nfactor a 2a/D = 1.70: {k:.2f} (quiebre)  {ks:.2f} (cuello suave)")
print(u[["poblacion", "tau_s", "beta_c_kPa_a_1mm_s", "beta_c_corregido_kPa_a_1mm_s",
         "beta_c_cuello_suave_kPa_a_1mm_s"]].round(3).to_string(index=False))
r = T.set_index("conf_nominal")
print(f"\nP3  cociente 1.70/1.42: publicado {0.160/0.007:.0f}x  corregido "
      f"{r.loc[1.70,'MCF7_corregido']/r.loc[1.42,'MCF7_corregido']:.1f}x  suave "
      f"{r.loc[1.70,'MCF7_corregido_cuello_suave']/r.loc[1.42,'MCF7_corregido_cuello_suave']:.1f}x")
print(f"P5  dP/E_app a 1.70: publicado 0.059-0.068  corregido {0.0591*k:.3f}-{0.0678*k:.3f}  "
      f"suave {0.0591*ks:.3f}-{0.0678*ks:.3f} kPa/kPa")
print(f"U* (arista inferior a beta = 10 kPa por mm/s): publicado 3.9 um/s  corregido ~{3.9*k:.0f}  suave ~{3.9*ks:.0f} um/s")
