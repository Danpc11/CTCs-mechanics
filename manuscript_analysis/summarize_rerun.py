"""Key manuscript numbers from a re-run (run inside the working copy's manuscript_analysis/)."""
import os, json, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "datos")
pd.set_option("display.width", 200)
print("kappa/mu =", os.environ.get("CTC_KAPPA_MU", "1e3"))
t = pd.read_csv(os.path.join(D, "tabla2_corregida_directa.csv"))
print("\nTable 2 (kPa):\n", t.pivot_table(index="conf_nominal", columns="estirpe", values="corregido").round(4).to_string())
g = os.path.join(D, "geometrias_170_corregido.csv")
if os.path.exists(g):
    print("\nEntry geometries at nominal 1.70:\n", pd.read_csv(g).round(4).to_string(index=False))
u = pd.read_csv(os.path.join(D, "umbral_dispositivo.csv"))
print("\nThresholds beta_c (kPa per mm/s):\n", u[["poblacion", "E_app_kPa", "dP_el_kPa", "tau_s", "beta_c_kPa_a_1mm_s"]].round(4).to_string(index=False))
b = dict(zip(u.poblacion, u.beta_c_kPa_a_1mm_s))
if "PBMC" in b:
    for k in ("MDA-MB-231", "MCF-7"):
        if k in b: print(f"  window {k}/PBMC = {b[k] / b['PBMC']:.0f}x")
c = pd.read_csv(os.path.join(D, "barrido_constriccion.csv"))
sig = c[c.dP_critica_kPa > 1e-3]
print(f"\nDesign rows with an appreciable barrier: {len(sig)} of {len(c)}; cortex limiting in {int((sig.compartimento_limitante == 'cortex').sum())}")
print("mean marginal split:", sig[[k for k in sig.columns if k.startswith("frac_marginal_")]].mean().round(3).to_dict())
s = os.path.join(D, "sobol_por_poblacion.csv")
if os.path.exists(s):
    print("\nSobol (dominant total-order indices):\n", pd.read_csv(s).head(20).to_string(index=False))
