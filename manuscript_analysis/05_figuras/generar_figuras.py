"""Regenera TODAS las figuras del manuscrito desde manuscript_analysis/datos/ (600 dpi, >= 500 que pide JMBBM).

    cd manuscript_analysis/05_figuras && python generar_figuras.py

Salida: <raiz del repo>/figures/Fig1.png ... Fig5.png y FigS1.png.
Fig. 2 y Fig. 3 usan la reduccion corregida (ver regenerar_fig2_fig3.py para el detalle).
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import pandas as pd
mpl.rcParams["savefig.dpi"] = 600
AQUI = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(AQUI, "..", "datos") + os.sep
OUT = os.path.abspath(os.path.join(AQUI, "..", "..", "figures")); os.makedirs(OUT, exist_ok=True)
sys.path[:0] = [os.path.join(AQUI, "..", "03_modalidades"), os.path.join(AQUI, "..", "01_solver")]
os.chdir(AQUI)

# fuente comun con las mismas correcciones que regenerar_fig2_fig3.py
K = 6.56
reg = open("regenerar_fig2_fig3.py", encoding="utf-8").read()
src = open("figuras_bmmb.py", encoding="utf-8").read()
for linea in reg.splitlines():
    if linea.startswith("src = src.replace("):
        exec(linea, {"K": K}, ns := {"src": src}); src = ns["src"]
src = src.replace("dpi=400", "dpi=600")
G = {"K": K}
exec(compile(src, "figuras_bmmb.py", "exec"), G)

hechas = []
def fig(nombre, fn, *args):
    p = os.path.join(OUT, nombre)
    fn(*args, p); hechas.append(nombre)

fig("Fig1.png", G["fig_size"], pd.read_csv(D + "dld_revisado.csv"), pd.read_csv(D + "dld_revisado_pureza.csv"))
t = pd.read_csv(D + "tabla2_corregida_directa.csv").dropna(subset=["corregido"])
t = t[~(((t.estirpe == "SK-BR-3") & (t.conf_nominal == 2.43)) | ((t.estirpe == "MDA-MB-231") & (t.conf_nominal == 1.90)))]
bc = pd.DataFrame(dict(geometria="referencia", estirpe=t.estirpe, confinamiento=t.conf_nominal, dP_critica_kPa=t.corregido))
fig("Fig2.png", G["fig_contact"], bc, pd.read_csv(D + "encuadre_mismo_dispositivo.csv"))
umb = pd.read_csv(D + "umbral_dispositivo_corregido.csv"); umb["beta_c_kPa_a_1mm_s"] = umb["beta_c_corregido_kPa_a_1mm_s"]
fig("Fig3.png", G["fig_rate"], umb, pd.read_csv(D + "pureza_rendimiento_caudal.csv"), pd.read_csv(D + "cota_ley_potencias_ventana.csv"))
fig("Fig4.png", G["fig_sobol"], pd.read_csv(D + "sobol_por_poblacion.csv"))
fig("Fig5.png", G["fig_design"], pd.read_csv(D + "diseno_optimizado.csv"))
fig("FigS1.png", G["fig_solver"], json.load(open(D + "validacion.json")))
print("figuras:", ", ".join(hechas), "->", OUT)
