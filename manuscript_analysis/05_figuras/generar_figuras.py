"""Regenera TODAS las figuras del manuscrito desde manuscript_analysis/datos/ (600 dpi, >= 500 que pide JMBBM).

    cd manuscript_analysis/05_figuras && python generar_figuras.py

Output: <repo>/figures/Fig1.png ... Fig7.png, FigS1.png and graphical_abstract.png.
Fig. 2 y Fig. 3 usan la reduccion corregida (ver regenerar_fig2_fig3.py para el detalle).
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import pandas as pd
mpl.rcParams["savefig.dpi"] = 600
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic", "mathtext.bf": "Liberation Sans:bold"})          # TrueType fonts embedded in PDF output
FMT = os.environ.get("CTC_FIG_FMT", "png")  # png or pdf (vector)
mpl.rcParams["font.family"] = "Liberation Sans"   # metric-compatible with Arial (JMBBM accepts Arial)
AQUI = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(AQUI, "..", "datos") + os.sep
OUT = os.path.abspath(os.path.join(AQUI, "..", "..", "figures")); os.makedirs(OUT, exist_ok=True)
sys.path[:0] = [os.path.join(AQUI, "..", "03_modalidades"), os.path.join(AQUI, "..", "01_solver")]
os.chdir(AQUI)

# fuente comun con las mismas correcciones que regenerar_fig2_fig3.py
K = 1.0 if os.environ.get("CTC_DIRECT") == "1" else 6.56   # 1.0: thresholds computed directly
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

fig(f"Fig2.{FMT}", G["fig_size"], pd.read_csv(D + "dld_revisado.csv"), pd.read_csv(D + "dld_revisado_pureza.csv"))
t = pd.read_csv(D + "tabla2_corregida_directa.csv").dropna(subset=["corregido"])
t = t[~(((t.estirpe == "SK-BR-3") & (t.conf_nominal == 2.43)) | ((t.estirpe == "MDA-MB-231") & (t.conf_nominal == 1.90)))]
bc = pd.DataFrame(dict(geometria="referencia", estirpe=t.estirpe, confinamiento=t.conf_nominal, dP_critica_kPa=t.corregido))
fig(f"Fig4.{FMT}", G["fig_contact"], bc, pd.read_csv(D + "encuadre_mismo_dispositivo.csv"))
umb = (pd.read_csv(D + "umbral_dispositivo.csv") if os.environ.get("CTC_DIRECT") == "1" else
       pd.read_csv(D + "umbral_dispositivo_corregido.csv").assign(beta_c_kPa_a_1mm_s=lambda u: u.beta_c_corregido_kPa_a_1mm_s))
fig(f"Fig5.{FMT}", G["fig_rate"], umb, pd.read_csv(D + "pureza_rendimiento_caudal.csv"), pd.read_csv(D + "cota_ley_potencias_ventana.csv"))
fig(f"Fig6.{FMT}", G["fig_sobol"], pd.read_csv(D + "sobol_por_poblacion.csv"))
fig(f"Fig7.{FMT}", G["fig_design"], pd.read_csv(D + "diseno_optimizado.csv"))
fig(f"FigS1.{FMT}", G["fig_solver"], json.load(open(D + "validacion.json")))
import runpy; runpy.run_path(os.path.join(AQUI, "figuras_concepto.py"))   # Fig1, Fig3, graphical abstract
print("figuras:", ", ".join(hechas), "->", OUT)
