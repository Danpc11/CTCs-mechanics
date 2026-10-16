"""Figura principal de la Fase 2: eje tumoral frente a fondo leucocitario."""
import sys, os, json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L
import paso2b_tumoral_leucocito as P2

try:
    apply_figure_style(sizes=(8, 7, 6))
except NameError:
    mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8,
                         "axes.labelsize": 8, "legend.fontsize": 7,
                         "xtick.labelsize": 6, "ytick.labelsize": 6,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 150, "savefig.dpi": 300,
                         "axes.linewidth": 0.6, "lines.linewidth": 1.3})

B = json.load(open(os.path.join(HERE, "paso2b.json")))
c = pd.read_csv(os.path.join(HERE, "pureza_rendimiento_caudal.csv"))
sob = pd.read_csv(os.path.join(HERE, "sensibilidad_sobol.csv"))

C_TUM, C_LEU = "#b03a1a", "#1b4f9c"
fig, axes = plt.subplots(1, 3, figsize=(7.3, 2.6))

# ---- (a) presion critica aparente frente a la velocidad de garganta
ax = axes[0]
U = np.logspace(-2, 5, 500)
for nm, v in P2.TUMORAL.items():
    tau = v["tau"] if v["tau"] is not None else 1.90
    dPc = v["dP"] * (1.0 + L.deborah(tau, U, P2.L_CARGA))
    ax.loglog(U, dPc, color=C_TUM, lw=1.4, zorder=4)

for nm, v in P2.LEUCO.items():
    dPc = v["dP"] * (1.0 + L.deborah(v["tau"], U, P2.L_CARGA))
    ax.loglog(U, dPc, color=C_LEU, lw=1.1, zorder=3)
ax.text(3.0, 2.2e3, "tumoral\n(MCF-7, SK-BR-3,\nMDA-MB-231)", fontsize=6,
        color=C_TUM, va="center", ha="left")
ax.text(3.0e2, 3.5e-2, "leucocito\n(PBMC, granulocito,\nHL-60)", fontsize=6,
        color=C_LEU, va="center", ha="left")
ax.loglog(U, 0.01 * U, color="#444444", lw=0.9, ls=(0, (3, 1.5)), zorder=2)
ax.text(2.2e-2, 1.1e-3, "presión motriz\n10 kPa por mm/s", fontsize=6,
        color="#444444", ha="left", va="bottom")
ax.set_xlabel("velocidad en la garganta $U$ (µm/s)")
ax.set_ylabel("presión crítica aparente (kPa)")
ax.set_title("El leucocito se endurece mucho\nmás tarde que la tumoral")
ax.set_xlim(1e-2, 1e5)
ax.set_ylim(5e-4, 1e4)

# ---- (b) rendimiento y arrastre leucocitario: la ventana de caudal
ax = axes[1]
for beta, lw, al in [(0.01, 1.6, 1.0), (0.30, 1.0, 0.5)]:
    g = c[(c.beta_kPa_s_um == beta) &
          (c.referencia_leucocitaria == "PBMC")].sort_values("U_um_s")
    ax.semilogx(g.U_um_s, g.rendimiento_tumoral, color=C_TUM, lw=lw, alpha=al,
                zorder=4)
    ax.semilogx(g.U_um_s, g.arrastre_leucocitario, color=C_LEU, lw=lw, alpha=al,
                zorder=3)
g = c[(c.beta_kPa_s_um == 0.01) &
      (c.referencia_leucocitaria == "PBMC")].sort_values("U_um_s")
ok = g[g.separacion > 0.9]
ax.axvspan(ok.U_um_s.min(), ok.U_um_s.max(), color="#cfcfcf", alpha=0.45,
           zorder=1, lw=0)
ax.text(np.sqrt(ok.U_um_s.min() * ok.U_um_s.max()), 0.47,
        "ventana útil\n4,4 décadas", fontsize=6, ha="center", va="center")
ax.text(2e-2, 0.88, "tumoral", fontsize=6, color=C_TUM, ha="left")
ax.text(2e-2, 0.06, "leucocito", fontsize=6, color=C_LEU, ha="left")
ax.set_xlabel("velocidad en la garganta $U$ (µm/s)")
ax.set_ylabel("fracción retenida")
ax.set_title("A 10 kPa por mm/s la ventana es\nancha; a 300 kPa se cierra")
ax.set_xlim(1e-2, 1e5)
ax.set_ylim(-0.05, 1.08)

# ---- (c) que sobrevive a la incertidumbre de los parametros transferidos
ax = axes[2]
k_abs = "dP_critica_tumoral_kPa"
k_sel = "selectividad_constriccion_tumoral_sobre_leucocito"
pars = ["E_app_tumoral", "E_app_leucocito", "a_tumoral",
        "Pi", "C", "an_a", "r_env", "r_nuc"]
ETQ = {"E_app_tumoral": "módulo aparente tumoral", "a_tumoral": "radio tumoral",
       "E_app_leucocito": "módulo aparente leucocito", "Pi": "$\\Pi$ citoplasma",
       "C": "$K_s/G_s$ córtex", "an_a": "radio núcleo/célula",
       "r_env": "envoltura nuclear", "r_nuc": "nucleoplasma"}
MED = {"E_app_tumoral": True, "E_app_leucocito": True, "a_tumoral": True,
       "Pi": False, "C": False, "an_a": False, "r_env": False, "r_nuc": False}
y = np.arange(len(pars))[::-1]
st_abs = [float(sob[(sob.salida == k_abs) & (sob.parametro == p)].ST.iloc[0])
          for p in pars]
st_sel = [float(sob[(sob.salida == k_sel) & (sob.parametro == p)].ST.iloc[0])
          for p in pars]
ax.barh(y + 0.19, st_abs, height=0.34, color="#9a9a9a", zorder=3,
        label="presión crítica absoluta")
ax.barh(y - 0.19, np.maximum(st_sel, 2e-4), height=0.34, color="#2a7f62",
        zorder=3, label="selectividad tumoral/leucocito")
ax.set_yticks(y)
ax.set_yticklabels([ETQ[p] + ("" if MED[p] else " ·") for p in pars])
for i, p in enumerate(pars):
    if not MED[p]:
        ax.get_yticklabels()[i].set_color("#8a4a00")
ax.set_xlabel("índice de Sobol total $S_T$")
ax.set_title("Ningún parámetro transferido (·) afecta\na la selectividad")
ax.set_xlim(0, 0.85)
ax.annotate("presión crítica absoluta", xy=(st_abs[3], y[3] + 0.19),
            xytext=(0.40, y[3] + 1.6), fontsize=6,
            color="#6a6a6a", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="#9a9a9a"))
ax.annotate("selectividad\ntumoral/leucocito", xy=(st_sel[1], y[1] - 0.19),
            xytext=(0.47, y[1] - 2.3), fontsize=6,
            color="#2a7f62", va="center", ha="left",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="#2a7f62"))

fig.tight_layout()
fig.subplots_adjust(top=0.76)
for i, axx in enumerate(axes):
    bb = axx.get_position()
    fig.text(bb.x0 - (0.055 if i < 2 else 0.145), 0.995, "abc"[i],
             fontweight="bold", fontsize=10, va="top", ha="left")
fig.savefig(os.path.join(HERE, "figura_caudal.png"), bbox_inches="tight")

r = fig.canvas.get_renderer()
txt = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
       if t.get_text().strip() and t.get_visible()]
ov = [(a.get_text()[:20], b.get_text()[:20])
      for i, (a, ba) in enumerate(txt) for b, bb in txt[i + 1:]
      if ba.overlaps(bb)]
print("solapes:", len(ov))
for o in ov[:10]:
    print("   ", o)
print("ventana util (um/s):", float(ok.U_um_s.min()), float(ok.U_um_s.max()))
print("ST selectividad de los transferidos (max):",
      max(st_sel[i] for i, p in enumerate(pars) if not MED[p]))
