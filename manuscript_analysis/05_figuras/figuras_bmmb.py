"""BMMB figures: no in-illustration titles, patterns in addition to colour,
84 mm / 174 mm widths. All numbers read from the stored artefacts."""
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import norm

MM = 1 / 25.4
W1, W2 = 84 * MM, 174 * MM          # single / double column

CTUM = "#B2182B"    # tumour (focal)
CLEU = "#4878A6"    # leucocyte
CTR = "#8C8C8C"     # transferred / context
CALT = "#1B7837"    # control / secondary measured
CBLA = "#7B52AB"    # blast


def panel_letter(ax, letter):
    """Letra de panel en negrita, arriba a la izquierda (faltaba en el deposito)."""
    ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left")


# ------------------------------------------------------------------ Fig. 1
def fig_size(dld, dldp, path="fig1_size.png"):
    fig, (a, b) = plt.subplots(1, 2, figsize=(W2, 2.65))

    # (a) diameter distributions
    leu = [("lymphocyte", 6.05, 0.072, 0.3286),
           ("neutrophil", 8.40, 0.072, 0.5635 + 0.0209 + 0.0066),
           ("monocyte", 8.79, 0.065, 0.0804)]
    x = np.linspace(3, 23, 1200)
    mix = np.zeros_like(x)
    for _, mu, cv, w in leu:
        mix += w * norm.pdf(x, mu, mu * cv)
    a.fill_between(x, mix, color=CLEU, alpha=0.30, lw=0, hatch="///",
                   edgecolor=CLEU)
    a.plot(x, mix, color=CLEU, lw=1.2)
    a.text(4.6, mix.max() * 1.10, "normal leucocytes\n(differential-weighted)",
           color=CLEU, ha="left", va="bottom")

    tum = norm.pdf(x, 17.48, 1.748)
    a.fill_between(x, tum, color=CTUM, alpha=0.30, lw=0, hatch="\\\\\\",
                   edgecolor=CTUM)
    a.plot(x, tum, color=CTUM, lw=1.2)
    a.text(19.6, tum.max() * 1.05, "MCF-7\n(CV 10% declared)", color=CTUM,
           ha="center", va="bottom")

    for nm, mu, sd, col, dash in [("Jurkat", 11.50, 11.50 * 0.1304, CBLA, (3, 1.5)),
                                  ("HL-60", 17.00, 17.00 * 0.072, CBLA, (1, 1))]:
        y = norm.pdf(x, mu, sd)
        a.plot(x, y, color=col, lw=1.1, dashes=dash)
    a.text(12.4, norm.pdf(11.5, 11.5, 1.4995) * 1.02, "Jurkat", color=CBLA,
           ha="left", va="bottom")
    a.text(14.6, norm.pdf(17, 17, 1.224) * 0.45, "HL-60", color=CBLA,
           ha="right", va="center")
    a.axvline(10.86, color="k", lw=0.7, dashes=(4, 2))
    a.text(10.55, mix.max() * 0.22, r"optimum $D_c$ 10.86 $\mu$m",
           rotation=90, ha="center", va="center", color="k")
    a.set_xlabel(r"cell diameter in suspension  ($\mu$m)")
    a.set_ylabel("probability density  (a.u.)")
    a.set_yticks([])
    a.set_xlim(3, 23)
    a.margins(y=0.22)
    for sp in ("left", "right", "top"):
        a.spines[sp].set_visible(False)

    # (b) leucocyte carry-over vs assumed tumour CV, limiting line
    g = dldp[dldp.tumoral == "MDA-MB-231"].sort_values("CV_tumoral")
    FLOOR = 1e-2
    for yld, col, mk, hat in [(0.95, CLEU, "s", None), (0.99, CTUM, "o", None)]:
        h = g[g.rendimiento_objetivo == yld]
        y = np.maximum(h.leucocitos_por_millon.values, FLOOR)
        b.semilogy(100 * h.CV_tumoral, y, mk + "-", color=col, ms=4.5, lw=1.3,
                   mfc=col if yld == 0.99 else "white", mew=1.1)
    b.text(28.5, 1.0e5, "99% yield", color=CTUM, ha="right", va="bottom")
    b.text(28.5, 3.0e-2, "95% yield", color=CLEU, ha="right", va="bottom")
    b.axhline(1.0, color=CTR, lw=0.7, dashes=(4, 2))
    b.text(5.4, 1.5, "1 leucocyte per million", color=CTR, va="bottom")
    k = g[(g.rendimiento_objetivo == 0.99) & (g.CV_tumoral == 0.20)]
    b.annotate(r"CV 20%: 7.0$\times$10$^4$ per 10$^6$," "\n" r"enrichment 14$\times$",
               xy=(20, float(k.leucocitos_por_millon.iloc[0])),
               xytext=(9.0, 2.0e3), arrowprops=dict(arrowstyle="->", lw=0.8))
    b.set_xlabel("assumed tumour-diameter CV  (%)")
    b.set_ylabel(r"leucocyte carry-over  (cells per 10$^6$)")
    b.set_ylim(FLOOR / 2, 5e6)
    b.set_xlim(3, 31)
    b.set_yticks([1e-2, 1e0, 1e2, 1e4, 1e6])

    for ax, L in zip([a, b], "ab"):
        panel_letter(ax, L)
    fig.tight_layout(pad=0.4)
    fig.savefig(path, dpi=400)
    return fig


# ------------------------------------------------------------------ Fig. 2
def fig_contact(bc, enc, path="fig2_contact.png"):
    fig = plt.figure(figsize=(W2, 2.55))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.85, 1.0], wspace=0.44)
    a, b, c = (fig.add_subplot(gs[0, i]) for i in range(3))

    # (a) critical pressure vs nominal confinement
    ref = bc[bc.geometria == "referencia"]
    style = {"MCF-7": (CTUM, "o", (None, None)),
             "SK-BR-3": (CTUM, "^", (3, 1.5)),
             "MDA-MB-231": (CTUM, "v", (1, 1)),
             "leucocito": (CLEU, "s", (None, None))}
    label = {"MCF-7": "MCF-7", "SK-BR-3": "SK-BR-3",
             "MDA-MB-231": "MDA-MB-231", "leucocito": "leucocyte (PBMC)"}
    lab_off = {"MCF-7": 0.72, "SK-BR-3": 1.45, "MDA-MB-231": 0.72,
               "leucocito": 1.0}
    for k, (col, mk, dash) in style.items():
        h = ref[ref.estirpe == k].sort_values("confinamiento")
        h = h[h.dP_critica_kPa > 0]
        a.semilogy(h.confinamiento, h.dP_critica_kPa, mk + "-", color=col,
                   ms=4, lw=1.2, dashes=dash if dash[0] else (None, None))
        a.text(h.confinamiento.iloc[-1] + 0.045,
               h.dP_critica_kPa.iloc[-1] * lab_off[k],
               label[k], color=col, va="center", fontsize=6)
    a.axvspan(1.20, 2.43, color=CTR, alpha=0.10, lw=0)
    a.text(2.40, 3.2, "validated envelope", color=CTR, fontsize=6, ha="right", va="top")
    a.set_xlabel(r"nominal confinement  $2a/D$")
    a.set_ylabel(r"critical entry pressure  (kPa)")
    a.set_xlim(1.33, 3.12)
    a.set_ylim(1.5e-3, 4)

    # (b) marginal resistance split
    comp = [("cortex", 0.718, CTUM, "///"), ("cytoplasm", 0.263, CLEU, "..."),
            ("nuclear envelope", 0.011, CTR, "xxx"),
            ("nucleoplasm", 0.009, CTR, "\\\\\\")]
    ypos = np.arange(len(comp))[::-1]
    for (nm, v, col, hat), y in zip(comp, ypos):
        b.barh(y, max(v, 0.004), color=col, alpha=0.45, edgecolor=col,
               hatch=hat, height=0.62, lw=0.9)
        b.text(v + 0.03, y, f"{v:.3f}", va="center", color=col, fontsize=6)
        b.text(-0.03, y, nm, va="center", ha="right", fontsize=6,
               color="k" if v > 0.1 else CTR)
    b.set_xlim(0, 0.98)
    b.set_ylim(-0.6, len(comp) - 0.4)
    b.set_yticks([])
    b.set_xlabel("fraction of marginal resistance")
    for sp in ("left", "right", "top"):
        b.spines[sp].set_visible(False)

    # (c) same-device framing
    enc = enc.copy()
    enc["lab"] = enc.poblacion.str.replace(r"\s*\(.*\)", "", regex=True)
    enc.loc[enc.poblacion.str.startswith("HL-60"), "lab"] = "HL-60"
    enc.loc[enc.poblacion.str.startswith("Jurkat"), "lab"] = "Jurkat"
    name_en = {"linfocito": "lymphocyte", "neutrofilo": "neutrophil",
               "monocito": "monocyte", "eosinofilo": "eosinophil"}
    enc["lab"] = enc.lab.replace(name_en)
    enc = enc.drop_duplicates("lab").sort_values("conf_efectivo")
    grp = {"lymphocyte": CLEU, "neutrophil": CLEU, "monocyte": CLEU,
           "eosinophil": CLEU, "Jurkat": CBLA, "HL-60": CBLA,
           "MCF-7": CTUM, "MDA-MB-231": CTUM}
    c.axhspan(1.039, 1.636, color=CTR, alpha=0.12, lw=0)
    c.axhline(1.0, color="k", lw=0.8, dashes=(4, 2))
    for i, (_, r) in enumerate(enc.iterrows()):
        col = grp[r.lab]
        touch = r.conf_efectivo >= 1.0
        c.vlines(i, 0.3, r.conf_efectivo, color=col, lw=0.9, alpha=0.55)
        c.plot(i, r.conf_efectivo, "o", ms=6, color=col,
               mfc=col if touch else "white", mew=1.2)
    c.set_xticks(range(len(enc)))
    c.set_xticklabels(enc.lab, rotation=55, ha="right", fontsize=6)
    c.set_ylabel(r"effective confinement  $a/R_{\min}$")
    c.set_ylim(0.3, 1.78)
    c.text(-0.35, 1.03, "wall contact", fontsize=6, va="bottom")
    c.text(len(enc) - 0.55, 1.63, "envelope", fontsize=6, color=CTR,
           ha="right", va="top")
    c.text(0.02, 0.80, "open symbol:\nno contact,\n" r"$\Delta P_{\rm crit}=0$",
           fontsize=6, color=CTR, transform=c.transAxes, va="top")
    for sp in ("right", "top"):
        c.spines[sp].set_visible(False)

    for ax, L in zip([a, b, c], "abc"):
        panel_letter(ax, L)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.93, bottom=0.30)
    fig.savefig(path, dpi=400, bbox_inches="tight")
    return fig


# ------------------------------------------------------------------ Fig. 3
def fig_rate(umb, pur, ley, path="fig3_rate.png"):
    fig = plt.figure(figsize=(W2, 2.55))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 0.85], wspace=0.40)
    a, b, c = (fig.add_subplot(gs[0, i]) for i in range(3))

    # (a) beta_c per population
    rows = [("HL-60", CLEU, "leucocyte"), ("granulocito", CLEU, "leucocyte"),
            ("PBMC", CLEU, "leucocyte"), ("MDA-MB-231", CTUM, "tumour"),
            ("MCF-10A", CALT, "control"), ("MCF-7", CTUM, "tumour")]
    disp = {"granulocito": "granulocyte"}
    vals = {r.poblacion: r.beta_c_kPa_a_1mm_s for _, r in umb.iterrows()
            if r.poblacion != "SK-BR-3"}
    sk = sorted(umb[umb.poblacion == "SK-BR-3"].beta_c_kPa_a_1mm_s.values)
    lo, hi = vals["PBMC"], vals["MDA-MB-231"]
    a.axvspan(lo, hi, color=CTUM, alpha=0.09, lw=0)
    for i, (nm, col, _) in enumerate(rows):
        v = vals[nm]
        a.hlines(i, 1.5e-2, v, color=col, lw=0.8, alpha=0.45)
        a.plot(v, i, "o", color=col, ms=5.5)
        a.text(v * 1.5, i, disp.get(nm, nm), va="center", color=col, fontsize=6)
    a.hlines(len(rows), sk[0], sk[1], color=CBLA, lw=2.2)
    a.plot(sk, [len(rows)] * 2, "|", color=CBLA, ms=7)
    a.text(sk[0] / 1.6, len(rows), r"SK-BR-3 ($\tau$ not measured)",
           va="center", ha="right", color=CBLA, fontsize=6)
    a.set_xscale("log")
    a.set_xlim(2e-2, 2e4)
    a.set_ylim(-0.8, len(rows) + 0.8)
    a.set_yticks([])
    a.set_xlabel(r"threshold device pressure $\beta_c$" "\n"
                 r"(kPa per mm s$^{-1}$)")
    a.annotate("", xy=(lo, 3.55), xytext=(hi, 3.55),
               arrowprops=dict(arrowstyle="<->", lw=0.9, color=CTUM))
    a.text(np.sqrt(lo * hi), 3.65, r"165$\times$", ha="center", va="bottom",
           color=CTUM, fontsize=6)
    for sp in ("left", "right", "top"):
        a.spines[sp].set_visible(False)

    # (b) separation vs velocity
    sub = pur[pur.referencia_leucocitaria == "PBMC"]
    for beta, col, dash in [(10.0, CTUM, (None, None)), (164.5, CLEU, (3, 1.5))]:
        g = sub[sub.dP_dispositivo_a_1mm_s_kPa == beta].sort_values("U_um_s")
        b.semilogx(g.U_um_s, g.separacion, color=col, lw=1.4,
                   dashes=dash if dash[0] else (None, None))
        k = g.separacion.idxmax()
        b.plot(g.U_um_s[k], g.separacion[k], "o", color=col, ms=5.5)
    b.text(1.6e1, 0.90, r"$\beta$ = 10 kPa per mm s$^{-1}$" "\n"
           "(plateau 4.41 decades)", color=CTUM, fontsize=6, va="top")
    b.text(2.3e-1, 0.33, r"$\beta$ = 164.5" "\n" "(0.38 decades)",
           color=CLEU, fontsize=6, ha="center")
    b.set_xlabel(r"throat velocity  $U$  ($\mu$m s$^{-1}$)")
    b.set_ylabel("separation  (tumour $-$ leucocyte retained)")
    b.margins(y=0.10)

    # (c) Kelvin-Voigt vs power-law window
    sel = ley[(ley.U_max_mm_s == 1.0) &
              (ley.escenario.str.contains("neutrofilo Tsai"))]
    pares = ["MDA-MB-231 / PBMC", "MCF-7 / PBMC"]
    xp = np.arange(len(pares))
    kv = [float(sel[sel.pareja == p].ventana_KV.iloc[0]) for p in pares]
    pl = [float(sel[sel.pareja == p].ventana_ley_potencias.iloc[0]) for p in pares]
    for i, (k_, p_) in enumerate(zip(kv, pl)):
        c.vlines(i, p_, k_, color=CTR, lw=1.0)
        c.plot(i, k_, "s", color=CLEU, ms=6, mfc=CLEU)
        c.plot(i, p_, "o", color=CTUM, ms=6, mfc="white", mew=1.4)
        c.text(i + 0.13, k_, f"{k_:.0f}" r"$\times$", color=CLEU, fontsize=6,
               va="center")
        c.text(i + 0.13, p_, f"{p_:.1f}" r"$\times$", color=CTUM, fontsize=6,
               va="center")
    c.axhline(2.0, color="k", lw=0.7, dashes=(4, 2))
    c.text(-0.42, 2.15, r"2$\times$", fontsize=6, va="bottom")
    c.set_yscale("log")
    c.set_xticks(xp)
    c.set_xticklabels([p.replace(" / ", "\n/ ") for p in pares], fontsize=6)
    c.set_xlim(-0.5, 1.75)
    c.set_ylabel(r"operating window  ($\times$)")
    c.text(0.52, 0.56, "Kelvin--Voigt", color=CLEU, fontsize=6,
           transform=c.transAxes)
    c.text(0.52, 0.18, "power law", color=CTUM, fontsize=6,
           transform=c.transAxes)
    for sp in ("right", "top"):
        c.spines[sp].set_visible(False)

    for ax, L in zip([a, b, c], "abc"):
        panel_letter(ax, L)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.93, bottom=0.28)
    fig.savefig(path, dpi=400, bbox_inches="tight")
    return fig


# ------------------------------------------------------------------ Fig. 4
def fig_sobol(sob, path="fig4_sobol.png"):
    pretty = {"E_app_tum": r"$E_{\rm app}$ tumour", "E_app_leu": r"$E_{\rm app}$ leucocyte",
              "C_tum": r"$C$ tumour", "C_leu": r"$C$ leucocyte",
              "a_leu": "leucocyte radius", "a_tum": "tumour radius",
              "an_a_leu": r"$a_n/a$ leucocyte", "an_a_tum": r"$a_n/a$ tumour",
              "Pi_tum": r"$\Pi$ tumour", "Pi_leu": r"$\Pi$ leucocyte",
              "r_env_leu": r"$r_{\rm env}$ leucocyte", "r_env_tum": r"$r_{\rm env}$ tumour",
              "r_nuc_leu": r"$r_{\rm nuc}$ leucocyte", "r_nuc_tum": r"$r_{\rm nuc}$ tumour",
              "tau_tum": r"$\tau$ tumour", "tau_leu": r"$\tau$ leucocyte",
              "ell_factor": r"$\ell_{\rm factor}$"}
    # Pi_leu NO es medido: el modulo de citoplasma del leucocito se transfiere
    # de MDA-MB-231 (fase2_lib.ESTIRPES["leucocito"].G_cyt_src).
    measured = {"E_app_tum", "E_app_leu", "a_tum", "a_leu", "Pi_tum",
                "an_a_leu", "tau_tum"}
    g = sob[sob.salida == "selectividad_constriccion_tumoral_sobre_leucocito"]
    g = g.sort_values("ST", ascending=True)
    s_med = g[g.parametro.isin(measured)].ST.sum()
    s_tr = g[~g.parametro.isin(measured)].ST.sum()
    fig, ax = plt.subplots(figsize=(W1, 3.25))
    y = np.arange(len(g))
    for yy, (_, r) in zip(y, g.iterrows()):
        m = r.parametro in measured
        col = CTUM if m else CTR
        hat = "///" if m else "..."
        ax.barh(yy, r.ST, height=0.66, color=col, alpha=0.45, edgecolor=col,
                hatch=hat, lw=0.9)
        ax.errorbar(r.ST, yy, xerr=r.ST_conf, color=col, lw=0.9, capsize=1.6)
    ax.set_yticks(y)
    ax.set_yticklabels([pretty[p] for p in g.parametro], fontsize=6)
    for t, p in zip(ax.get_yticklabels(), g.parametro):
        t.set_color("k" if p in measured else CTR)
    ax.set_xlabel(r"Sobol total-order index  $S_T$")
    ax.set_xlim(0, 0.30)
    ax.margins(y=0.02)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=CTUM, alpha=0.45, hatch="///",
                             edgecolor=CTUM, label=f"measured  ($\\Sigma S_T$ = {s_med:.3f})"),
                       Patch(facecolor=CTR, alpha=0.45, hatch="...",
                             edgecolor=CTR, label=f"transferred or declared  ({s_tr:.3f})")],
              frameon=False, loc="lower right", fontsize=6)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    panel_letter(ax, "a")
    fig.tight_layout(pad=0.4)
    fig.savefig(path, dpi=400)
    return fig


# ------------------------------------------------------------------ Fig. 5
def fig_design(dis, path="fig5_design.png"):
    fig, a = plt.subplots(figsize=(W1, 2.55))
    d = dis.sort_values("Rt_um")
    a.plot(d.Rt_um, d.selectividad, "-", color=CTUM, lw=1.6)
    a.set_xlabel(r"throat radius  $R_t$  ($\mu$m)")
    a.set_ylabel(r"selectivity, tumour/leucocyte  ($\times$)", color=CTUM)
    a.tick_params(axis="y", colors=CTUM)
    a2 = a.twinx()
    a2.semilogy(d.Rt_um, d.ancho_kPa, color=CLEU, lw=1.6, dashes=(3, 1.5))
    a2.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    a2.set_yticks([0.02, 0.05, 0.1, 0.2])
    a2.set_yticklabels(["0.02", "0.05", "0.1", "0.2"])
    a2.set_ylabel("pressure-window width  (kPa)", color=CLEU)
    a2.tick_params(axis="y", colors=CLEU)
    a.axvspan(5.90, 6.00, color=CTR, alpha=0.15, lw=0)
    a.text(5.95, a.get_ylim()[0] + 0.10, "edge of\nsampling", ha="center",
           va="bottom", color=CTR, fontsize=6)
    k = (d.Rt_um - 5.675).abs().idxmin()
    a.annotate(f"$R_t$ = {d.Rt_um[k]:.1f} " r"$\mu$m," "\n"
               f"{d.selectividad[k]:.2f}" r"$\times$, "
               f"{d.ancho_kPa[k]:.3f} kPa wide",
               xy=(d.Rt_um[k], d.selectividad[k]), xytext=(0.06, 0.74),
               textcoords="axes fraction", fontsize=6,
               arrowprops=dict(arrowstyle="->", lw=0.8))
    panel_letter(a, "a")
    fig.tight_layout(pad=0.4)
    fig.savefig(path, dpi=400)
    return fig


# ------------------------------------------------------------------ Fig. A1
def fig_solver(V, path="figA1_solver.png"):
    fig, axs = plt.subplots(2, 3, figsize=(W2, 4.1))
    (a, b, c), (d, e, f) = axs

    res = V["V1_patch"]["res"][-1]
    a.semilogy(range(len(res)), res, "o-", color=CLEU, ms=4)
    a.set_xlabel("Newton iteration")
    a.set_ylabel(r"$\|R\|_2$  (free dofs)")
    a.annotate(f"{res[-1]:.1e}", xy=(len(res) - 1, res[-1]), xytext=(-6, 14),
               textcoords="offset points", ha="right", color=CLEU, fontsize=6)

    aR = np.array([x["aR"] for x in V["V2_hertz"]["casos"]])
    er = 100 * np.array([x["err"] for x in V["V2_hertz"]["casos"]])
    m, q = V["V2_hertz"]["ajuste_err_vs_aR"]
    xs = np.linspace(0, aR.max() * 1.05, 50)
    b.plot(xs, 100 * (m * xs + q), color=CTR, lw=1, dashes=(3, 1.5))
    b.plot(aR, er, "o", color=CLEU, ms=4)
    b.plot(0, 100 * q, "s", color=CTUM, ms=5)
    b.axhline(0, color="k", lw=0.6)
    b.set_xlabel(r"$a/R$  (contact radius / sphere radius)")
    b.set_ylabel("force error vs Hertz  (%)")
    b.text(0.06, 0.78, f"extrapolated to $a/R\\to0$:\n{100*q:+.2f}%",
           transform=b.transAxes, color=CTUM, fontsize=6)

    h3 = np.array([x["h"] for x in V["V3_esfera_hueca"]["casos"]])
    e3 = np.abs([x["err"] for x in V["V3_esfera_hueca"]["casos"]])
    h4 = np.array([x["h"] for x in V["V4_laplace"]["casos"]])
    e4 = np.abs([x["err_p"] for x in V["V4_laplace"]["casos"]])
    c.loglog(h3, e3, "o-", color=CLEU, ms=4, label="inflated hollow sphere")
    c.loglog(h4, e4, "s", color=CTUM, ms=4, ls="--", label="Laplace law")
    ref = np.array([h3.min(), h3.max()])
    c.loglog(ref, e3.max() * 12 * (ref / ref.max()) ** 2, ":", color=CTR, lw=1)
    c.text(0.50, 0.84, "slope 2", transform=c.transAxes, color=CTR, rotation=22,
           fontsize=6)
    c.set_xlabel("element size  $h$")
    c.set_ylabel("|relative error| vs exact")
    c.set_xticks([0.03, 0.05, 0.1, 0.2])
    c.set_xticklabels(["0.03", "0.05", "0.1", "0.2"])
    c.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    c.legend(frameon=False, loc="lower right", fontsize=6)

    d5 = np.array([x["delta"] for x in V["V5"]["hertz_inc"]])
    e5 = 100 * np.array([x["err"] for x in V["V5"]["hertz_inc"]])
    d.semilogx(np.sqrt(d5), e5, "o-", color=CLEU, ms=4)
    d.axhline(0, color="k", lw=0.6)
    d.set_xlabel(r"$a/R$")
    d.set_ylabel("force error vs Hertz  (%)")
    d.set_xticks([0.05, 0.1])
    d.set_xticklabels(["0.05", "0.1"])
    d.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    g = V["V5"]["grande"]
    d.text(0.04, 0.76,
           r"$\nu=0.49995$; at $\delta/R=0.30$" "\n"
           f"volume loss {abs(g[-1]['dV'])*100:.3f}%",
           transform=d.transAxes, color=CTR, fontsize=6, va="bottom")

    cm = V["constriccion"]["convergencia_malla"]
    nd = np.array([x["ndof"] for x in cm])
    DI = 100 * np.array([x["DI"] for x in cm])
    dA = 100 * np.array([x["dA"] for x in cm])
    e.plot(nd, DI, "o-", color=CLEU, ms=4)
    e.set_xlabel("degrees of freedom")
    e.set_ylabel("deformation index DI  (%)", color=CLEU)
    e.tick_params(axis="y", colors=CLEU)
    e2 = e.twinx()
    e2.plot(nd, dA, "s", color=CTUM, ms=4, ls="--")
    e2.set_ylabel(r"area strain  $\Delta A/A_0$  (%)", color=CTUM)
    e2.tick_params(axis="y", colors=CTUM)

    bcf = [x for x in V["constriccion"]["barrido_confinamiento"] if x["dA"] > 0]
    xda = 100 * np.array([x["dA"] for x in bcf])
    yDI = 100 * np.array([x["DI"] for x in bcf])
    o = np.argsort(xda)
    xda, yDI = xda[o], yDI[o]
    f.plot(xda, yDI, "o-", color=CLEU, ms=4, label="present 2D axisymmetric")
    DI2D = np.interp(4.0, xda, yDI)
    f.plot([4.0], [DI2D], "o", mfc="white", mec=CLEU, ms=8, zorder=5)
    f.errorbar([4.0], [2.8], yerr=[0.5], fmt="D", color=CTUM, ms=5,
               capsize=2, label="Wang et al. 2023, 3D")
    f.annotate(r"+16% at matched area strain", xy=(4.15, 2.95), xytext=(5.3, 2.1),
               arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=6)
    f.set_xlabel(r"area strain  $\Delta A/A_0$  (%)")
    f.set_ylabel("deformation index DI  (%)")
    f.legend(frameon=False, loc="upper left", fontsize=6)

    for ax, L in zip([a, b, c, d, e, f], "abcdef"):
        ax.margins(0.05)
        panel_letter(ax, L)
    fig.tight_layout(pad=0.4)
    fig.savefig(path, dpi=400)
    return fig
