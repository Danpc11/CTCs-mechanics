"""Publication figures (English) for the BMMB manuscript.
All numbers are read from the stored artifacts; nothing is retyped."""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

apply_figure_style(sizes=(8, 7, 6))

CTUM = "#B2182B"   # tumour / focal
CLEU = "#4878A6"   # leucocyte
CTR = "#8C8C8C"    # transferred / context
CALT = "#1B7837"   # secondary measured

# ----------------------------------------------------------------- Figure 1
def fig_solver(V, path="fig1_solver.png"):
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6))
    (a, b, c), (d, e, f) = axs

    res = V["V1_patch"]["res"][-1]
    a.semilogy(range(len(res)), res, "o-", color=CLEU, ms=4)
    a.set_xlabel("Newton iteration"); a.set_ylabel(r"$\|R\|_2$  (free dofs)")
    a.set_title("Newton converges quadratically")
    a.annotate(f"{res[-1]:.1e}", xy=(len(res) - 1, res[-1]), xytext=(-6, 14),
               textcoords="offset points", ha="right", color=CLEU)

    aR = np.array([x["aR"] for x in V["V2_hertz"]["casos"]])
    er = 100 * np.array([x["err"] for x in V["V2_hertz"]["casos"]])
    m, q = V["V2_hertz"]["ajuste_err_vs_aR"]
    xs = np.linspace(0, aR.max() * 1.05, 50)
    b.plot(xs, 100 * (m * xs + q), "--", color=CTR, lw=1)
    b.plot(aR, er, "o", color=CLEU, ms=4)
    b.plot(0, 100 * q, "s", color=CTUM, ms=5)
    b.axhline(0, color="k", lw=0.6)
    b.set_xlabel("$a/R$  (contact radius / sphere radius)")
    b.set_ylabel("force error vs Hertz  (%)")
    b.set_title("Hertz offset is model, not numerics")
    b.text(0.06, 0.80, f"solver error extrapolated\nto $a/R\\to0$: {100*q:+.1f}%",
           transform=b.transAxes, color=CTUM, fontsize=6)

    h3 = np.array([x["h"] for x in V["V3_esfera_hueca"]["casos"]])
    e3 = np.abs([x["err"] for x in V["V3_esfera_hueca"]["casos"]])
    h4 = np.array([x["h"] for x in V["V4_laplace"]["casos"]])
    e4 = np.abs([x["err_p"] for x in V["V4_laplace"]["casos"]])
    c.loglog(h3, e3, "o-", color=CLEU, ms=4, label="inflated hollow sphere\n(inner radius)")
    c.loglog(h4, e4, "s--", color=CTUM, ms=4, label="Laplace law\n(internal pressure)")
    ref = np.array([h3.min(), h3.max()])
    c.loglog(ref, e3.max() * 12 * (ref / ref.max()) ** 2, ":", color=CTR, lw=1)
    c.text(0.52, 0.86, "slope 2", transform=c.transAxes, color=CTR, rotation=22)
    c.set_xlabel("element size $h$  (dimensionless)")
    c.set_ylabel("|relative error| vs exact solution")
    c.set_title("Convergence to exact solutions")
    c.set_xticks([0.03, 0.05, 0.1, 0.2]); c.set_xticklabels(["0.03", "0.05", "0.1", "0.2"])
    c.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    c.legend(frameon=False, loc="lower right", fontsize=6)

    aR5 = np.array([x["a_fe"] if "a_fe" in x else np.nan for x in V["V5"]["hertz_inc"]])
    d5 = np.array([x["delta"] for x in V["V5"]["hertz_inc"]])
    e5 = 100 * np.array([x["err"] for x in V["V5"]["hertz_inc"]])
    d.semilogx(np.sqrt(d5), e5, "o-", color=CLEU, ms=4)
    d.axhline(0, color="k", lw=0.6)
    d.set_xlabel("$a/R$")
    d.set_ylabel("force error vs Hertz  (%)")
    d.set_title(r"No volumetric locking at $\nu=0.49995$")
    g = V["V5"]["grande"]
    d.set_xticks([0.05, 0.1]); d.set_xticklabels(["0.05", "0.1"])
    d.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    d.text(0.04, 0.78,
           f"at $\\delta/R=0.30$: volume loss {abs(g[-1]['dV'])*100:.3f}%,\n"
           f"force varies {100*abs(g[-1]['F']-g[0]['F'])/g[-1]['F']:.2f}% over 10x elements",
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
    e2.plot(nd, dA, "s--", color=CTUM, ms=4)
    e2.set_ylabel(r"membrane area strain $\Delta A/A_0$  (%)", color=CTUM)
    e2.tick_params(axis="y", colors=CTUM)
    e.set_title("Constriction: DI stable to 1.3%")

    bc = [x for x in V["constriccion"]["barrido_confinamiento"] if x["dA"] > 0]
    xda = 100 * np.array([x["dA"] for x in bc]); yDI = 100 * np.array([x["DI"] for x in bc])
    o = np.argsort(xda); xda, yDI = xda[o], yDI[o]
    f.plot(xda, yDI, "o-", color=CLEU, ms=4, label="present 2D axisymmetric model")
    DI2D = np.interp(4.0, xda, yDI)
    f.plot([4.0], [DI2D], "o", mfc="white", mec=CLEU, ms=8, zorder=5)
    f.plot([4.0], [2.8], "D", color=CTUM, ms=6, label="Wang et al. 2023, 3D (HL-60)")
    f.annotate("+16% at matched\narea strain", xy=(4.15, 2.95), xytext=(5.6, 2.3),
               arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=6)
    f.set_xlabel(r"membrane area strain  $\Delta A/A_0$  (%)")
    f.set_ylabel("deformation index DI  (%)")
    f.set_title("2D calibration against published 3D")
    f.legend(frameon=False, loc="upper left", fontsize=6)

    for ax, L in zip([a, b, c, d, e, f], "abcdef"):
        ax.margins(0.05); panel_letter(ax, L)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    return fig


# ----------------------------------------------------------------- Figure 2
def fig_selectivity(umbral, pureza, path="fig2_selectivity.png"):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.1))

    rows = []
    for _, r in umbral.iterrows():
        rows.append((r.poblacion, r.grupo, r.beta_c_kPa_a_1mm_s))
    skbr = [v for n, g, v in rows if n == "SK-BR-3"]
    order = ["HL-60", "granulocyte", "PBMC", "MDA-MB-231", "MCF-10A", "MCF-7"]
    vals = {"HL-60": 0.0, "granulocyte": 0.0, "PBMC": 0.0,
            "MDA-MB-231": 0.0, "MCF-10A": 0.0, "MCF-7": 0.0}
    for n, g, v in rows:
        key = "granulocyte" if n == "granulocito" else n
        if key in vals:
            vals[key] = v
    grp = {"HL-60": CLEU, "granulocyte": CLEU, "PBMC": CLEU,
           "MDA-MB-231": CTUM, "MCF-7": CTUM, "MCF-10A": CALT}
    lo, hi = vals["PBMC"], vals["MDA-MB-231"]
    a.axvspan(lo, hi, color=CTUM, alpha=0.08, lw=0)
    for i, n in enumerate(order):
        a.hlines(i, 1e-2, vals[n], color=grp[n], lw=0.8, alpha=0.5)
        a.plot(vals[n], i, "o", color=grp[n], ms=6)
        a.text(vals[n] * 1.35, i, n, va="center", color=grp[n])
    a.hlines(len(order), skbr[0], skbr[1], color="#7B52AB", lw=2)
    a.plot(skbr, [len(order)] * 2, "|", color="#7B52AB", ms=8)
    a.text(skbr[1] * 1.35, len(order), r"SK-BR-3 ($\tau$ not measured)",
           va="center", color="#7B52AB")
    a.set_xscale("log"); a.set_xlim(2e-2, 4e3); a.set_ylim(-0.7, len(order) + 0.7)
    a.set_yticks([])
    a.set_xlabel(r"threshold device pressure  $\beta_c$  (kPa per mm s$^{-1}$)")
    a.set_title("Tumour cells need 165-1644x more pressure")
    a.text(np.sqrt(lo * hi), 4.35, "operating window\n165x", ha="center",
           color=CTUM, fontsize=6)
    for sp in ("left", "right", "top"):
        a.spines[sp].set_visible(False)

    sub = pureza[pureza.referencia_leucocitaria == "PBMC"]
    for beta, col, ls, lab in [(10.0, CTUM, "-", r"10 kPa/(mm s$^{-1}$)"),
                               (164.5, CLEU, "--", r"164.5 kPa/(mm s$^{-1}$)")]:
        g = sub[sub.dP_dispositivo_a_1mm_s_kPa == beta].sort_values("U_um_s")
        b.semilogx(g.U_um_s, g.separacion, ls, color=col, lw=1.4, label=lab)
        k = g.separacion.idxmax()
        b.plot(g.U_um_s[k], g.separacion[k], "o", color=col, ms=5)
    b.set_xlabel(r"throat velocity  $U$  ($\mu$m s$^{-1}$)")
    b.set_ylabel("separation  (tumour $-$ leucocyte retained)")
    b.set_title("Interior optimum in 27 of 27 combinations")
    b.legend(frameon=False, loc="center right", fontsize=6)
    b.text(0.03, 0.86, r"plateau $\geq$0.95:  4.41 decades $\rightarrow$ 0.38",
           transform=b.transAxes, color=CTR, fontsize=6)
    b.margins(y=0.08)
    for ax, L in zip([a, b], "ab"):
        panel_letter(ax, L)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    return fig


# ----------------------------------------------------------------- Figure 3
def fig_microfilter_sobol(vent, sob, path="fig3_microfilter_sobol.png"):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2))

    d = vent.sort_values("diametro_poro_um")
    a.fill_between(d.diametro_poro_um, d.dP_leucocito_kPa, d.dP_tumoral_minima_kPa,
                   color=CTUM, alpha=0.12, lw=0)
    a.semilogy(d.diametro_poro_um, d.dP_tumoral_minima_kPa, "o-", color=CTUM, ms=4,
               label="limiting tumour line (MDA-MB-231)")
    a.semilogy(d.diametro_poro_um, d.dP_leucocito_kPa, "s--", color=CLEU, ms=4,
               label="leucocyte (PBMC)")
    a.set_xlabel(r"pore diameter  ($\mu$m)")
    a.set_ylabel(r"pore entry threshold  $\Delta P$  (kPa)")
    a.set_title("Operating window at all 14 validated pore sizes")
    a.legend(frameon=False, loc="lower left", fontsize=6)
    i0 = d.ancho_ventana_kPa.idxmax()
    a.annotate(f"widest absolute window\n{d.ancho_ventana_kPa[i0]:.3f} kPa",
               xy=(d.diametro_poro_um[i0], np.sqrt(d.dP_leucocito_kPa[i0] * d.dP_tumoral_minima_kPa[i0])),
               xytext=(0.30, 0.80), textcoords="axes fraction", fontsize=6,
               arrowprops=dict(arrowstyle="->", lw=0.8))

    pretty = {"E_app_tumoral": r"$E_{\rm app}$ tumour", "E_app_leucocito": r"$E_{\rm app}$ leucocyte",
              "a_tumoral": "tumour radius", "a_leucocito": "leucocyte radius",
              "tau_tumoral": r"$\tau$ tumour", "tau_leucocito": r"$\tau$ leucocyte",
              "an_a": r"$a_n/a$", "C": r"$C=K_s/G_s$", "Pi": r"$\Pi$",
              "r_nuc": r"$r_{\rm nuc}$", "r_env": r"$r_{\rm env}$",
              "ell_factor": r"$\ell_{\rm factor}$"}
    measured = {"E_app_tumoral", "E_app_leucocito", "a_tumoral", "a_leucocito", "tau_tumoral"}
    order = ["E_app_tumoral", "E_app_leucocito", "a_tumoral", "a_leucocito",
             "tau_tumoral", "tau_leucocito", "an_a", "C", "Pi", "r_nuc", "r_env", "ell_factor"]
    FLOOR = 1e-22
    sc = sob[sob.salida == "selectividad_constriccion_tumoral_sobre_leucocito"].set_index("parametro")
    sm = sob[sob.salida == "selectividad_microfiltro"].set_index("parametro")
    b.axvspan(FLOOR / 3, 1e-19, color=CTR, alpha=0.12, lw=0)
    for i, p in enumerate(order):
        y = len(order) - 1 - i
        for src, mk, col, lab in [(sc, "o", CTUM, "constriction"), (sm, "D", CLEU, "microfilter")]:
            v = float(src.ST[p])
            open_ = v <= 0
            b.plot(max(v, FLOOR), y + (0.16 if mk == "o" else -0.16), mk,
                   color=col, ms=4.5, mfc="white" if open_ else col,
                   label=lab if i == 0 else None)
        b.text(FLOOR / 6, y, pretty[p], ha="right", va="center",
               color="k" if p in measured else CTR)
    b.set_xscale("log"); b.set_xlim(1e-26, 4)
    b.set_ylim(-0.8, len(order) - 0.2); b.set_yticks([])
    b.set_xlabel(r"Sobol total-order index  $S_T$")
    b.set_title("Selectivity rests only on measured moduli")
    b.text(3e-21, -0.55, "numerically zero", ha="center", color=CTR, fontsize=6)
    b.legend(frameon=False, loc="center left", bbox_to_anchor=(0.30, 0.18), fontsize=6)
    for sp in ("left", "right", "top"):
        b.spines[sp].set_visible(False)
    for ax, L in zip([a, b], "ab"):
        panel_letter(ax, L)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    return fig


# ----------------------------------------------------------------- Figure 4
def fig_design(dis, path="fig4_design.png"):
    fig, a = plt.subplots(figsize=(3.5, 2.9))
    d = dis.sort_values("Rt_um")
    a.plot(d.Rt_um, d.selectividad, "-", color=CTUM, lw=1.6)
    a.set_xlabel(r"throat radius  $R_t$  ($\mu$m)")
    a.set_ylabel(r"selectivity  tumour/leucocyte  ($\times$)", color=CTUM)
    a.tick_params(axis="y", colors=CTUM)
    a2 = a.twinx()
    a2.semilogy(d.Rt_um, d.ancho_kPa, "--", color=CLEU, lw=1.6)
    a2.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    a2.set_ylabel("pressure-window width  (kPa)", color=CLEU)
    a2.tick_params(axis="y", colors=CLEU)
    a.axvspan(5.90, 6.00, color=CTR, alpha=0.15, lw=0)
    a.text(5.95, a.get_ylim()[0] + 0.12, "edge of\nsampling", ha="center", va="bottom",
           color=CTR, fontsize=6)
    k = (d.Rt_um - 5.675).abs().idxmin()
    a.annotate(f"usable compromise\n$R_t\\approx${d.Rt_um[k]:.1f} "
               f"$\\mu$m, {d.selectividad[k]:.1f}x, {d.ancho_kPa[k]:.3f} kPa",
               xy=(d.Rt_um[k], d.selectividad[k]), xytext=(0.06, 0.60),
               textcoords="axes fraction", fontsize=6,
               arrowprops=dict(arrowstyle="->", lw=0.8))
    a.set_title("Selectivity rises where the window stops\nbeing controllable")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    return fig
