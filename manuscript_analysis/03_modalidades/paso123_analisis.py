"""Pasos 1, 2 y 3: post-proceso del paisaje energetico FEM.

  paso 1  presion critica de paso y tiempo de transito en constriccion
  paso 2  barrido de caudal y numero tipo Deborah  (nucleo del articulo)
  paso 3  umbral de paso por poro de microfiltro y ventana de operacion
"""
import sys, os, json, glob, itertools
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L

COMPARTIMENTOS = ["cortex", "citoplasma", "nucleo", "envoltura"]
MU_F = 1e-6            # viscosidad del medio 1e-3 Pa s  -> nN s/um^2
BETA_NOM = 0.01        # dP_motriz = beta * U ; 10 kPa a U = 1 mm/s = 1e3 um/s
DP_GRID = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]        # kPa


def cargar():
    regs = []
    for f in sorted(glob.glob(os.path.join(HERE, "paso1_fem_s*.json"))):
        regs += json.load(open(f))["registros"]
    df = pd.DataFrame(regs)
    return df


def splines(df):
    """W_i(R) por (estirpe, geometria) con PCHIP sobre los puntos CONVERGIDOS,
    anclado en W=0 a R=a (sin confinamiento).  Monotono por construccion."""
    S = {}
    for (nm, geo), g in df[df.converge].groupby(["estirpe", "geometria"]):
        g = g.sort_values("Rt")
        a = float(g.a.iloc[0])
        R = np.append(g.Rt.values, a * 1.02)
        out = {}
        for c in COMPARTIMENTOS + ["total"]:
            col = "W_" + c if c != "total" else "W_total"
            if col not in g:
                continue
            W = np.append(g[col].values, 0.0)
            o = np.argsort(R)
            out[c] = PchipInterpolator(R[o], np.maximum(W[o], 0.0),
                                       extrapolate=False)
        S[(nm, geo)] = (out, a)
    return S


# ---------------------------------------------------------------------
if __name__ == "__main__":
    df = cargar()
    print("registros FEM: %d ; convergidos: %d ; fuera de envolvente: %d"
          % (len(df), int(df.converge.sum()), int((~df.dentro_envolvente).sum())))
    no_conv = df[~df.converge][["estirpe", "geometria", "conf_nominal",
                                "lam_alcanzado"]]
    print("no convergidos:\n", no_conv.to_string(index=False) if len(no_conv)
          else "  ninguno")
    S = splines(df)

    # =================================================================
    # PASO 1 -- presion critica y tiempo de transito
    # =================================================================
    filas = []
    for (nm, geo), (sp, a) in S.items():
        zt = float(df[(df.estirpe == nm) & (df.geometria == geo)].zt.iloc[0])
        Ltap0 = float(df[(df.estirpe == nm) & (df.geometria == geo)].Ltap.iloc[0])
        Rwide = 1.6 * a
        est = L.ESTIRPES[nm]
        for conf_d in [1.20, 1.42, 1.70, 1.90, 2.10, 2.43]:
            Rt = a / conf_d
            theta = L.angulo_de_Ltap(Ltap0, Rt, Rwide)
            pc = L.presion_critica(a, Rt, Rwide, zt, Ltap0,
                                   sp["total"], sp["total"].derivative())
            if not np.isfinite(pc["dP_crit"]):
                continue
            # compartimento limitante: mayor |dW_i/dR| en el punto critico
            marg = {}
            for c in COMPARTIMENTOS:
                if c in sp:
                    d = sp[c].derivative()(pc["R_crit"])
                    marg[c] = abs(float(d)) if np.isfinite(d) else 0.0
            tot = sum(marg.values()) or 1.0
            lim = max(marg, key=marg.get)
            fila = dict(estirpe=nm, geometria=geo, zt_um=zt, longitud_garganta_um=2 * zt,
                        Ltap_um=Ltap0, theta_entrada_deg=round(theta, 2),
                        a_um=a, Rt_um=round(Rt, 4), confinamiento=conf_d,
                        Gs_mN_m=est.Gs if est.Gs else None,
                        dP_critica_kPa=pc["dP_crit"],
                        L_entrada_um=pc["L_entrada"],
                        conf_en_punto_critico=pc["conf_crit"],
                        compartimento_limitante=lim,
                        tau_s=est.tau, tau_disponible=est.tau is not None)
            for c in COMPARTIMENTOS:
                fila["frac_marginal_" + c] = marg.get(c, 0.0) / tot
            for dP in DP_GRID:
                if est.tau is None:
                    fila["t_transito_s_dP%g" % dP] = None
                    continue
                t = L.tiempo_transito(a, Rt, Rwide, zt, Ltap0, sp["total"],
                                      sp["total"].derivative(), dP, est.tau)
                fila["t_transito_s_dP%g" % dP] = t
            filas.append(fila)
    c1 = pd.DataFrame(filas)
    c1.to_csv(os.path.join(HERE, "barrido_constriccion.csv"), index=False)
    print("\n== PASO 1: presion critica (kPa), geometria de referencia ==")
    piv = c1[c1.geometria == "referencia"].pivot_table(
        index="confinamiento", columns="estirpe", values="dP_critica_kPa")
    print(piv.round(3).to_string())
    print("\ncompartimento limitante (geometria de referencia):")
    print(c1[c1.geometria == "referencia"].pivot_table(
        index="confinamiento", columns="estirpe",
        values="compartimento_limitante", aggfunc="first").to_string())
    # el compartimento limitante solo tiene sentido donde la barrera es
    # apreciable: a confinamiento nominal 1.20 el EFECTIVO es 1.04 (apenas hay
    # contacto) y dP_crit ~ 0, con lo que el reparto es ruido numerico.
    sig = c1[c1.dP_critica_kPa > 1e-3]
    print("\nfraccion marginal media por compartimento (dP_crit > 1e-3 kPa, "
          "%d de %d filas):" % (len(sig), len(c1)))
    print(sig[["frac_marginal_" + c for c in COMPARTIMENTOS]].mean().round(4).to_string())
    print("compartimento limitante, recuento:")
    print(sig.compartimento_limitante.value_counts().to_string())

    # =================================================================
    # PASO 2 -- barrido de caudal / numero de Deborah
    # =================================================================
    GEO2, CONF2 = "referencia", 1.70
    U = np.logspace(-2.0, 5.5, 240)       # um/s  (1e-5 - 316 mm/s)
    lineas_tau = [k for k in S if k[1] == GEO2 and
                  L.ESTIRPES[k[0]].tau is not None]
    base = {}
    for (nm, geo) in sorted(lineas_tau):
        sp, a = S[(nm, geo)]
        g = df[(df.estirpe == nm) & (df.geometria == geo)]
        zt, Lt = float(g.zt.iloc[0]), float(g.Ltap.iloc[0])
        Rt = a / CONF2
        pc = L.presion_critica(a, Rt, 1.6 * a, zt, Lt, sp["total"],
                               sp["total"].derivative())
        base[nm] = dict(dP_el=pc["dP_crit"], L_carga=pc["L_entrada"], a=a,
                        Rt=Rt, tau=L.ESTIRPES[nm].tau,
                        beta_c=pc["dP_crit"] * L.ESTIRPES[nm].tau
                        / pc["L_entrada"])
    print("\n== PASO 2: base a confinamiento nominal %.2f (efectivo %.3f) =="
          % (CONF2, float(df[(df.geometria == GEO2) &
                            (df.conf_nominal == CONF2)].conf_efectivo.iloc[0])))
    print("  beta_c = dP_el*tau/L_carga es el umbral de presion por unidad de")
    print("  velocidad por encima del cual la estirpe puede pasar a ALGUN caudal.")
    for nm, b in base.items():
        print("  %-11s dP_el=%8.4f kPa  L_carga=%5.2f um  tau=%4.2f s  "
              "beta_c=%8.4f kPa s/um  (= %7.1f kPa a 1 mm/s)"
              % (nm, b["dP_el"], b["L_carga"], b["tau"], b["beta_c"],
                 b["beta_c"] * 1000.0))

    # --- barrido de beta: el regimen de operacion del dispositivo
    bmin = min(b["beta_c"] for b in base.values())
    bmax = max(b["beta_c"] for b in base.values())
    BETAS = np.unique(np.concatenate([
        np.logspace(np.log10(bmin) - 1.2, np.log10(bmax) + 1.2, 25),
        [BETA_NOM]]))

    rng = np.random.default_rng(20251014)
    NMC = 4000
    pobl = {}
    for nm, b in base.items():
        e = L.ESTIRPES[nm]
        sa = (e.a_rng[1] - e.a_rng[0]) / 2.0 if e.a_rng[0] else 0.0
        st = (e.tau_rng[1] - e.tau_rng[0]) / 2.0
        a_s = rng.normal(e.a, max(sa, 1e-9), NMC).clip(0.5 * e.a, 2 * e.a)
        t_s = rng.normal(e.tau, st, NMC).clip(0.05, 20.0)
        # dP_crit propto Gs/a = E_app/K  ->  NO depende de a a confinamiento
        # fijo; la dispersion de E_app no esta publicada para las lineas
        # tumorales, asi que el ancho de estas curvas es una SUBESTIMACION.
        dP_s = np.full(NMC, b["dP_el"])
        L_s = b["L_carga"] * (a_s / e.a)
        pobl[nm] = (dP_s, t_s, L_s)

    fil2, sep_rows, resumen2 = [], [], {}
    pares = list(itertools.combinations(sorted(base), 2))
    for beta in BETAS:
        cap = {}
        for nm, b in base.items():
            dP_s, t_s, L_s = pobl[nm]
            fc = []
            for u in U:
                m, _, _ = L.margen_caudal(u, dP_s, t_s, L_s, beta)
                fc.append(float(np.mean(m <= 1.0)))
            cap[nm] = np.array(fc)
            uc = L.U_critica(b["dP_el"], b["tau"], b["L_carga"], beta)
            for u, f in zip(U, cap[nm]):
                fil2.append(dict(
                    estirpe=nm, beta_kPa_s_um=beta,
                    dP_dispositivo_a_1mm_s_kPa=beta * 1000.0,
                    U_um_s=u, U_mm_s=u / 1000.0, confinamiento_nominal=CONF2,
                    geometria=GEO2, tau_s=b["tau"], dP_el_kPa=b["dP_el"],
                    L_carga_um=b["L_carga"],
                    De=L.deborah(b["tau"], u, b["L_carga"]),
                    dP_motriz_kPa=beta * u,
                    dP_crit_aparente_kPa=b["dP_el"] * (
                        1 + L.deborah(b["tau"], u, b["L_carga"])),
                    m_nominal=L.margen_caudal(u, b["dP_el"], b["tau"],
                                              b["L_carga"], beta)[0],
                    m_inf_nominal=beta * b["L_carga"] / (b["dP_el"] * b["tau"]),
                    U_critica_um_s=uc, fraccion_capturada=f))
        for p in pares:
            d = np.abs(cap[p[0]] - cap[p[1]])
            i = int(np.argmax(d))
            interior = bool(0 < i < len(d) - 1 and d[i] > d[-1] + 1e-3
                            and d[i] > d[0] + 1e-3)
            sep_rows.append(dict(
                par="%s vs %s" % p, beta_kPa_s_um=beta,
                dP_dispositivo_a_1mm_s_kPa=beta * 1000.0,
                separacion_max=float(d[i]), U_optima_um_s=float(U[i]),
                U_optima_mm_s=float(U[i]) / 1000.0, maximo_interior=interior,
                sep_en_U_minima=float(d[0]), sep_en_U_maxima=float(d[-1]),
                caida_tras_el_maximo=float(d[i] - d[-1])))
    c2 = pd.DataFrame(fil2)
    c2.to_csv(os.path.join(HERE, "barrido_caudal.csv"), index=False)
    sp2 = pd.DataFrame(sep_rows)
    sp2.to_csv(os.path.join(HERE, "separacion_caudal.csv"), index=False)

    print("\n== PASO 2: regimenes de operacion ==")
    print("  estirpes que pueden pasar a algun caudal, por presion de dispositivo:")
    for beta in BETAS[::4]:
        pasan = [nm for nm, b in base.items() if beta > b["beta_c"]]
        print("    %8.1f kPa a 1 mm/s : %s"
              % (beta * 1000.0, ", ".join(pasan) if pasan else "ninguna"))
    print("\n  velocidad critica U* por estirpe (um/s), por presion de dispositivo:")
    print("    %10s " % "kPa@1mm/s" + " ".join("%12s" % n for n in base))
    for beta in BETAS:
        us = [L.U_critica(b["dP_el"], b["tau"], b["L_carga"], beta)
              for b in base.values()]
        print("    %10.1f " % (beta * 1000) +
              " ".join("%12s" % ("inf" if not np.isfinite(u) else "%.4g" % u)
                       for u in us))
    print("\n  separacion maxima y su caudal optimo (par MCF-7 / MDA-MB-231):")
    sub = sp2[sp2.par == "MCF-7 vs MDA-MB-231"].sort_values("beta_kPa_s_um")
    print(sub[["dP_dispositivo_a_1mm_s_kPa", "separacion_max", "U_optima_mm_s",
               "maximo_interior", "sep_en_U_maxima"]].round(4).to_string(index=False))

    nm_any = bool(sp2.maximo_interior.any())
    resumen2 = dict(
        no_monotonia_en_algun_beta=nm_any,
        beta_c_por_estirpe={k: v["beta_c"] for k, v in base.items()},
        dP_dispositivo_umbral_kPa_a_1mm_s={
            k: v["beta_c"] * 1000.0 for k, v in base.items()},
        pares_con_maximo_interior=sorted(
            sp2[sp2.maximo_interior].par.unique().tolist()),
        separacion_max_global=float(sp2.separacion_max.max()),
        U_optima_en_el_maximo_global_mm_s=float(
            sp2.loc[sp2.separacion_max.idxmax(), "U_optima_mm_s"]),
        beta_en_el_maximo_global=float(
            sp2.loc[sp2.separacion_max.idxmax(), "beta_kPa_s_um"]))
    # ventana de beta en la que SOLO pasa la de tau mas corto
    orden = sorted(base, key=lambda k: base[k]["beta_c"])
    resumen2["ventana_beta_selectiva"] = dict(
        estirpe_que_pasa=orden[0],
        beta_inferior=base[orden[0]]["beta_c"],
        beta_superior=base[orden[1]]["beta_c"],
        factor=base[orden[1]]["beta_c"] / base[orden[0]]["beta_c"],
        kPa_a_1mm_s=[base[orden[0]]["beta_c"] * 1000.0,
                     base[orden[1]]["beta_c"] * 1000.0])
    print("\n  ventana de presion en la que SOLO pasa %s: %.1f - %.1f kPa a "
          "1 mm/s (factor %.2f)"
          % (orden[0], base[orden[0]]["beta_c"] * 1000.0,
             base[orden[1]]["beta_c"] * 1000.0,
             base[orden[1]]["beta_c"] / base[orden[0]]["beta_c"]))

    # =================================================================
    # PASO 3 -- microfiltro
    # =================================================================
    # Poro de microfiltro: geometria de garganta larga (la mas proxima a un
    # poro recto).  El canto vivo no define pendiente de pared, asi que la
    # barrera se reparte sobre una longitud de deglucion ell = ell_factor * a.
    # REJILLA COMUN de diametros de poro para que las estirpes sean comparables
    # (cada una tiene su propio radio, asi que el confinamiento difiere).
    GEO3 = "garganta_larga"
    DIAM_PORO = np.arange(7.0, 15.01, 0.5)
    fil3 = []
    for (nm, geo), (sp, a) in S.items():
        if geo != GEO3:
            continue
        est = L.ESTIRPES[nm]
        g = df[(df.estirpe == nm) & (df.geometria == geo) & df.converge]
        Rt_ok = np.sort(g.Rt.values)
        for Dp in DIAM_PORO:
            Rp = Dp / 2.0
            if Rp < Rt_ok.min() or Rp > a:
                continue              # fuera del rango FEM cubierto
            W = float(sp["total"](min(Rp, Rt_ok.max())))
            if not np.isfinite(W):
                continue
            conf = a / Rp
            for ell in (0.5, 1.0, 2.0):
                fil3.append(dict(
                    estirpe=nm, a_um=a, diametro_poro_um=float(Dp),
                    Rp_um=Rp, confinamiento_nominal=conf,
                    dentro_envolvente=bool(L.dentro_envolvente(conf)),
                    ell_factor=ell, W_nN_um=W, Gs_mN_m=est.Gs,
                    E_app_kPa=est.E_app,
                    dP_umbral_kPa=L.presion_poro(W, Rp, a, ell),
                    tau_disponible=est.tau is not None))
    c3 = pd.DataFrame(fil3)
    c3.to_csv(os.path.join(HERE, "barrido_microfiltro.csv"), index=False)
    print("\n== PASO 3: umbral de paso por poro (kPa, ell = a) ==")
    c3v = c3[(c3.ell_factor == 1.0) & c3.dentro_envolvente]
    p3 = c3v.pivot_table(index="diametro_poro_um", columns="estirpe",
                         values="dP_umbral_kPa")
    print(p3.round(4).to_string())
    print("\nventana de operacion (el leucocito pasa, la tumoral se retiene):")
    vent = []
    for d_, row in p3.iterrows():
        if "leucocito" not in row or not np.isfinite(row.get("leucocito",
                                                            np.nan)):
            continue
        tum = row[[c for c in row.index
                   if c not in ("leucocito", "MCF-10A")]].dropna()
        if not len(tum):
            continue
        lo, hi = float(row["leucocito"]), float(tum.min())
        vent.append(dict(diametro_poro_um=float(d_), dP_leucocito_kPa=lo,
                         dP_tumoral_minima_kPa=hi,
                         tumoral_limitante=str(tum.idxmin()),
                         ancho_ventana_kPa=hi - lo,
                         factor_selectividad=hi / lo if lo else np.nan,
                         ventana_existe=bool(hi > lo)))
        print("  poro %5.2f um: %8.4f - %8.4f kPa  (ancho %8.4f, factor %.2f)"
              "  limita %s" % (d_, lo, hi, hi - lo, hi / lo if lo else np.nan,
                               tum.idxmin()))
    pd.DataFrame(vent).to_csv(os.path.join(HERE, "ventana_microfiltro.csv"),
                              index=False)

    json.dump({"paso2_base": {k: {kk: float(vv) for kk, vv in v.items()}
                              for k, v in base.items()},
               "paso2_resumen": resumen2,
               "beta_nominal_kPa_s_um": BETA_NOM,
               "confinamiento_paso2": CONF2,
               "geometria_paso3": GEO3,
               "ventana": vent},
              open(os.path.join(HERE, "paso123.json"), "w"), indent=1, default=str)
