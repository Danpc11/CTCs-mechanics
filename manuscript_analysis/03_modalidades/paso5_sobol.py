"""Paso 5-B: sustituto + analisis de sensibilidad global de Sobol.

Que conclusiones sobreviven a la incertidumbre de los parametros TRANSFERIDOS
(todos los nucleares y de membrana) y cuales no.
"""
import sys, os, json, glob, itertools
import numpy as np
import pandas as pd
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L

COMP = ["cortex", "citoplasma", "nucleo", "envoltura"]
CONF_DIS = 1.70          # confinamiento nominal de diseno
ZT, LTAP = 4.0, 5.0
A_REF = 8.74


# ---------------------------------------------------------------- sustituto
def _z(X):
    """(c, an_a, C, Pi, r_nuc, r_env) -> coordenadas del ajuste."""
    c, an_a, C, Pi, r_nuc, r_env = [np.asarray(v, dtype=float) for v in X]
    return np.vstack([c, an_a, np.log(C), np.log(Pi), np.log(r_nuc),
                      np.log(r_env)])


def _base(Z):
    """Base polinomica: grado 2 completo en las 6 variables, mas grados 3 y 4
    en el confinamiento c y los acoplos c^2 * z_j.  f recorre tres ordenes de
    magnitud con c, asi que el grado 2 no basta (error mediano 22%)."""
    n = Z.shape[1]
    cols = [np.ones(n)] + [Z[i] for i in range(6)]
    for i in range(6):
        for j in range(i, 6):
            cols.append(Z[i] * Z[j])
    cols += [Z[0]**3, Z[0]**4, Z[0]**5]
    cols += [Z[0]**2 * Z[j] for j in range(1, 6)]
    cols += [Z[0]**3 * Z[j] for j in range(1, 6)]
    return np.vstack(cols).T


def ajustar(df, col):
    Z = _z([df.c, df.an_a, df.C, df.Pi, df.r_nuc, df.r_env])
    A = _base(Z)
    y = np.log(np.maximum(df[col].values, 1e-12))
    rng = np.random.default_rng(3)
    idx = rng.permutation(len(y))
    ntr = int(0.8 * len(y))
    tr, te = idx[:ntr], idx[ntr:]
    coef, *_ = np.linalg.lstsq(A[tr], y[tr], rcond=None)
    pred_te = A[te] @ coef
    err = np.exp(np.abs(pred_te - y[te])) - 1.0
    coef_all, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef_all, dict(n=len(y), err_mediano_holdout=float(np.median(err)),
                          err_p90_holdout=float(np.percentile(err, 90)),
                          err_max_holdout=float(err.max()),
                          R2=float(1 - np.var(A @ coef_all - y) / np.var(y)))


def f_sustituto(coef, c, an_a, C, Pi, r_nuc, r_env):
    Z = _z([np.atleast_1d(c), np.atleast_1d(an_a), np.atleast_1d(C),
            np.atleast_1d(Pi), np.atleast_1d(r_nuc), np.atleast_1d(r_env)])
    return np.exp(_base(Z) @ coef)


def dP_critica_sustituto(coef, a, Gs, conf_nom, an_a, C, Pi, r_nuc, r_env,
                         zt=ZT, Ltap=LTAP, nz=1501):
    """Presion critica (kPa) a partir del sustituto, con la MISMA definicion
    que presion_critica(): maximo de |dW/ds|/(pi R_w^2) dentro de la banda
    valida de la envolvente.  Vectorizado sobre las muestras de Sobol."""
    a = np.atleast_1d(a) * 1.0
    Rwide = 1.6 * a
    Rt = a / conf_nom
    Rmin_val = L.perfil_Rw(0.0, a / L.CONF_MAX, Rwide, zt, Ltap)
    Rmax_val = L.perfil_Rw(0.0, a / L.CONF_MIN, Rwide, zt, Ltap)
    z = np.linspace(0.0, zt + 6.0 * Ltap, nz)[1:, None]          # z = 0 excluded (sign(0) = 0)
    Rw = L.perfil_Rw(z, Rt[None, :], Rwide[None, :], zt, Ltap)
    dR = np.abs(L.dRw_dz(z, Rt[None, :], Rwide[None, :], zt, Ltap))
    # f was fitted against the NOMINAL confinement of each state; the state whose neck equals the local wall
    # radius R_w has nominal throat Rt_eq = (R_w - g0 R_wide)/(1 - g0)
    g0 = 0.5 * (1.0 + np.tanh(-zt / Ltap))
    Rt_eq = np.maximum((Rw - g0 * Rwide[None, :]) / (1.0 - g0), 1e-9)
    cc = a[None, :] / Rt_eq
    dc_dR = a[None, :] / (Rt_eq ** 2 * (1.0 - g0))
    # dW/dR = Gs a^2 * df/dc * dc/dR = -Gs a^3 f'(c) / R^2
    eps = 1e-4
    # cc tiene forma (nz, n) y .ravel() en orden C recorre z-mayor: para cada
    # indice de z, TODAS las muestras.  Los parametros deben replicarse con
    # np.tile (no np.repeat, que los desalinearia con el camino).
    rep = [np.tile(v, nz) for v in (an_a, C, Pi, r_nuc, r_env)]
    fp = (f_sustituto(coef, (cc + eps).ravel(), *rep)
          - f_sustituto(coef, (cc - eps).ravel(), *rep)) / (2 * eps)
    fp = fp.reshape(cc.shape)
    dWdR = Gs[None, :] * a[None, :]**2 * np.abs(fp) * dc_dR
    q = dWdR * dR / (np.pi * Rw**2)
    val = (Rw >= Rmin_val[None, :]) & (Rw <= Rmax_val[None, :])
    q = np.where(val, q, -np.inf)
    return q.max(axis=0)


# ---------------------------------------------------------------- principal
if __name__ == "__main__":
    regs = []
    for f in sorted(glob.glob(os.path.join(HERE, "paso5b_diseno_s*.json"))):
        regs += json.load(open(f))["registros"]
    d = pd.DataFrame(regs)
    print("diseno: %d puntos, %d convergen (%.0f%%)"
          % (len(d), int(d.converge.sum()), 100 * d.converge.mean()))
    d = d[d.converge & (d.f > 0)].copy()

    coef, diag = ajustar(d, "f")
    print("sustituto de f = W/(Gs a^2): R2=%.5f  error mediano hold-out=%.3f%%"
          "  p90=%.3f%%  max=%.3f%%"
          % (diag["R2"], 100 * diag["err_mediano_holdout"],
             100 * diag["err_p90_holdout"], 100 * diag["err_max_holdout"]))
    coefs_c, diags_c = {}, {}
    for c_ in COMP:
        col = "f_" + c_
        if col in d and (d[col] > 0).sum() > 60:
            coefs_c[c_], diags_c[c_] = ajustar(d[d[col] > 0], col)
            print("   %-11s R2=%.4f  err mediano=%.2f%%"
                  % (c_, diags_c[c_]["R2"],
                     100 * diags_c[c_]["err_mediano_holdout"]))

    # ---- problema de Sobol
    RG = L.RANGOS if hasattr(L, "RANGOS") else None
    nombres = ["an_a", "C", "Pi", "r_nuc", "r_env",
               "E_app_tumoral", "E_app_leucocito", "a_tumoral", "a_leucocito",
               "tau_tumoral", "tau_leucocito", "ell_factor"]
    limites = [[0.40, 0.60],        # an/a            TRANSFERIDO
               [2.40, 37.0],        # C de Skalak: Ks = (1+2C) Gs ; extremos no emparejados de PC-3 (Lu et al. 2025)  TRANSFERIDO
               [0.03, 1.00],        # Pi              REGULARIZACION (Fase 1)
               [1.40, 4.00],        # G_nuc/G_cyt     TRANSFERIDO
               [2.00, 8.00],        # Gs_n/Gs         TRANSFERIDO
               [1.527, 2.772],      # E_app tumoral   MEDIDO (rango del panel)
               [0.400, 0.660],      # E_app leucocito MEDIDO
               [8.51, 10.00],       # a tumoral       MEDIDO/TRANSFERIDO
               [8.50, 8.50 * 1.0001],   # a leucocito TRANSFERIDO (sin rango)
               [0.37, 3.76],        # tau tumoral  MEDIDO (horquilla del panel
                                    #   MDA-MB-231 a MCF-7, +/-1 SD)
               [8.575e-3, 3.23],    # tau LEUCOCITO: el parametro critico SIN
                                    #   medir.  Minimo = anclaje de celula
                                    #   entera de HL-60 (eta/E = 3.43/400);
                                    #   maximo = el mayor tau de mama medido,
                                    #   para dejar que ROMPA la conclusion.
               [0.50, 2.00]]        # ell de deglucion  DECLARADO
    problema = {"num_vars": len(nombres), "names": nombres, "bounds": limites}
    Xs = sobol_sample.sample(problema, 2048, calc_second_order=False)
    print("muestras de Sobol: %d" % len(Xs))

    an_a, C, Pi, r_nuc, r_env = Xs[:, 0], Xs[:, 1], Xs[:, 2], Xs[:, 3], Xs[:, 4]
    Et, El = Xs[:, 5], Xs[:, 6]
    at, al = Xs[:, 7], Xs[:, 8]
    tT, tL = Xs[:, 9], Xs[:, 10]
    ell = Xs[:, 11]
    Gst = L.Gs_de_Eapp(Et, at)
    Gsl = L.Gs_de_Eapp(El, al)

    dP_t = dP_critica_sustituto(coef, at, Gst, CONF_DIS, an_a, C, Pi, r_nuc, r_env)
    dP_l = dP_critica_sustituto(coef, al, Gsl, CONF_DIS, an_a, C, Pi, r_nuc, r_env)
    selectividad = dP_t / dP_l

    # beta_c = dP_el * tau / L_carga   (umbral de presion del dispositivo por
    # encima del cual la poblacion escapa a ALGUN caudal).  El EJE PRINCIPAL es
    # tumoral frente a leucocito: la ventana de operacion existe si y solo si
    # beta_c(leucocito) < beta_c(tumoral).
    Lc_t = np.array([L.longitud_de_carga(ai, ai / CONF_DIS, 1.6 * ai, ZT, LTAP)
                     for ai in at])
    Lc_l = np.array([L.longitud_de_carga(ai, ai / CONF_DIS, 1.6 * ai, ZT, LTAP)
                     for ai in al])
    beta_t = dP_t * tT / Lc_t
    beta_l = dP_l * tL / Lc_l
    razon_beta = beta_t / beta_l

    # compartimento limitante: fraccion marginal del cortex
    frac_cortex = None
    if set(COMP) <= set(coefs_c):
        marg = {}
        for c_ in COMP:
            eps = 1e-4
            cu = np.full_like(an_a, CONF_DIS)
            fp = (f_sustituto(coefs_c[c_], cu + eps, an_a, C, Pi, r_nuc, r_env)
                  - f_sustituto(coefs_c[c_], cu - eps, an_a, C, Pi,
                                r_nuc, r_env)) / (2 * eps)
            marg[c_] = np.abs(fp)
        tot = sum(marg.values())
        frac_cortex = marg["cortex"] / tot
        frac_nuc = (marg["nucleo"] + marg["envoltura"]) / tot

    # umbral de poro del microfiltro
    Rp_t, Rp_l = at / CONF_DIS, al / CONF_DIS
    cu2 = np.full_like(an_a, CONF_DIS)
    Wt = Gst * at**2 * f_sustituto(coef, cu2, an_a, C, Pi, r_nuc, r_env)
    Wl = Gsl * al**2 * f_sustituto(coef, cu2, an_a, C, Pi, r_nuc, r_env)
    dPp_t = Wt / (ell * at * np.pi * Rp_t**2)
    dPp_l = Wl / (ell * al * np.pi * Rp_l**2)
    sel_poro = dPp_t / dPp_l

    salidas = {"dP_critica_tumoral_kPa": dP_t,
               "selectividad_constriccion_tumoral_sobre_leucocito": selectividad,
               "razon_beta_c_tumoral_sobre_leucocito": razon_beta,
               "dP_umbral_poro_tumoral_kPa": dPp_t,
               "selectividad_microfiltro": sel_poro}
    if frac_cortex is not None:
        salidas["fraccion_marginal_cortex"] = frac_cortex
        salidas["fraccion_marginal_nuclear"] = frac_nuc

    filas = []
    for k, y in salidas.items():
        Si = sobol_analyze.analyze(problema, np.asarray(y, dtype=float),
                                   calc_second_order=False, print_to_console=False)
        for i, nm in enumerate(nombres):
            filas.append(dict(salida=k, parametro=nm,
                              S1=float(Si["S1"][i]), S1_conf=float(Si["S1_conf"][i]),
                              ST=float(Si["ST"][i]), ST_conf=float(Si["ST_conf"][i]),
                              media=float(np.mean(y)),
                              p5=float(np.percentile(y, 5)),
                              p50=float(np.percentile(y, 50)),
                              p95=float(np.percentile(y, 95))))
    sob = pd.DataFrame(filas)
    sob.to_csv(os.path.join(HERE, "sensibilidad_sobol.csv"), index=False)

    print("\n== indices totales ST (los 5 mayores por salida) ==")
    for k in salidas:
        s = sob[sob.salida == k].nlargest(5, "ST")
        print("  %s   [p5=%.4g p50=%.4g p95=%.4g]"
              % (k, s.p5.iloc[0], s.p50.iloc[0], s.p95.iloc[0]))
        for _, r in s.iterrows():
            print("      %-18s ST=%.3f  S1=%.3f" % (r.parametro, r.ST, r.S1))

    # ---- que conclusiones sobreviven
    concl = {
        "C1 la celula tumoral exige mas presion que el leucocito en constriccion":
            float(np.mean(selectividad > 1.0)),
        "C1b selectividad de constriccion > 2x":
            float(np.mean(selectividad > 2.0)),
        "C2 el cortex es el compartimento limitante (>50% del margen)":
            float(np.mean(frac_cortex > 0.5)) if frac_cortex is not None else None,
        "C2b los compartimentos nucleares aportan <10% del margen":
            float(np.mean(frac_nuc < 0.10)) if frac_cortex is not None else None,
        "C3 existe ventana de operacion: beta_c(tumoral) > beta_c(leucocito)":
            float(np.mean(razon_beta > 1.0)),
        "C3b la ventana de presion de dispositivo supera 10x":
            float(np.mean(razon_beta > 10.0)),
        "C3c la ventana supera 100x":
            float(np.mean(razon_beta > 100.0)),
        "C4 el microfiltro retiene tumoral y deja pasar leucocito":
            float(np.mean(sel_poro > 1.0)),
        "C4b ventana de microfiltro > 2x":
            float(np.mean(sel_poro > 2.0)),
    }
    print("\n== fraccion de muestras en que se cumple cada conclusion ==")
    for k, v in concl.items():
        print("  %-62s %s" % (k, "n/a" if v is None else "%.4f" % v))

    # cota sobre tau del leucocito: hasta que valor sobrevive la ventana
    cota = []
    for thr in (0.5, 0.9, 0.99):
        lo, hi = 8.575e-3, 3.23
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            frac = float(np.mean(razon_beta[tL <= mid] > 1.0)) \
                if (tL <= mid).sum() > 30 else 1.0
            if frac >= thr:
                lo = mid
            else:
                hi = mid
        cota.append(dict(fraccion_exigida=thr, tau_leucocito_max_s=lo))
    print("\n== cota sobre tau del leucocito (Sobol) ==")
    for c_ in cota:
        print("   la ventana se cumple en >=%.0f%% de las muestras con "
              "tau_leu <= %.4g s" % (100 * c_["fraccion_exigida"],
                                     c_["tau_leucocito_max_s"]))
    concl["C3d cota de tau leucocitario (s) al 90% de las muestras"] = \
        [c_["tau_leucocito_max_s"] for c_ in cota if c_["fraccion_exigida"] == 0.9][0]

    json.dump({"cota_tau_leucocito": cota,
               "diagnostico_sustituto": diag,
               "diagnostico_por_compartimento": diags_c,
               "problema": problema, "n_muestras": int(len(Xs)),
               "conclusiones": concl,
               "percentiles": {k: dict(p5=float(np.percentile(v, 5)),
                                       p50=float(np.percentile(v, 50)),
                                       p95=float(np.percentile(v, 95)))
                               for k, v in salidas.items()}},
              open(os.path.join(HERE, "paso5_sobol.json"), "w"), indent=1)
