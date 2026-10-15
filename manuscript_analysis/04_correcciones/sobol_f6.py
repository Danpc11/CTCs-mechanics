"""Fase 6: Sobol con MUESTREO INDEPENDIENTE por poblacion.

El Sobol de la Fase 2 compartia an_a, C, Pi, r_nuc y r_env entre la poblacion
tumoral y la leucocitaria.  Con W = Gs a^2 f(c) y la MISMA f para las dos, esos
parametros se cancelan EXACTAMENTE en el cociente: los indices totales de 1e-20
y 1e-31 que se reportaron son la huella en coma flotante de esa cancelacion, no
una medida de insensibilidad.

Aqui cada poblacion recibe sus propios valores, con los rangos MEDIDOS:
  - Pi leucocito 0.162 (banda declarada +/-30%) frente a 0.032-0.078 tumorales
  - an_a leucocito 0.625-0.763 (neutrofilo-linfocito) frente a 0.40-0.60 tumoral
  - radio leucocito 3.03-4.40 um (linfocito-monocito) frente a 8.51-10.0 um
"""
import sys, os, json, glob
import numpy as np
import pandas as pd
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L
from paso5_sobol import ajustar, f_sustituto, dP_critica_sustituto, COMP, CONF_DIS, ZT, LTAP

N_BASE = 2048

NOMBRES = ["an_a_tum", "an_a_leu", "C_tum", "C_leu", "Pi_tum", "Pi_leu",
           "r_nuc_tum", "r_nuc_leu", "r_env_tum", "r_env_leu",
           "E_app_tum", "E_app_leu", "a_tum", "a_leu",
           "tau_tum", "tau_leu", "ell_factor"]
LIMITES = [[0.40, 0.60],      # an_a tumoral      TRANSFERIDO (Wang 2023/2025)
           [0.625, 0.763],    # an_a leucocito    MEDIDO (neutrofilo 0.625 - linfocito 0.763)
           [2.40, 37.0],      # C tumoral         TRANSFERIDO
           [2.40, 37.0],      # C leucocito       TRANSFERIDO (rango igual, sorteo independiente)
           [0.032, 0.078],    # Pi tumoral        MEDIDO (MCF-7 0.032 - MCF-10A 0.078)
           [0.113, 0.211],    # Pi leucocito      MEDIDO 0.162, banda DECLARADA de +/-30%
           [1.40, 4.00],      # r_nuc tumoral     TRANSFERIDO
           [1.40, 4.00],      # r_nuc leucocito   TRANSFERIDO
           [2.00, 8.00],      # r_env tumoral     TRANSFERIDO
           [2.00, 8.00],      # r_env leucocito   TRANSFERIDO
           [1.527, 2.772],    # E_app tumoral     MEDIDO (SBAT)
           [0.380, 0.660],    # E_app leucocito   MEDIDO (granulocito-PBMC, dRT-DC)
           [8.51, 10.00],     # a tumoral         MEDIDO/TRANSFERIDO
           [3.025, 4.395],    # a leucocito       MEDIDO (linfocito 6.05/2 - monocito 8.79/2)
           [0.37, 3.76],      # tau tumoral       MEDIDO
           [8.575e-3, 3.23],  # tau leucocito     NO MEDIDO (anclaje a horquilla completa)
           [0.50, 2.00]]      # ell_factor        DECLARADO


if __name__ == "__main__":
    d = pd.DataFrame(json.load(open(os.path.join(HERE, "diseno_f6.json")))["registros"])
    print("diseno ampliado: %d puntos, %d convergen (%.0f%%)"
          % (len(d), int(d.converge.sum()), 100 * d.converge.mean()), flush=True)
    d = d[d.converge & (d.f > 0)].copy()
    print("   cobertura an_a: %.3f-%.3f   Pi: %.4f-%.3f"
          % (d.an_a.min(), d.an_a.max(), d.Pi.min(), d.Pi.max()), flush=True)

    coef, diag = ajustar(d, "f")
    print("sustituto f: R2=%.5f  err mediano hold-out=%.2f%%  p90=%.2f%%  max=%.2f%%"
          % (diag["R2"], 100 * diag["err_mediano_holdout"],
             100 * diag["err_p90_holdout"], 100 * diag["err_max_holdout"]), flush=True)
    coefs_c, diags_c = {}, {}
    for c_ in COMP:
        col = "f_" + c_
        if col in d and (d[col] > 0).sum() > 60:
            coefs_c[c_], diags_c[c_] = ajustar(d[d[col] > 0], col)

    problema = {"num_vars": len(NOMBRES), "names": NOMBRES, "bounds": LIMITES}
    X = sobol_sample.sample(problema, N_BASE, calc_second_order=False)
    print("muestras: %d x %d" % X.shape, flush=True)
    c = {n: X[:, i] for i, n in enumerate(NOMBRES)}

    Gst = L.Gs_de_Eapp(c["E_app_tum"], c["a_tum"])
    Gsl = L.Gs_de_Eapp(c["E_app_leu"], c["a_leu"])

    # --- ENCUADRE A: confinamiento nominal igualado (convencion de la Fase 2)
    dP_t = dP_critica_sustituto(coef, c["a_tum"], Gst, CONF_DIS,
                                c["an_a_tum"], c["C_tum"], c["Pi_tum"],
                                c["r_nuc_tum"], c["r_env_tum"])
    dP_l = dP_critica_sustituto(coef, c["a_leu"], Gsl, CONF_DIS,
                                c["an_a_leu"], c["C_leu"], c["Pi_leu"],
                                c["r_nuc_leu"], c["r_env_leu"])
    selectividad = dP_t / dP_l

    Lc_t = np.array([L.longitud_de_carga(ai, ai / CONF_DIS, 1.6 * ai, ZT, LTAP)
                     for ai in c["a_tum"]])
    Lc_l = np.array([L.longitud_de_carga(ai, ai / CONF_DIS, 1.6 * ai, ZT, LTAP)
                     for ai in c["a_leu"]])
    razon_beta = (dP_t * c["tau_tum"] / Lc_t) / (dP_l * c["tau_leu"] / Lc_l)

    Rp_t, Rp_l = c["a_tum"] / CONF_DIS, c["a_leu"] / CONF_DIS
    cu = np.full_like(c["an_a_tum"], CONF_DIS)
    Wt = Gst * c["a_tum"]**2 * f_sustituto(coef, cu, c["an_a_tum"], c["C_tum"],
                                           c["Pi_tum"], c["r_nuc_tum"], c["r_env_tum"])
    Wl = Gsl * c["a_leu"]**2 * f_sustituto(coef, cu, c["an_a_leu"], c["C_leu"],
                                           c["Pi_leu"], c["r_nuc_leu"], c["r_env_leu"])
    sel_poro = (Wt / (c["ell_factor"] * c["a_tum"] * np.pi * Rp_t**2)) / \
               (Wl / (c["ell_factor"] * c["a_leu"] * np.pi * Rp_l**2))

    # --- fraccion marginal por compartimento (poblacion tumoral)
    frac_cortex = frac_nuc = None
    if set(COMP) <= set(coefs_c):
        marg = {}
        for c_ in COMP:
            eps = 1e-4
            fp = (f_sustituto(coefs_c[c_], cu + eps, c["an_a_tum"], c["C_tum"],
                              c["Pi_tum"], c["r_nuc_tum"], c["r_env_tum"])
                  - f_sustituto(coefs_c[c_], cu - eps, c["an_a_tum"], c["C_tum"],
                                c["Pi_tum"], c["r_nuc_tum"], c["r_env_tum"])) / (2 * eps)
            marg[c_] = np.abs(fp)
        tot = sum(marg.values())
        frac_cortex = marg["cortex"] / tot
        frac_nuc = (marg["nucleo"] + marg["envoltura"]) / tot

    salidas = {"dP_critica_tumoral_kPa": dP_t,
               "dP_critica_leucocito_kPa": dP_l,
               "selectividad_constriccion_tumoral_sobre_leucocito": selectividad,
               "razon_beta_c_tumoral_sobre_leucocito": razon_beta,
               "selectividad_microfiltro": sel_poro}
    if frac_cortex is not None:
        salidas["fraccion_marginal_cortex"] = frac_cortex
        salidas["fraccion_marginal_nuclear"] = frac_nuc

    filas = []
    for k, y in salidas.items():
        Si = sobol_analyze.analyze(problema, np.asarray(y, float),
                                   calc_second_order=False, print_to_console=False)
        for i, nm in enumerate(NOMBRES):
            filas.append(dict(salida=k, parametro=nm,
                              S1=float(Si["S1"][i]), S1_conf=float(Si["S1_conf"][i]),
                              ST=float(Si["ST"][i]), ST_conf=float(Si["ST_conf"][i]),
                              media=float(np.mean(y)), p5=float(np.percentile(y, 5)),
                              p50=float(np.percentile(y, 50)),
                              p95=float(np.percentile(y, 95))))
    sob = pd.DataFrame(filas)
    sob.to_csv(os.path.join(HERE, "sobol_por_poblacion.csv"), index=False)

    print("\n== indices totales ST (los 6 mayores por salida) ==", flush=True)
    for k in salidas:
        s = sob[sob.salida == k].nlargest(6, "ST")
        print("  %s  [p5=%.4g p50=%.4g p95=%.4g]"
              % (k, s.p5.iloc[0], s.p50.iloc[0], s.p95.iloc[0]))
        for _, r in s.iterrows():
            print("      %-12s ST=%.4f +/- %.4f   S1=%.4f" % (r.parametro, r.ST, r.ST_conf, r.S1))

    concl = {
        "la tumoral exige mas presion que el leucocito": float(np.mean(selectividad > 1.0)),
        "selectividad > 2x": float(np.mean(selectividad > 2.0)),
        "selectividad > 4x": float(np.mean(selectividad > 4.0)),
        "microfiltro retiene tumoral y deja pasar leucocito": float(np.mean(sel_poro > 1.0)),
        "ventana de microfiltro > 2x": float(np.mean(sel_poro > 2.0)),
        "existe ventana de operacion beta_c(tum) > beta_c(leu)": float(np.mean(razon_beta > 1.0)),
        "ventana de presion de dispositivo > 10x": float(np.mean(razon_beta > 10.0)),
    }
    if frac_cortex is not None:
        concl["el cortex es el limitante (>50% del margen)"] = float(np.mean(frac_cortex > 0.5))
        concl["los compartimentos nucleares aportan <10%"] = float(np.mean(frac_nuc < 0.10))
    print("\n== fraccion de muestras en que se cumple ==", flush=True)
    for k, v in concl.items():
        print("  %-55s %.4f" % (k, v))

    json.dump({"problema": problema, "n_muestras": int(len(X)),
               "diagnostico_sustituto": diag, "conclusiones": concl,
               "percentiles": {k: dict(p5=float(np.percentile(v, 5)),
                                       p50=float(np.percentile(v, 50)),
                                       p95=float(np.percentile(v, 95)),
                                       media=float(np.mean(v)))
                               for k, v in salidas.items()}},
              open(os.path.join(HERE, "sobol_f6.json"), "w"), indent=1)
    np.save(os.path.join(HERE, "coef_sustituto_f6.npy"), coef)
    print("\nlisto", flush=True)
