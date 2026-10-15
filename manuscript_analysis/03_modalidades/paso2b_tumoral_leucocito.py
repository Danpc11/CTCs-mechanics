"""Paso 2-B: EJE PRINCIPAL -- tumoral frente al fondo leucocitario.

Probabilidad de captura ANALITICA (no Monte Carlo) para poder resolver colas
de 1e-7, que es lo que exige una pureza con 1 CTC por millon de leucocitos.

  captura  <=>  m(U) <= 1  <=>  tau >= tau_umbral(U) = (L/dP_el) (beta - dP_el/U)
  con tau ~ N(mu, sigma) truncada en tau > 0:
      P(captura) = P(tau >= tau_umbral)

El umbral de dispositivo por encima del cual una estirpe escapa a algun caudal es
      beta_c = dP_el * tau / L_carga .
"""
import sys, os, json
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L

K_EAPP = L.K_EAPP
CONF = 1.70

# dP_el calculado por FEM a confinamiento nominal 1.70 (paso123.json) y
# L_carga geometrica.  dP_el es EXACTAMENTE proporcional a E_app (= Gs/a * cte),
# asi que se escala a cada referencia leucocitaria con su E_app medido.
J = json.load(open(os.path.join(HERE, "paso123.json")))
c1 = pd.read_csv(os.path.join(HERE, "barrido_constriccion.csv"))
ref = c1[(c1.geometria == "referencia") & (c1.confinamiento == CONF)]
DP_EL = dict(zip(ref.estirpe, ref.dP_critica_kPa))
L_CARGA = float(J["paso2_base"]["MCF-7"]["L_carga"])

# ---- grupo TUMORAL (las tres lineas con modulo aparente medido)
TUMORAL = {
    "MCF-7":      dict(E=2.650, dP=DP_EL["MCF-7"], tau=3.23, s_tau=0.53,
                       tau_src="medido"),
    "SK-BR-3":    dict(E=2.772, dP=DP_EL["SK-BR-3"], tau=None, s_tau=None,
                       tau_src="NO MEDIDO"),
    "MDA-MB-231": dict(E=1.527, dP=DP_EL["MDA-MB-231"], tau=0.57, s_tau=0.20,
                       tau_src="medido"),
}
# MCF-10A es la linea epitelial NO TUMORIGENICA del panel y ademas la unica sin
# modulo aparente medido (su Gs esta transferido de MCF-7).  Agruparla dentro de
# 'tumoral' seria incorrecto desde el punto de vista biologico, asi que el grupo
# principal son las tres lineas tumorigenicas; pero se reporta TAMBIEN el
# rendimiento del grupo de las cuatro lineas de mama, que es la agregacion
# literalmente pedida, para que la diferencia sea visible y no una omision.
CONTROL = {"MCF-10A": dict(E=None, dP=DP_EL["MCF-10A"], tau=0.77, s_tau=0.37,
                           tau_src="medido; E_app NO medido, Gs transferido")}
MAMA_4 = dict(TUMORAL); MAMA_4.update(CONTROL)

# ---- referencias LEUCOCITARIAS (modulo aparente medido, dRT-DC)
# tau del leucocito: NO hay microrreologia de citoplasma.  Unico anclaje:
# viscosidad aparente de CELULA ENTERA de HL-60 por Kelvin-Voigt,
# eta = 3.43 +/- 0.16 Pa s con E = 0.400 +/- 0.02 kPa (10.1038/s41467-019-08370-3)
# -> tau = eta/E = 8.575e-3 s.  Es de CELULA ENTERA, no de citoplasma, y por eso
# no es homologo a los tau de microrreologia de las lineas de mama.
TAU_LEU = 3.43 / 400.0
S_TAU_LEU = TAU_LEU * np.hypot(0.16 / 3.43, 0.02 / 0.400)
LEUCO = {
    "PBMC":        dict(E=0.620, dP=DP_EL["leucocito"] * 0.620 / 0.620),
    "granulocito": dict(E=0.430, dP=DP_EL["leucocito"] * 0.430 / 0.620),
    "HL-60":       dict(E=0.400, dP=DP_EL["leucocito"] * 0.400 / 0.620),
}
for v in LEUCO.values():
    v.update(tau=TAU_LEU, s_tau=S_TAU_LEU, tau_src="anclaje de celula entera")

RATIO_CTC = 1e-6          # 1 CTC por millon de leucocitos (declarado)


def p_captura(U, beta, dP_el, tau_mu, tau_sd):
    """P(tau >= tau_umbral) con tau normal truncada en 0."""
    U = np.asarray(U, dtype=float)
    thr = (L_CARGA / dP_el) * (beta - dP_el / U)
    z0 = -tau_mu / max(tau_sd, 1e-12)
    norma = 1.0 - norm.cdf(z0)
    p = (1.0 - norm.cdf((thr - tau_mu) / max(tau_sd, 1e-12))) / norma
    return np.clip(np.where(thr <= 0.0, 1.0, p), 0.0, 1.0)


def beta_c(dP_el, tau):
    return dP_el * tau / L_CARGA


if __name__ == "__main__":
    out = {"CONF_nominal": CONF, "L_carga_um": L_CARGA,
           "tau_leucocito_s": TAU_LEU, "tau_leucocito_sd_s": float(S_TAU_LEU),
           "tau_leucocito_fuente": "eta/E de celula entera de HL-60 "
                                   "(10.1038/s41467-019-08370-3); NO es "
                                   "microrreologia de citoplasma",
           "ratio_CTC_leucocito": RATIO_CTC}

    print("== umbral de dispositivo beta_c = dP_el*tau/L_carga ==")
    print("%-14s %8s %10s %8s %14s" % ("poblacion", "E_app", "dP_el", "tau",
                                       "kPa a 1 mm/s"))
    filas_bc = []
    for grupo, dic in (("tumoral", TUMORAL), ("control", CONTROL),
                       ("leucocito", LEUCO)):
        for nm, v in dic.items():
            taus = ([v["tau"]] if v["tau"] is not None
                    else [0.57, 3.23])       # SK-BR-3: bracket declarado
            for t in taus:
                bc = beta_c(v["dP"], t)
                filas_bc.append(dict(grupo=grupo, poblacion=nm,
                                     E_app_kPa=v["E"], dP_el_kPa=v["dP"],
                                     tau_s=t, tau_fuente=v["tau_src"],
                                     beta_c_kPa_s_um=bc,
                                     beta_c_kPa_a_1mm_s=bc * 1000.0))
                print("%-14s %8s %10.5f %8.4g %14.3f"
                      % (nm + ("" if v["tau"] is not None else " (bracket)"),
                         "%.3f" % v["E"] if v["E"] else "  --  ",
                         v["dP"], t, bc * 1000.0))
    pd.DataFrame(filas_bc).to_csv(
        os.path.join(HERE, "umbral_dispositivo.csv"), index=False)

    # ---- cota sobre el tau leucocitario que preserva la conclusion
    dP_tum_min = min(v["dP"] for v in TUMORAL.values())
    nm_min = min(TUMORAL, key=lambda k: TUMORAL[k]["dP"])
    tau_tum_min = TUMORAL["MDA-MB-231"]["tau"]
    cotas = {}
    for nm, v in LEUCO.items():
        cotas[nm] = {
            "frente_a_la_tumoral_mas_blanda": float(
                (dP_tum_min / v["dP"]) * tau_tum_min),
            "frente_a_MCF-7": float((TUMORAL["MCF-7"]["dP"] / v["dP"])
                                    * TUMORAL["MCF-7"]["tau"]),
            "margen_sobre_el_anclaje": float(
                (dP_tum_min / v["dP"]) * tau_tum_min / TAU_LEU)}
    out["cota_tau_leucocito_s"] = cotas
    out["tumoral_limitante"] = nm_min
    print("\n== cota sobre tau del leucocito (s) que mantiene la ventana ==")
    print("  (la ventana existe mientras beta_c(leucocito) < beta_c(tumoral))")
    for nm, c in cotas.items():
        print("  %-12s tau_leu < %.3f s frente a %s  (= %.0f veces el anclaje "
              "de %.4g s)" % (nm, c["frente_a_la_tumoral_mas_blanda"], nm_min,
                              c["margen_sobre_el_anclaje"], TAU_LEU))

    # ---- barrido de caudal: rendimiento, arrastre leucocitario y pureza
    U = np.logspace(-2, 5, 220)
    BET = np.array([0.01, 0.02, 0.05, 0.10, 0.1645, 0.20, 0.30, 0.60, 1.0])
    filas, res = [], {}
    for beta in BET:
        cap_t = {}
        for nm, v in TUMORAL.items():
            if v["tau"] is not None:
                cap_t[nm] = p_captura(U, beta, v["dP"], v["tau"], v["s_tau"])
            else:    # SK-BR-3: banda con el bracket de tau
                lo = p_captura(U, beta, v["dP"], 0.57, 0.20)
                hi = p_captura(U, beta, v["dP"], 3.23, 0.53)
                cap_t[nm] = 0.5 * (lo + hi)
                cap_t[nm + "_lo"], cap_t[nm + "_hi"] = lo, hi
        rend = np.mean([cap_t[nm] for nm in TUMORAL], axis=0)
        cap10 = p_captura(U, beta, CONTROL["MCF-10A"]["dP"],
                          CONTROL["MCF-10A"]["tau"], CONTROL["MCF-10A"]["s_tau"])
        rend4 = np.mean([cap_t[nm] for nm in TUMORAL] + [cap10], axis=0)
        for nml, vl in LEUCO.items():
            arr = p_captura(U, beta, vl["dP"], vl["tau"], vl["s_tau"])
            pureza = (rend * RATIO_CTC) / np.maximum(
                rend * RATIO_CTC + arr, 1e-300)
            enriq = np.where(arr > 0, rend / np.maximum(arr, 1e-300), np.inf)
            sep = rend - arr
            for i, u in enumerate(U):
                filas.append(dict(
                    beta_kPa_s_um=beta, dP_dispositivo_a_1mm_s_kPa=beta * 1000,
                    U_um_s=u, U_mm_s=u / 1000.0, referencia_leucocitaria=nml,
                    rendimiento_tumoral=rend[i],
                    arrastre_leucocitario=arr[i],
                    separacion=sep[i], pureza=pureza[i],
                    factor_enriquecimiento=min(enriq[i], 1e300),
                    captura_MCF7=cap_t["MCF-7"][i],
                    captura_SKBR3=cap_t["SK-BR-3"][i],
                    captura_MDAMB231=cap_t["MDA-MB-231"][i],
                    captura_MCF10A=cap10[i],
                    rendimiento_mama_4lineas=rend4[i],
                    separacion_mama_4lineas=rend4[i] - arr[i]))
            j = int(np.argmax(sep))
            res["%g|%s" % (beta, nml)] = dict(
                separacion_max=float(sep[j]), U_optima_um_s=float(U[j]),
                maximo_interior=bool(0 < j < len(U) - 1
                                     and sep[j] > sep[-1] + 1e-6
                                     and sep[j] > sep[0] + 1e-6),
                sep_U_min=float(sep[0]), sep_U_max=float(sep[-1]),
                pureza_en_el_optimo=float(pureza[j]),
                rendimiento_en_el_optimo=float(rend[j]))
    c = pd.DataFrame(filas)
    c.to_csv(os.path.join(HERE, "pureza_rendimiento_caudal.csv"), index=False)
    out["optimos"] = res

    print("\n== separacion tumoral - leucocito frente al caudal ==")
    print("%10s %-12s %10s %12s %9s %9s %9s"
          % ("kPa@1mm/s", "referencia", "sep_max", "U_opt(um/s)", "interior",
             "sep(Umin)", "sep(Umax)"))
    for k, v in res.items():
        b, nml = k.split("|")
        print("%10.1f %-12s %10.4f %12.4g %9s %9.4f %9.4f"
              % (float(b) * 1000, nml, v["separacion_max"], v["U_optima_um_s"],
                 v["maximo_interior"], v["sep_U_min"], v["sep_U_max"]))

    # mejor punto de operacion global por pureza x rendimiento
    c["score"] = c.pureza * c.rendimiento_tumoral
    best = c.loc[c.score.idxmax()]
    out["mejor_punto_operacion"] = {k: (None if pd.isna(best[k]) else
                                        float(best[k]) if not isinstance(best[k], str)
                                        else best[k])
                                    for k in ("beta_kPa_s_um",
                                              "dP_dispositivo_a_1mm_s_kPa",
                                              "U_um_s", "rendimiento_tumoral",
                                              "arrastre_leucocitario", "pureza",
                                              "factor_enriquecimiento")}
    out["mejor_punto_operacion"]["referencia_leucocitaria"] = best.referencia_leucocitaria
    print("\n== mejor punto de operacion (pureza x rendimiento, 1 CTC / 1e6 leucocitos) ==")
    for k, v in out["mejor_punto_operacion"].items():
        print("   %-32s %s" % (k, v))
    # ---- ANCHO del intervalo de caudal util y SENSIBILIDAD al ancho de la
    #      distribucion de tau leucocitario.  La dispersion publicada de
    #      tau_leu procede del ERROR ESTANDAR de eta y E (CV 3.3%), no de la
    #      variabilidad celula a celula: la cola gaussiana resultante es
    #      irrealmente estrecha y la pureza que sale de ella es una COTA
    #      SUPERIOR.  Se reporta el arrastre frente al CV supuesto.
    print("\n== ancho del intervalo de caudal con separacion > 0.9 (decadas) ==")
    anchos = []
    for beta in BET:
        g = c[(c.beta_kPa_s_um == beta) &
              (c.referencia_leucocitaria == "PBMC")].sort_values("U_um_s")
        ok = g[g.separacion > 0.9]
        if len(ok) > 1:
            dec = float(np.log10(ok.U_um_s.max() / ok.U_um_s.min()))
            anchos.append(dict(beta_kPa_s_um=float(beta),
                               dP_a_1mm_s_kPa=float(beta * 1000),
                               U_min_um_s=float(ok.U_um_s.min()),
                               U_max_um_s=float(ok.U_um_s.max()),
                               decadas=dec))
            print("   %8.1f kPa@1mm/s : U = %.4g - %.4g um/s  (%.2f decadas)"
                  % (beta * 1000, ok.U_um_s.min(), ok.U_um_s.max(), dec))
    out["ancho_intervalo_caudal"] = anchos

    print("\n== arrastre leucocitario en el optimo frente al CV supuesto de tau_leu ==")
    cvs = [0.033, 0.10, 0.25, 0.50, 1.00]
    sens = []
    b_op, u_op = 0.01, float(res["0.01|PBMC"]["U_optima_um_s"])
    for cv in cvs:
        fila = dict(CV_tau_leucocito=cv)
        for nml, vl in LEUCO.items():
            arr = float(p_captura(np.array([u_op]), b_op, vl["dP"], TAU_LEU,
                                  cv * TAU_LEU)[0])
            fila["arrastre_" + nml] = arr
            fila["pureza_" + nml] = (1.0 * RATIO_CTC) / (1.0 * RATIO_CTC + arr) \
                if arr > 0 else 1.0
        sens.append(fila)
        print("   CV=%.3f  arrastre PBMC=%.3e  pureza=%.4f"
              % (cv, fila["arrastre_PBMC"], fila["pureza_PBMC"]))
    out["sensibilidad_CV_tau_leucocito"] = sens
    out["punto_operacion_evaluado"] = dict(beta_kPa_s_um=b_op, U_um_s=u_op)
    pd.DataFrame(sens).to_csv(os.path.join(HERE, "sensibilidad_tau_leucocito.csv"),
                              index=False)

    # comparacion de las dos agregaciones
    g10 = c[(c.beta_kPa_s_um == 0.01) & (c.referencia_leucocitaria == "PBMC")]
    j = int(np.argmax(g10.separacion.values))
    out["agregacion"] = dict(
        grupo_principal="3 lineas tumorigenicas (MCF-7, SK-BR-3, MDA-MB-231)",
        grupo_alternativo="4 lineas de mama (incluye MCF-10A, no tumorigenica "
                          "y sin E_app medido)",
        rendimiento_3lineas_en_el_optimo=float(g10.rendimiento_tumoral.values[j]),
        rendimiento_4lineas_en_el_optimo=float(
            g10.rendimiento_mama_4lineas.values[j]),
        separacion_3lineas=float(g10.separacion.values[j]),
        separacion_4lineas=float(g10.separacion_mama_4lineas.values[j]),
        nota="beta_c de MCF-10A = 40,4 kPa por mm/s queda entre MDA-MB-231 "
             "(17,0) y MCF-7 (169,4), asi que incluirla NO cambia ni la linea "
             "limitante del grupo ni la ventana de presion de dispositivo.")
    print("\n== agregacion: 3 lineas tumorigenicas frente a las 4 de mama ==")
    for k_, v_ in out["agregacion"].items():
        print("   %-34s %s" % (k_, v_))

    n_int = sum(1 for v in res.values() if v["maximo_interior"])
    out["n_combinaciones_con_maximo_interior"] = n_int
    out["n_combinaciones"] = len(res)
    print("\nmaximo interior en %d de %d combinaciones (beta x referencia)"
          % (n_int, len(res)))

    json.dump(out, open(os.path.join(HERE, "paso2b.json"), "w"), indent=1,
              default=str)
