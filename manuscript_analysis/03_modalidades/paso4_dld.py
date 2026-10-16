"""Paso 4: descarte cuantitativo de DLD por solapamiento de tamanos.

No se modela DLD: el solver solo es valido en regimen de CONTACTO (2a/D >= 1.2)
y DLD opera sin contacto.  El descarte se hace con los radios medidos en
suspension y con el propio criterio de diametro critico de DLD.

CORRELACION DE DLD (Davis 2006, empirica):   Dc = 1.4 * g * eps^0.48
con g el hueco entre postes y eps la fraccion de desplazamiento de fila.
Su exactitud publicada es del orden del 10% de Dc.
"""
import sys, os, json
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L

# diametro (um) y semidispersion PUBLICADA como ERROR ESTANDAR, no desviacion
# tipica.  Usar SE como si fuera SD da el solapamiento MINIMO compatible con
# los datos: SD >= SE siempre.  Es por tanto una cota inferior del solapamiento
# y una cota superior de la separabilidad.
DIAM = {
    "MCF-7":      dict(D=17.48, se=0.31, fuente="10.1039/c7ra10750a", transferido=False),
    "MDA-MB-231": dict(D=17.24, se=0.21, fuente="10.1039/c7ra10750a", transferido=False),
    "SK-BR-3":    dict(D=20.00, se=None, fuente="10.1038/s41378-020-0150-6",
                       transferido=True),
    "MCF-10A":    dict(D=17.48, se=None, fuente="10.1039/c7ra10750a",
                       transferido=True),
    "leucocito":  dict(D=17.00, se=None, fuente="10.1007/s10237-023-01705-6",
                       transferido=True),
}
RES_DLD_REL = 0.10          # exactitud relativa de la correlacion de Dc


def dc_davis(g, eps):
    return 1.4 * g * eps**0.48


def solape_gauss(m1, s1, m2, s2, n=200001):
    """Coeficiente de solapamiento OVL = int min(f1, f2), y coeficiente de
    Bhattacharyya."""
    lo = min(m1 - 6 * s1, m2 - 6 * s2)
    hi = max(m1 + 6 * s1, m2 + 6 * s2)
    x = np.linspace(lo, hi, n)
    f1, f2 = norm.pdf(x, m1, s1), norm.pdf(x, m2, s2)
    ovl = float(np.trapezoid(np.minimum(f1, f2), x))
    bc = float(np.trapezoid(np.sqrt(f1 * f2), x))
    return ovl, bc


def youden_max(m_t, s_t, m_l, s_l):
    """Maximo indice de Youden (sens+esp-1) sobre todos los Dc posibles,
    tratando 'tumoral' como positivo y reteniendo los de D > Dc."""
    dc = np.linspace(min(m_t, m_l) - 6 * max(s_t, s_l),
                     max(m_t, m_l) + 6 * max(s_t, s_l), 40001)
    sens = 1.0 - norm.cdf(dc, m_t, s_t)       # tumoral desplazada
    esp = norm.cdf(dc, m_l, s_l)              # leucocito en zigzag
    J = sens + esp - 1.0
    i = int(np.argmax(J))
    return float(J[i]), float(dc[i]), float(sens[i]), float(esp[i])


if __name__ == "__main__":
    out = {"correlacion": "Dc = 1.4 g eps^0.48 (Davis 2006)",
           "exactitud_Dc_relativa": RES_DLD_REL, "filas": [], "notas": []}

    # ---- 1. argumento que NO necesita ninguna dispersion: brecha de medias
    #         frente a la resolucion propia del metodo
    Dl = DIAM["leucocito"]["D"]
    print("ARGUMENTO 1 - brecha de medias frente a la resolucion de Dc")
    print("%-12s %8s %10s %10s %12s %s"
          % ("estirpe", "D(um)", "brecha", "brecha/D", "res.Dc(um)", "separable"))
    for nm, d in DIAM.items():
        if nm == "leucocito":
            continue
        brecha = d["D"] - Dl
        res = RES_DLD_REL * d["D"]
        sep = abs(brecha) > res
        out["filas"].append(dict(
            estirpe=nm, D_um=d["D"], D_transferido=d["transferido"],
            D_leucocito_um=Dl, brecha_um=round(brecha, 3),
            brecha_relativa=round(brecha / d["D"], 5),
            resolucion_Dc_um=round(res, 3), separable_por_tamano=bool(sep),
            veces_bajo_resolucion=round(res / abs(brecha), 2) if brecha else None,
            fuente=d["fuente"]))
        print("%-12s %8.2f %10.2f %10.4f %12.2f  %s%s"
              % (nm, d["D"], brecha, brecha / d["D"], res,
                 "SI" if sep else "NO",
                 "  [D TRANSFERIDO]" if d["transferido"] else ""))

    # geometrias de DLD que darian ese Dc
    geos = []
    for eps in (0.05, 0.10, 0.20):
        g = brentq(lambda gg: dc_davis(gg, eps) - Dl, 1.0, 500.0)
        geos.append(dict(eps=eps, hueco_um=round(g, 2), Dc_um=round(Dl, 2),
                         Dc_incertidumbre_um=round(RES_DLD_REL * Dl, 2)))
    out["geometrias_para_Dc_leucocito"] = geos
    print("\n  geometrias de array que centran Dc en el leucocito (%.1f um):" % Dl)
    for g in geos:
        print("    eps=%.2f -> hueco g=%.1f um ; Dc = %.1f +/- %.1f um"
              % (g["eps"], g["hueco_um"], g["Dc_um"], g["Dc_incertidumbre_um"]))

    # ---- 2. solapamiento de distribuciones con SE-como-SD (cota inferior)
    print("\nARGUMENTO 2 - solapamiento con SE tratado como SD (solapamiento MINIMO)")
    sl = DIAM["MCF-7"]["se"]            # el leucocito no tiene dispersion
    out["notas"].append(
        "El leucocito NO tiene dispersion de tamano publicada en la tabla; "
        "para el solapamiento se le asigna la de MCF-7 (SE = %.2f um) y se "
        "declara como transferida." % sl)
    print("%-12s %9s %9s %9s %9s %9s" % ("estirpe", "OVL", "Bhatt.", "J_max",
                                         "Dc_opt", "sens/esp"))
    sol = []
    for nm in ("MCF-7", "MDA-MB-231"):
        d = DIAM[nm]
        ovl, bc = solape_gauss(d["D"], d["se"], Dl, sl)
        J, dco, se_, sp_ = youden_max(d["D"], d["se"], Dl, sl)
        sol.append(dict(estirpe=nm, sigma_um=d["se"], sigma_tipo="SE como SD",
                        OVL=round(ovl, 4), bhattacharyya=round(bc, 4),
                        J_max=round(J, 4), Dc_optimo_um=round(dco, 3),
                        sensibilidad=round(se_, 4), especificidad=round(sp_, 4)))
        print("%-12s %9.4f %9.4f %9.4f %9.3f  %.3f/%.3f"
              % (nm, ovl, bc, J, dco, se_, sp_))
    out["solapamiento_SE_como_SD"] = sol

    # ---- 3. solapamiento frente al CV supuesto (sin inventar un CV)
    print("\nARGUMENTO 3 - solapamiento frente al coeficiente de variacion supuesto")
    cvs = np.array([0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20])
    curva = []
    print("   CV      sigma(um)   OVL(MCF-7)  J_max(MCF-7)  OVL(MDA)  J_max(MDA)")
    for cv in cvs:
        fila = dict(CV=float(cv))
        for nm in ("MCF-7", "MDA-MB-231"):
            d = DIAM[nm]
            s_t, s_l = cv * d["D"], cv * Dl
            ovl, bc = solape_gauss(d["D"], s_t, Dl, s_l)
            J, dco, _, _ = youden_max(d["D"], s_t, Dl, s_l)
            fila["OVL_" + nm] = round(ovl, 4)
            fila["J_" + nm] = round(J, 4)
        curva.append(fila)
        print("  %.3f   %8.3f   %10.4f  %12.4f  %8.4f  %9.4f"
              % (cv, cv * 17.48, fila["OVL_MCF-7"], fila["J_MCF-7"],
                 fila["OVL_MDA-MB-231"], fila["J_MDA-MB-231"]))
    out["solapamiento_vs_CV"] = curva

    # CV al que el solapamiento supera el 50% y el 80%
    umbrales = {}
    for nm in ("MCF-7", "MDA-MB-231"):
        d = DIAM[nm]
        for thr in (0.5, 0.8):
            def fn(cv, nm=nm, d=d, thr=thr):
                o, _ = solape_gauss(d["D"], cv * d["D"], Dl, cv * Dl)
                return o - thr
            try:
                umbrales["CV_OVL%d_%s" % (thr * 100, nm)] = round(
                    brentq(fn, 1e-4, 0.5), 5)
            except ValueError:
                umbrales["CV_OVL%d_%s" % (thr * 100, nm)] = None
    out["CV_umbral"] = umbrales
    print("\n  CV al que el solapamiento alcanza 50% / 80%:")
    for k, v in umbrales.items():
        print("    %-24s %s" % (k, "%.4f" % v if v else "n/a"))

    json.dump(out, open(os.path.join(HERE, "paso4_dld.json"), "w"), indent=1)
