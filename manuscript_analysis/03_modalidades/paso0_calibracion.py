"""Paso 0: calibracion del cortex contra el modulo aparente MEDIDO.

Es el paso que sustituye la transferencia de Gs de la Fase 1 (PC-3 para las
lineas de mama, HL-60 para el leucocito) por una calibracion contra el unico
observable mecanico medido linea a linea.

Sonda: traccion normal seguidora de patron cuadrupolar sigma_rr = p0 P2(cos th)
sobre la celula LIBRE (sin contacto), que es el modo que impone una pinza
acustica de haz unico; E_app se define operativamente como el modulo de la
esfera homogenea casi incompresible que mostraria la misma amplitud bajo la
misma traccion.  Se usa la respuesta linealizada.

Con la regularizacion de citoplasma de produccion (Pi = G_cyt a / Gs = 0.10,
heredada de la Fase 1) TODO escala con Gs, de modo que la relacion es
exactamente afin,  E_app = K_EAPP * Gs / a,  y la calibracion se reduce a
invertirla:  Gs = E_app a / K_EAPP.  Este script (i) calibra la constante
adimensional de la sonda, (ii) comprueba que K_EAPP es constante al variar Gs,
el radio y la malla, y (iii) emite el Gs de cada estirpe.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fase2_lib as L
HERE = os.path.dirname(os.path.abspath(__file__))
HREL = 20.0

K = L.calibrar_sonda_P2(hrel=HREL)
print("constante adimensional de la sonda P2:  K = %.5f  (hrel = %.0f)" % (K, HREL))

# --- (ii) K_EAPP = E_app * a / Gs debe ser constante bajo Pi fijo
print("\ncomprobacion de que K_EAPP = E_app*a/Gs es constante (Pi = %.2f):"
      % L.PI_REG)
print("   %6s %6s %6s %11s %10s" % ("Gs", "a", "hrel", "E_app(kPa)", "K_EAPP"))
chk = []
for Gs, a_, hr in [(2.0, 8.74, HREL), (6.0, 8.74, HREL), (12.0, 8.74, HREL),
                   (12.0, 8.74, 14.0), (12.0, 8.74, 26.0),
                   (12.0, 10.0, HREL), (12.0, 8.50, HREL)]:
    e = L.modulo_aparente(a_, 0.5 * a_, Gs, 10.0, L.PI_REG * Gs / a_,
                          L.ESTIRPES["MCF-7"].r_nuc,
                          L.ESTIRPES["MCF-7"].r_env, K=K, hrel=hr)
    k_e = e * a_ / Gs
    chk.append(dict(Gs=Gs, a=a_, hrel=hr, E_app=e, K_EAPP=k_e))
    print("   %6.1f %6.2f %6.0f %11.5f %10.5f" % (Gs, a_, hr, e, k_e))
ks = np.array([c["K_EAPP"] for c in chk])
print("   K_EAPP = %.4f ; dispersion %.2f %% ; valor en uso %.3f"
      % (ks.mean(), 100 * (ks.max() / ks.min() - 1), L.K_EAPP))

# --- (iii) Gs por estirpe
out = {"K_P2": K, "hrel": HREL, "Pi": L.PI_REG, "K_EAPP_en_uso": L.K_EAPP,
       "K_EAPP_comprobacion": chk, "lineas": []}
print("\n%-12s %10s %10s %10s %8s" % ("estirpe", "E_app medido", "Gs calibrado",
                                     "Gs Fase 1", "cociente"))
for nm, e in L.ESTIRPES.items():
    if e.E_app is None:
        gs = L.Gs_de_Eapp(L.ESTIRPES["MCF-7"].E_app, L.ESTIRPES["MCF-7"].a)
        out["lineas"].append(dict(
            estirpe=nm, E_app_medido=None, Gs=gs, Gs_literatura=e.Gs_lit,
            nota="E_app NO MEDIDO: Gs transferido de MCF-7 y declarado"))
        print("%-12s %10s %10.3f %10.3f %8s   [E_app NO MEDIDO: Gs de MCF-7]"
              % (nm, "--", gs, e.Gs_lit, "--"))
        continue
    gs = L.Gs_de_Eapp(e.E_app, e.a)
    out["lineas"].append(dict(
        estirpe=nm, E_app_medido=e.E_app, E_app_rango=list(e.E_app_rng),
        E_app_fuente=e.E_app_src, a_um=e.a, Gs=gs, Gs_literatura=e.Gs_lit,
        Gs_lit_fuente=e.Gs_lit_src, cociente_cal_lit=gs / e.Gs_lit,
        Gs_en_los_extremos_de_E_app=[L.Gs_de_Eapp(v, e.a) if v else None
                                     for v in e.E_app_rng]))
    print("%-12s %10.3f %10.3f %10.3f %8.2f"
          % (nm, e.E_app, gs, e.Gs_lit, gs / e.Gs_lit))

json.dump(out, open(os.path.join(HERE, "calibracion_Gs.json"), "w"), indent=1)
