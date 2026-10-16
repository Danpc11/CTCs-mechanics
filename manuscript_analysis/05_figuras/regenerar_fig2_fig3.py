"""Regenera Fig. 2 y Fig. 3 con la reduccion corregida (07_correccion_reduccion).
Fig. 2a: presiones criticas de tabla2_corregida_directa.csv (FEM re-ejecutado, geometria de referencia).
Fig. 3a: beta_c corregidos (umbral_dispositivo_corregido.csv).
Fig. 3b: mismas curvas; las etiquetas de beta se multiplican por el factor comun k = 6.56, porque la
separacion depende de beta solo via beta/dP_el y beta/beta_c (exacto con un factor comun)."""
import os
import pandas as pd
os.makedirs('salida', exist_ok=True)
D = '../datos/'
K = 6.56
src = open('figuras_bmmb.py', encoding='utf-8').read()
src = src.replace('a.text(2.40, 3.2, "validated envelope", color=CTR, fontsize=6, ha="right", va="top")', 'a.text(1.36, 7.5, "assessed range", color=CTR, fontsize=6, ha="left", va="top")')
src = src.replace('"envelope"', '"assessed range"')
src = src.replace('a.text(h.confinamiento.iloc[-1] + 0.045,', 'a.text(2.48,')
src = src.replace('lab_off = {"MCF-7": 0.72, "SK-BR-3": 1.45, "MDA-MB-231": 0.72,', 'lab_off = {"MCF-7": 1.3, "SK-BR-3": 1.0, "MDA-MB-231": 0.72,').replace('a.set_ylim(1.5e-3, 4)', 'a.set_ylim(1.5e-3, 10)')
src = src.replace('r"$\\beta$ = 10 kPa per mm s$^{-1}$"', 'rf"$\\beta$ = {10.0 * K:.0f} kPa per mm s$^{{-1}}$"')
src = src.replace('r"$\\beta$ = 164.5"', 'rf"$\\beta$ = {164.5 * K:.0f}"')
src = src.replace('r"165$\\times$"', 'rf"{vals[\'MDA-MB-231\'] / vals[\'PBMC\']:.0f}$\\times$"')
src = src.replace('"Kelvin--Voigt"', '"Kelvin–Voigt"')
G = {'K': K}
exec(compile(src, 'figuras_bmmb.py', 'exec'), G)
t = pd.read_csv(D + 'tabla2_corregida_directa.csv').dropna(subset=['corregido'])
# puntos FEM no convergidos (n.c. en la Tabla 2): no se dibujan ni se interpolan
t = t[~(((t.estirpe == 'SK-BR-3') & (t.conf_nominal == 2.43)) | ((t.estirpe == 'MDA-MB-231') & (t.conf_nominal == 1.90)))]
bc = pd.DataFrame(dict(geometria='referencia', estirpe=t.estirpe, confinamiento=t.conf_nominal, dP_critica_kPa=t.corregido))
G['fig_contact'](bc, pd.read_csv(D + 'encuadre_mismo_dispositivo.csv'), 'salida/fig2_contact.png')
umb = pd.read_csv(D + 'umbral_dispositivo_corregido.csv')
umb['beta_c_kPa_a_1mm_s'] = umb['beta_c_corregido_kPa_a_1mm_s']
G['fig_rate'](umb, pd.read_csv(D + 'pureza_rendimiento_caudal.csv'), pd.read_csv(D + 'cota_ley_potencias_ventana.csv'), 'salida/fig3_rate.png')
print('ok')
