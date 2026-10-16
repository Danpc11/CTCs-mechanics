"""Ensamblado de la tabla de parametros mecanicos por compartimento y subtipo.
Cada fila lleva fuente (DOI), tecnica, tipo celular original y bandera de transferencia.
Ningun valor se inventa: 'nominal' es o bien el valor publicado, o bien el punto medio
de un rango publicado (marcado en nominal_es_punto_medio)."""
import pandas as pd, numpy as np

R = []
def add(cel, comp, prop, nom, mn, mx, unit, doi, tec, orig, transf, just="", verb="", mid=False):
    R.append(dict(tipo_celular=cel, compartimento=comp, propiedad=prop,
                  valor_nominal=nom, minimo=mn, maximo=mx, unidad=unit,
                  fuente_doi=doi, tecnica=tec, tipo_celular_original=orig,
                  es_transferido=transf, nominal_es_punto_medio=mid,
                  justificacion_transferencia=just, cita_verbatim=verb))

# --- Constantes de Wang et al. 2023 (HL-60), tal como se publican -------------
W23 = "10.1007/s10237-023-01705-6"; W25 = "10.1007/s10237-024-01920-9"
mu0_w23 = 32.5e-3          # Pa s, viscosidad media del fluido del canal
a_w23   = 8.5e-6           # m, radio de la celula HL-60
V_w23   = 0.89e-2          # m/s
Ca_w23, C_w23 = 0.75, 10.0
Gs_w23  = 0.39e-3          # N/m, valor publicado
Ks_w23  = 8.10e-3          # N/m, valor publicado
# comprobacion de consistencia Ca = mu0*V/Gs
Gs_chk  = mu0_w23*V_w23/Ca_w23
assert abs(Gs_chk - Gs_w23)/Gs_w23 < 0.02, (Gs_chk, Gs_w23)
# viscosidades de membrana adimensionales publicadas -> dimensionales
eta_s_shear = 20.0  * mu0_w23 * a_w23      # N s / m
eta_s_elong = 100.0 * mu0_w23 * a_w23
ETA_S_SHEAR_uN = eta_s_shear*1e6           # uN s / m
ETA_S_ELONG_uN = eta_s_elong*1e6

# ============================ GEOMETRIA ======================================
add("MCF-7","geometria","radio celular (suspension)", 8.74, 8.59, 8.89, "um",
    "10.1039/c7ra10750a","analisis de imagen en microcanal confinante","MCF-7", False,
    "", "cell diameter 17.48 +/- SE 0.31 um")
add("MDA-MB-231","geometria","radio celular (suspension)", 8.62, 8.51, 8.73, "um",
    "10.1039/c7ra10750a","analisis de imagen en microcanal confinante","MDA-MB-231", False,
    "", "cell diameter 17.24 +/- SE 0.21 um")
add("SK-BR-3","geometria","radio celular (suspension)", 10.0, 8.6, 10.0, "um",
    "10.1038/s41378-020-0150-6","microscopia optica (pinzas acusticas, valor agregado del panel)",
    "MCF-7 / MDA-MB-231 / SK-BR-3 (valor conjunto)", True,
    "Unico dato de tamano que incluye SK-BR-3; se publica como media del panel completo, no por linea. "
    "El minimo se fija en el radio medido de MDA-MB-231 en el mismo rango de lineas.",
    "The average diameter of the cells and AHSs is ~20 um")
add("MCF-10A","geometria","radio celular (suspension)", 8.74, 8.51, 10.0, "um",
    "10.1039/c7ra10750a","analisis de imagen en microcanal confinante","MCF-7 (mismo dispositivo y panel)", True,
    "MCF-10A se mide en el mismo panel y dispositivo que MCF-7 y MDA-MB-231 pero su diametro "
    "no se publica numericamente (solo en figura). Se transfiere el valor de MCF-7 y el rango "
    "cubre el intervalo medido del panel.", "")
add("leucocito","geometria","radio celular (suspension)", 8.5, 8.5, 8.5, "um",
    W23,"analisis de imagen (configuracion experimental de Fregin et al. 2019)","HL-60", True,
    "Leucocito de referencia representado por HL-60, la linea sobre la que Wang et al. calibran el modelo.",
    "The cell has a radius of a = 8.5 um")

add("leucocito","geometria","radio nuclear", 4.25, 3.5, 5.0, "um",
    "10.1074/jbc.M112.441535","citado por Wang et al. 2023 (10.1007/s10237-023-01705-6)","HL-60", False,
    "", "according to Rowat et al (2013), the nucleus radius is typically 3.5-5 um", mid=True)
for cl in ["MCF-10A","MCF-7","SK-BR-3","MDA-MB-231"]:
    add(cl,"geometria","relacion de radios nucleo/celula (an/a)", 0.5, 0.4, 0.6, "-",
        W23,"ajuste del modelo multicompartimento 3D a experimento (IB-LBM)","HL-60 / PC-3 / K-562", True,
        "No se encontro ningun valor publicado de radio nuclear ni de relacion N/C para lineas de mama "
        "en suspension. Se transfiere el valor de mejor ajuste de Wang et al. (0.5), identico en HL-60 (2023) "
        "y PC-3/K-562 (2025), con el rango derivado del radio nuclear de HL-60 (3.5-5 um sobre a = 8.5 um).",
        "the nucleus has a size ratio of an/a = 0.5 ... in the range of 0.4 to 0.6")
add("leucocito","geometria","relacion de radios nucleo/celula (an/a)", 0.5, 0.4, 0.6, "-",
    W23,"ajuste del modelo multicompartimento 3D a experimento (IB-LBM)","HL-60", False, "",
    "an/a = 0.5 ... in the range of 0.4 to 0.6")

# ===================== MEMBRANA / CORTEX =====================================
add("leucocito","membrana_cortex","modulo de corte de membrana Gs", 0.39, 0.39, 0.39, "mN/m",
    W23,"ajuste IB-LBM 3D a transito en microcanal constrenido (Ca = 0.75, C = 10)","HL-60", False, "",
    "Ca = 0.75 and C = 10 correspond to cell membrane elastic moduli values of Gs = 0.39 mN m-1 and Ks = 8.10 mN m-1")
add("leucocito","membrana_cortex","modulo de dilatacion de area Ks", 8.10, 8.10, 8.10, "mN/m",
    W23,"ajuste IB-LBM 3D (ley Skalak, C = 10 con Ks = (1+2C) Gs)","HL-60", False, "",
    "Gs = 0.39 mN m-1 and Ks = 8.10 mN m-1")
add("leucocito","membrana_cortex","viscosidad de membrana en corte eta_s^s", round(ETA_S_SHEAR_uN,2),
    round(ETA_S_SHEAR_uN,2), 14.0, "uN s/m", W23,
    "ajuste IB-LBM 3D; valor adimensional eta_s^s/(mu0 a) = 20 redimensionado con mu0 = 32.5 mPa s y a = 8.5 um",
    "HL-60", False,
    "Maximo tomado del valor medido en granulocitos (8.5 uN s/m, Evans & Yeung 1989) y del rango de K-562 "
    "(7.8-14 uN s/m, Wang et al. 2025), que acotan por arriba la dispersion real entre celulas leucocitarias.",
    "the cell membrane viscosity follows Eq. (14) with eta_s^e/mu0 a = 100, eta_s^s/mu0 a = 20")
add("leucocito","membrana_cortex","viscosidad de membrana en elongacion eta_s^e", round(ETA_S_ELONG_uN,2),
    round(ETA_S_ELONG_uN,2), round(ETA_S_ELONG_uN,2), "uN s/m", W23,
    "ajuste IB-LBM 3D; valor adimensional eta_s^e/(mu0 a) = 100 redimensionado","HL-60", False,
    "", "eta_s^e/mu0 a = 100")
add("leucocito","membrana_cortex","viscosidad de membrana eta_s (granulocito)", 8.5, 8.5, 8.5, "uN s/m",
    "10.1016/S0006-3495(89)82660-8","aspiracion por micropipeta (citado por Wang et al. 2025)","granulocito humano", False,
    "", "comparable to the membrane of granulocytes that has eta_s = 8.5 uN s m-1 (Evans and Yeung 1989)")
add("leucocito","membrana_cortex","tension cortical persistente", 0.035, 0.035, 0.035, "mN/m",
    "10.1016/S0006-3495(89)82660-8","aspiracion por micropipeta","granulocito humano", False, "",
    "the outer cortex of the cell maintains a small persistent tension of approximately 0.035 dyn/cm")
add("leucocito","membrana_cortex","reserva de area de membrana (A_max/A_0)", 2.15, 2.1, 2.2, "-",
    "10.1016/S0006-3495(89)82660-8","aspiracion por micropipeta","granulocito humano", False, "",
    "final areas were measured to be 2.1 to 2.2 times the area of the initial spherical shape", mid=True)

# Lineas de mama: transferencia obligada desde PC-3 (unica linea de tumor solido calibrada
# en el mismo modelo) -- Wang et al. 2025, Tabla 1.
PC3_Gs = (5.8, 20.7); PC3_Ks = (120.3, 434.0); PC3_eta = (8.7, 22.5)
for cl in ["MCF-10A","MCF-7","SK-BR-3","MDA-MB-231"]:
    just_pc3 = ("No existe ninguna medida de modulo de membrana/cortex para lineas de mama. PC-3 (prostata) "
                "es la unica linea de TUMOR SOLIDO epitelial parametrizada dentro de este mismo modelo "
                "multicompartimento, por lo que es la transferencia mas razonable; el rango une las dos "
                "geometrias de canal publicadas (constrenido y cross-slot).")
    add(cl,"membrana_cortex","modulo de corte de membrana Gs", round(np.mean(PC3_Gs),2), PC3_Gs[0], PC3_Gs[1],
        "mN/m", W25,"inferencia por ajuste IB-LBM 3D a transito en dos geometrias de microcanal","PC-3", True,
        just_pc3, "Gs (mN/m) PC-3: 10.4-20.7 (constricted) / 5.8-14.4 (cross-slot)", mid=True)
    add(cl,"membrana_cortex","modulo de dilatacion de area Ks", round(np.mean(PC3_Ks),1), PC3_Ks[0], PC3_Ks[1],
        "mN/m", W25,"inferencia por ajuste IB-LBM 3D","PC-3", True, just_pc3,
        "Ks (mN/m) PC-3: 217-434 / 120.3-302.6", mid=True)
    add(cl,"membrana_cortex","viscosidad de membrana eta_s", round(np.mean(PC3_eta),2), PC3_eta[0], PC3_eta[1],
        "uN s/m", W25,"inferencia por ajuste IB-LBM 3D","PC-3", True, just_pc3,
        "eta_s (uN s m-1) PC-3: 9.7-22.5 / 8.7-16.1", mid=True)
# K-562 como cota de linea hematologica (contexto para el leucocito)
add("leucocito","membrana_cortex","modulo de corte de membrana Gs (linea hematologica K-562)",
    4.4, 3.6, 5.2, "mN/m", W25,"inferencia por ajuste IB-LBM 3D","K-562", True,
    "K-562 es la unica linea hematologica con Gs inferido en este modelo; se incluye como cota superior "
    "independiente para el leucocito frente al valor de HL-60 (0.39 mN/m), que es 10 veces menor.",
    "Gs (mN/m) K-562: 4.4-5.2 / 3.6-4.8", mid=True)

# ============================ CITOPLASMA =====================================
CYT = "10.1039/d4na00003j"
cyt_tec = "espectroscopia rotacional magnetica (MRS) con microhilos intracelulares"
for cl, eta, deta, G, dG, tau, dtau, eta1 in [
    ("MCF-10A", 41.6, 8.7, 79.3, 7.3, 0.77, 0.37, 36.3),
    ("MCF-7",   56.4,16.6, 32.9, 6.0, 3.23, 0.53, 65.9),
    ("MDA-MB-231",10.7,5.4, 38.6, 5.8, 0.57, 0.20, 12.0)]:
    add(cl,"citoplasma","viscosidad estatica del citoplasma", eta, round(eta-deta,1), round(eta+deta,1),
        "Pa s", CYT, cyt_tec, cl, False,
        f"Segunda determinacion (Tabla 1). La primera determinacion del mismo trabajo da {eta1} Pa s; "
        "la dispersion real con el tamano de sonda es mucho mayor (ver fila de rango de sensibilidad).",
        f"eta_{cl} = {eta} +/- {deta} Pa s")
    add(cl,"citoplasma","modulo elastico del citoplasma", G, round(G-dG,1), round(G+dG,1), "Pa",
        CYT, cyt_tec, cl, False, "", f"G = {G} +/- {dG} Pa")
    add(cl,"citoplasma","tiempo de relajacion del citoplasma", tau, round(tau-dtau,2), round(tau+dtau,2), "s",
        CYT, cyt_tec, cl, False, "", f"tau = {tau} +/- {dtau} s")
# rango de sensibilidad por dependencia con el tamano de sonda (dispersion real medida)
add("MCF-10A","citoplasma","viscosidad del citoplasma (rango por escala de sonda)", 41.6, 15.0, 160.0,
    "Pa s", CYT, cyt_tec+"; barrido de longitud de hilo L = 2-6 um","MCF-10A", False,
    "", "as L increases from 2 um to 6 um, the MCF-10A viscosity varies from 15 to 160 Pa s")
add("MDA-MB-231","citoplasma","viscosidad del citoplasma (rango por escala de sonda)", 10.7, 3.0, 25.0,
    "Pa s", CYT, cyt_tec+"; barrido de longitud de hilo L = 2-6 um","MDA-MB-231", False,
    "", "for MDA-MB-231 the variation is from 3 Pa s to 25 Pa s over the same interval")
add("SK-BR-3","citoplasma","viscosidad estatica del citoplasma", 41.6, 10.7, 56.4, "Pa s",
    CYT, cyt_tec, "MCF-10A (nominal) / panel MCF-10A-MCF-7-MDA-MB-231 (rango)", True,
    "SK-BR-3 no aparece en ningun estudio de microrreologia intracelular. Se transfiere el valor de MCF-10A "
    "(epitelial, no mesenquimal, como SK-BR-3) y el rango abarca las tres lineas medidas en el mismo trabajo.", "")
add("SK-BR-3","citoplasma","modulo elastico del citoplasma", 79.3, 32.9, 79.3, "Pa",
    CYT, cyt_tec, "MCF-10A (nominal) / panel de tres lineas (rango)", True,
    "Misma transferencia que la viscosidad: no hay dato de SK-BR-3.", "")
add("leucocito","citoplasma","viscosidad aparente del interior celular (37 C)", 100.0, 100.0, 200.0,
    "Pa s", "10.1016/S0006-3495(89)82660-8","aspiracion por micropipeta, modelo de nucleo liquido-corteza",
    "granulocito humano", False,
    "Minimo = valor a 37 C, maximo = valor a 23 C (dependencia termica publicada).",
    "values on the order of 2 x 10(3) poise at 23 degrees C, lower values of 1 x 10(3) poise at 37 degrees C")
add("leucocito","citoplasma","viscosidad efectiva de citoplasma usada en el modelo continuo", 32.5, 10.0, 100.0,
    "mPa s", W23,"ajuste IB-LBM 3D con razon de viscosidades lambda = 1 (rango de literatura citado)",
    "HL-60", False,
    "ATENCION: esta viscosidad 'efectiva' de modelo continuo es 3 ordenes de magnitud menor que la "
    "viscosidad microrreologica local medida en las lineas de mama (10-60 Pa s). No son la misma magnitud "
    "fisica y no deben mezclarse.",
    "the effective cytoplasmic viscosity of biological cells ranges roughly from 10 to 100 mPa s")
# Referencia externa de celulas tumorales por micropipeta (cuartiles publicados)
add("referencia_tumoral_no_mama","citoplasma","viscosidad citoplasmatica (cuartiles)", 42.1, 16.7, 110.3,
    "Pa s","10.1098/rsos.181707","aspiracion continua por micropipeta a -10 kPa (n = 652)","H1299 (pulmon)", False,
    "Se incluye como contexto de dispersion entre lineas tumorales, no como parametro de simulacion.",
    "16.7 Pa s, 42.1 Pa s, 110.3 Pa s and 74% for H1299 cells at -10 kPa")
add("referencia_tumoral_no_mama","citoplasma","viscosidad citoplasmatica (cuartiles)", 489.8, 144.8, 1390.7,
    "Pa s","10.1098/rsos.181707","aspiracion continua por micropipeta a -10 kPa (n = 785)","A549 (pulmon)", False,
    "Idem: muestra que la dispersion entre lineas tumorales abarca mas de un orden de magnitud.",
    "144.8 Pa s, 489.8 Pa s, 1390.7 Pa s, and 81% for A549 cells at -10 kPa")

# ============================== NUCLEO =======================================
for cl in ["MCF-10A","MCF-7","SK-BR-3","MDA-MB-231","leucocito"]:
    add(cl,"nucleo","razon de modulos capsula nuclear / membrana celular (Gsn/Gs)", 2.0, 2.0, 8.0, "-",
        W23,"ajuste IB-LBM 3D (hipotesis del modelo, verificada por sensibilidad)","HL-60 / PC-3 / K-562", True,
        "Valor de referencia del modelo de Wang et al., identico en 2023 y 2025. El maximo (8) es el valor "
        "maximo explorado por los autores, que no reduce significativamente la deformacion celular.",
        "The elastic moduli of the cell nucleus membrane are also assumed to be twice those of the cell membrane")
add("MCF-10A","nucleo","modulo elastico de la red de lamina nuclear", 25.0, 25.0, 25.0, "mN/m",
    "10.1242/jcs.01357","aspiracion por micropipeta de envolturas nucleares aisladas","ovocito de Xenopus", True,
    "Unica medida directa del modulo elastico de la envoltura nuclear como red bidimensional, en la misma "
    "magnitud (N/m) que los modulos de membrana del modelo. Transferencia obligada: no existe equivalente "
    "en celulas de mama ni en celulas humanas.",
    "yields a network elastic modulus, unaffected by nucleoplasm, that averages 25 mN/m")
for cl in ["MCF-7","SK-BR-3","MDA-MB-231","leucocito"]:
    add(cl,"nucleo","modulo elastico de la red de lamina nuclear", 25.0, 25.0, 25.0, "mN/m",
        "10.1242/jcs.01357","aspiracion por micropipeta de envolturas nucleares aisladas","ovocito de Xenopus", True,
        "Igual que para MCF-10A: unica medida directa disponible de la lamina como red 2D.",
        "averages 25 mN/m")
add("MCF-10A","nucleo","constante elastica nuclear (extension corta, <3 um)", 0.65, 0.41, 0.70, "nN/um",
    "10.1091/mbc.e16-09-0653","micromanipulacion con micropipetas de nucleos aislados","HeLa y MEF vimentina-nula", True,
    "Escala de rigidez nuclear a deformacion pequena medida sobre nucleos aislados humanos (HeLa). "
    "El minimo corresponde a nucleos con cromatina descompactada en el mismo trabajo.",
    "from 0.60 to 0.41 nN/um in MEF V-/- and from 0.70 to 0.47 nN/um in HeLa")
for cl in ["MCF-7","SK-BR-3","MDA-MB-231","leucocito"]:
    add(cl,"nucleo","constante elastica nuclear (extension corta, <3 um)", 0.65, 0.41, 0.70, "nN/um",
        "10.1091/mbc.e16-09-0653","micromanipulacion de nucleos aislados","HeLa y MEF vimentina-nula", True,
        "Misma transferencia que para MCF-10A.", "0.60 to 0.41 nN/um ... 0.70 to 0.47 nN/um")
for cl in ["MCF-10A","MCF-7","SK-BR-3","MDA-MB-231","leucocito"]:
    add(cl,"nucleo","endurecimiento por deformacion (factor, extension >3 um)", 2.0, 1.5, 2.5, "-",
        "10.1091/mbc.e16-09-0653","micromanipulacion de nucleos aislados","HeLa, BJ5-ta, MEF vimentina-nula", True,
        "Justifica usar una ley hiperelastica (tipo Skalak) para la capsula nuclear en vez de un solido lineal. "
        "Transferido: no medido en lineas de mama.",
        "Strain stiffening was consistently observed as a 1.5- to 2.5-fold increase", mid=True)
add("MCF-10A","nucleo","razon de rigidez nucleo/citoplasma (in situ)", 1.4, 1.4, 4.0, "-",
    "10.1016/j.jbiomech.2009.01.024","MEF de celula con inclusion nuclear ajustado a compresion no confinada",
    "condrocito articular", True,
    "Minimo = razon in situ (1.4, E_cito = 3.75 kPa y E_nuc = 5.25 kPa); maximo = razon sobre nucleos AISLADOS "
    "del mismo tipo celular (3-4 veces mas rigidos, 10.1006/bbrc.2000.2360). La horquilla 1.4-4 es la "
    "discrepancia real publicada entre nucleo in situ y nucleo aislado.",
    "the ratio of cytoplasmic to nuclear stiffness was ... 1.4 (E_Cyto = 3.5 kPa, E_Nuc = 5 kPa)")
for cl in ["MCF-7","SK-BR-3","MDA-MB-231","leucocito"]:
    add(cl,"nucleo","razon de rigidez nucleo/citoplasma (in situ)", 1.4, 1.4, 4.0, "-",
        "10.1016/j.jbiomech.2009.01.024","MEF ajustado a compresion no confinada; cota superior de 10.1006/bbrc.2000.2360",
        "condrocito articular", True, "Misma transferencia que para MCF-10A.",
        "isolated nuclei ... were 3-4 times stiffer and nearly twice as viscous as the cytoplasm")
add("MCF-10A","nucleo","razon de viscosidad nucleo/citoplasma", 2.0, 2.0, 2.0, "-",
    "10.1006/bbrc.2000.2360","aspiracion por micropipeta de nucleos aislados","condrocito articular", True,
    "Unico dato encontrado que compara viscosidad nuclear y citoplasmatica del mismo tipo celular.",
    "nearly twice as viscous as the cytoplasm")
for cl in ["MCF-7","SK-BR-3","MDA-MB-231","leucocito"]:
    add(cl,"nucleo","razon de viscosidad nucleo/citoplasma", 2.0, 2.0, 2.0, "-",
        "10.1006/bbrc.2000.2360","aspiracion por micropipeta de nucleos aislados","condrocito articular", True,
        "Misma transferencia que para MCF-10A.", "nearly twice as viscous as the cytoplasm")
add("MDA-MB-231","nucleo","rigidez nuclear por AFM (fenotipo mesenquimal vs epitelial)", 682.6, 682.6, 1278.0,
    "Pa","10.1038/s41598-026-56141-0","mapeo de fuerza por AFM","U251 (glioblastoma)", True,
    "Dato direccional, no absoluto: cuantifica cuanto BAJA la rigidez nuclear en una transicion a fenotipo "
    "mesenquimal (1278 -> 682.6 Pa, factor 0.53). Se usa para justificar que MDA-MB-231 (mesenquimal) tenga "
    "un nucleo mas blando que MCF-7/SK-BR-3 (epiteliales) en el analisis de sensibilidad, no como valor de entrada.",
    "a reduction in nuclear stiffness, especially in U251 cells (1278 Pa -> 682.6 Pa)")

# ======================= CELULA ENTERA (calibracion) =========================
add("MDA-MB-231","celula_entera","modulo de Young aparente", 1.527, 1.527, 1.527, "kPa",
    "10.1038/s41378-020-0150-6","pinzas acusticas de haz unico (SBAT) con esferas de agarosa patron","MDA-MB-231", False,
    "", "1.527 kPa for MDA-MB-231")
add("MCF-7","celula_entera","modulo de Young aparente", 2.650, 2.650, 2.650, "kPa",
    "10.1038/s41378-020-0150-6","pinzas acusticas de haz unico (SBAT)","MCF-7", False, "", "2.650 kPa for MCF-7")
add("SK-BR-3","celula_entera","modulo de Young aparente", 2.772, 2.772, 2.772, "kPa",
    "10.1038/s41378-020-0150-6","pinzas acusticas de haz unico (SBAT)","SK-BR-3", False,
    "UNICO parametro mecanico de SK-BR-3 encontrado en toda la busqueda.", "2.772 kPa for SKBR-3")
add("leucocito","celula_entera","modulo de Young aparente (granulocito)", 0.43, 0.38, 0.48, "kPa",
    "10.1038/s41467-019-08370-3","citometria de deformabilidad en tiempo real dinamica (dRT-DC)","granulocito humano", False,
    "", "granulocytes (E = 0.43 +/- 0.05 kPa)")
add("leucocito","celula_entera","modulo de Young aparente (PBMC)", 0.62, 0.58, 0.66, "kPa",
    "10.1038/s41467-019-08370-3","dRT-DC","celula mononuclear de sangre periferica", False, "",
    "PBMCs (E = 0.62 +/- 0.04 kPa)")
add("leucocito","celula_entera","modulo de Young aparente (HL-60)", 0.40, 0.38, 0.44, "kPa",
    "10.1038/s41467-019-08370-3","dRT-DC","HL-60", False,
    "Maximo = control con DMSO del mismo trabajo (0.43 +/- 0.01 kPa).",
    "an apparent Young's modulus of E = 0.40 +/- 0.02 kPa and E = 0.43 +/- 0.01 kPa")
add("leucocito","celula_entera","viscosidad aparente (HL-60, Kelvin-Voigt)", 3.43, 3.27, 3.72, "Pa s",
    "10.1038/s41467-019-08370-3","dRT-DC, tiempo caracteristico tau = eta/E","HL-60", False,
    "Maximo = control con DMSO (3.67 +/- 0.05 Pa s).",
    "eta = 3.43 +/- 0.16 Pa s for the wildtype to eta = 3.67 +/- 0.05 Pa s for the DMSO control")

df = pd.DataFrame(R)
df.to_csv("parametros_compartimento.csv", index=False)
print("filas:", len(df), "| transferidas:", int(df.es_transferido.sum()))
print(df.groupby(["tipo_celular","compartimento"]).size().to_string())
print("\neta_s^s HL-60 =", round(ETA_S_SHEAR_uN,3), "uN s/m ; eta_s^e =", round(ETA_S_ELONG_uN,3),
      "uN s/m ; Gs comprobado =", round(Gs_chk*1e3,3), "mN/m")
