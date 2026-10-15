# Método

## 1. Población celular

Cada subpoblación tiene seis parámetros con distribución marginal truncada (`data/poblaciones.yaml`):
diámetro en suspensión D, razón núcleo/célula a_n/a, módulo aparente E_app, tensión cortical γ,
viscosidad η y tiempo de relajación τ. Se muestrean por hipercubo latino con una semilla derivada de
(semilla global, nombre de población). Los parámetros son **independientes** (declarado): no hay
medidas conjuntas de tamaño y rigidez por subpoblación.

Derivadas por célula (manuscrito, Sect. 2.5): G_s = E_app·a/K_E con K_E = 2.608, y
G_cyt = Π·G_s/a con Π = 0.10 (regularización declarada).

| subpoblación | D (µm) | a_n/a | E_app (kPa) | γ (pN/µm) | estado principal |
|---|---|---|---|---|---|
| neutrófilo | 8.40 ± 0.60 | 0.625 | 0.43 (granulocito) | 37.0 ± 4.8 | medido |
| linfocito | 6.05 | 0.763 | 0.62 (PBMC) | 74.7 ± 9.8 | medido (sd de D transferida) |
| monocito | 8.79 ± 0.57 | 0.638 | 0.62 (PBMC) | transferida | E_app y γ transferidos |
| eosinófilo, basófilo | = neutrófilo | | | | transferido |

Fuentes: Downey et al. 1990 (Coulter); Rosenbluth et al. 2006; Schmid-Schönbein et al. 1980;
Fregin et al. 2019 (dRT-DC); Guo et al. 2012, Lab Chip (tensión cortical de neutrófilo y linfocito con
la misma técnica); Evans & Yeung 1989; Chávez-González et al. 2023 (diferencial). Bufi et al. 2015
(microplacas) indican que el monocito podría ser más rígido que el linfocito T; no se mezcla con
dRT-DC por ser otra técnica.

**Huecos declarados:** dispersión célula a célula de E_app (CV 25 % declarado), tensión cortical y
viscosidad del monocito, viscosidad del linfocito en reposo, CV del diámetro tumoral (10 %).

## 2. Dispositivos

Perfil común `R_w(z) = R_n + (R_ext − R_n)·h(z)` con apertura exacta R_n (ver `geometria.py`).

| dispositivo | z_t (µm) | L (µm) | R_ext (µm) | E pared |
|---|---|---|---|---|
| membrana PC (track-etched) | 5 (espesor 10) | 1.0 (borde) | 25 | 2.0–2.6 GPa |
| membrana PET (track-etched) | 6 (espesor 12) | 1.0 | 25 | 2.0–2.7 GPa |
| membrana PDMS (moldeada) | 10 (espesor 20) | 2.0 | 25 | 1.32–2.97 MPa |
| canal PDMS | 5 | 5 | 15 | 1.32–2.97 MPa |
| canal PMMA | 5 | 10 | 15 | 1.8–3.4 GPa |

Las membranas de PC track-etched existen comercialmente con poros de 5, 8, 10 y 12 µm (Whatman
Nuclepore). La densidad de poros de PET 8 µm es 6×10⁴ /cm² (Millicell). Espesores y radios de borde
son **declarados**: el fabricante no los publica en las fichas recuperadas.

**Efecto del material:** entra por la geometría (espesor, borde, transición) y por la rigidez de
pared (Winkler en serie con la penalización). Comprobado con MCF-7 en canal a 12 µm: PMMA, PC y PET
son indistinguibles de una pared rígida (< 0.01 %); el PDMS reduce W un 3.5 % y la presión crítica
un 4 %. No se modelan adhesión, fricción ni deformación del canal bajo presión.

## 3. FEM

Solver del manuscrito sin cambios (Apéndice A). Una continuación por (célula, dispositivo) cierra la
apertura desde ρ = 1.02a hasta la apertura más pequeña pedida con a/ρ ≤ `conf_max_intento`,
registrando W(ρ) en cada paso aceptado y métricas completas en cada apertura objetivo.

## 4. Presión crítica

* `dP_local`: reducción axial con W indexado por el radio real del cuello (canales).
* `dP_trabajo` = W(R_n)/(πR_n²·l_factor·a) (membranas).
* `dP_gota` = 2γ(1/R_n − 1/a): cota de gota líquida, independiente del FEM, para leucocitos.

## 5. Análisis

Ventana de presión por dispositivo y apertura: [P99 de la mezcla leucocitaria ponderada por el
diferencial, P1 de cada línea tumoral]; presión de operación en la media geométrica; rendimiento
tumoral y arrastre leucocitario por millón. Umbral de dispositivo Kelvin–Voigt β_c = ΔP·τ/L_carga.

## 6. Verificación incluida

`make pruebas`: patch test del solver (error < 10⁻¹⁰), apertura exacta y derivada nula del perfil,
reproducibilidad del muestreo, calibración del córtex del manuscrito (G_s = 8.881 mN/m), barrera
finita con contacto leve, marca de cota inferior, cota de Laplace.
