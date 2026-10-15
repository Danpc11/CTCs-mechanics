# DLD rehecho con radios de leucocito primarios — de resultado negativo a criterio complementario

**Qué cambia.** La Fase 2 descartó DLD usando un leucocito de **17,00 µm de diámetro transferido de
HL-60**, que es una línea de leucemia promielocítica, no un leucocito circulante. Con los diámetros
**medidos en suspensión** de las tres subpoblaciones reales el descarte se cae por completo: la
brecha de medias pasa del 1,4–2,8 % al **49–65 %** y el índice de Youden del 0,38–0,56 al **0,9999**.

**Qué NO cambia.** El criterio mecánico sigue teniendo sentido, pero por una razón distinta y mejor:
DLD separa **población** (tumoral frente a fondo leucocitario) y lo hace bien; lo que no separa es
**fenotipo dentro de la población tumoral**, ni resuelve las dos subpoblaciones que sí solapan en
tamaño (blastos leucémicos circulantes y el extremo grande del monocito). El artículo deja de
competir con ISET/Parsortix y se sitúa **aguas abajo** de ellos.

---

## 1. Entradas

### 1.1 Diámetros de leucocito adoptados (en suspensión, medidos)

| subpoblación | D (µm) | dispersión usada | técnica | fuente |
|---|---|---|---|---|
| linfocito  | **6,05** | CV 7,2 % *(transferido del neutrófilo)* | citometría de flujo, extinción axial de luz | 10.1152/jappl.1990.69.5.1767 |
| neutrófilo | **8,40** | CV 7,2 % (SD 0,6 µm sobre 8,3 µm) | contador Coulter / microscopía de fase | 10.1152/jappl.1990.69.5.1767 · 10.1529/biophysj.105.067496 |
| monocito   | **8,79** | CV 6,5 % (rango Coulter 7,72–9,99 → SD ≈ rango/4) | contador Coulter | 10.1152/jappl.1990.69.5.1767 |
| eosinófilo | 8,40 *(transferido del neutrófilo)* | CV 7,2 % | — | — |
| basófilo   | 8,40 *(transferido del neutrófilo)* | CV 7,2 % | — | — |

Tres técnicas independientes concuerdan en el neutrófilo: Coulter 8,4 µm (rango 7,3–9,7),
citometría de flujo 8,3 µm (rango 7,14–10,25) y microscopía de fase 8,3 ± 0,6 µm. El volumen por
micropipeta de 299 ± 32 µm³ equivale a 8,30 ± 0,30 µm, un cuarto acuerdo. **Ninguna de estas
medidas se acerca a los 17 µm de HL-60.**

*Huecos declarados:* no se recuperó dispersión publicada del linfocito (se le transfiere el CV del
neutrófilo, que es la hipótesis **conservadora** porque maximiza el solapamiento); las cifras de
Downey 1990, Schmid-Schönbein 1980 y Ting-Beall 1993 se recuperaron a través de BioNumbers citando
la tabla primaria, porque los PDF originales están de pago y no son accesibles en este entorno.

### 1.2 Diferencial leucocitario

Intervalos de referencia de donantes sanos (varones, N = 156, percentiles 2,5–97,5;
DOI 10.11613/BM.2023.020701): neutrófilos 38,9–68,0 %, linfocitos 21–48 %, monocitos 5–11 %,
eosinófilos 0–4 %, basófilos 0–1 %. La fracción central usada para ponderar se deriva del punto
medio del intervalo de recuentos **absolutos** normalizado al total: **neutrófilo 0,5635,
linfocito 0,3286, monocito 0,0804, eosinófilo 0,0209, basófilo 0,0066**. Es una derivación, no un
valor publicado como tal, y se declara así.

### 1.3 Diámetros tumorales

MCF-7 17,48 µm y MDA-MB-231 17,24 µm (medidos, 10.1039/c7ra10750a); SK-BR-3 20,00 µm y MCF-10A
17,48 µm (transferidos). La semidispersión publicada es un **error estándar**, no una desviación
típica, así que no sirve como σ poblacional: se adopta **CV = 10 % declarado** y se barre de 5 % a
30 %.

---

## 2. Vía 1 — brecha de medias frente a la exactitud de D_c

Correlación de Davis `D_c = 1,4·g·ε^0,48`, exactitud publicada 10 % de `D_c`.

| tumoral | leucocito | brecha (µm) | brecha relativa | resolución de D_c (µm) | veces **sobre** resolución | ¿separable? |
|---|---|---|---|---|---|---|
| MCF-7 | linfocito | 11.43 | 65.4 % | 0.83 | **13.7×** | **sí** |
| MCF-7 | neutrofilo | 9.08 | 51.9 % | 1.07 | **8.5×** | **sí** |
| MCF-7 | monocito | 8.69 | 49.7 % | 1.09 | **8.0×** | **sí** |
| MCF-7 | eosinofilo | 9.08 | 51.9 % | 1.07 | **8.5×** | **sí** |
| MCF-7 | basofilo | 9.08 | 51.9 % | 1.07 | **8.5×** | **sí** |
| MDA-MB-231 | linfocito | 11.19 | 64.9 % | 0.83 | **13.5×** | **sí** |
| MDA-MB-231 | neutrofilo | 8.84 | 51.3 % | 1.07 | **8.3×** | **sí** |
| MDA-MB-231 | monocito | 8.45 | 49.0 % | 1.09 | **7.8×** | **sí** |
| MDA-MB-231 | eosinofilo | 8.84 | 51.3 % | 1.07 | **8.3×** | **sí** |
| MDA-MB-231 | basofilo | 8.84 | 51.3 % | 1.07 | **8.3×** | **sí** |

Compárese con la tabla de la Fase 2, donde las mismas dos líneas quedaban **3,6 y 7,2 veces
POR DEBAJO** de la resolución. La inversión es de un factor 30–100 y su única causa es haber
sustituido HL-60 por el leucocito real.

## 3. Vía 2 — solapamiento de distribuciones y Youden, ponderado por el diferencial

| tumoral | D medio de la mezcla leucocitaria (µm) | solapamiento OVL | Youden máximo | D_c óptimo (µm) | sensibilidad | especificidad |
|---|---|---|---|---|---|---|
| MCF-7 | 7.66 | 1.02e-04 | **0.99990** | 10.86 | 0.99992 | 0.99997 |
| MDA-MB-231 | 7.66 | 1.32e-04 | **0.99987** | 10.82 | 0.99990 | 0.99997 |
| SK-BR-3 *(D transferido)* | 7.66 | 7.00e-06 | **0.99999** | 11.23 | 0.99999 | 1.00000 |
| MCF-10A *(D transferido)* | 7.66 | 1.02e-04 | **0.99990** | 10.86 | 0.99992 | 0.99997 |

Con los CV medidos el solapamiento de la mezcla leucocitaria completa con MCF-7 es de
**1,0 × 10⁻⁴**, y haría falta **multiplicar por 3,10 el CV leucocitario** (es decir, llevar el CV
del neutrófilo del 7,2 % al 22,3 %) para que el solapamiento llegue siquiera al **1 %**, y por 7,56
(CV 54,4 %) para que llegue al 10 %. El 50 % es inalcanzable con cualquier dispersión gaussiana.
La Fase 2 afirmaba lo contrario —que bastaba un CV del 1–2 %— porque medía el solapamiento de dos
poblaciones **cuyas medias distaban 0,24–0,48 µm**.

## 4. Lo que la separación por tamaño cuesta en pureza, y por qué hace falta el criterio mecánico

El índice de Youden no es la magnitud de diseño: con 1 CTC por 10⁶ leucocitos lo que decide es el
**arrastre leucocitario** en el valor de `D_c` que da el rendimiento tumoral exigido.

| CV tumoral supuesto | rendimiento objetivo | D_c (µm) | arrastre leucocitario | leucocitos por 10⁶ | pureza | enriquecimiento |
|---|---|---|---|---|---|---|
| 5 % | 95 % | 16.04 | 0.00e+00 | 0.0 | 1 | ∞ |
| 5 % | 99 % | 15.45 | 0.00e+00 | 0.0 | 1 | ∞ |
| 10 % | 95 % | 14.61 | 0.00e+00 | 0.0 | 1 | ∞ |
| 10 % | 99 % | 13.41 | 0.00e+00 | 0.0 | 1 | ∞ |
| 15 % | 95 % | 13.17 | 1.78e-15 | 0.0 | 1 | 1e+06× |
| 15 % | 99 % | 11.38 | 4.79e-07 | 0.5 | 0.674 | 6.74e+05× |
| 20 % | 95 % | 11.73 | 2.16e-08 | 0.0 | 0.978 | 9.78e+05× |
| 20 % | 99 % | 9.35 | 4.79e-02 | 47927.4 | 2.07e-05 | 20.7× |
| 30 % | 95 % | 8.85 | 1.70e-01 | 170298.1 | 5.58e-06 | 5.6× |
| 30 % | 99 % | 5.28 | 9.87e-01 | 987290.0 | 1e-06 | 1× |

**Lectura.** Mientras la dispersión de tamaño de la población tumoral se mantenga en un CV ≤ 15 %,
DLD por sí solo deja el arrastre leucocitario por debajo de 2 células por millón y la pureza por
encima del 36 %. A **CV = 20 % con 99 % de rendimiento** el arrastre salta a 7,0 × 10⁴ por millón
y la pureza se desploma a 1,4 × 10⁻⁵ (enriquecimiento 14×). Es decir: **el punto de fallo de DLD no
es el tamaño del leucocito sino la cola pequeña de la propia CTC**, y el criterio mecánico opera
justamente ahí, porque el contraste de módulo aparente no depende del radio (§ Sobol).
La cifra de CV tumoral es un parámetro **declarado**, no medido: es el hueco que esta vía deja
abierto y el que habría que cerrar experimentalmente.

## 5. Qué subpoblación sigue solapando — la que justifica el criterio mecánico

| tumoral | subpoblación | D (µm) | brecha relativa | solapamiento OVL | Youden máximo | ¿separable por D_c? |
|---|---|---|---|---|---|---|
| MCF-7 | Jurkat (blasto linfoide T) | 11.50 | 34.2 % | **0.065** | **0.9346** | sí |
| MDA-MB-231 | Jurkat (blasto linfoide T) | 11.50 | 33.3 % | **0.075** | **0.9252** | sí |
| MCF-7 | HL-60 (blasto mieloide) | 17.00 | 2.8 % | **0.796** | **0.1771** | **NO** |
| MDA-MB-231 | HL-60 (blasto mieloide) | 17.00 | 1.4 % | **0.827** | **0.1254** | **NO** |
| MCF-7 | monocito grande (p97.5 del rango Coulter) | 9.99 | 42.9 % | **0.002** | **0.9984** | sí |
| MDA-MB-231 | monocito grande (p97.5 del rango Coulter) | 9.99 | 42.0 % | **0.002** | **0.9980** | sí |

- **El monocito grande no es el problema.** Incluso en el percentil 97,5 del rango Coulter
  (9,99 µm) el solapamiento con MCF-7 es del 0,16 % y Youden 0,998.
- **El blasto linfoide sí empieza a serlo**: Jurkat, 11,5 ± 1,5 µm, da 6,5 % de solapamiento con
  MCF-7 y Youden 0,935 — todavía separable, pero ya no con holgura.
- **El blasto mieloide es exactamente el caso que la Fase 2 analizó sin saberlo.** HL-60, 17,0 µm,
  solapa el **79,6 %** con MCF-7 y hunde Youden a **0,177**: no es separable por `D_c`. Es decir,
  la sección 3.5 del manuscrito no demostraba que DLD no separa CTC de leucocitos; demostraba que
  **DLD no separa CTC de blastos de leucemia mieloide**, que es un enunciado verdadero,
  interesante y completamente distinto.

Esto da la frase exacta para el artículo: la selección por tamaño resuelve el fondo leucocitario
**normal** con holgura, y falla en (i) pacientes con blastos circulantes y (ii) la cola de CTC
pequeñas. El criterio mecánico es complementario porque su contraste —el cociente de módulos
aparentes— es, según el análisis de Sobol, **independiente del radio**.

---

## 6. Reposicionamiento: de competidor a complemento

| | antes (Fase 2) | ahora |
|---|---|---|
| afirmación | «DLD no puede discriminar CTC de leucocitos» | «DLD discrimina la población con Youden 0,9999; no discrimina fenotipo ni blastos» |
| leucocito de referencia | HL-60, 17,0 µm (línea de leucemia) | linfocito 6,05 / neutrófilo 8,40 / monocito 8,79 µm (medidos) |
| brecha de medias | 1,4–2,8 % | 49–65 % |
| Youden máximo | 0,381–0,561 | 0,9999 |
| CV para 50 % de solapamiento | 1,0–2,1 % | inalcanzable (1 % exige CV leucocitario del 22 %) |
| posición frente a ISET/Parsortix | rival al que se refuta | etapa previa sobre la que se construye |
| qué aporta el criterio mecánico | sustituir al tamaño | separar fenotipo dentro de la población ya enriquecida, y cubrir blastos y CTC pequeñas |

**Cifras del manuscrito que quedan invalidadas:** las de la sección 3.5 y las dos frases del
resumen y del título que se apoyan en ellas (ver `correcciones_fase6.md`).
