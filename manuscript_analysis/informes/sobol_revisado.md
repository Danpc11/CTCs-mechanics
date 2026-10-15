# Sobol remuestreado por población — la robustez anterior era cancelación algebraica

## 0. Declaración obligatoria

**El análisis de Sobol de la Fase 2 medía una cancelación estructural, no una insensibilidad.**
Con `W = Gs·a²·f(c, a_n/a, C, Π, r_nuc, r_env)` y **la misma `f` evaluada para las dos poblaciones**,
los parámetros de forma entran idénticos en numerador y denominador del cociente
`ΔP_t/ΔP_l` y **se cancelan exactamente**. Los índices totales de **1,3 × 10⁻²⁰** (constricción) y
**2,4 × 10⁻³¹** (microfiltro) que se reportaron no son una medida de robustez: son el residuo en
coma flotante de una resta de cantidades iguales. Cualquier conclusión del tipo «la selectividad no
depende de ningún parámetro transferido» obtenida con ese montaje es **vacua por construcción**.

La cancelación se rompe en cuanto las dos poblaciones reciben valores distintos, y los datos exigen
que los reciban: la `Π` medida del leucocito es **0,162** frente a **0,032–0,078** de las líneas de
mama (factor 2,1–5,1), y su razón núcleo/célula medida es **0,625–0,763** frente al **0,50**
transferido que se usaba para todas.

## 1. Montaje nuevo

Sustituto `f = W/(Gs·a²)` reajustado sobre un diseño FEM **ampliado**: hipercubo latino con
`a_n/a ∈ [0,40, 0,90]` (el original llegaba a 0,60 y no cubría al linfocito), unión de dos barridos,
**640 puntos de los que 473 convergen (74 %)**. Calidad del sustituto:
**R² = 0,9979, error mediano en hold-out 3,5 %, p90 15,6 %, máximo 27,3 %** — comparable al original.
Cobertura efectiva en la esquina leucocitaria: 88 puntos con `a_n/a ∈ [0,62, 0,77]` y 103 con
`Π ∈ [0,11, 0,22]`.

**17 parámetros con sorteo independiente por población** (los de forma, por duplicado), 2048 muestras
base → **38 912 evaluaciones**. Rangos:

| parámetro | rango | procedencia |
|---|---|---|
| `an_a_tum` | 0,40 – 0,60 | transferido (Wang 2023/2025) |
| **`an_a_leu`** | **0,625 – 0,763** | **MEDIDO**: neutrófilo 0,625 (10.1529/biophysj.105.067496), linfocito 0,763 (10.1182/blood.V56.5.866.866) |
| `C_tum`, `C_leu` | 2,40 – 37,0 cada uno | transferido, sorteo independiente |
| **`Pi_tum`** | **0,032 – 0,078** | **MEDIDO** (microrreología rotacional, 10.1039/d4na00003j) |
| **`Pi_leu`** | **0,113 – 0,211** | **MEDIDO** 0,162, banda **declarada** de ±30 % |
| `r_nuc_tum/leu`, `r_env_tum/leu` | 1,40–4,00 y 2,00–8,00 cada uno | transferido, sorteo independiente |
| `E_app_tum` | 1,527 – 2,772 kPa | medido (SBAT) |
| `E_app_leu` | 0,380 – 0,660 kPa | medido (granulocito–PBMC, dRT-DC) |
| `a_tum` | 8,51 – 10,00 µm | medido/transferido |
| **`a_leu`** | **3,025 – 4,395 µm** | **MEDIDO** (linfocito 6,05 µm ÷ 2 – monocito 8,79 µm ÷ 2). El Sobol anterior lo fijaba en 8,50 µm, el radio de HL-60 |
| `tau_tum` | 0,37 – 3,76 s | medido |
| `tau_leu` | 8,575 × 10⁻³ – 3,23 s | **no medido** (horquilla completa, para dejar que rompa la conclusión) |
| `ell_factor` | 0,50 – 2,00 | declarado |

## 2. Índices totales reales de la selectividad

Selectividad de constricción, **p5 = 3.06 · p50 = 5.40 · p95 = 9.67**
(el análisis anterior daba 2,98 / 4,40 / 6,42):

| parámetro | S_T | ± | S1 | clase |
|---|---|---|---|---|
| `E_app_tum` | **0.2454** | 0.0180 | 0.2153 | medido |
| `E_app_leu` | **0.2244** | 0.0148 | 0.2029 | medido |
| `C_tum` | **0.1330** | 0.0118 | 0.1215 | transferido/declarado |
| `a_leu` | **0.1036** | 0.0081 | 0.0934 | medido |
| `C_leu` | **0.1013** | 0.0095 | 0.0839 | transferido/declarado |
| `an_a_leu` | **0.0836** | 0.0061 | 0.0756 | medido |
| `an_a_tum` | **0.0422** | 0.0042 | 0.0414 | transferido/declarado |
| `Pi_tum` | **0.0352** | 0.0028 | 0.0309 | medido |
| `r_env_leu` | **0.0296** | 0.0025 | 0.0287 | transferido/declarado |
| `Pi_leu` | **0.0274** | 0.0021 | 0.0233 | medido |
| `a_tum` | **0.0194** | 0.0013 | 0.0184 | medido |
| `r_nuc_leu` | **0.0103** | 0.0010 | 0.0086 | transferido/declarado |
| `r_env_tum` | **0.0041** | 0.0005 | 0.0018 | transferido/declarado |
| `r_nuc_tum` | **0.0013** | 0.0002 | 0.0002 | transferido/declarado |
| `tau_tum` | **0.0000** | 0.0000 | 0.0000 | medido |
| `tau_leu` | **0.0000** | 0.0000 | 0.0000 | transferido/declarado |

**Reparto de `S_T` (suma = 1,061, es decir interacciones modestas):**

| clase de parámetro | suma de `S_T` | antes (Fase 2) |
|---|---|---|
| medidos (`E_app` ×2, `a` ×2, `Π` ×2, `a_n/a` leucocito) | **0,739** | ≈ 0,64 (los dos `E_app` más el radio) |
| **transferidos y declarados** (`C` ×2, `a_n/a` tumoral, `r_nuc` ×2, `r_env` ×2) | **0,322** | **≈ 10⁻²⁰ (cancelación)** |

> **El resultado real.** La selectividad **sí** depende de los parámetros transferidos: explican el
> **32 % de su varianza**, no el 10⁻¹⁸ %. Lo que sobrevive —y ahora es un resultado legítimo— es que
> **los dos módulos aparentes medidos siguen siendo los dos parámetros dominantes**
> (`E_app_tum` S_T = 0,245 ± 0,018 y `E_app_leu` S_T = 0,224 ± 0,015, juntos el **47 %**), y que el
> conjunto de los parámetros medidos explica casi tres cuartas partes de la varianza. Es una
> afirmación **mucho más fuerte** que la anterior porque ya no es una identidad algebraica: hay un
> contraste de información real entre lo medido y lo transferido, y lo medido gana.

El tercer parámetro en importancia es `C = K_s/G_s` de la población tumoral (S_T = 0,133), que es
**transferido**: es el principal candidato a medir si se quiere estrechar la incertidumbre de la
selectividad. Le siguen el radio del leucocito (0,104), `C` del leucocito (0,101) y la razón
núcleo/célula del leucocito (0,084) — esta última **medida en esta fase**, y que antes valía 10⁻²⁰.

`tau_tum`, `tau_leu` y `ell_factor` dan S_T = 0 **exactamente** sobre la selectividad de
constricción: eso sí es una independencia estructural genuina y correcta, porque la presión crítica
no contiene ninguno de los tres.

## 3. Las demás salidas

**Presión crítica tumoral** (p5 0.0886 / p50 0.1358 / p95 0.1956 kPa):

| parámetro | S_T | ± | S1 | clase |
|---|---|---|---|---|
| `E_app_tum` | **0.5213** | 0.0330 | 0.5096 | medido |
| `C_tum` | **0.2857** | 0.0200 | 0.2746 | transferido/declarado |
| `an_a_tum` | **0.0907** | 0.0064 | 0.0955 | transferido/declarado |
| `Pi_tum` | **0.0747** | 0.0055 | 0.0711 | medido |
| `a_tum` | **0.0414** | 0.0030 | 0.0398 | medido |
| `r_env_tum` | **0.0087** | 0.0011 | 0.0019 | transferido/declarado |

Con `Π` acotada a su **rango medido** (0,032–0,078 en vez de 0,03–1,0), su índice total cae de
**0,530 a 0,075** y el dominio pasa al módulo aparente medido (**S_T = 0,521**) y a `C`
transferido (0,286). **Esto corrige una de las conclusiones negativas más fuertes de la Fase 2**:
el valor absoluto de la presión crítica **no** está dominado por la regularización de citoplasma
una vez que `Π` se restringe a lo que se ha medido; está dominado por `E_app` medido, con `C`
transferido en segundo lugar.

**Presión crítica del leucocito** (p5 0.0163 / p50 0.0250 / p95 0.0373 kPa): dominada por `E_app_leu` (0,396), el radio
(0,187), la razón núcleo/célula (0,158) y `C` (0,155). El radio del leucocito, que antes era una
constante, es ahora el segundo parámetro.

**Selectividad de microfiltro** (p5 1.22 / p50 2.14 / p95 3.71; antes 2,77 / 4,06 / 5,85):
mismo reparto cualitativo, `E_app_tum` 0,260 y `E_app_leu` 0,238, con `C` de las dos poblaciones
aportando 0,300 en conjunto. **La mediana baja de 4,06 a 2.14.**

**Razón `β_c`** (p5 1.44 / p50 7.03 / p95 70.78; antes 1,21 / 5,69 / 52,2): sigue dominada por
`tau_leu`, con S_T = 1,0014 ± 0,1934. Como `S_T ≤ 1` por construcción, el exceso es ruido del
estimador sobre una salida de cola pesada; la lectura correcta es que `tau_leu` explica
esencialmente toda la varianza. **Nada cambia aquí: el τ del leucocito sigue siendo el hueco
dominante.**

**Fracción marginal del córtex** (p50 0.590) y **nuclear** (p50 0.0125): ahora dominadas
por `C_tum` (0,767) y `an_a_tum` (0,879) respectivamente, los dos **transferidos**. Honestidad
obligada: el resultado «el córtex limita» descansa sobre un parámetro transferido.

## 4. Qué conclusiones sobreviven al muestreo independiente

| conclusión | fracción de muestras | antes (Fase 2) |
|---|---|---|
| la tumoral exige mas presion que el leucocito | **1.0000** | 1,0000 |
| selectividad > 2x | **0.9984** | 1,0000 |
| selectividad > 4x | **0.8064** | — |
| microfiltro retiene tumoral y deja pasar leucocito | **0.9880** | 1,0000 |
| ventana de microfiltro > 2x | **0.5770** | 1,0000 |
| existe ventana de operacion beta_c(tum) > beta_c(leu) | **0.9783** | 0,9684 |
| ventana de presion de dispositivo > 10x | **0.3631** | 0,2844 |
| el cortex es el limitante (>50% del margen) | **0.8478** | 0,9303 |
| los compartimentos nucleares aportan <10% | **1.0000** | 0,9242 |

**Sobrevive con holgura:** que la tumoral exige más presión que el leucocito (1,0000), que la
selectividad supera 2× (0,9984) e incluso 4× (0,8064), que los compartimentos nucleares aportan
menos del 10 % del margen (**1,0000**, mejor que el 0,9242 anterior) y que existe ventana de
operación `β_c` (0,9783).

**Se degrada y hay que reescribirlo:**
- La ventana de microfiltro superior a 2× cae de **1,0000 a 0,5770**. El microfiltro deja de tener
  ventana garantizada: la tiene en el 98,8 % de las muestras, pero solo supera 2× en poco más de la
  mitad.
- «El córtex es el compartimento limitante» cae de 0,9303 a **0,8478**.
- La ventana de presión de dispositivo superior a 10× sube de 0,2844 a 0,3631, pero sigue siendo
  minoritaria y sigue dependiendo del `τ` leucocitario no medido.

## 5. Un tercer hallazgo que el cambio de radios obliga a reportar

Con los radios reales aparece un problema de encuadre que la Fase 2 no podía ver. El Sobol —el viejo
y el nuevo— compara las dos poblaciones **a confinamiento nominal igualado (2a/D = 1,70)**, es decir,
cada población en *su propia* garganta. Eso sólo era físicamente inocuo porque el radio del leucocito
valía 8,50 µm, casi el de la célula tumoral. **En un mismo dispositivo no lo es.**

Fijando la garganta para retener la CTC (`R_t = a_t/1,70 = 5,14 µm`, radio efectivo en el cuello
6,63 µm por el perfil tanh):

| población | radio (µm) | 2a/D nominal en ESE dispositivo | confinamiento efectivo | ¿contacto? | `ΔP_crít` del modelo |
|---|---|---|---|---|---|
| linfocito | 3,03 | 0,588 | 0,456 | **no** | **= 0** |
| neutrófilo | 4,20 | 0,817 | 0,634 | **no** | **= 0** |
| monocito | 4,40 | 0,855 | 0,663 | **no** | **= 0** |
| Jurkat (blasto linfoide) | 5,75 | 1,118 | 0,868 | no | = 0 |
| **HL-60 (lo que usaba la Fase 2)** | 8,50 | **1,653** | 1,283 | **sí** | calculable |
| MDA-MB-231 | 8,62 | 1,677 | 1,301 | sí | calculable |
| MCF-7 | 8,74 | 1,700 | 1,319 | sí | calculable |

**En una constricción dimensionada para retener una CTC, el leucocito circulante normal no llega
siquiera a tocar la pared.** La presión crítica que el modelo le asigna es cero y la selectividad es
formalmente **no acotada**, no 4,40×. La cifra de 4,40× describe un experimento que compara cada
población en un dispositivo distinto, y sólo parecía físico porque el radio del leucocito estaba
tomado de una línea de leucemia.

No existe además **ninguna** garganta que meta a las dos poblaciones dentro de la envolvente
validada: poner al neutrófilo en el umbral de contacto exige `R_t = 3,50 µm`, y ahí MCF-7 queda en
2a/D = 2,497, por encima del límite superior de 2,43.

Las dos correcciones de esta fase convergen en el mismo enunciado: el criterio mecánico de
constricción **sólo opera sobre células de diámetro ≳ 12,3 µm** (`a ≥ 1,2·R_t`), que es exactamente
la población —blastos circulantes, monocitos grandes, CTC pequeñas— en la que DLD también falla
(`dld_revisado.md` §5). Ésa, y no la competencia con el tamaño, es la justificación del criterio
mecánico.
