# Correcciones de la Fase 6 — qué cifra del manuscrito queda invalidada y por cuál se sustituye

Esta vía es sólo cálculo y literatura. No se ha tocado el manuscrito. Lo que sigue es la lista de
cifras que hay que cambiar, con el valor nuevo y la razón.

---

## 1. Cifras INVALIDADAS y su sustituto

### 1.1 Sección 3.5 y el descarte de DLD — invalidada entera

| cifra del manuscrito | valor nuevo | por qué |
|---|---|---|
| leucocito de referencia: **17,00 µm** de diámetro | **linfocito 6,05 · neutrófilo 8,40 · monocito 8,79 µm** | 17,00 µm está transferido de HL-60, una línea de leucemia promielocítica. Cuatro medidas independientes (Coulter, extinción axial de luz, microscopía de fase, volumen por micropipeta) sitúan al neutrófilo en 8,3–8,4 µm |
| brecha de medias **0,48 µm (2,75 %)** MCF-7 y **0,24 µm (1,39 %)** MDA-MB-231 | **9,08 µm (51,9 %)** y **8,84 µm (51,3 %)** frente al neutrófilo; **11,43 µm (65,4 %)** y **11,19 µm (64,9 %)** frente al linfocito | mismo cálculo, radio corregido |
| «**3,6 y 7,2 veces por debajo** de la resolución de `D_c`» | «**8,5 y 8,3 veces POR ENCIMA**» (frente al neutrófilo) | inversión completa del signo del argumento |
| solapamiento mínimo **0,439 / 0,618** | **1,0 × 10⁻⁴** (mezcla leucocitaria completa ponderada por el diferencial) | — |
| índice de Youden máximo **0,561 / 0,381** | **0,99990 / 0,99987** | — |
| «basta un CV de tamaño del **1,04–2,06 %** para que el solapamiento supere el 50 %» | **hace falta multiplicar por 3,10 el CV leucocitario medido (del 7,2 % al 22,3 %) sólo para llegar al 1 % de solapamiento**; el 50 % es inalcanzable | — |
| conclusión: «el criterio de `D_c` **no puede discriminar** CTC de leucocitos» | «**discrimina la población con Youden 0,9999**; no discrimina fenotipo, ni blastos circulantes, ni la cola de CTC pequeñas» | — |

**El resultado negativo no desaparece, cambia de sujeto.** Lo que la Fase 2 demostró, sin saberlo,
es que **DLD no separa CTC de blastos de leucemia mieloide**: HL-60 (17,0 µm) solapa el **79,6 %**
con MCF-7 y hunde Youden a **0,177**. Ese enunciado es verdadero, acotado y publicable; el otro no.

### 1.2 El resultado de Sobol sobre la selectividad — invalidado como «robustez»

| cifra del manuscrito | valor nuevo | por qué |
|---|---|---|
| «la selectividad **no depende de NINGÚN parámetro transferido**; sus índices son del orden de 10⁻²⁰, cero a precisión de máquina» | **los transferidos explican el 32 % de la varianza** (`C_tum` 0,133, `C_leu` 0,101, `an_a_tum` 0,042, `r_env` 0,034, `r_nuc` 0,012) | con `W = Gs·a²·f` y la MISMA `f` para las dos poblaciones, esos parámetros entran idénticos arriba y abajo del cociente y se cancelan **por construcción**. Verificado en el CSV original: S_T = 1,14–1,38 × 10⁻²⁰ (constricción) y 2,37–7,47 × 10⁻³¹ (microfiltro) — residuo en coma flotante de una resta de iguales |
| selectividad de constricción **2,98 / 4,40 / 6,42** (p5/p50/p95) | **3,06 / 5,40 / 9,67** | muestreo independiente con `Π`, `a_n/a` y radio del leucocito medidos |
| selectividad de microfiltro **2,77 / 4,06 / 5,85** | **1,22 / 2,14 / 3,72** | idem |
| «ventana de microfiltro > 2× en el **100 %** de las muestras» | **57,7 %** | degradación real; hay que reescribir la afirmación |
| «el córtex es el limitante en el **93,0 %** de las muestras» | **84,8 %** | — |
| «la presión crítica absoluta la domina `Π` (S_T = **0,530**)» | **S_T(`Π`) = 0,075**; domina `E_app` medido (**0,521**), con `C` transferido segundo (0,286) | el 0,530 salía de barrer `Π` sobre [0,03, 1,0]; restringida a su **rango medido** (0,032–0,078) su peso cae un factor 7. **Esto rehabilita el valor absoluto de la presión crítica**, que el manuscrito declara no-resultado |
| «los compartimentos nucleares aportan < 10 % en el 92,4 % de las muestras» | **100,0 %** | se refuerza |

### 1.3 Contraste de córtex tumoral/leucocito

| cifra del manuscrito | valor nuevo |
|---|---|
| «4,4× (MCF-7/PBMC) a 5,3× (SK-BR-3/PBMC)» | **2,50× (MDA-MB-231/PBMC) a 5,26× (SK-BR-3/PBMC)** — el manuscrito omite la línea más blanda, que es la limitante |
| «el ajuste de HL-60 de Wang es **5,2 veces** más blando que el `E_app` medido de un PBMC» | contra el `E_app` medido del **propio HL-60** (0,400 kPa → `Gs` = 1,304 mN/m), el ajuste de Wang (0,39) queda **3,34×** más blando. El 5,2× comparaba tipos celulares distintos |

*(Ambas ya estaban identificadas en el triaje; se recogen aquí para que la lista sea completa.)*

### 1.4 El marco «confinamiento nominal igualado» deja de ser inocuo

Con el radio de HL-60 (8,50 µm) el leucocito quedaba a 2a/D = 1,653 en la misma garganta que retiene
a la CTC. Con los radios reales queda a **0,59–0,86**, es decir **sin contacto con la pared**, y la
presión crítica que el modelo le asigna es **exactamente cero**. La selectividad en un mismo
dispositivo no es 4,40× sino **no acotada**. Hay que declarar que las cifras de selectividad
describen cada población en su propia garganta, y añadir que **no existe ninguna garganta que meta a
las dos poblaciones dentro de la envolvente validada** (el neutrófilo en el umbral de contacto exige
`R_t = 3,50 µm`, donde MCF-7 queda en 2a/D = 2,497 > 2,43).

---

## 2. Cifras que se MANTIENEN y ahora están mejor respaldadas

| afirmación | estado |
|---|---|
| la célula tumoral en suspensión es **más rígida** que el leucocito | **confirmada con una sola técnica**: estirador óptico, mismo instrumento y protocolo, 2,33× (MDA-MB-231/PBMC) frente al 2,46× del cruce de técnicas del manuscrito — **5,3 % de discrepancia** (10.3390/cancers13051119) |
| la selectividad la dominan los dos módulos aparentes **medidos** | **confirmada y ahora no trivial**: juntos explican el 47 % de la varianza, y el conjunto de los parámetros medidos el 74 % |
| los compartimentos nucleares no limitan el paso dentro de la envolvente | **reforzada**: < 10 % del margen en el 100 % de las muestras |
| la tumoral exige más presión que el leucocito | 1,0000 de las muestras |
| selectividad de constricción > 2× | 0,9984 |
| el 104 / 83 del barrido de constricción | **verificado exactamente** en el CSV: 104 = 83 (barrera > 10⁻³ kPa, córtex limitante en las 83) + 20 (`ΔP` = 0, todas en 2a/D = 1,20) + 1 (leucocito, garganta larga, 2a/D = 1,42, 0,000918 kPa) |

---

## 3. Lo que el FEM aporta al cociente de `E_app` — la pregunta del revisor, contestada

| línea | cociente crudo de `E_app` | FEM con `Π` = 0,10 común | FEM con `Π` medido por población | aporte del FEM |
|---|---|---|---|---|
| MCF-7 / PBMC | 4,274 | 4,395 | **2,863** | **−33,0 %** |
| SK-BR-3 / PBMC | 4,471 | 5,260 | **4,233** | **−5,3 %** |
| MDA-MB-231 / PBMC | 2,463 | 2,498 | **1,944** | **−21,1 %** |

Con `Π` común el FEM devolvía el cociente de módulos con un **+1,4 % a +17,7 %** de corrección: el
revisor tenía razón, no aportaba casi nada. Con la `Π` **medida** por población el FEM corrige el
cociente crudo **a la baja entre un 5 y un 33 %**, y lo hace por una razón física que un cociente de
módulos no puede capturar: la `Π` medida del leucocito (0,162) está **por encima** del valor común y
la de las líneas tumorales (0,032–0,078) **por debajo**, de modo que el contraste de citoplasma
**juega en contra** de la selectividad. La selectividad mediana de **4,40× baja a 2,86×** por este
solo efecto.

Dicho de otro modo: **el FEM no infla la selectividad, la desinfla**, y ése es precisamente el tipo
de aportación que lo justifica. Con los radios y las razones núcleo/célula medidas además, vuelve a
subir hasta una mediana de **5,47×** sobre las nueve parejas línea × subpoblación — pero por un
efecto geométrico distinto y declarado.

---

## 4. Huecos declarados que esta fase NO cierra

1. **`τ` del leucocito circulante.** Sigue anclado en una medida de célula entera de HL-60
   (8,58 ms) y sigue dominando la razón `β_c` (S_T = 1,00 ± 0,19). Es el hueco número uno.
2. **`τ` de SK-BR-3.** Sin medir.
3. **Exponente de ley de potencias de PBMC y granulocito humanos en suspensión.** Se usa THP-1
   (α = 0,25 por AFM, 0,44 por SF-DC) y el neutrófilo de Tsai 1993 (α = 0,48 en gran deformación).
4. **CV del diámetro de las líneas tumorales.** La semidispersión publicada es un error estándar.
   El CV del 10 % usado en el análisis de DLD es **declarado**, y es el parámetro que decide si DLD
   basta o no (pureza > 36 % con CV ≤ 15 %; enriquecimiento 14× con CV = 20 % y 99 % de rendimiento).
5. **Dispersión de tamaño del linfocito.** No recuperada; se le transfiere el CV del neutrófilo.
6. **Contraste de una sola técnica para MCF-7 y SK-BR-3.** Sólo existe para MDA-MB-231.
7. **Valores numéricos del mechano-NPS para Jurkat frente a MCF-7/MCF-10A medidos lado a lado**
   (10.1038/micronano.2017.91): sólo existen como figura. Es la segunda confirmación de una sola
   técnica más barata de conseguir.
8. **`C = K_s/G_s`**, tercer parámetro en importancia sobre la selectividad (S_T = 0,133) y
   enteramente transferido. Es el candidato número uno a medir si se quiere estrechar la
   incertidumbre de la selectividad.
9. **Ruta de recuperación de las cifras clásicas.** Los diámetros de Downey 1990,
   Schmid-Schönbein 1980 y Ting-Beall 1993 se obtuvieron a través de BioNumbers citando las tablas
   primarias; los PDF originales están de pago y no fueron accesibles en este entorno. Conviene
   verificarlos contra el original antes de enviar.
10. **Gossett 2012**: sólo se recuperó el resumen; la descripción del régimen se apoya en el resumen
    y en descripciones secundarias, no en la sección de métodos original.

---

## 5. Efecto neto sobre el encuadre del artículo

El título y el resumen tenían que apoyarse en «no es el tamaño». Con los radios reales **sí es el
tamaño**, y con holgura. Pero las dos correcciones apuntan al mismo sitio y dan un encuadre mejor:

- La separación por tamaño resuelve el **fondo leucocitario normal** (Youden 0,9999) y deja de ser
  el rival del artículo: pasa a ser la etapa previa.
- El criterio mecánico opera donde el tamaño falla, y ese dominio está ahora **cuantificado**: el
  criterio de constricción sólo actúa sobre células de diámetro ≳ 12,3 µm, que es exactamente la
  población —blastos circulantes, monocitos grandes, CTC pequeñas— en la que el solapamiento de
  tamaños sobrevive (Youden 0,177 frente a HL-60, 0,935 frente a Jurkat).
- La afirmación central —tumoral más rígida que leucocito en suspensión— tiene ahora **confirmación
  de una sola técnica** con un 5,3 % de discrepancia, que es lo que el trabajo necesitaba al no
  haber validación experimental propia.
