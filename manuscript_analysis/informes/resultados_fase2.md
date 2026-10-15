# Fase 2 — Física de las modalidades con el eje temporal como discriminante

**Eje de resultado (reorientado durante la ejecución): tumoral frente al fondo leucocitario.**
La comparación entre subtipos de mama pasa a ser secundaria. Todos los números de este informe
salen de los artefactos guardados, releídos para escribirlo.

---

## 0. Lo que hay que saber antes de leer los números

### 0.1 La envolvente de la Fase 1 estaba expresada en confinamiento NOMINAL

El perfil de pared tanh de la Fase 1,
`R_w(z) = R_t + (R_wide − R_t)·½[1+tanh((|z|−z_t)/L_tap)]`, **no alcanza el radio nominal
`R_t` en el cuello**: en `z = 0` vale `R_t + g₀(R_wide − R_t)` con
`g₀ = ½[1+tanh(−z_t/L_tap)] = 0,168` para `z_t/L_tap = 0,8`. El «confinamiento 2a/D = a/R_t»
de la Fase 1 es por tanto **nominal**, y el efectivo (`a/R_min`) es siempre menor:

| 2a/D nominal | 1,05 | 1,12 | **1,20** | 1,30 | 1,42 | 1,55 | 1,70 | 1,90 | 2,10 | 2,30 | **2,43** | 2,65 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| efectivo `a/R_min` | 0,942 | 0,989 | **1,039** | 1,100 | 1,170 | 1,241 | 1,319 | 1,415 | 1,504 | 1,586 | **1,636** | 1,716 |

Esto explica y confirma el hallazgo de la Fase 1 de que «no hay contacto por debajo de
2a/D = 1,2»: el confinamiento efectivo ahí es 1,04, es decir, la célula apenas roza la pared.
**Se ha conservado la geometría exacta de la Fase 1** (es aquella en la que se validaron la
convergencia de malla y la calibración 2D–3D de +16 %) y se mantiene la envolvente en su
convención original, 1,20 ≤ 2a/D nominal ≤ 2,43. Se probó renormalizar el perfil para que
`R_w(0) = R_t`, y se descartó: cambia la geometría validada y lleva el confinamiento efectivo
hasta 2,43, régimen en el que la deformación de área de membrana superaría el 60 % y el
problema deja de converger.

### 0.2 El córtex se ha CALIBRADO contra el módulo aparente medido

La Fase 1 dejó `Gs` enteramente transferido de PC-3 (13,25 mN/m) para las cuatro líneas de mama
y de HL-60 (0,39 mN/m) para el leucocito. Aquí se sustituye esa transferencia por una
calibración contra el único observable mecánico medido por línea, el **módulo de Young aparente
en suspensión** (pinzas acústicas, 10.1038/s41378-020-0150-6).

Sonda de calibración: tracción normal seguidora de patrón cuadrupolar
`σ_rr = p₀·P₂(cos θ)` sobre la célula libre — el modo que impone una pinza acústica de haz
único — y `E_app` definido operativamente como el módulo de la esfera homogénea casi
incompresible que mostraría la misma amplitud de deformación bajo la misma tracción. La
constante adimensional de la sonda, calibrada sobre esfera homogénea de módulo conocido, es
**K = 1,1856** (invariante en radio, módulo y amplitud; 0,02 % de variación al refinar la malla
de h = a/14 a a/20). Se usa la respuesta **linealizada** (primer paso de Newton desde el estado
sin deformar), que es exacta en el régimen de pequeña deformación en que se mide `E_app` y
evita la bifurcación de arrugamiento de la membrana sin rigidez de flexión.

Con la regularización de citoplasma de §0.3 todo escala con `Gs`, de modo que la relación
resulta **exactamente afín**, `E_app = K_Eapp·Gs/a`. Comprobado variando `Gs` (2 a 12 mN/m), el
radio (8,50 a 10,0 µm) y la malla (h = a/14 a a/26): **K_Eapp = 2,611 con una dispersión del
1,69 %**. Se usa 2,608, a 0,13 % de esa media, con lo que la calibración se reduce a
`Gs = E_app·a/2,608`:

| línea | `E_app` medido (kPa) | `Gs` calibrado (mN/m) | `Gs` transferido en Fase 1 | cociente |
|---|---|---|---|---|
| MCF-7 | 2,650 | **8,881** | 13,25 (PC-3) | 0,67 |
| SK-BR-3 | 2,772 | **10,629** | 13,25 (PC-3) | 0,80 |
| MDA-MB-231 | 1,527 | **5,047** | 13,25 (PC-3) | 0,38 |
| leucocito (PBMC) | 0,620 | **2,021** | 0,39 (HL-60) | 5,18 |
| MCF-10A | no medido | 8,881 *(transferido de MCF-7)* | 13,25 (PC-3) | — |

**Consecuencia cuantitativa importante:** el contraste de córtex entre tumoral y leucocito que
se deduce de los datos medidos es de **4,4× (MCF-7/PBMC) a 5,3× (SK-BR-3/PBMC)**, no el 34× que
implicaba la combinación de valores transferidos PC-3 frente a HL-60. El ajuste de HL-60 de
Wang et al. es 5,2 veces más blando que lo que exige el módulo aparente medido de un PBMC.

### 0.3 El citoplasma: por qué se usa la regularización de la Fase 1 y no el módulo medido

La Fase 1 fijó `μ_cit = 0,1·Gs/a`, es decir `Π ≡ G_cit·a/Gs = 0,1`, y verificó que `DI` cambia
solo +4,5 % y `ΔA/A₀` +22,5 % sobre un barrido de 33× en `μ_cit` (Π de 0,03 a 1,0). Los módulos
elásticos de citoplasma **medidos** por microrreología rotacional implican
Π = 0,032 (MCF-7), 0,066 (MDA-MB-231), 0,075 (SK-BR-3), 0,078 (MCF-10A), 0,162 (leucocito).

Con Π por debajo de ≈ 0,03 el problema **deja de estar bien puesto**: el córtex rígido rodea un
interior casi vacío, la membrana no tiene rigidez de flexión y el conjunto admite
configuraciones arrugadas de igual energía. Se comprobó de forma reproducible: la convergencia
falla en el mismo λ con cualquier penalización (kn·f de 20 a 200) y cualquier malla (h = a/10 o
a/14), y el fallo es monótono en Π — el leucocito (Π = 0,20) converge en todos los puntos con
cero reducciones de paso, MDA-MB-231 (Π = 0,066) en casi todos, MCF-7 (Π = 0,032) en casi
ninguno.

Se adopta por tanto **Π = 0,10** para todas las estirpes y se declara: el módulo elástico de
citoplasma de microrreología es una magnitud **local, dependiente de la escala de sonda** (el
propio trabajo publica 15–160 Pa·s al variar la longitud del microhilo de 2 a 6 µm), no el
módulo continuo de la gota. Lo que **sí** es robusto a la escala es `τ = η/G`, cociente de dos
medidas de la misma sonda — y `τ` es exactamente lo que fija el eje temporal de este trabajo.

### 0.4 Estructura dimensional (verificada, no supuesta)

Con `Π`, `C`, `a_n/a`, `r_nuc` y `r_env` fijos, la energía almacenada se factoriza como
`W = Gs·a²·f(c)`. El colapso se comprueba sobre las cinco estirpes, que difieren en `a` y `Gs`:
la dispersión relativa de `f = W/(Gs·a²)` a confinamiento nominal fijo es de **1,0 % a 5,9 %**
(1,0–1,0 % a 2a/D ≥ 2,3; 2,6–5,9 % entre 1,2 y 2,1), compatible con la incertidumbre de malla
que la Fase 1 cuantificó (1,3 % en DI, 2,3 % en ΔA). De ahí se sigue que
**ΔP_crítica ∝ Gs/a ∝ E_app**, y en efecto `ΔP_crit/E_app` vale 0,0591–0,0678 kPa/kPa sobre las
cuatro líneas con módulo medido (±7 %).

---

## 1. Constricción: presión crítica y tiempo de tránsito

### 1.1 Qué se calculó

El FEM entrega el paisaje energético `W(R)` de confinamiento radial (célula centrada en una
garganta de radio `R`). La entrada se reduce a una coordenada axial `s` con la aproximación
declarada `W(s) ≈ W_FEM(R_w(s))`, de donde la fuerza axial resistente es
`dW/ds = |dW/dR|·|dR_w/ds|` y el avance exige `ΔP·πR_w² ≥ dW/ds`. La presión crítica es el
máximo de `dW/ds /(πR_w²)` **restringido a la banda válida de la envolvente**.

**Reparametrización del camino de carga.** Con el camino original de la Fase 1 (la pared cierra
linealmente de `R_wide` a `R_t`) la célula no toca la pared durante la mayor parte de la
continuación: el primer contacto ocurre en λ = 0,941 / 0,713 / 0,607 a confinamiento nominal
1,20 / 1,70 / 2,43, es decir **el 61–94 % del camino transcurre sin contacto** y todo el
aplastamiento se concentra al final. Se reparametriza para arrancar en el primer contacto
(`R_min(λ=0) = 1,02·a`), lo que reparte los pasos de continuación sobre el aplastamiento. Es una
mejora **puramente numérica**: el estado en λ = 1 es idéntico y la geometría no cambia.

El máximo se alcanza en el propio cuello (`z → 0⁺`): el perfil `tanh(|z|)` tiene ahí un
**quiebre**, de modo que `|dR_w/dz|` salta de 0 a su valor finito justo donde el radio es
mínimo. Eso obliga a rejilla fina: de `nz = 301` a `nz = 12001` el valor sube de 0,1368 a
0,1609 kPa (MCF-7, 2a/D = 1,70). Se usa `nz = 6001`, con un residuo del orden del 1 % frente al
valor convergido — **la presión crítica reportada es, por este motivo, una ligera subestimación**.

### 1.2 Presión crítica (kPa), geometría de referencia (`z_t` = 4 µm, `L_tap` = 5 µm, θ = 41,5°)

| 2a/D nominal | MCF-7 | SK-BR-3 | MDA-MB-231 | MCF-10A | leucocito |
|---|---|---|---|---|---|
| 1,20 | 0,000 | 0,000 | 0,000 | 0,000 | 0,000 |
| 1,42 | 0,007 | 0,008 | 0,004 | 0,007 | 0,002 |
| 1,70 | **0,160** | **0,188** | **0,091** | 0,160 | **0,037** |
| 1,90 | 0,363 | 0,423 | 0,207 | 0,363 | 0,083 |
| 2,10 | 0,712 | 0,823 | 0,406 | 0,712 | 0,163 |
| 2,43 | 1,368 | 1,594 | 0,784 | 1,368 | 0,317 |

El cero a 2a/D = 1,20 no es un artefacto: el confinamiento efectivo ahí es 1,039 y el modelo
cuasi-estático de sólido, sin acoplamiento fluido, predice deformación nula sin contacto. Es la
cota inferior dura de la envolvente, y por ello **todas las presiones críticas de este informe
son cotas inferiores**: el tramo del camino con 2a/D < 1,2, donde una célula real sí se deforma
por el campo hidrodinámico, aporta cero a la barrera en este modelo.

### 1.3 Geometría: el ángulo de entrada manda

A confinamiento fijo (2a/D = 1,70) y para MCF-7, la presión crítica recorre un factor **42**
según la geometría de entrada:

| geometría | `z_t` (µm) | `L_tap` (µm) | θ efectivo | ΔP_crít (kPa) |
|---|---|---|---|---|
| entrada suave / garganta larga | 8,0 | 10,0 | 23,9° | **0,096** |
| referencia | 4,0 | 5,0 | 41,5° | 0,160 |
| garganta corta | 2,0 | 2,5 | 60,5° | 0,289 |
| entrada abrupta | 4,0 | 2,0 | 65,7° | **4,037** |

El rango completo del barrido es 0,0009–4,037 kPa (104 filas de diseño en
`barrido_constriccion.csv`).

### 1.4 Qué compartimento limita el paso

Sobre las 83 filas con barrera apreciable (ΔP_crít > 10⁻³ kPa), **el córtex limita el paso en
las 83**, sin una sola excepción. Reparto medio de la resistencia marginal `|dW_i/dR|`:

| compartimento | fracción del margen |
|---|---|
| córtex | **0,718** |
| citoplasma | 0,263 |
| envoltura nuclear | 0,011 |
| nucleoplasma | 0,009 |

**Esto contradice la expectativa de que el núcleo limite el paso en confinamiento fuerte.**
Dentro de la envolvente validada (confinamiento efectivo ≤ 1,64, es decir radio de garganta
≥ 5,34 µm frente a un radio nuclear de 4,37 µm) el núcleo todavía no es el cuello de botella:
con `a_n/a = 0,5` y los parámetros nucleares transferidos, los compartimentos nucleares aportan
el 2,0 % de la resistencia marginal. El régimen en que el núcleo dominaría exige radio de
garganta por debajo del radio nuclear, lo que cae fuera de la envolvente del solver. Se reporta
así, no como se esperaba.

### 1.5 Tiempo de tránsito (etapa de entrada)

Con respuesta de Kelvin–Voigt `σ = G(ε + τ·dε/dt)`, y como la energía es lineal en los módulos a
campo de deformación fijo, la fuerza generalizada se multiplica por `(1 + τ·ε̇/ε)`; de ahí
`ṡ = (m−1)·ε/(τ·dε/ds)` con `m = ΔP·πR_w²/(dW/ds)`. Tiempos a 2a/D = 1,70, geometría de
referencia:

| ΔP (kPa) | MDA-MB-231 (τ=0,57 s) | MCF-10A (τ=0,77 s) | MCF-7 (τ=3,23 s) |
|---|---|---|---|
| 0,5 | 26,0 ms | 68,2 ms | 286,2 ms |
| 1 | 12,3 ms | 30,4 ms | 127,6 ms |
| 10 | 1,17 ms | 2,78 ms | 11,66 ms |
| 50 | 0,23 ms | 0,55 ms | 2,31 ms |

El tiempo escala con `τ` a presión fija: el cociente MCF-7/MDA-MB-231 vale **9,96–11,00 sobre
las siete presiones del barrido** (11,00 a 0,5 kPa, 10,40 a 1 kPa, 9,96 a 50 kPa) — no es
exactamente constante, decrece levemente al subir la presión porque la barrera elástica pesa
menos frente al exceso de presión. Compárese con un cociente de `τ` de 5,7× y de presión crítica
de 1,76×: el tiempo de tránsito amplifica el contraste de `τ`. **SK-BR-3 y el
leucocito no aparecen: su `τ` no está medido y no se ha inventado.**

### 1.6 Convergencia

De 300 puntos FEM, **244 convergen y 56 no** (45 dentro de la envolvente, 11 fuera). Los fallos
se concentran en confinamiento nominal ≥ 1,90 (54 de 56) y en la geometría de entrada abrupta
(26 de 56). **75 puntos quedan fuera de la envolvente** por construcción del barrido (2a/D =
1,05, 1,12 y 2,65 en las cinco estirpes y cinco geometrías): se calcularon, se marcaron con
`dentro_envolvente = False` y se excluyen de todas las conclusiones.

---

## 2. Barrido de caudal — eje tumoral frente a leucocito

### 2.0 Cómo se agrega el grupo tumoral

El grupo principal son las **tres líneas tumorigénicas** con módulo aparente medido (MCF-7,
SK-BR-3, MDA-MB-231). **MCF-10A se reporta aparte**, por dos razones: es la línea epitelial *no
tumorigénica* del panel, de modo que agruparla dentro de «tumoral» sería incorrecto; y es además
la única sin `E_app` medido, con lo que su `Gs` está transferido de MCF-7. Como la agregación
literal pedida eran las cuatro líneas de mama, se calcula **también** esa variante y se compara:
en el punto de operación el rendimiento es **1,000 con tres líneas y 1,000 con cuatro**, y
`β_c`(MCF-10A) = 40,4 kPa por mm/s queda entre MDA-MB-231 (17,0) y MCF-7 (169,4), de manera que
incluirla **no cambia ni la línea limitante del grupo ni la ventana de presión de dispositivo**.
La columna `rendimiento_mama_4lineas` de `pureza_rendimiento_caudal.csv` lleva la variante
completa.

### 2.1 El grupo adimensional

Bajo una rampa de duración `t_carga = L_carga/U` un sólido de Kelvin–Voigt responde con rigidez
aparente `G(1 + τ/t_carga) = G(1 + De)`, con `De = τ·U/L_carga`. Luego
`ΔP_crít(U) = ΔP_el·(1 + De)`. En operación a caudal constante la presión motriz crece
linealmente con el caudal, `ΔP_motriz = β·U`, donde `β` es la resistencia hidráulica del
dispositivo por el área de garganta — un parámetro de diseño declarado. El margen de paso es

> `m(U) = m_∞·De/(1+De)`, con `m_∞ = β·L_carga/(ΔP_el·τ)`.

La célula escapa si `m > 1`. **Consecuencia: si `m_∞ ≤ 1` la célula no escapa a NINGÚN caudal**,
porque al aumentar el caudal se endurece exactamente al mismo ritmo al que crece la presión que
la empuja. El umbral de diseño es

> `β_c = ΔP_el·τ/L_carga`  (presión de dispositivo por unidad de velocidad).

### 2.2 El τ del leucocito: el parámetro crítico sin medir

No existe microrreología de citoplasma de ningún leucocito. El único anclaje utilizable es la
**viscosidad aparente de célula entera** de HL-60 por Kelvin–Voigt,
η = 3,43 ± 0,16 Pa·s con E = 0,400 ± 0,02 kPa (dRT-DC, 10.1038/s41467-019-08370-3), que da

> **τ_leucocito = η/E = 8,58 × 10⁻³ s**

Se declara expresamente: es **de célula entera, no de citoplasma**, y por tanto *no es homólogo*
a los τ de microrreología rotacional de las líneas de mama (0,57–3,23 s). Toda conclusión que
dependa de él va acotada en §2.5.

### 2.3 Umbrales de dispositivo (`z_t` = 4 µm, 2a/D nominal 1,70, `L_carga` = 3,06 µm)

| población | `E_app` (kPa) | `ΔP_el` (kPa) | τ (s) | `β_c` (kPa por mm/s) |
|---|---|---|---|---|
| SK-BR-3 | 2,772 | 0,188 | *bracket* 0,57–3,23 | **35,0 – 198,5** |
| MCF-7 | 2,650 | 0,160 | 3,23 | **169,4** |
| MCF-10A | — | 0,160 | 0,77 | 40,4 |
| MDA-MB-231 | 1,527 | 0,091 | 0,57 | **17,0** |
| leucocito PBMC | 0,620 | 0,0367 | 8,58e-3 | **0,103** |
| leucocito granulocito | 0,430 | 0,0254 | 8,58e-3 | **0,071** |
| leucocito HL-60 | 0,400 | 0,0237 | 8,58e-3 | **0,066** |

SK-BR-3 entra con su módulo aparente medido y un **bracket declarado** de τ (el intervalo de las
líneas de mama medidas), no con un valor inventado.

> **La ventana de presión de dispositivo entre el leucocito y la célula tumoral más blanda
> (MDA-MB-231) es de 0,103 a 17,0 kPa por mm/s: un factor 165. Frente a MCF-7 el factor es
> 1 644.**

### 2.4 ¿Hay un caudal óptimo? Sí, y es no monótono

**La no monotonía se confirma: hay un máximo interior en las 27 combinaciones de presión de
dispositivo × referencia leucocitaria examinadas.** La separación (fracción tumoral retenida
menos arrastre leucocitario) vale 0 a caudal muy bajo (se retiene todo), alcanza 1,000 en un
óptimo, y decae al seguir subiendo el caudal cuando también escapa la célula tumoral. Las dos
aristas tienen origen distinto:

- **arista inferior** (empieza a escapar el leucocito): la fija el contraste **elástico**,
  `U ≈ ΔP_el,leu/β`, y es prácticamente independiente de `τ_leu`;
- **arista superior** (empieza a escapar la tumoral): la fija `τ` a través de `β_c`, y
  **no existe** si `β < β_c(tumoral)`, en cuyo caso la tumoral se retiene a cualquier caudal.

Ancho del intervalo de caudal con separación > 0,9:

| presión de dispositivo (kPa por mm/s) | intervalo útil de `U` | décadas |
|---|---|---|
| 10 | 3,9 µm/s – >100 mm/s | **4,41** |
| 20 | 1,86 – 14,6 µm/s | 0,89 |
| 50 | 0,77 – 2,50 µm/s | 0,51 |
| 100 | 0,37 – 1,03 µm/s | 0,45 |
| 300 | 0,13 – 0,32 µm/s | 0,38 |
| 1000 | 0,038 – 0,091 µm/s | 0,38 |

**Lectura honesta:** el óptimo existe y es interior, pero a presión de dispositivo baja
(10 kPa por mm/s) la meseta útil abarca 4,4 décadas de caudal, de modo que el caudal **no es un
knob fino** en ese régimen; lo que selecciona el régimen de operación es la razón
presión/velocidad `β`. El caudal óptimo sí depende del subtipo a través de su velocidad crítica
`U*` — a 192 kPa por mm/s, `U*` vale 0,51 µm/s (MDA-MB-231), 1,02 (MCF-10A) y 5,63 (MCF-7), un
rango de 11× — pero esa dependencia se vuelve operativamente relevante solo por encima de
`β_c(tumoral)`, donde la ventana se estrecha a ~0,4 décadas.

Mejor punto de operación por pureza × rendimiento, con un ratio declarado de 1 CTC por 10⁶
leucocitos: **10 kPa por mm/s y U = 3,9 µm/s**, con rendimiento tumoral 1,000 y arrastre
leucocitario por debajo de la precisión de la aritmética.

### 2.5 Cuánto aguanta la conclusión si el τ del leucocito es mayor

La ventana existe mientras `β_c(leucocito) < β_c(tumoral)`, lo que acota el τ leucocitario:

| referencia | cota de τ_leu frente a la tumoral más blanda | margen sobre el anclaje |
|---|---|---|
| PBMC | **1,420 s** | 166× |
| granulocito | 2,047 s | 239× |
| HL-60 | 2,201 s | 257× |

Es decir: **el τ del leucocito tendría que ser 166 veces mayor que el valor anclado — y
comparable al de MCF-7, la línea de mama de relajación más lenta — para que la ventana se
cerrase.** Además, el arrastre leucocitario en el punto de operación se mantiene por debajo de
1,8 × 10⁻¹³ incluso suponiendo un coeficiente de variación del 100 % en τ_leu
(`sensibilidad_tau_leucocito.csv`), de modo que la conclusión no depende de la forma de la cola.
Se advierte no obstante que la dispersión publicada de τ_leu (CV 3,3 %) procede del **error
estándar** de η y E, no de la variabilidad célula a célula: las purezas calculadas son por eso
**cotas superiores**, y lo robusto aquí es la separación de umbrales `β_c`, no el valor absoluto
de la pureza.

---

## 3. Microfiltro: umbral de paso por poro

Poro de canto vivo: no define pendiente de pared, así que la barrera se reparte sobre una
longitud de deglución `ℓ = ℓ_factor·a` (la célula debe avanzar del orden de un radio para quedar
engullida); `ΔP_poro = W(R_p)/(ℓ·πR_p²)`. `ℓ_factor` es un parámetro **declarado** y entra en la
sensibilidad: a poro de 10 µm y MCF-7 el umbral vale 0,935 / 0,468 / 0,234 kPa para
`ℓ_factor` = 0,5 / 1 / 2 (escala exactamente como 1/ℓ).

Se usa la geometría de garganta larga (`z_t` = 8 µm), la más próxima a un poro recto disponible;
se declara que un poro real recto confinaría la célula en toda su longitud y daría una barrera
mayor, por lo que estos umbrales son **cotas inferiores**.

**Existe ventana de operación en los 14 diámetros de poro dentro de envolvente (7,5–14,0 µm),
sin excepción** (`ventana_microfiltro.csv`). Con `ℓ_factor` = 1:

| diámetro de poro (µm) | leucocito (kPa) | tumoral mínima (kPa) | factor de selectividad |
|---|---|---|---|
| 7,5 | 0,578 | 1,521 | 2,63 |
| 9,0 | 0,195 | 0,525 | 2,69 |
| 10,0 | 0,089 | 0,244 | 2,74 |
| 12,0 | 0,013 | 0,038 | 2,94 |
| 14,0 | 0,0003 | 0,0016 | 4,62 |

La tumoral limitante es **MDA-MB-231 en los 14 diámetros** (es la de menor módulo aparente). El
ancho absoluto de la ventana es máximo en el poro más estrecho (0,943 kPa a 7,5 µm) y cae a
0,0012 kPa a 14 µm, mientras el factor relativo crece: el poro estrecho da margen operativo
absoluto, el ancho da contraste relativo pero a presiones impracticablemente pequeñas.

**Dirección del efecto, que conviene no confundir:** la célula tumoral en suspensión es *más
rígida* que el leucocito, luego exige *más* presión para pasar. La ventana es
`ΔP_leucocito < ΔP_aplicada < ΔP_tumoral`: a esa presión el leucocito pasa y la tumoral queda
retenida. Es lo contrario del relato habitual de «la célula cancerosa es más blanda», que
procede de medidas en célula adherida.

---

## 4. DLD: descarte cuantitativo (resultado negativo)

No se modela. El solver solo es válido en régimen de contacto (2a/D ≥ 1,2) y DLD opera sin
contacto. El descarte se hace con los radios medidos en suspensión y con el propio criterio de
DLD, por tres vías independientes.

### 4.1 La brecha de medias está por debajo de la resolución del método

Correlación empírica de Davis, `D_c = 1,4·g·ε^0,48`, con exactitud publicada del orden del 10 %
de `D_c`. Frente a un leucocito de 17,0 µm de diámetro:

| línea | diámetro (µm) | brecha (µm) | brecha relativa | resolución de `D_c` (µm) | ¿separable? |
|---|---|---|---|---|---|
| MCF-7 | 17,48 | 0,48 | 2,75 % | 1,75 | **no** — 3,6× bajo resolución |
| MDA-MB-231 | 17,24 | 0,24 | 1,39 % | 1,72 | **no** — 7,2× bajo resolución |
| MCF-10A | 17,48 *(transferido)* | 0,48 | 2,75 % | 1,75 | no |
| SK-BR-3 | 20,00 *(transferido)* | 3,00 | 15,0 % | 2,00 | sí |

**Este argumento no supone ninguna dispersión:** compara la brecha de medias con la exactitud
intrínseca de la correlación de `D_c`. Las geometrías de array que centrarían `D_c` en el
leucocito son g = 51,1 / 36,7 / 26,3 µm para ε = 0,05 / 0,10 / 0,20, todas con `D_c` = 17,0 ±
1,7 µm. Nótese que **la única línea que sería separable por tamaño es precisamente aquella cuyo
diámetro no está medido individualmente** (SK-BR-3 se publica como media agregada del panel).

### 4.2 Solapamiento de distribuciones, con la dispersión mínima compatible con los datos

Tratando el error estándar publicado como si fuera desviación típica — lo que da el
solapamiento **mínimo** posible, porque SD ≥ SE siempre:

| línea | σ usada (µm) | coeficiente de solapamiento | Bhattacharyya | índice de Youden máximo | sens./espec. |
|---|---|---|---|---|---|
| MCF-7 | 0,31 | **0,439** | 0,741 | **0,561** | 0,781 / 0,781 |
| MDA-MB-231 | 0,21 | **0,618** | 0,870 | **0,381** | 0,815 / 0,566 |

Al leucocito no se le publica dispersión de tamaño; se le asigna la de MCF-7 y se declara como
transferida.

### 4.3 Y crece deprisa con cualquier dispersión realista

El solapamiento supera el 50 % en cuanto el coeficiente de variación del tamaño llega a
**2,06 %** (MCF-7) o **1,04 %** (MDA-MB-231), y el 80 % a 5,50 % y 2,77 % respectivamente — muy
por debajo de cualquier dispersión de tamaño plausible en una población celular. A un CV del
10 %, el índice de Youden máximo cae a 0,111 (MCF-7) y 0,056 (MDA-MB-231).

**Conclusión:** el criterio de diámetro crítico de DLD no puede discriminar células tumorales de
mama frente a leucocitos. Es un resultado negativo defendible y es el contraejemplo que
justifica usar un criterio mecánico en lugar de uno geométrico.

---

## 5. Sensibilidad global de Sobol

### 5.1 Montaje

Sustituto de `f = W/(Gs·a²)` sobre seis entradas adimensionales (`c`, `a_n/a`, `C = K_s/G_s`,
`Π`, `r_nuc`, `r_env`), ajustado sobre un hipercubo latino de 192 puntos FEM de los que
**141 convergen (73 %)**. Base polinómica de grado 2 completo más grados 3–5 en el confinamiento
(f recorre tres órdenes de magnitud con `c`): **R² = 0,9994, error mediano en hold-out del
2,9 %, p90 12,3 %, máximo 42,1 %**. Sobol de primer y total orden con 28 672 muestras sobre
12 parámetros, con los nucleares y de membrana — todos transferidos — incluidos explícitamente.

*Desviación declarada:* el rango de `Π` se acotó a [0,03, 1,0], el verificado por la Fase 1, en
vez de [0,005, 1,0]; por debajo de 0,03 el problema no está bien puesto (§0.3) y el diseño se
atascaba sin producir información.

### 5.2 Qué domina cada salida (índice total `S_T`)

| salida | p5 / p50 / p95 | parámetros dominantes |
|---|---|---|
| ΔP crítica tumoral (kPa) | 0,128 / 0,250 / 0,407 | **Π 0,530 ± 0,036**, E_app tumoral 0,278 ± 0,021, C 0,160 ± 0,012, a_n/a 0,044 |
| **Selectividad constricción tumoral/leucocito** | **2,98 / 4,40 / 6,42** | **E_app tumoral 0,551 ± 0,034, E_app leucocito 0,418 ± 0,025**, radio tumoral 0,044; **todos los transferidos ≈ 10⁻²⁰** |
| Razón `β_c` tumoral/leucocito | 1,21 / 5,69 / 52,2 | **τ leucocito 1,057 ± 0,386**, τ tumoral 0,117 ± 0,057 |
| Umbral de poro tumoral (kPa) | 0,150 / 0,345 / 0,836 | ℓ_factor 0,629 ± 0,049, Π 0,204 ± 0,020, E_app tumoral 0,118, C 0,102 |
| Selectividad microfiltro | 2,77 / 4,06 / 5,85 | E_app tumoral 0,575 ± 0,039, E_app leucocito 0,437 ± 0,026 |
| Fracción marginal del córtex | 0,489 / 0,602 / 0,685 | Π 0,622 ± 0,069, a_n/a 0,346 ± 0,022, C 0,295 ± 0,034 |
| Fracción marginal nuclear | 0,010 / 0,037 / 0,110 | a_n/a 0,603 ± 0,045, Π 0,422 ± 0,035, r_env 0,042 |

El índice total de `τ_leucocito` sobre la razón `β_c` sale **1,057 ± 0,386**, por encima de la
unidad. Como `S_T ≤ 1` por construcción, ese exceso es **ruido del estimador**, no un valor con
significado: con 2048 muestras base y una salida de cola muy pesada (p95/p50 = 9,2) el intervalo
de confianza es del ±37 %. La lectura correcta es que `τ_leucocito` explica *esencialmente toda*
la varianza de esa salida, no que la explique «más del 100 %». Los demás índices tienen
intervalos del 7–12 % y se leen sin esa reserva.

El resultado estructural más importante: **la presión crítica absoluta la domina la
regularización de citoplasma `Π` (53 % de la varianza), pero la selectividad tumoral/leucocito
no depende de NINGÚN parámetro transferido** — sus índices son del orden de 10⁻²⁰, cero a
precisión de máquina — porque los parámetros nucleares y de membrana son comunes a las dos
poblaciones y se cancelan exactamente en el cociente. La selectividad queda fijada por los dos
módulos aparentes **medidos** y por el radio.

### 5.3 Qué conclusiones sobreviven

| conclusión | fracción de muestras en que se cumple |
|---|---|
| La tumoral exige más presión que el leucocito en constricción | **1,0000** |
| …y la selectividad supera 2× | **1,0000** |
| El microfiltro retiene tumoral y deja pasar leucocito | **1,0000** |
| …con ventana > 2× | **1,0000** |
| El córtex es el compartimento limitante (> 50 % del margen) | **0,9303** |
| Los compartimentos nucleares aportan < 10 % del margen | **0,9242** |
| Existe ventana de operación: `β_c`(tumoral) > `β_c`(leucocito) | **0,9684** |
| …y la ventana de presión supera 10× | 0,2844 |
| …y supera 100× | 0,0258 |

**Sobrevive con holgura:** que la célula tumoral se retiene donde el leucocito pasa, tanto en
constricción como en microfiltro, con un factor mediano de 4,1–4,4× y un p5 de 2,8–3,0×; y que
el córtex, no el núcleo, es el compartimento limitante dentro de la envolvente.

**No sobrevive:** la magnitud de la ventana de presión de dispositivo. Supera 10× solo en el
28 % de las muestras y 100× solo en el 2,6 %, porque la domina el τ leucocitario, que no está
medido: la ventana se cumple en ≥ 99 % de las muestras solo si τ_leu ≤ 2,28 s. Tampoco sobrevive
ningún valor absoluto de presión crítica ni de umbral de poro: están dominados por `Π` y por
`ℓ_factor`, que son parámetros declarados del modelo, no datos.

---

## 6. Resumen de lo que salió y lo que no

**Salió:**
1. La selectividad tumoral/leucocito es robusta y está determinada enteramente por datos
   medidos: factor 4,4× mediano (p5 2,98, p95 6,42) en constricción y 4,1× en microfiltro, con
   cero sensibilidad a los parámetros transferidos.
2. Existe un caudal óptimo y la separación es **no monótona** en caudal (27/27 combinaciones),
   con las dos aristas de origen distinto: la inferior elástica, la superior viscosa (fijada por τ).
3. La ventana de presión de dispositivo entre leucocito y tumoral es de 165× a 1 644×, y
   sobrevive a que el τ del leucocito sea hasta 166 veces mayor que su anclaje.
4. El microfiltro tiene ventana de operación en todos los diámetros de poro validados.
5. El descarte de DLD es firme por tres vías independientes.
6. La calibración del córtex contra el módulo aparente medido reduce el contraste
   tumoral/leucocito de 34× (valores transferidos) a 4,4–5,3× (dato medido).

**No salió como se esperaba, y se reporta así:**
1. **El núcleo no limita el paso** en ninguna parte de la envolvente validada: lo hace el córtex
   en las 83 filas con barrera apreciable, con los compartimentos nucleares aportando el 2,0 %.
   El régimen en que el núcleo dominaría exige radio de garganta por debajo del radio nuclear y
   cae fuera de la envolvente.
2. **El caudal no es un knob fino** en el régimen de baja presión de dispositivo: la meseta útil
   abarca 4,4 décadas a 10 kPa por mm/s. El parámetro de control que selecciona el régimen es la
   razón presión/velocidad, no el caudal.
3. Ningún valor absoluto de presión es un resultado: la presión crítica está dominada por la
   regularización `Π` y el umbral de poro por `ℓ_factor`, ambos parámetros declarados.

**Huecos declarados, sin rellenar:** τ de SK-BR-3 y del leucocito (este último anclado en una
medida de célula entera, no homóloga); dispersión de `E_app` de las líneas tumorales; módulo de
citoplasma a escala continua; y la parte del camino de entrada con 2a/D < 1,2, que este modelo
no puede describir y que hace de todas las presiones críticas cotas inferiores.

---

## 7. Artefactos

| archivo | contenido |
|---|---|
| `barrido_constriccion.csv` | 104 filas: presión crítica, tiempo de tránsito a 7 presiones, compartimento limitante y reparto marginal, por estirpe × geometría × confinamiento |
| `barrido_caudal.csv` | 18 720 filas: De, presión crítica aparente, margen, `U*` y fracción capturada por estirpe × 26 presiones de dispositivo × 240 caudales |
| `separacion_caudal.csv` | separación por pares de subtipo y presión de dispositivo, con test de máximo interior |
| `pureza_rendimiento_caudal.csv` | **eje principal**: rendimiento tumoral, arrastre leucocitario, pureza y enriquecimiento frente al caudal, por referencia leucocitaria |
| `umbral_dispositivo.csv` | `β_c` por población |
| `sensibilidad_tau_leucocito.csv` | arrastre y pureza frente al CV supuesto de τ leucocitario |
| `barrido_microfiltro.csv` | 246 filas: umbral de poro frente a diámetro, estirpe y `ℓ_factor` |
| `ventana_microfiltro.csv` | ventana de operación por diámetro de poro |
| `dld_solapamiento.md` | argumento cuantitativo del descarte de DLD |
| `sensibilidad_sobol.csv` | `S1` y `S_T` con intervalos, 7 salidas × 12 parámetros |
| `figura_caudal.png` | figura principal (eje tumoral/leucocito y qué sobrevive al Sobol) |
| `paso0_calibracion.py` | calibración del córtex y verificación de que `K_Eapp` es constante |
| `fase2_lib.py`, `paso*.py` | código reproducible |
