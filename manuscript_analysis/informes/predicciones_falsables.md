# Predicciones falsables del modelo

El trabajo es puramente computacional. Estas son las predicciones que el modelo hace de forma
comprometida, con el experimento concreto que las refutaría, la magnitud a medir, el umbral que
discrimina y la resolución necesaria. Todas las cifras de la columna «predicción» salen de los
artefactos de las Fases 1–2 y 6, no de estimaciones.

---

## P1 — Existe una ventana de presión de dispositivo entre leucocito y célula tumoral

- **Predicción.** El umbral `β_c = ΔP_el·τ/L_carga` separa las dos poblaciones:
  **PBMC 0,103 · granulocito 0,071 · HL-60 0,066** frente a **MDA-MB-231 17,0 · MCF-10A 40,4 ·
  MCF-7 169,4 kPa por mm/s**. La ventana es de **165× (frente a la tumoral más blanda) a 1 644×
  (frente a MCF-7)**. Con reología de ley de potencias se comprime a **3,3–17,5×**
  (`cota_ley_potencias.md`) pero no desaparece.
- **Qué la refutaría.** Medir `β_c(leucocito) ≥ β_c(MDA-MB-231)`, es decir encontrar que el
  leucocito también queda retenido a cualquier caudal en un dispositivo en el que la tumoral queda
  retenida. Equivale a `τ_leucocito ≥ 1,420 s`, **166 veces** el anclaje actual de célula entera de
  HL-60 (8,58 ms) y comparable al τ de MCF-7.
- **Magnitud a medir.** `τ` del leucocito circulante (PBMC y granulocito por separado) **en
  suspensión**, como cociente η/E de la misma sonda.
- **Umbral discriminante.** `τ_leu = 1,42 s`. Por debajo, la predicción sobrevive; por encima, cae.
  Es un umbral cómodo: hay dos órdenes y medio de margen.
- **Técnica y resolución.** dRT-DC (10.1038/s41467-019-08370-3) o citometría de deformabilidad por
  flujo cortante (10.7554/eLife.78823), que entregan η y E de la misma medida. Basta resolver `τ` a
  **un factor 2**: el margen es de 166×, así que incluso un error de un orden de magnitud no decide
  nada. **Es la predicción más barata de comprobar y la que más sostiene el artículo.**

## P2 — La separación es NO MONÓTONA en caudal, con máximo interior

- **Predicción.** La separación (fracción tumoral retenida menos arrastre leucocitario) vale 0 a
  caudal muy bajo, alcanza **1,000** en un óptimo interior y decae al seguir subiendo el caudal.
  Verificado en **27 de 27** combinaciones de presión de dispositivo × referencia leucocitaria. Las
  dos aristas tienen origen distinto: la inferior es **elástica** (`U ≈ ΔP_el,leu/β`, casi
  independiente de τ_leu) y la superior es **viscosa** (fijada por τ a través de `β_c`).
- **Qué la refutaría.** Una curva de separación **monótona creciente** en caudal sobre el rango
  accesible, o un óptimo situado fuera del intervalo predicho.
- **Magnitud a medir.** Fracción retenida de cada población frente a la velocidad media en la
  garganta, con las dos poblaciones en la misma corrida (mezcla marcada).
- **Umbral discriminante.** El ancho del intervalo útil depende fuertemente de β: **4,41 décadas a
  10 kPa por mm/s**, y luego **0,895 / 0,511 / 0,447 / 0,384 / 0,384 décadas** a 20 / 50 / 100 / 300
  / 1000. Para ver la arista superior hay que operar **por encima de `β_c`(tumoral)**; por debajo la
  curva es monótona *por predicción del propio modelo* y el experimento no discrimina.
- **Técnica y resolución.** Microscopía de alta velocidad en chip con conteo por población
  (dos fluoróforos) y barrido de caudal con **≥ 3 puntos por década** sobre ≥ 1,5 décadas,
  operando a β ≈ 50–300 kPa por mm/s donde la ventana mide 0,4–0,5 décadas. Precisión de conteo
  necesaria: resolver una caída de retención del 10 % → **≥ 400 células por punto** por población.
  **Es el experimento más exigente de los cinco** y el que más depende de fijar β con precisión.

## P3 — El margen de caudal se recupera ensanchando la garganta, y lo hace de forma muy abrupta

- **Predicción.** A geometría de referencia la presión crítica de MCF-7 cae
  **0,160 → 0,007 → 0,000 kPa** al pasar de 2a/D nominal 1,70 → 1,42 → 1,20: un factor **23×** en
  un solo escalón de confinamiento y cero por debajo de 1,20, donde el confinamiento efectivo es
  1,039 y no hay contacto. Para el leucocito, 0,037 → 0,002 → 0,000 kPa.
- **Qué la refutaría.** Medir una dependencia **suave** de la presión crítica con el radio de
  garganta (por ejemplo un factor < 3 entre 2a/D = 1,70 y 1,42), o barrera medible por debajo de
  2a/D = 1,20.
- **Magnitud a medir.** Presión mínima de paso frente al radio de garganta, a caudal cuasi-estático.
- **Umbral discriminante.** **Factor 10×** entre 2a/D = 1,70 y 1,42. El modelo predice 23×; un
  factor medido por debajo de 10 lo refuta.
- **Técnica y resolución.** Chip de constricción única con control de presión y detección de paso
  por impedancia o vídeo; serie de al menos **4 radios de garganta** entre 4,5 y 7,5 µm. Resolución
  de presión: **± 0,002 kPa** (2 Pa) para medir el escalón de 0,007 kPa — está al límite de lo
  rutinario, de modo que conviene medir el escalón 1,90 → 1,70 (0,363 → 0,160 kPa), donde basta
  ± 0,02 kPa.
- **Nota honesta.** Toda presión crítica del modelo es **cota inferior**: el tramo del camino con
  2a/D < 1,2 aporta cero porque el modelo cuasi-estático de sólido no describe deformación sin
  contacto. El experimento medirá valores iguales o mayores, nunca menores.

## P4 — Dentro de la envolvente, el córtex y no el núcleo limita el paso

- **Predicción.** En las **83 filas de diseño con barrera apreciable** el córtex es el compartimento
  limitante **sin una sola excepción**, con un reparto medio de la resistencia marginal de
  **córtex 0,718 · citoplasma 0,263 · envoltura nuclear 0,011 · nucleoplasma 0,009**: los
  compartimentos nucleares aportan el **2,0 %**.
- **Qué la refutaría.** Que una intervención sobre el núcleo mueva la presión crítica tanto o más
  que una intervención sobre el córtex de actina.
- **Magnitud a medir.** Presión crítica de paso antes y después de dos intervenciones separadas:
  (i) despolimerización de actina (latrunculina B o citocalasina D) y (ii) descondensación de
  cromatina (tricostatina A) o degradación de lámina A/C.
- **Umbral discriminante.** El modelo predice que la intervención nuclear cambia `ΔP_crít` **menos
  del 5 %** y la cortical **del orden del 70 %**. **La predicción cae si la intervención nuclear
  mueve `ΔP_crít` más de un 20 %.**
- **Técnica y resolución.** La misma plataforma de P3, o RT-DC con dosis crecientes de LatB —
  existe curva dosis-respuesta publicada para comparar (EC50 ≈ 26 nM en HL-60,
  10.7554/eLife.78823). Resolución de presión crítica: **± 10 %**, holgada frente al contraste
  70 % / 5 %.
- **Alcance declarado.** La predicción vale **solo dentro de la envolvente** (confinamiento efectivo
  ≤ 1,64, radio de garganta ≥ 5,34 µm frente a un radio nuclear de 4,37 µm). Con gargantas por
  debajo del radio nuclear el modelo no dice nada y la expectativa de la literatura (núcleo
  limitante) puede perfectamente ser cierta. **El experimento debe declarar el radio de garganta**,
  o no refuta nada.

## P5 — La presión crítica es proporcional al módulo aparente medido

- **Predicción.** `ΔP_crít ∝ Gs/a ∝ E_app`. Medido sobre las cuatro líneas con módulo publicado,
  `ΔP_crít/E_app` = **0,0591–0,0678 kPa/kPa**, es decir constante dentro del **± 7 %**.
- **Qué la refutaría.** Que el cociente `ΔP_crít/E_app` varíe más de **2×** entre líneas, lo que
  indicaría que algún parámetro no capturado (reserva de área, adhesión, forma nuclear) domina.
- **Magnitud a medir.** Presión crítica en la misma garganta para ≥ 3 líneas con `E_app` medido en
  la misma técnica.
- **Umbral discriminante.** Factor **2** en el cociente. El modelo predice ± 7 %.
- **Técnica y resolución.** Chip de constricción única + estirador óptico o SBAT sobre las mismas
  suspensiones. Resolución: `ΔP_crít` a ± 15 % y `E_app` a ± 15 %, que propagan a ± 21 % — bastante
  por debajo del umbral de 2×.
- **Por qué importa.** Es la predicción que convierte el FEM en algo más que un cociente de
  `E_app`: fija la **constante de proporcionalidad** y con ella el valor absoluto de la presión de
  operación, que es lo que un diseñador necesita.

## P6 — El tiempo de tránsito amplifica el contraste de τ

- **Predicción.** A presión fija, el cociente de tiempos de tránsito MCF-7/MDA-MB-231 vale
  **9,96–11,00** sobre las siete presiones del barrido (11,00 a 0,5 kPa; 10,40 a 1 kPa; 9,96 a
  50 kPa), frente a un cociente de τ de **5,7×** y de presión crítica de **1,76×**. El tiempo de
  tránsito no es un proxy de τ: lo amplifica aproximadamente por el factor de presión crítica.
- **Qué la refutaría.** Medir un cociente de tiempos de tránsito **igual al cociente de τ** (≈ 5,7)
  o menor.
- **Magnitud a medir.** Tiempo de tránsito por la constricción de MCF-7 y MDA-MB-231 a presión
  aplicada igualada, y τ de las dos por una técnica que dé η y E.
- **Umbral discriminante.** Cociente de tiempos **≥ 8**. Entre 5,7 y 8 la predicción queda indecisa.
- **Técnica y resolución.** Citometría de deformabilidad por constricción con vídeo a
  ≥ 2000 fotogramas/s (los tiempos predichos van de 0,23 ms a 286 ms según presión; a 10 kPa son
  1,2–11,7 ms, que es el rango cómodo). Resolución temporal **± 0,2 ms**.

## P7 — Kelvin–Voigt frente a ley de potencias: la forma de `ΔP_crít(U)`

- **Predicción.** Kelvin–Voigt predice `ΔP_crít(U) = ΔP_el·(1 + τU/L_carga)`, es decir **lineal** en
  el caudal; la ley de potencias predice `ΔP_crít ∝ U^α` con α medido de **0,147 (MDA-MB-231) ·
  0,186 (MCF-10A) · 0,234 (MCF-7)** (10.1038/s42003-020-01330-4). Sobre dos décadas de caudal los
  dos modelos difieren en un factor **100^(1−α) = 34 a 51** (50,8 para α = 0,147; 42,5 para
  α = 0,186; 34,0 para α = 0,234), es decir el cuadrado del factor por década de
  `cota_ley_potencias.md` (7,13 / 6,52 / 5,83).
- **Qué la refutaría.** Cualquiera de los dos: es una predicción **discriminante**, no un respaldo.
  Si `ΔP_crít` crece linealmente, el `β_c` del artículo es real; si crece como `U^α`, `β_c` no
  existe y lo único defendible es la ventana relativa comprimida.
- **Magnitud a medir.** Presión crítica de paso frente a la velocidad media en la garganta.
- **Umbral discriminante.** El exponente logarítmico local `d ln ΔP_crít / d ln U`: **≈ 1** para
  Kelvin–Voigt (a De ≫ 1), **0,15–0,25** para ley de potencias. Son inconfundibles.
- **Técnica y resolución.** Chip de constricción con barrido de presión motriz sobre **≥ 2 décadas**
  y ≥ 5 puntos; se necesita medir `ΔP_crít` con un error del **± 20 %** para separar pendiente 1 de
  pendiente 0,2 — exigencia muy modesta.
- **Por qué está aquí.** Es la predicción que decide entre dos versiones del propio artículo, y es
  la más barata después de P1. Debería ser el primer experimento de cualquier grupo que quiera
  continuar este trabajo.

---

## Resumen operativo

| predicción | coste experimental | qué decide | umbral |
|---|---|---|---|
| **P1** ventana de `β_c` | bajo (una corrida de dRT-DC) | la afirmación central de diseño | `τ_leu` < 1,42 s |
| **P7** forma de `ΔP_crít(U)` | bajo | qué versión del artículo es la correcta | pendiente log-log 1 vs 0,2 |
| **P5** `ΔP_crít ∝ E_app` | medio | que el FEM aporta la constante, no solo el cociente | cociente constante dentro de 2× |
| **P4** córtex vs núcleo | medio | el resultado negativo sobre el núcleo | intervención nuclear < 20 % |
| **P3** dependencia con el radio de garganta | medio | la forma de la barrera | factor ≥ 10× entre 1,70 y 1,42 |
| **P6** amplificación del tiempo de tránsito | medio-alto | el eje temporal como discriminante | cociente de tiempos ≥ 8 |
| **P2** no monotonía en caudal | **alto** | el régimen de operación óptimo | óptimo interior en 0,4–0,5 décadas |
