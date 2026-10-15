# Contraste de una sola técnica — el flanco del cruce pinzas-acústicas / RT-DC

**Resultado principal.** Sí existe un conjunto de datos en el que la línea de mama y el leucocito
están medidos con **la misma técnica, el mismo instrumento y el mismo protocolo**, y **confirma el
signo y casi exactamente la magnitud** del contraste que sostiene el artículo. No es RT-DC, es
**estirador óptico**.

---

## 1. Lo que se buscó y lo que hay

| familia de técnica | líneas de mama | leucocitos | ¿las dos en el mismo trabajo? |
|---|---|---|---|
| pinzas acústicas (SBAT) | MCF-7, SK-BR-3, MDA-MB-231 — 10.1038/s41378-020-0150-6 | **ninguno** | no |
| RT-DC / dRT-DC | **ninguna** encontrada | PBMC, granulocito, HL-60 (10.1038/s41467-019-08370-3); sangre completa, 21 donantes (10.7554/eLife.29213) | **no** |
| citometría de deformabilidad por flujo cortante (SF-DC) | ninguna | THP-1, K562, HL-60 (10.7554/eLife.78823) | no |
| detección por nanoporo (mechano-NPS) | MCF-7, MCF-10A | **Jurkat**, medido *lado a lado* en el mismo dispositivo (10.1038/micronano.2017.91) | **sí**, pero los valores solo existen como figura |
| mz-visco-NPS (G′, G″ resueltos en frecuencia) | MCF-10A, MCF-7, MDA-MB-231 (10.1002/advs.202406013) | ninguno | no |
| **estirador óptico (OS)** | **MDA-MB-231** + CTC de paciente | **PBMC** + HL-60 | **SÍ, con valores numéricos publicados en el texto** (10.3390/cancers13051119) |
| estiramiento hidrodinámico (xDC) | — | leucocitos y células malignas nativas de derrame pleural (10.1073/pnas.1200107109) | sí, pero otra procedencia de muestra |

**Declaración honesta: no se encontró ningún trabajo de RT-DC que mida líneas de carcinoma de mama.**
La cobertura de RT-DC sobre leucocitos es amplia y la de las líneas de mama es nula en lo recuperable
aquí; el camino para cerrar el flanco no es RT-DC sino el estirador óptico.

---

## 2. El dato que blinda la afirmación central

**Hayer *et al.*, *Cancers* 13:1119 (2021), DOI 10.3390/cancers13051119** — estirador óptico de haz
doble, pasos de tensión de 0,38 / 0,76 / 1,14 Pa (fuerzas de pico 80 / 160 / 240 pN), células en
suspensión y sin contacto, ajuste de Kelvin–Voigt. Midieron en el **mismo instrumento**: PBMC de
donante sano, PBMC de paciente con cáncer de mama, la línea MDA-MB-231, la línea HL-60 y 3641
células CD45-negativas de pacientes con carcinoma de mama.

Deformación relativa mediana publicada: **MDA-MB-231 = 0,012 · PBMC = 0,028 · HL-60 = 0,023**.

A paso de tensión igualado la deformación es inversa al módulo, de modo que el contraste de rigidez
MDA-MB-231/PBMC medido con una sola técnica es **0,028/0,012 = 2,33×**.

| | contraste MDA-MB-231 / PBMC |
|---|---|
| cruce de técnicas del manuscrito (SBAT 1,527 kPa ÷ dRT-DC 0,620 kPa) | **2,46×** |
| una sola técnica (estirador óptico, deformación relativa) | **2,33×** |
| discrepancia | **5,3 %** |

**El signo coincide y la magnitud coincide dentro del 5,3 %.** En la muestra clínica, además, los
PBMC se deformaron **el doble** que las CTC candidatas CD45-negativas de pacientes con carcinoma de
mama, en el mismo instrumento y las mismas condiciones — el mismo signo otra vez, ahora con células
de paciente y no con línea celular.

Esto es lo que había que encontrar: el cruce de técnicas podía estar inventando el contraste, y
resulta que no lo inventa.

---

## 3. Lo que el mismo dato SÍ invalida: HL-60 como leucocito de referencia

La pareja HL-60/PBMC da **signos opuestos** según la técnica:

| | contraste HL-60 / PBMC |
|---|---|
| dRT-DC (misma técnica, mismo trabajo) | 0,400/0,620 = **0,65×** — HL-60 **más blanda** |
| estirador óptico (misma técnica, mismo trabajo) | 0,028/0,023 = **1,22×** — HL-60 **más rígida** |
| factor de discrepancia | **1,89, con inversión de signo** |

Es decir: la única pareja en la que las dos técnicas se contradicen es precisamente la que involucra
**HL-60**, la línea de leucemia que el fallo 1 ya obligó a retirar como sustituto del leucocito
circulante. Las dos correcciones de esta fase apuntan en la misma dirección y se refuerzan: HL-60
no es un leucocito, ni por tamaño ni por mecánica.

---

## 4. Gossett *et al.* 2012 — qué midió exactamente y en qué régimen

**Gossett, Tse, Lee, Ying, Lindgren, Yang, Rao, Clark & Di Carlo, *PNAS* 109(20):7630–7635 (2012),
DOI 10.1073/pnas.1200107109.**

- **Técnica**: citometría de deformabilidad por **flujo extensional** (*hydrodynamic stretching*).
  Enfoque inercial lleva cada célula al punto de estancamiento de un cruce de canales, donde el
  flujo extensional la deforma **sin contacto con la pared**; cámara de alta velocidad y análisis
  automático. Caudal ~2000 células/s.
- **Observable**: *deformabilidad* = razón de aspecto en el instante de deformación máxima (métrica
  de forma, **no** un módulo). No hay ajuste constitutivo: no se publica E ni τ.
- **Régimen**: inercial, número de Reynolds alto, deformación en escala de **microsegundos** y
  deformaciones grandes. Es el extremo opuesto del nuestro: cuasi-estático, dominado por el
  contacto, escala de milisegundos a segundos.
- **Muestra**: células **nativas de derrame pleural** — leucocitos y células malignas de la misma
  muestra de paciente. No son líneas cultivadas y no proceden de sangre.
- **Resultado**: las células malignas del derrame resultaron **más deformables** que los leucocitos
  de la misma muestra, y los derrames malignos mostraron además células notablemente más grandes.
  Clasificación de estado de enfermedad con sensibilidad 91 % y especificidad 86 %.

**Qué significa para nosotros.** Gossett es un contraejemplo real al signo de nuestro contraste, y
la frase del manuscrito de que «el relato de células más blandas viene solo de células adherentes»
**no resiste** ese contraejemplo: Gossett mide células en suspensión. Lo defendible no es elegir
bando sino declarar las dos diferencias que separan los dos resultados:

1. **Procedencia de la muestra.** Derrame pleural frente a sangre. Una célula maligna que lleva
   semanas en una cavidad serosa no es una CTC recién desprendida en el torrente; el propio trabajo
   señala que las células malignas del derrame son más grandes. El dato de estirador óptico que sí
   es homólogo al nuestro —CTC de sangre de pacientes con carcinoma de mama frente a PBMC— da el
   signo **contrario al de Gossett y el mismo que el nuestro**.
2. **Régimen de deformación.** Microsegundos e inercia frente a milisegundos y contacto. Con
   reología de ley de potencias el módulo aparente escala como la frecuencia elevada a α, y α es
   distinto en las dos poblaciones (§`cota_ley_potencias.md`): **el orden de rigidez puede
   invertirse al cambiar la escala temporal de sondeo**, y eso es exactamente lo que se observa.

Esa es la frase que debe ir en la discusión: el signo del contraste depende de la escala temporal de
la sonda y de la procedencia de la muestra; nuestra conclusión está condicionada al régimen
cuasi-estático de contacto y a células de origen sanguíneo, y en ese régimen tiene **dos**
confirmaciones independientes de una sola técnica (línea celular y célula de paciente).

---

## 5. Lo que el cruce de técnicas sí le debe al revisor

Wu *et al.*, *Nature Methods* 15:491 (2018), DOI 10.1038/s41592-018-0015-1, midieron el mismo tipo
celular —MCF-7 entre otros— con AFM, microrreología de partículas, estirador óptico, reología de
placas paralelas y citometría de torsión magnética: los módulos medios difieren **en al menos dos
órdenes de magnitud** entre métodos, mientras que el error sistemático relativo **dentro** de cada
método se mantiene por debajo del 20 %. La conclusión operativa es la que conviene al artículo: los
**valores absolutos** de módulo no son comparables entre técnicas, pero los **cocientes medidos
dentro de una misma técnica** sí lo son. El artículo ya construye todo sobre cocientes (decisión de
diseño de la Fase 2); lo que faltaba era demostrar que el cociente sobrevive a medirse dentro de una
sola técnica, y §2 lo demuestra con un 5,3 % de discrepancia.

---

## 6. Huecos que quedan

- **No hay medida de una sola técnica para MCF-7 ni SK-BR-3 frente a leucocito.** El contraste de
  una sola técnica está anclado únicamente en MDA-MB-231, que es además la línea tumoral *más
  blanda* del panel: es el caso **más desfavorable**, lo que juega a favor, pero el factor 4,27×
  que el manuscrito atribuye a MCF-7/PBMC no tiene confirmación de una sola técnica.
- **El dato de una sola técnica es una métrica de forma, no un módulo.** La conversión
  «deformación ∝ 1/E» vale a paso de tensión igualado y en respuesta lineal; ambas condiciones se
  cumplen en el protocolo publicado, pero no es una medida directa de E.
- **Los valores de mechano-NPS para Jurkat frente a MCF-7/MCF-10A medidos lado a lado solo existen
  como figura** y no se han podido extraer en este entorno. Es la vía más directa que queda para
  una segunda confirmación de una sola técnica, y es barata: basta pedir los datos de origen.
