# Triaje de tu revisión del manuscrito

Verifiqué contra los artefactos cada punto que toca una cifra. Resultado: **los dos problemas de
fondo más graves son correctos**, y ambos dañan resultados que yo te había presentado como sólidos.

---

## A. Entra tal cual — verificado, tienes razón

### A1. El tamaño del leucocito invalida el descarte de DLD (tu punto 1) — CRÍTICO
Usamos a = 8,5 µm (diámetro 17 µm) transferido de **HL-60, que es una línea de leucemia, no un
leucocito circulante**. Recalculado con tamaños reales:

| Leucocito | diámetro | brecha frente a MCF-7 | ¿separable por D_c (±10%)? |
|---|---|---|---|
| HL-60 (lo que usamos) | 17,0 µm | 0,48 µm = 2,8 % | no |
| Neutrófilo | 11,0 µm | 6,48 µm = **58,9 %** | **sí** |
| Linfocito | 8,0 µm | 9,48 µm = **118,5 %** | **sí** |

Con tamaños reales **DLD sí separa, y con holgura**. Eso tumba la sección 3.5 entera, el "not size"
del título, y la frase del resumen. Es además el fundamento de ISET y Parsortix, que citamos.

### A2. La robustez del Sobol es cancelación algebraica, no un hallazgo (tu punto 3) — CRÍTICO
Confirmado en el CSV: los índices totales de a_n/a, C, Π, r_nuc y r_env sobre la selectividad valen
**1,3·10⁻²⁰** (constricción) y **2,4·10⁻³¹** (microfiltro). Eso no es insensibilidad, es
**cancelación exacta manifestándose como ruido de coma flotante**: con W = Gs·a²·f(c) y la misma f
para las dos poblaciones, esos parámetros se van en el cociente por construcción.

Y la cancelación se rompe exactamente donde tú dices: la **Π medida del leucocito es 0,162 frente a
0,032–0,078 de las tumorales** (2,5–5× de diferencia) y el Sobol usó **Π común de 0,10**; la razón
núcleo/célula del linfocito está muy por encima de 0,50. Hay que remuestrear por población.

Esto desmonta lo que yo te vendí como "el resultado que blinda el artículo". No blinda nada: con
los parámetros cancelados, la selectividad **es** el cociente de E_app (≈4,3 para MCF-7/PBMC), y la
pregunta de qué aporta el FEM a esa cifra es legítima.

### A3. Contraste de córtex mal reportado
MCF-7/PBMC = 4,394× · SK-BR-3/PBMC = 5,259× · **MDA-MB-231/PBMC = 2,497×**. El rango correcto es
**2,50–5,26×**, no 4,4–5,3×: el manuscrito omite la línea más blanda, que es justo la limitante.

### A4. "HL-60 5,2× más blando que PBMC" compara tipos distintos
Contra el E_app medido del propio HL-60 (0,400 kPa) → Gs = 1,304 mN/m, y el ajuste de Wang (0,39)
queda **3,34×** más blando, no 5,2×.

### A5. Las nueve inconsistencias numéricas restantes
Todas son de redacción o de trazabilidad y hay que arreglarlas una a una: las dos constantes de
sonda sin relacionar, el umbral de meseta (>0,9 en texto, ≥0,95 en figura), la convergencia de malla
(0,3 % vs 1,3 %), las tres definiciones de DI, el límite de lisis frente a la envolvente y la
reserva de área, el pie de la tabla 6, y la relación entre 104 filas / 300 puntos / 83 filas / 75
fuera de envolvente.

### A6. Forma de la revista
Todo correcto y todo pendiente: resumen de 150–250 palabras (el nuestro tiene ~550), plantilla
sn-jnl con `[iicol]`, **sin títulos dentro de las ilustraciones** (las cuatro figuras los llevan
—esto choca de frente con la convención de figura que apliqué—), pies con "Fig." en negrita y sin
punto final, patrones además de color en 2a y 3b, anchos de 84/174 mm, tensores en negrita y "d"
redonda, Statements and Declarations obligatorio, nomenclatura ESM_1, y uso de IA documentado en
Métodos.

### A7. Las frases de bitácora
Tienes razón en todas. "The reviewer's obvious question", "verified twice independently", "we did
not invent one", "renormalisation was tested and discarded". Y **"We take it"** es temerario si Sui
o Feng arbitran: lo escribí yo en el brief y es un error de tono.

### A8. El cribado de 157 registros sin método, y Paul et al. vol. 17 (2017)
Correctos los dos.

---

## B. Entra, pero con matiz

### B1. Gossett 2012 (tu punto 2)
De acuerdo en que hay que discutirlo explícitamente y en que la frase de que "la narrativa de más
blandas viene solo de células adherentes" no resiste ese contraejemplo. El matiz: Gossett midió
**derrames pleurales primarios por estiramiento hidrodinámico sin contacto**, que es otro régimen y
otro tipo de muestra. Lo defendible no es elegir bando sino declarar que **el signo del contraste
depende de técnica y de régimen**, y que nuestra conclusión está condicionada a los datos de pinzas
acústicas. Buscar RT-DC con ambas poblaciones medidas igual es la acción correcta y es viable.

### B2. Reología de ley de potencias (tu punto 4)
Sí a la frase en limitaciones, y se puede acotar analíticamente sin rehacer el solver. El estudio
de sensibilidad completo con ley de potencias es otro proyecto.

### B3. U* ≈ 3,9 µm/s demasiado lento
Correcto como objeción, pero el golpe no es entero: ese es el óptimo a β = 10 kPa/(mm/s), y ahí la
**meseta útil abarca 4,41 décadas de caudal**. Se puede operar mucho más rápido perdiendo poca
separación. Hay que presentarlo así en lugar de conceder el punto completo.

### B4. Validación experimental (tu punto 7)
De acuerdo en que es el punto más débil para BMMB. Tu chip de poro de 8 µm está dentro de la
envolvente y un ensayo mínimo MCF-7 frente a PBMC es exactamente lo que lo cierra.

---

## C. No entra, o no como está

### C1. El título que propones es seguro pero plano
*"A calibrated multi-compartment finite-element cell model for comparing microfluidic capture
modalities for breast circulating tumour cells"* describe el método, no la afirmación. Si A1 tumba
"not size", el título debería apoyarse en lo único que sobrevive intacto a toda tu revisión: **la
inversión de la regla de operación** — se fija la presión para que pase el leucocito y quede
retenida la tumoral. Eso no depende ni del tamaño ni del Sobol.

### C2. Tu resumen propuesto arrastra la afirmación que tu propio punto 1 derriba
Dice que DLD *"cannot resolve the size difference between the two populations"*. Con tamaños reales
de leucocito, sí la resuelve. Tú mismo lo señalas al final ("si el punto 1 cambia el resultado de
DLD, esa frase se ajusta"), pero entonces el resumen hay que rehacerlo **después** de resolver A1,
no antes. Las 215 palabras y el registro sí son el objetivo correcto.
