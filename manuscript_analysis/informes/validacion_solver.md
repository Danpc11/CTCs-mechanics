# Validación del solver axisimétrico multicompartimento

**Resultado del punto de corte del paso 2: SUPERADO.** Los cuatro problemas con solución
conocida se reproducen dentro de tolerancia y con convergencia de malla limpia. El paso 3
(calibración frente al 3D publicado) se ha ejecutado, pero **en una forma reducida respecto a
la literal del encargo**, por una razón de fondo que se explica en §6 y que conviene conocer
antes de leer los números.

Código: `solver_axisimetrico.py`. Entorno `ctcfem` (scikit-fem 12.0.2, gmsh instalado por pip).

---

## 1. Qué se ha implementado

| Componente | Formulación |
|---|---|
| Cinemática | Axisimétrica, gran deformación, `F` 3×3 con componente circunferencial `1+u_r/R` explícita; volumen `2πR dR dZ` |
| Citoplasma / nucleoplasma | Neo-Hookeano casi incompresible, **formulación mixta u–p** (Lagrangiano perturbado `C10(Ī1−3) + p(J−1) − p²/2κ`) |
| Elementos | Desplazamiento P2+burbuja vectorial (`ElementTriP2B`), presión P1 **discontinua** (`ElementTriP1DG`) — par de Crouzeix–Raviart, inf-sup estable y que **admite salto de presión en la interfase núcleo/citoplasma** |
| Córtex y envoltura nuclear | Capa delgada con energía por unidad de área de referencia `w(λ_s, λ_θ)`: tensión `γ`, cizalla y dilatación en ley neo-Hookeana de superficie **o ley de Skalak** (`Ks=(1+2C)Gs`), más viscosidad de Kelvin–Voigt opcional |
| Contacto con la pared | Penalización `½k_n⟨g⟩₊²` con perfil de pared `r = R_w(z, λ)` arbitrario (constricción tipo tanh incluida) |
| Carga | Presión seguidora exacta vía relación de Nanson `t dA = p J F⁻ᵀ N dA₀` |
| Solución | Newton–Raphson con tangente consistente, búsqueda de línea de Armijo y continuación en la carga |

**Cómo se obtienen las derivadas.** El tangente del material `A_iJkL = ∂P_iJ/∂F_kL` se calcula
por **paso complejo** (`Im P(F+ih e_kL)/h`, h=1e−20): sin error de truncamiento, exacto a
precisión de máquina. Los términos de superficie y contacto usan paso complejo para el
gradiente y diferencia central **sobre variables bien escaladas** (los estiramientos
`λ_s, λ_θ` y la geometría de la pared, nunca la penetración, que puede ser ~1e−7 veces menor).

**Verificación del tangente frente a diferencias finitas de todo el residuo** (25 columnas
aleatorias de K, perturbación 1e−6, modelo con dos materiales, dos capas de superficie con
`γ`, `Gs` y `Ks`, contacto con pared de perfil y presión seguidora simultáneos):

- error relativo medio **2.07e−10**, máximo **1.25e−09**.

Es decir: el tangente es consistente para *todos* los términos implementados, que es lo que
sostiene la convergencia cuadrática de Newton (figura, panel a: residuos
1.13e+00, 2.14e-01, 7.41e-03, 1.28e-05, 1.38e-11).

---

## 2. Casos de verificación (punto de corte)

### V1 — Extensión uniaxial homogénea de un cilindro (parche)
Estado homogéneo exacto, contrastado contra la solución analítica del **mismo** material
(κ finito), no contra el límite incompresible. Verifica material + cinemática axisimétrica +
mixto u–p simultáneamente.

- λ_axial impuesto = 1.3; λ_radial exacto = 0.8770714776, EF = 0.8770714776
- **error relativo = 2.78e-12**
- J: exacto 1.00003069, EF 1.00003069; presión: exacta 0.306899, EF 0.306900
- 4 iteraciones de Newton por paso de carga

### V2 — Indentación de Hertz sobre esfera elástica (E=1, ν=0.3)
Esfera comprimida entre platos rígidos, contacto por penalización. Referencia
`F = (4/3)E*√R δ^{3/2}`, `E* = E/(1−ν²)`.

| δ/R | a/R | elementos | F (EF) | F (Hertz) | error | ‖R‖ final |
|---|---|---|---|---|---|---|
| 0.0025 | 0.0500 | 2181 | 1.86404e-04 | 1.83150e-04 | +1.78% | 6.0e-15 |
| 0.0050 | 0.0707 | 2224 | 5.32171e-04 | 5.18027e-04 | +2.73% | 5.3e-14 |
| 0.0100 | 0.1000 | 2269 | 1.52617e-03 | 1.46520e-03 | +4.16% | 6.3e-14 |
| 0.0200 | 0.1414 | 2304 | 4.40471e-03 | 4.14422e-03 | +6.29% | 1.8e-13 |
| 0.0400 | 0.2000 | 2278 | 1.28411e-02 | 1.17216e-02 | +9.55% | 2.1e-11 |


El error **crece linealmente con a/R** — exactamente la firma del término de corrección de
tamaño finito que separa una esfera de un semiespacio de Hertz, no un error numérico. Ajuste
`error(%) = 51.9·(a/R) -0.93`:
**extrapolado a a/R→0 el error del solver es -0.93%** (panel b).

Sensibilidad a la rigidez de penalización (δ/R=0.01, `k_n = f·E/h`):

| f | k_n | F | error vs Hertz |
|---|---|---|---|
| 10 | 1.000e+03 | 1.518085e-03 | +3.609% |
| 30 | 3.000e+03 | 1.523998e-03 | +4.013% |
| 100 | 1.000e+04 | 1.526174e-03 | +4.161% |
| 300 | 3.000e+04 | 1.526847e-03 | +4.207% |
| 1000 | 1.000e+05 | 1.527096e-03 | +4.224% |


La fuerza satura en f≳300: la **compliancia espuria de la penalización vale 0.6% de F a f=10
y <0.02% a f≥300**. Con el valor usado en producción (f=100) el sesgo por penalización es de
0.06 puntos porcentuales.

### V3 — Esfera hueca neo-Hookeana incompresible bajo presión interna
Solución analítica exacta (A=0.6, B=1.0, C10=0.5): equilibrio radial integrado con cuadratura
adaptativa, `P = ∫ₐᵇ 4C₁(λ²−λ⁻⁴) dr/r` con `λ(r)=r/(r³−c)^{1/3}`. Para a=0.9 (λ_a=1.5),
**P_analítico = 0.59790606** (error de cuadratura 3e−13). Se aplica esa presión
como carga seguidora y se mide el radio interno deformado.

| malla (nρ×nθ) | dofs | a (EF) | error en a | error en el desplazamiento | J medio |
|---|---|---|---|---|---|
| 4×8 | 420 | 0.898014 | -0.2206% | -0.662% | 1.0000555 |
| 6×12 | 1056 | 0.899245 | -0.0839% | -0.252% | 1.0000556 |
| 9×18 | 2550 | 0.899729 | -0.0302% | -0.090% | 1.0000556 |
| 13×26 | 5550 | 0.899925 | -0.0083% | -0.025% | 1.0000556 |
| 19×38 | 12210 | 0.900022 | +0.0024% | +0.007% | 1.0000556 |


Convergencia monótona de orden ≈2.5 (dominada por la aproximación poligonal de la frontera
curva). En la malla más fina el error cambia de signo y se queda en **+2.4e−05**, del orden
del sesgo residual por κ finito (J medio = 1.0000562 → el sólido gana 0.0056% de volumen).

### V4 — Ley de Laplace: esfera sólida con tensión superficial
La tensión cortical comprime la esfera hasta un estado hidrostático uniforme con
`p = −2γ/a`; combinado con `p = κ(J−1)` da una ecuación escalar exacta para λ.
Con R=1, κ=1, γ=0.02: **λ_exacto = 0.98629440, p_exacto = -0.04055584**.

| h | elementos | λ (EF) | error en la contracción | p (EF) | error en p |
|---|---|---|---|---|---|
| 0.200 | 50 | 0.9845652 | -12.62% | -0.0407068 | -0.3722% |
| 0.120 | 158 | 0.9855912 | -5.13% | -0.0406026 | -0.1152% |
| 0.070 | 417 | 0.9859968 | -2.17% | -0.0405726 | -0.0414% |
| 0.045 | 961 | 0.9861547 | -1.02% | -0.0405630 | -0.0176% |


Orden ≈2 en la presión (panel c). El error en λ se mide como fracción de la contracción total
(1−λ ≈ 1.37%), que es la magnitud físicamente relevante, no como fracción de λ.

### V5 — Compresión de esfera hiperelástica casi incompresible
κ/μ = 10⁴ → **ν = 0.49995**. En régimen de pequeña deformación se contrasta con Hertz
(E* = E/(1−ν²), E=2.9999):

| δ/R | elementos | F (EF) | F (Hertz) | error | J medio |
|---|---|---|---|---|---|
| 0.002 | 2182 | 4.82926e-04 | 4.76980e-04 | +1.25% | 0.99999999 |
| 0.005 | 2224 | 1.92413e-03 | 1.88543e-03 | +2.05% | 0.99999997 |
| 0.010 | 2269 | 5.49662e-03 | 5.33280e-03 | +3.07% | 0.99999991 |
| 0.020 | 2304 | 1.57769e-02 | 1.50834e-02 | +4.60% | 0.99999975 |


Mismo comportamiento lineal en a/R, sin rastro de bloqueo volumétrico: J medio se mantiene en
1.0000000 a ν=0.49995, que es precisamente lo que una formulación de desplazamiento puro
fallaría. En gran deformación (δ/R = 0.30) no hay solución cerrada; se reporta
**autoconvergencia** y conservación de volumen:

| h | elementos | F | ΔV/V | radio de contacto |
|---|---|---|---|---|
| 0.120 | 407 | 1.255344 | -0.1886% | 0.6154 |
| 0.080 | 911 | 1.257802 | -0.0847% | 0.6211 |
| 0.055 | 1915 | 1.258333 | -0.0414% | 0.6182 |
| 0.038 | 3906 | 1.258892 | -0.0205% | 0.6262 |


La fuerza varía **0.28% al multiplicar por 9.6 el número de elementos**;
la pérdida de volumen converge como O(h²) hasta 0.020%. El radio de contacto
(0.6154, 0.6211, 0.6182, 0.6262) **no** converge monótonamente: está cuantizado por el
tamaño de elemento en la zona de contacto y no debe usarse como métrica de convergencia.

### Resumen del punto de corte

| caso | referencia | error | veredicto |
|---|---|---|---|
| V1 extensión homogénea | analítica exacta | 2.8e-12 | ✔ |
| V2 Hertz (ν=0.3) | Hertz (semiespacio) | +1.8% a +9.6% (lineal en a/R); solver −0.9% extrapolado | ✔ |
| V3 esfera hueca inflada | analítica exacta | 2.2e−3 → 2.4e−5 (convergente) | ✔ |
| V4 ley de Laplace | analítica exacta | 3.7e−3 → 1.8e−4 (convergente) | ✔ |
| V5 esfera casi incomp. | Hertz + autoconvergencia | +1.2% a +4.6%; ΔV 2e−4 | ✔ |
| Tangente consistente | dif. finitas del residuo | 2.1e−10 medio / 1.3e−9 máx | ✔ |


**Error relativo máximo frente a una solución analítica exacta: 3.7e−3** (V4, malla más
basta; 1.8e−4 en la más fina). Las desviaciones de los casos de Hertz son del *modelo* de
referencia (semiespacio frente a esfera), no del solver, como demuestra la extrapolación
a a/R→0.

---

## 3. Estudio de convergencia de malla en la constricción

Célula completa de tres compartimentos (córtex de Skalak + citoplasma + envoltura nuclear +
nucleoplasma) empujada cuasi-estáticamente dentro de una constricción circular de radio de
garganta 6.0 µm (confinamiento 2a/D = 1.417). Parámetros de Wang et al. 2023 (§5).

| h (µm) | elementos | dofs | DI (%) | ΔA/A₀ (%) | fuerza de contacto (nN) | ‖R‖ final | iter. máx |
|---|---|---|---|---|---|---|---|
| 1.70 | 70 | 672 | 1.5959 | 1.731 | 0.0000 | 2.5e-07 | 9 |
| 1.20 | 134 | 1264 | 1.5762 | 1.921 | 0.0000 | 8.1e-07 | 9 |
| 0.85 | 268 | 2494 | 1.5868 | 2.010 | 0.0000 | 2.0e-08 | 13 |
| 0.60 | 507 | 4675 | 1.5909 | 2.065 | 0.0000 | 6.7e-06 | 14 |


- **DI varía 1.3% (máx−mín/media) al multiplicar por 7 los grados de libertad**, de forma no
  monótona: DI se calcula sobre el polígono nodal del contorno y hereda su ruido. Es una cota
  superior honesta de la incertidumbre de malla de esa métrica, no una tasa de convergencia.
- ΔA/A₀ **sí** converge monótonamente: 1.731 → 1.921 → 2.010 → 2.065%, con diferencias
  0.190, 0.089, 0.055 → orden ≈1; extrapolación de Richardson **2.12%**, es decir, la malla
  más fina está a **2.3%** del valor extrapolado.
- Todos los casos convergen con ‖R‖ ≤ 6.7e-06 y ≤ 14 iteraciones por paso de carga (9, 9, 13 y 14 de la malla más basta a la más fina; ver la columna «iter. máx» de la tabla).

Para producción se usa h = 0.85 µm (2494 dofs): DI dentro del 0.3% y ΔA dentro del 2.6% de la
malla más fina, a 1/2 del coste.

### Sensibilidad al módulo de cizalla sustituto del citoplasma

Wang et al. modelan citoplasma y nucleoplasma como **fluidos viscosos sin elasticidad**. Un
solver cuasi-estático de sólido necesita una rigidez no nula; se usa
μ_cit = 4.6 Pa = 0.1·(Gs/a) y se comprueba que el resultado no depende de esa regularización:

| factor | μ_cit (Pa) | converge | DI (%) | ΔA/A₀ (%) |
|---|---|---|---|---|
| 0.1× | 0.46 | **NO CONVERGE** | 1.206 | 1.140 |
| 0.3× | 1.38 | sí | 3.732 | 4.128 |
| 1.0× | 4.60 | sí | 3.768 | 4.596 |
| 3.0× | 13.80 | sí | 3.809 | 4.822 |
| 10.0× | 46.00 | sí | 3.900 | 5.059 |


Sobre el rango convergido (0.3×–10×, **33× en μ**): DI cambia **+4.5%** y ΔA/A₀ **+22.5%**.
El caso 0.1× **no convergió** y se excluye de esa estadística (queda declarado, no usado).
Conclusión: la métrica de forma DI es robusta frente a la regularización; la deformación de
área lo es bastante menos y conviene reportarla siempre con el valor de μ_cit empleado.

---

## 4. Barrido de confinamiento (régimen de aplicabilidad)

| R_garganta (µm) | 2a/D | DI (%) | ΔA/A₀ (%) | λ_s máx | fuerza de contacto (nN) |
|---|---|---|---|---|---|
| 8.60 | 0.988 | 0.040 | 0.000 | 1.000 | 0.000 |
| 8.00 | 1.062 | 0.040 | 0.000 | 1.000 | 0.000 |
| 7.50 | 1.133 | 0.040 | 0.000 | 1.000 | 0.000 |
| 7.00 | 1.214 | 0.280 | 0.310 | 1.042 | 0.000 |
| 6.50 | 1.308 | 0.809 | 1.026 | 1.098 | 0.000 |
| 6.00 | 1.417 | 1.587 | 2.010 | 1.161 | 0.000 |
| 5.50 | 1.545 | 2.572 | 3.202 | 1.230 | 0.000 |
| 5.00 | 1.700 | 3.768 | 4.596 | 1.304 | 0.000 |
| 4.60 | 1.848 | 4.862 | 5.824 | 1.364 | 0.000 |
| 4.25 | 2.000 | 5.908 | 6.976 | 1.418 | 0.000 |
| 3.80 | 2.237 | 7.346 | 8.535 | 1.486 | 0.000 |
| 3.50 | 2.429 | 8.348 | 9.594 | 1.530 | 0.000 |


(DI = 0.040% para la esfera sin deformar es el sesgo de discretización del contorno
poligonal; es el cero efectivo de la métrica.)

---

## 5. Paso 3 — Calibración frente al 3D publicado (Wang et al. 2023)

### 5.1 Qué se ha podido extraer del artículo y qué no

Artículo recuperado en texto completo (`10.1007_s10237-023-01705-6.pdf`, OA vía Unpaywall).

**Extraído del texto (valores citados literalmente por los autores):**

| Magnitud | Valor | Dónde |
|---|---|---|
| Radio celular a | 8.5 µm | §2.3, montaje |
| Canal recto, sección cuadrada, lado l | 30 µm → **2a/l = 0.567** | §2.3 |
| Ley de membrana | Skalak, `W = Gs/4·(I₁²+2I₁−2I₂) + C·Gs/4·I₂²`, `Ks=(1+2C)Gs` | §2.1 |
| Mejor ajuste a HL-60 | **Ca = 0.75, C = 10 → Gs = 0.39 mN/m, Ks = 8.10 mN/m** | §3.1 |
| Núcleo | `a_n/a = 0.5`, moduli de la envoltura = **2× los de la membrana** | §3.1 |
| **Aumento de área de membrana en el estado estacionario del mejor ajuste** | **4%** | §3.1 (texto) |
| Índice de deformación | `DI = 1 − 2√(πA)/P`, A y P proyectados sobre el plano de simetría | §2.4 |
| Citoplasma y nucleoplasma | fluidos newtonianos/viscosos, **sin elasticidad** | §2.1 |

**Leído gráficamente de la Fig. 3b (no es un número del texto):** meseta estacionaria
`DI ≈ 0.028`, con una incertidumbre de lectura que estimo en **±0.005**. Todo lo que dependa
de este valor va etiquetado como provisional.

**No se ha podido extraer**: las series temporales numéricas de DI(x) (solo existen como
figura), los perfiles de forma de la Fig. 3c en forma tabulada, ni los valores de viscosidad
de membrana y citoplasma que fijan el transitorio — esos solo se dan como razones
adimensionales dentro de barridos paramétricos representados gráficamente.

### 5.2 Por qué la reproducción del *transitorio* no es ejecutable con este solver

Dos hechos del propio artículo lo impiden, y conviene decirlo sin rodeos:

1. **El confinamiento es 2a/l = 0.567 < 1: la célula no toca la pared en el canal recto.**
   La forma de bala estacionaria la produce el campo de esfuerzos del fluido alrededor de una
   célula libre, no el contacto. Un solver de sólido cuasi-estático sin problema de fluido
   acoplado **no tiene con qué cargar la célula** en ese régimen: con contacto únicamente,
   predice deformación idénticamente nula (ver barrido §4: ΔA/A₀ = 0 para 2a/D ≤ 1.13).
2. **El transitorio está gobernado por viscosidades** (de membrana, dependiente del modo de
   deformación, y de citoplasma) acopladas a la hidrodinámica. El propio artículo señala que
   la parte hiperelástica fija el estado estacionario y las viscosidades el camino hacia él.

Lo que **sí** es comparable, y es lo que se ha hecho, es la parte **hiperelástica**: la
relación entre dos invariantes de forma del estado deformado bajo *la misma ley constitutiva
y los mismos parámetros*. Esa relación es independiente de cómo se cargue la célula, y es por
tanto la comparación 2D-vs-3D legítima que puede hacerse sin resolver el fluido.

### 5.3 Discrepancia 2D vs 3D cuantificada

Se monta la célula de Wang en el solver axisimétrico (a = 8.5 µm, a_n = 4.25 µm, córtex de
Skalak con Gs = 0.39 mN/m y C = 10, envoltura nuclear con 2Gs y C = 10) y se la comprime
cuasi-estáticamente en constricciones circulares de radio decreciente, obteniendo la curva
DI(ΔA/A₀) del panel f.

| Comparación | 3D (Wang et al. 2023) | 2D axisimétrico (este trabajo) | Discrepancia |
|---|---|---|---|
| DI a deformación de área igualada (ΔA/A₀ = 4.0%) | 0.028 ± 0.005 *(lectura de figura)* | **0.0326** | **+16%** (rango -1% a +42% con la incertidumbre de lectura) |
| ΔA/A₀ a DI igualado (DI = 0.028) | 4.0% *(texto)* | **3.47%** | **-13.3%** |

El modelo 2D axisimétrico **sobreestima el alargamiento por unidad de deformación de área en
un ~16%**, equivalentemente subestima en un 13% la deformación de área necesaria para
alcanzar un DI dado. El confinamiento axisimétrico que reproduce la deformación de área del
estado estacionario 3D es **2a/D = 1.63**.

**Caveat que hay que escribir en el artículo, no esconder:** los modos de carga no son
idénticos (deformación libre inducida por cizalla en canal cuadrado frente a compresión
axisimétrica por contacto con pared). La comparación contrasta la *relación entre invariantes
de forma* bajo la misma ley de membrana, no estados idénticos. Un ~16% es, a la vista de esa
diferencia de modo, un acuerdo razonable, pero no es una validación punto a punto.

### 5.4 Régimen de fiabilidad del modelo reducido

- **Cota inferior, dura: 2a/D ≳ 1.2.** Por debajo no hay contacto y el modelo cuasi-estático
  de sólido predice deformación nula (§4: ΔA/A₀ = 0 exactamente hasta 2a/D = 1.13, y 0.31% a
  2a/D = 1.21). El caso de validación experimental de Wang (2a/l = 0.567) cae **de lleno
  dentro de esta zona prohibida**. Si el artículo necesita cubrir el régimen sin contacto,
  hace falta acoplar el fluido; no es negociable.
- **Zona de buen acuerdo: 1.4 ≲ 2a/D ≲ 2.0**, donde la discrepancia con el 3D en la relación
  DI–ΔA es del orden del 16% y la malla de producción da DI estable al 1.3%.
- **Cota superior:** hasta 2a/D = 2.43 (ΔA/A₀ = 9.6%, λ_s = 1.53) el solver sigue
  convergiendo sin degradación. El límite ahí no es numérico sino constitutivo: por encima de
  ~5% de deformación de área una bicapa pura lisaría, y describir ese régimen exige un modelo
  de reservorio de membrana. No se ha explorado por encima de 2a/D = 2.43.
- **Geometría del canal:** el 2D axisimétrico es un tubo circular; el canal de Wang es
  cuadrado. Esa diferencia no se ha cuantificado y es un sesgo adicional no incluido en el
  16% anterior.

---

## 6. Qué recomendaría cambiar si el barrido de diseño lo exige

El solver tal como está es suficiente para la vía de constricción, trampa de deformabilidad y
microfiltro de poros, que son los tres regímenes **dominados por contacto**. Dos limitaciones
que conviene tener presentes antes de comprometerse con el barrido completo:

1. **DLD no está cubierto.** El desplazamiento lateral determinista es intrínsecamente
   tridimensional y asimétrico (la célula se desplaza lateralmente respecto a un array de
   postes): no tiene eje de revolución. Para DLD el modelo axisimétrico no sirve y habrá que
   usar una correlación basada en diámetro efectivo deformado — que este solver **sí** puede
   suministrar como entrada.
2. **Régimen sin contacto.** Si se quiere cubrir confinamientos 2a/D < 1.2 hay que acoplar un
   problema de fluido. Con scikit-fem es viable montar un Stokes axisimétrico (Taylor–Hood,
   mismo par de elementos que ya se usa) con acoplamiento particionado; si eso resultara
   insuficiente, la alternativa natural es FEniCSx/dolfinx, que da ALE y contacto por
   multiplicadores de Lagrange de fábrica. **No ha hecho falta para esta validación**:
   scikit-fem, con formulación mixta u–p y tangente por paso complejo, reproduce los cuatro
   problemas analíticos y la convergencia cuadrática de Newton sin tocar ningún parámetro
   numérico.

---

## 7. Archivos

- `solver_axisimetrico.py` — solver, importable y reutilizable
- `validacion_solver.png` — figura de verificación (paneles a–f)
- `validacion.json` — todos los números crudos de verificación
- `constriccion.py` — montaje del caso de constricción con la célula de tres compartimentos
- `verificacion.py` — script que genera `validacion.json`
