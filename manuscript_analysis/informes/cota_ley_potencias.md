# Cota analítica para reología de ley de potencias

**Conclusión en una línea.** Con reología de ley de potencias **el umbral de diseño `β_c` deja de
existir como umbral**: formalmente colapsa a cero y toda célula acaba pasando a algún caudal. Lo que
sobrevive es la **ventana relativa** entre leucocito y tumoral, que se comprime de los
**165×–1646× de Kelvin–Voigt a 3,3×–17,2×** con los exponentes medidos. La conclusión «existe
ventana de operación» aguanta; la conclusión «la ventana es de dos a tres órdenes de magnitud» no.

---

## 1. De dónde sale `β_c` y por qué depende de Kelvin–Voigt

Con `σ = G(ε + τ·dε/dt)` y una rampa de carga de duración `t_L = L_carga/U`, la razón
`dε/dt / ε = 1/t_L` es constante, así que la rigidez aparente es `G(1 + De)` con
`De = τ/t_L = τ·U/L_carga`. Luego

> `ΔP_crít(U) = ΔP_el·(1 + De)`  (Kelvin–Voigt)

y, frente a una presión motriz `ΔP = β·U`, el margen de paso es
`m(U) = β·U / [ΔP_el(1+De)] = m_∞·De/(1+De)` con `m_∞ = β·L_carga/(ΔP_el·τ)`.

El punto clave es que **`ΔP_crít` crece LINEALMENTE en `U`**: numerador y denominador del margen
crecen al mismo ritmo asintótico, el margen **satura** en `m_∞`, y por eso existe un umbral finito

> `β_c = ΔP_el·τ / L_carga`,  tal que si `β ≤ β_c` la célula **no escapa a ningún caudal**.

Esa saturación es un artefacto del término viscoso lineal. Es el único ingrediente del que depende
la existencia de `β_c`.

## 2. La misma cuenta con ley de potencias

Sustitúyase el elemento viscoso por un material de ley de potencias,
`G*(ω) = G_0·(ω/ω_0)^α` con `0 < α < 1`, equivalentemente módulo de relajación
`G(t) = G_0·(t/t_0)^(−α)/Γ(1−α)`. Bajo la MISMA rampa de deformación `ε(t) = ε̇·t`:

> `σ(t_L) = ε̇·∫₀^{t_L} G(t_L−s)·ds = ε(t_L)·G_0·(t_0/t_L)^α / Γ(2−α)`

es decir, rigidez aparente `G_app ∝ t_L^(−α) ∝ U^α`. Normalizando en el tiempo de carga de
referencia `t_ref` en el que se recupera el valor cuasi-estático (`U_ref = L_carga/t_ref`):

> `ΔP_crít(U) = ΔP_el·(U/U_ref)^α`   (ley de potencias)

La condición de escape `β·U ≥ ΔP_crít(U)` da ahora

> `U_esc = [ ΔP_el / (β·U_ref^α) ]^{1/(1−α)}`

que es **finita para todo `β > 0` siempre que `α < 1`**. El margen
`m(U) = (β/ΔP_el)·U_ref^α·U^{1−α}` **crece sin cota** en vez de saturar.

> **Cota 1 (existencia).** Con `α < 1` no hay umbral: `β_c → 0`. El enunciado «si `β < β_c` la
> célula no pasa a ninguna presión» es **específico de Kelvin–Voigt** y debe declararse como tal.

## 3. Cuánto se desplaza el umbral si se acota el caudal máximo del dispositivo

Un dispositivo real tiene un caudal máximo `U_max`. Exigiendo `U_esc ≤ U_max` se recupera un umbral
**efectivo**:

> `β_c^eff = ΔP_el·U_ref^(−α)·U_max^(α−1)`,  y por tanto
> **`β_c^eff / β_c^KV = Λ^{−(1−α)}` con `Λ = U_max/U_ref`**

(con `β_c^KV = ΔP_el·τ/L_carga` y `τ ≡ t_ref`). El límite `α → 1` devuelve exactamente `β_c^KV`, y
`α → 0` devuelve `ΔP_el/U_max`, el umbral puramente elástico. La cota es por tanto continua y
contiene los dos casos límite.

**Desplazamiento por década de rango dinámico de caudal:**

| α | factor por el que BAJA `β_c` por cada década de `Λ` |
|---|---|
| 0,147 (MDA-MB-231, medido) | **7,13×** |
| 0,186 (MCF-10A, medido) | 6,52× |
| 0,234 (MCF-7, medido) | 5,83× |
| 0,25 (THP-1 por AFM) | 5,62× |
| 0,44 (THP-1 por SF-DC) | 3,63× |
| 0,48 (neutrófilo, Tsai 1993) | **3,31×** |

Con tres décadas de rango de caudal y α = 0,2, `β_c^eff` cae a
`Λ^{−(1−α)} = 1000^{−0,8} = 3,98 × 10⁻³` de su valor de Kelvin–Voigt: un factor **≈ 251**
(= 6,31³, con 6,31 el factor por década de la tabla anterior interpolado a α = 0,2).

## 4. Lo que de verdad importa: la ventana entre poblaciones

La afirmación del artículo no es sobre `β_c` absoluto sino sobre la **razón**
`β_c(tumoral)/β_c(leucocito)`. Con exponentes propios de cada población:

> `β_c^eff,t / β_c^eff,l = (ΔP_el,t/ΔP_el,l) · (τ_t^{α_t} / τ_l^{α_l}) · (U_max/L_carga)^{α_t−α_l}`

que con `α_t = α_l = α` se reduce a

> **`ventana(α) = (ΔP_el,t/ΔP_el,l) · (τ_t/τ_l)^α`**  frente a  `ventana(KV) = (ΔP_el,t/ΔP_el,l) · (τ_t/τ_l)`

Es decir: **el contraste de `τ` entra elevado a α en vez de linealmente.** El factor de compresión es
exactamente `(τ_t/τ_l)^{1−α}`.

Con `L_carga = 3,06 µm`, `U_max = 1 mm/s` y los parámetros de la Fase 2:

| pareja | α tumoral | α leucocito | ventana Kelvin–Voigt | ventana ley de potencias | factor de compresión |
|---|---|---|---|---|---|
| MDA-MB-231 / PBMC | 0,147 | 0,147 | 165,6× | **4,62×** | 35,9 |
| MDA-MB-231 / PBMC | 0,147 | 0,25 (AFM THP-1) | 165,6× | **4,15×** | 39,9 |
| MDA-MB-231 / PBMC | 0,147 | 0,44 (SF-DC THP-1) | 165,6× | **3,41×** | 48,5 |
| MDA-MB-231 / PBMC | 0,147 | 0,48 (neutrófilo) | 165,6× | **3,28×** | 50,6 |
| MCF-7 / PBMC | 0,234 | 0,234 | 1646,4× | **17,51×** | 94,0 |
| MCF-7 / PBMC | 0,234 | 0,25 | 1646,4× | **17,23×** | 95,6 |
| MCF-7 / PBMC | 0,234 | 0,44 | 1646,4× | **14,16×** | 116,2 |
| MCF-7 / PBMC | 0,234 | 0,48 | 1646,4× | **13,59×** | 121,1 |

> **Cota 2 (magnitud).** Para la línea tumoral más blanda, la ventana cae de **165× a 3,3–4,6×**;
> para MCF-7, de **1646× a 13,6–17,5×**. La ventana **sigue existiendo y sigue superando 2×** en
> todos los escenarios de exponente examinados, pero pierde entre 1,5 y 2 órdenes de magnitud.
>
> Nótese además el signo del tercer factor: como el exponente medido del leucocito
> (α = 0,25–0,48) es **mayor** que el de la línea tumoral (α = 0,147–0,234), el leucocito se endurece
> **más deprisa** con el caudal, y la ventana se estrecha *un poco más* al subir `U_max`. El efecto
> es modesto (factor 1,1–1,4 sobre el caso de α común) pero va en contra, no a favor.

## 5. Exponentes usados y su procedencia

| población | α | técnica | fuente |
|---|---|---|---|
| MCF-10A | **0,186 ± 0,049** | AFM, modelo de Ting de ley de potencias, 0,1–10 µm/s | 10.1038/s42003-020-01330-4 |
| MCF-7 | **0,234 ± 0,060** | idem | 10.1038/s42003-020-01330-4 |
| MDA-MB-231 | **0,147 ± 0,061** | idem | 10.1038/s42003-020-01330-4 |
| THP-1 (monocítica), **en suspensión** | **0,44** | citometría de deformabilidad por flujo cortante (SF-DC) | 10.7554/eLife.78823 |
| THP-1, adherida | **0,25** | AFM a ~1 s⁻¹, mismo trabajo | 10.7554/eLife.78823 |
| neutrófilo humano | **0,48** (= 1 − b, con b = 0,52 ± 0,09) | aspiración por micropipeta, viscosidad citoplásmica frente a cizalla media | 10.1016/S0006-3495(93)81238-4 |
| rango general de células | 0,1–0,3 (G′ y G″) | microrreología oscilatoria, múltiples técnicas | 10.1126/sciadv.abn6093 |

Los tres exponentes de las líneas de mama están medidos **sobre nuestras propias líneas**, lo que
hace esta cota mucho menos especulativa de lo que sería con un α genérico.

## 6. Lo que esta cota NO hace

- No es un estudio de sensibilidad: el solver sigue siendo Kelvin–Voigt y no se ha rehecho.
- La cota supone que la **geometría de la barrera elástica no cambia**, es decir que `ΔP_el` sigue
  siendo el mismo y solo cambia el factor de endurecimiento por velocidad. Con ley de potencias
  estricta no existe un límite elástico verdadero (`G → ∞` cuando `t_L → 0` y `G → 0` cuando
  `t_L → ∞`); `ΔP_el` se reinterpreta como el valor en `t_ref`, y `t_ref` queda declarado, no medido.
- **No hay α medido de PBMC ni de granulocito humano en suspensión.** Se usa THP-1 (línea monocítica)
  y el neutrófilo de Tsai 1993 en gran deformación. Es el mismo hueco que ya tenía `τ_leucocito`.
