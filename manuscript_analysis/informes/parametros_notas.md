# Parámetros mecánicos por compartimento y por subtipo — procedencia, transferencias y huecos

Documento de respaldo de `parametros_compartimento.csv` (86 filas; 53 marcadas como transferidas).
Regla aplicada sin excepción: **ningún valor numérico se inventa ni se interpola**. Cuando una fuente
publica un rango y la simulación necesita un escalar, el `valor_nominal` es el punto medio de **ese rango
publicado** y la fila lleva `nominal_es_punto_medio = True`; esto es una convención declarada, no un dato.

---

## 1. Cómo se obtuvo la literatura

Dos barridos independientes, con los vocabularios deliberadamente separados (el vocabulario
«circulating tumor cell» + «capture» no alcanza esta literatura, como estaba previsto):

- **PubMed** (`search_articles`): 24 consultas en dos rondas — tránsito por microcanal constreñido,
  deformación inducida por flujo, núcleos aislados, lamina A/C, aspiración por micropipeta del núcleo,
  AFM en líneas de mama, mecánica de la envoltura nuclear, tensión cortical, viscosidad citoplasmática,
  relación N/C, tamaño celular, viscosidad de membrana, tiempos de relajación. 810 PMID únicos,
  784 con metadatos recuperables.
- **OpenAlex** (`openalex_search_works`): 11 consultas equivalentes, 480 registros, 203 con resumen
  reutilizable; 87 aportaban registros no cubiertos por PubMed.
- Cribado combinado de **632 resúmenes** por extracción determinista de patrones
  «número + unidad mecánica» (kPa, Pa, mN/m, µN·s/m, Pa·s, dyn/cm, µm) cruzada con menciones de línea
  celular y de compartimento. 113 resúmenes contenían al menos un valor con unidad mecánica;
  34 combinaban línea celular del panel con terminología mecánica.
- De ahí se seleccionaron 16 artículos para texto completo; **12 se obtuvieron** (Unpaywall/PMC) y
  4 resultaron inaccesibles (ver §6).

---

## 2. Parametrización de Wang et al. — transcripción literal

Es la referencia obligada y el punto de calibración del solver 2D. Todo lo que sigue está transcrito
tal como lo publican los autores.

### 2.1 Wang et al. 2023, *Biomech. Model. Mechanobiol.* 22:1129–1143 — doi 10.1007/s10237-023-01705-6

Célula HL-60 (leucemia humana), geometría y flujo tomados del experimento de Fregin et al. 2019.

**Geometría y flujo**
- «The cell has a radius of `a = 8.5 μm` and is released with zero velocity at `x = −190 μm`.»
- Constricción con ángulo convergente constante de 45°, seguida de canal recto de sección cuadrada
  de anchura `l = 30 μm` y longitud 300 μm, conectado a una sección divergente de 45°.
- Razón de confinamiento `2a/l = 0.57`.
- Fluido del canal: PBS con 1 % (p/v) metilcelulosa, ley de potencia con `m = 0.60 Pa s`,
  `α = 0.64`, `γ̇₀ = 1 s⁻¹`. Densidad del fluido `ρ = 1065 kg m⁻³`.
- «A typical flow rate is `8 nL s⁻¹`, corresponding to an average flow speed in the narrow straight
  channel with `V = 0.89 cm s⁻¹`. The average fluid viscosity `μ₀`, calculated in the channel cross
  section at x = 150 μm, is `32.5 mPa s`. These lead to a flow Reynolds number of about `0.01`.»

**Membrana (bicapa + córtex, viscoelástica)**
- Números adimensionales de control: número capilar `Ca = μ₀V/Gs`; parámetro de dilatación de área
  de la ley de Skalak `C`, con `Ks = (1 + 2C) Gs` (Wang 2023: C = 10, Gs = 0.39, Ks = 8.10 mN/m lo confirma); razón de viscosidades citoplasma/fluido `λ`; confinamiento `2a/l`;
  razón de tamaños núcleo/célula `aₙ/a`.
- Mejor ajuste: «the parameters of `Ca = 0.75` and `C = 10` correspond to cell membrane elastic
  moduli values of `Gs = 0.39 mN m⁻¹` and `Ks = 8.10 mN m⁻¹`.»
  *(Comprobación propia: Gs = μ₀V/Ca = 0.0325 × 0.0089 / 0.75 = 0.386 mN m⁻¹ — consistente.)*
- Viscosidad de membrana dependiente del modo de deformación:
  `η_s^e/(μ₀a) = 100`, `η_s^s/(μ₀a) = 20`. Redimensionando con μ₀ = 32.5 mPa·s y a = 8.5 µm:
  **η_s^s = 5.53 µN·s·m⁻¹** y **η_s^e = 27.6 µN·s·m⁻¹**. (En la Fig. 7 usan además `η* = 40`.)
- Flexión: «A small bending resistance `k_c = 0.001 Gs a²` has been used in the present study to
  prevent membrane wrinkles.» Se desprecia la viscosidad de dilatación de área.

**Citoplasma**
- «A unity viscosity ratio seems to give the best fit. In the present setup, this corresponds to a
  cytoplasm viscosity of `32.5 mPa s`. It has been found in previous studies that the effective
  cytoplasmic viscosity of biological cells ranges roughly from `10 to 100 mPa s`
  (Luby-Phelps 1999; Mogilner and Manhart 2018).»

**Núcleo (microcápsula menor y más rígida)**
- «the nucleus has a size ratio of `aₙ/a = 0.5`. The elastic moduli of the cell nucleus membrane are
  also assumed to be twice those of the cell membrane» → `Gsn/Gs = 2`.
- «according to Rowat et al (2013), the nucleus radius is typically `3.5–5 μm`, leading to a size
  ratio of aₙ/a in the range of `0.4 to 0.6`.»

### 2.2 Wang et al. 2025, *Biomech. Model. Mechanobiol.* 24:487–505 — doi 10.1007/s10237-024-01920-9

Tabla 1, «Inferred mechanical properties of PC-3 and K-562 cells», transcrita literalmente:

| | Canal constreñido PC-3 | Canal constreñido K-562 | Cross-slot PC-3 | Cross-slot K-562 |
|---|---|---|---|---|
| `a` (μm) | 7.4−8.6 | 5.4–6 | 6.8−9.8 | 6.2−6.8 |
| `Gs` (mN/m) | 10.4−20.7 | 4.4−5.2 | 5.8−14.4 | 3.6−4.8 |
| `Ks` (mN/m) | 217–434 | 93–109 | 120.3−302.6 | 75.6−100.8 |
| `η_s` (μN s m⁻¹) | 9.7−22.5 | 7.8–14 | 8.7−16.1 | 6.9−8.4 |

- Incrementos del barrido paramétrico: `ΔCa = 0.1`, `ΔC = 1`, `Δη = 5`. Incertidumbre relativa media
  resultante de Gs, Ks y η_s: 10.6 %, 11.8 % y 4.9 % (PC-3); 3.8 %, 6.1 % y 5.6 % (K-562).
- Hipótesis mantenidas: `λ` = 1 (dependencia «mild» en el rango `0.2 < λ < 2`), `aₙ/a = 0.5`
  (mejor ajuste), `Gsn/Gs = 2` — y comprueban que una cápsula nuclear 8 veces más rígida «does not
  significantly reduce the cell deformation».
- Contexto que los propios autores citan para situar sus valores: módulo de Young de PC-3
  `1.4−1.97 kPa` (Faria et al. 2008; Lekka 2016) frente a K-562 `0.42−0.60 kPa`
  (Lange et al. 2016; Islam et al. 2018); viscosidad de membrana de granulocitos
  `η_s = 8.5 μN s m⁻¹` (Evans y Yeung 1989) y de eritrocito `0.1–1 μN s m⁻¹`.
- Diferencia metodológica relevante para nuestro 2D: en 2025 **no** necesitan viscosidad de membrana
  dependiente del modo de deformación, porque cada una de las dos geometrías genera por sí sola modos
  de corte y de extensión.

---

## 3. De dónde sale cada bloque de la tabla

| Bloque | Fuente primaria | Técnica | Estado |
|---|---|---|---|
| Radio celular MCF-7, MDA-MB-231 | 10.1039/c7ra10750a | análisis de imagen en microcanal confinante | medido |
| Radio celular SK-BR-3 | 10.1038/s41378-020-0150-6 | valor agregado del panel («~20 µm») | parcial |
| Radio celular MCF-10A | 10.1039/c7ra10750a | — | transferido desde MCF-7 |
| Radio celular y nuclear leucocito | 10.1007/s10237-023-01705-6 / 10.1074/jbc.M112.441535 | imagen | medido |
| aₙ/a y Gsn/Gs (todas las líneas de mama) | Wang 2023/2025 | ajuste del modelo | transferido |
| Gs, Ks, η_s de membrana (líneas de mama) | 10.1007/s10237-024-01920-9 (PC-3) | inferencia por ajuste IB-LBM | transferido |
| Gs, Ks, η_s de membrana (leucocito) | Wang 2023 (HL-60) + Wang 2025 (K-562) | ajuste IB-LBM | medido/inferido |
| Tensión cortical y reserva de área | 10.1016/S0006-3495(89)82660-8 | aspiración por micropipeta | medido (granulocito) |
| Viscosidad, módulo y τ de citoplasma (MCF-10A, MCF-7, MDA-MB-231) | 10.1039/d4na00003j | espectroscopia rotacional magnética | **medido por subtipo** |
| Citoplasma SK-BR-3 | 10.1039/d4na00003j | — | transferido desde MCF-10A |
| Viscosidad de citoplasma leucocitario | 10.1016/S0006-3495(89)82660-8 | micropipeta | medido |
| Lamina nuclear (25 mN/m) | 10.1242/jcs.01357 | micropipeta de envolturas aisladas | transferido (Xenopus) |
| Constante elástica nuclear y endurecimiento | 10.1091/mbc.e16-09-0653 | micromanipulación de núcleos aislados | transferido (HeLa/MEF) |
| Razón núcleo/citoplasma de rigidez y viscosidad | 10.1016/j.jbiomech.2009.01.024 y 10.1006/bbrc.2000.2360 | MEF + micropipeta | transferido (condrocito) |
| Dirección del efecto mesenquimal sobre el núcleo | 10.1038/s41598-026-56141-0 | mapeo de fuerza por AFM | transferido (U251) |
| Módulo de Young de célula entera del panel de mama | 10.1038/s41378-020-0150-6 | pinzas acústicas de haz único | **medido por subtipo** |
| Célula entera leucocitaria (HL-60, granulocito, PBMC) | 10.1038/s41467-019-08370-3 | dRT-DC | medido |

---

## 4. Transferencias: qué se transfiere, desde dónde y por qué

1. **Membrana/córtex de las cuatro líneas de mama ← PC-3.** No existe ninguna medida de módulo de
   membrana o de córtex para líneas de mama expresada en N/m. PC-3 es la **única línea de tumor sólido
   epitelial parametrizada dentro de este mismo modelo multicompartimento**, lo que hace que la
   transferencia sea compatible por construcción con el solver. El rango une las dos geometrías de
   canal publicadas (constreñido y cross-slot): Gs 5.8–20.7 mN/m, Ks 120.3–434 mN/m, η_s 8.7–22.5 µN·s/m.
2. **aₙ/a = 0.5 para las cuatro líneas de mama ← HL-60/PC-3/K-562.** Es el valor de mejor ajuste,
   idéntico en 2023 y 2025 para tres líneas de origen distinto; el rango 0.4–0.6 procede del radio
   nuclear medido de HL-60 (3.5–5 µm sobre a = 8.5 µm).
3. **Lamina nuclear ← ovocito de Xenopus (25 mN/m).** Es la única medida directa del módulo elástico
   de la envoltura nuclear tratada como **red bidimensional**, es decir, en la misma magnitud física
   que los módulos de membrana del modelo. No hay equivalente humano.
4. **Rigidez nuclear a deformación pequeña ← HeLa y MEF vimentina-nula (0.41–0.70 nN/µm).** Núcleos
   aislados humanos y murinos; es la escala de rigidez nuclear mejor caracterizada a deformaciones
   fisiológicas (<3 µm, ~30 % de deformación).
5. **Razones núcleo/citoplasma ← condrocito articular.** El único par de fuentes que mide los dos
   compartimentos del **mismo** tipo celular: 1.4 in situ (E_cito = 3.75 kPa, E_núc = 5.25 kPa) frente
   a 3–4× sobre núcleos aislados. La horquilla 1.4–4 **es** la discrepancia publicada, no una
   incertidumbre inventada.
6. **Citoplasma de SK-BR-3 ← MCF-10A.** SK-BR-3 es epitelial y no mesenquimal, igual que MCF-10A y
   MCF-7; el rango abarca las tres líneas medidas en el mismo trabajo (10.7–56.4 Pa·s).
7. **Efecto mesenquimal sobre el núcleo ← U251 (glioblastoma).** Se transfiere la **dirección y la
   magnitud relativa** (1278 → 682.6 Pa, factor 0.53), no el valor absoluto. Sirve para justificar que
   MDA-MB-231 reciba un núcleo más blando que MCF-7/SK-BR-3 en el análisis de sensibilidad.
8. **Leucocito ← HL-60 + granulocito + PBMC.** El leucocito de referencia se construye con HL-60
   (línea sobre la que Wang et al. calibran) y se acota con granulocito primario y PBMC medidos por
   dRT-DC en el mismo trabajo.

### SK-BR-3 (HER2+): prácticamente todo transferido

Se buscó explícitamente («SK-BR-3»/«SKBR3» cruzado con rigidez, deformabilidad, módulo, biomecánica,
AFM, reología, viscoelasticidad, tamaño celular y núcleo: 123 registros). **El barrido completo
devolvió un único parámetro mecánico de SK-BR-3**: módulo de Young aparente de 2.772 kPa por pinzas
acústicas de haz único (10.1038/s41378-020-0150-6), más el diámetro agregado del panel (~20 µm) del
mismo trabajo. Todo lo demás de SK-BR-3 en la tabla está transferido y marcado como tal. No se
rellenó ninguna celda con un valor intermedio entre MCF-7 y MDA-MB-231.

---

## 5. Rangos para el análisis de sensibilidad

Los `minimo`/`maximo` no son un porcentaje alrededor del nominal. Proceden de:

- **Dispersión intra-estudio publicada** (±SD o ±SE): citoplasma de MCF-10A/MCF-7/MDA-MB-231,
  módulos de célula entera de HL-60/granulocito/PBMC.
- **Dispersión entre determinaciones del mismo trabajo**: el estudio de microrreología da dos
  determinaciones independientes (η_MCF-10A = 36.3 y 41.6 Pa·s; η_MCF-7 = 65.9 y 56.4; η_MDA-MB-231 =
  12.0 y 10.7).
- **Dependencia con la escala de sonda, medida**: MCF-10A 15→160 Pa·s y MDA-MB-231 3→25 Pa·s al variar
  la longitud del microhilo de 2 a 6 µm. Este es el rango honesto para la viscosidad citoplasmática y
  es mucho mayor que la barra de error nominal.
- **Dispersión entre geometrías de dispositivo**: Gs, Ks y η_s de PC-3 y K-562 difieren entre canal
  constreñido y cross-slot (p. ej. Gs de PC-3: 10.4–20.7 frente a 5.8–14.4 mN/m). Esto importa
  directamente: nuestro trabajo compara tres familias de dispositivo, y los propios datos de calibración
  muestran que el módulo inferido **depende de la geometría de medida**.
- **Discrepancia entre técnicas, no promediada**: la razón de rigidez núcleo/citoplasma vale 1.4
  in situ y 3–4 en núcleos aislados; la tabla recoge 1.4–4 y el nominal se deja en el extremo in situ.
- **Rango explorado por los autores del modelo**: Gsn/Gs entre 2 y 8; λ entre 0.2 y 2.

---

## 6. Huecos abiertos (ningún valor utilizable encontrado)

1. **Radio nuclear y relación N/C de las cuatro líneas de mama en suspensión.** Ninguna de las
   consultas (PubMed y OpenAlex, incluidas búsquedas específicas de volumen/área/diámetro nuclear)
   devolvió un valor publicado. Es el hueco más caro del proyecto, porque aₙ/a es uno de los dos ejes
   del mapa de fases propuesto. Hoy aₙ/a entra como parámetro transferido de Wang et al. (0.5, rango
   0.4–0.6) y debe tratarse como variable de barrido, no como dato.
2. **Módulo de corte o de dilatación de área de membrana de cualquier línea de mama.** Cero fuentes.
   Todo viene de PC-3.
3. **Viscosidad de membrana (η_s) de cualquier línea de mama.** Cero fuentes.
4. **Tensión cortical de cualquier línea de mama.** Cero fuentes utilizables; el único valor de
   tensión cortical recuperado es de granulocito (0.035 mN/m).
5. **Cualquier propiedad de citoplasma o núcleo de SK-BR-3.** Cero fuentes (ver §4).
6. **Tiempo de relajación nuclear.** No se encontró ninguna medida directa; solo el endurecimiento por
   deformación (factor 1.5–2.5) y la razón de viscosidad nuclear/citoplasmática (~2) del condrocito.
7. **Incompatibilidad de magnitudes en la viscosidad citoplasmática.** La viscosidad «efectiva» del
   citoplasma en el modelo continuo de Wang et al. es 32.5 mPa·s (rango de literatura citado
   10–100 mPa·s), mientras que la microrreología intracelular en las líneas de mama da 10–60 Pa·s:
   **tres órdenes de magnitud de diferencia**. No son la misma magnitud física (viscosidad efectiva de
   gota continua frente a viscosidad local sondeada a escala de micras) y la tabla las mantiene en
   filas separadas con unidades distintas. Resolver qué valor alimenta el solver 2D es una decisión de
   modelado pendiente, no un problema de datos.

### Fuentes identificadas pero inaccesibles (texto completo cerrado)

- 10.1039/c7lc00599g — mechanophenotyping en microcanal constreñido, MCF-10A y MDA-MB-231; el resumen
  publica 182 ± 34.74, 260 ± 52, 360 ± 75 y 763 ± 93 Pa pero **no permite asignar cada valor a su línea
  celular**, por lo que no se incorporó ninguna fila. Candidato prioritario a petición de PDF.
- 10.1016/s0021-9290(99)00175-x — revisión de aspiración por micropipeta (Hochmuth).
- 10.1006/bbrc.2000.2360 — Guilak et al., núcleos aislados; solo se usó la razón cualitativa del
  resumen («3-4 times stiffer and nearly twice as viscous as the cytoplasm»).
- 10.1016/j.jbiomech.2006.09.023 — mecánica del núcleo en aspiración por micropipeta.
- 10.1242/jcs.01357 — el valor de 25 mN/m procede del resumen; el texto completo no es accesible.

### Advertencia de versión

El texto completo recuperado para 10.1039/d4na00003j (viscosidad citoplasmática por subtipo) es la
**preprint de bioRxiv** depositada por los autores, no la versión final de *Nanoscale Advances*. Los
valores de la Tabla 1 deben reverificarse contra la versión publicada antes de enviar el manuscrito.
Igualmente, 10.1038/s41598-026-56141-0 se recuperó como *accepted manuscript / article in press*.
