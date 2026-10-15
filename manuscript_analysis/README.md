# Código del artículo — captura microfluídica de CTC de mama

Modelo multicompartimento 2D axisimétrico por elementos finitos y su aplicación a la comparación
de modalidades de captura. Todo el código que produjo los resultados del manuscrito.

## Entorno

```bash
conda create -n ctcfem -c conda-forge python=3.12 numpy scipy matplotlib pandas \
      scikit-fem meshio gmsh shapely scikit-learn sympy
conda activate ctcfem
pip install SALib          # solo para los análisis de Sobol
```
`pygmsh` NO está en conda-forge y rompe la resolución del entorno; no hace falta.

## Estructura

| carpeta | qué contiene |
|---|---|
| `01_solver/` | El solver. `solver_axisimetrico.py` es el núcleo: formulación mixta u-p con par Crouzeix-Raviart (P2+burbuja / P1DG, que admite salto de presión en la interfase núcleo/citoplasma), Newton-Raphson con continuación en la carga, contacto por penalización y tangente del material por paso complejo. `verificacion.py` reproduce las cuatro soluciones analíticas. `constriccion.py` es el barrido geométrico. |
| `02_parametros/` | `build_params.py` ensambla la tabla de parámetros por compartimento desde las fuentes de literatura, con DOI, técnica y bandera de transferencia. |
| `03_modalidades/` | La física de las tres modalidades. `fase2_lib.py` es la librería común; los `paso*.py` son las etapas en orden: calibración del córtex contra el módulo aparente medido, barrido FEM, análisis, eje tumoral/leucocito, DLD y diseño/Sobol. |
| `04_correcciones/` | Reparaciones posteriores a la revisión: `diseno_f6.py` rehace el barrido con Π por población y `sobol_f6.py` el análisis de sensibilidad con muestreo independiente por población. |
| `05_figuras/` | `figuras_bmmb.py` genera el juego final de figuras del manuscrito (sin títulos dentro de los paneles, con tramas además de color). |
| `06_hidraulica_chip/` | `hidraulica_membrana.py` — red hidráulica del chip de dos canales separados por membrana de policarbonato. Independiente del resto. |
| `datos/` | Todas las tablas de entrada y de resultados en CSV/JSON. |
| `informes/` | Los informes técnicos que documentan cada fase, incluida la auditoría que invalidó dos resultados y el triaje de la revisión. |

## Orden de ejecución

```
01_solver/verificacion.py            # punto de corte: 4 soluciones analíticas
02_parametros/build_params.py        # tabla de parámetros
03_modalidades/paso0_calibracion.py  # córtex contra E_app medido
03_modalidades/paso1_fem.py          # barrido FEM (reparte en lotes; ~2 h en 10 núcleos)
03_modalidades/paso123_analisis.py
03_modalidades/paso2b_tumoral_leucocito.py
03_modalidades/paso4_dld.py
03_modalidades/paso5_diseno.py  y  paso5_sobol.py
04_correcciones/diseno_f6.py  y  sobol_f6.py
05_figuras/figuras_bmmb.py
```
`06_hidraulica_chip/hidraulica_membrana.py` se ejecuta suelto y no depende de nada anterior.

## Límites del código, que son los del artículo

- **Envolvente de validez: 1,20 ≤ 2a/D ≤ 2,43** en confinamiento nominal (efectivo 1,04–1,64).
  Por debajo no hay contacto con la pared y el solver predice deformación idénticamente nula.
  No extrapolar fuera de ese rango.
- **No hay fluido resuelto.** La presión se impone; no sale de Navier–Stokes. Por eso DLD no se
  modela con este solver.
- **Convergencia ~74 %.** Los puntos que no convergen no se interpolan ni se rellenan: se declaran.
- **Las presiones críticas son cotas inferiores** (el tramo con 2a/D < 1,2 aporta cero) y sus
  valores absolutos dependen de la regularización de citoplasma Π, que es un parámetro declarado.
- **El problema deja de estar bien puesto con Π < 0,03**: córtex rígido sobre interior casi vacío
  admite configuraciones arrugadas de igual energía. El fallo es reproducible y monótono en Π.
- **53 de 86 parámetros son transferidos.** La bandera está en `datos/parametros_compartimento.csv`.
  El análisis de `04_correcciones/sobol_f6.py` cuantifica cuánto pesan: 32 % de la varianza de la
  selectividad.

## Advertencia sobre dos resultados corregidos

`informes/correcciones_fase6.md` lista las cifras que una revisión posterior invalidó. Las dos
importantes: el descarte de DLD descansaba en un radio de leucocito tomado de una línea de leucemia,
y la aparente insensibilidad del Sobol a los parámetros transferidos era cancelación algebraica, no
robustez. Los scripts de `04_correcciones/` son los que reparan ambas. **Usa esos resultados, no los
de `03_modalidades/` donde se solapen.**

## Correcciones de método (octubre 2026) — `07_correccion_reduccion/`

Revisión del método con re-ejecución del solver de producción. **Usar estos resultados para toda
presión absoluta.**

| defecto | dónde | efecto | corrección |
|---|---|---|---|
| W tabulado contra R_t **nominal** y evaluado en el radio **geométrico** R_w(z) | `paso123_analisis.splines` + `fase2_lib.presion_critica` | presiones críticas ×3.8–37 demasiado bajas; cambia la forma de ΔP(c) (P3: 23× → 4.2×) | `reduccion.spline_por_cuello`: W indexado por el radio real del cuello R_min |
| "cero exacto" a 2a/D = 1.20 | banda válida en radio nominal + `np.sign(0) = 0` | ΔP = 0 para cualquier W aunque haya contacto (conf. efectivo 1.039) | `z > 0` estricto y banda en confinamiento **efectivo** |
| máximo en el quiebre del perfil tanh(\|z\|) | `perfil_Rw` | la "convergencia en nz" era dz → 0 sobre una discontinuidad | perfil de cuello suave como análisis de sensibilidad (×1.6–2.1 más bajo, cocientes ±5 %) |
| Π_leucocito tratado como medido | `figuras_bmmb.py` (Fig. 4) y manuscrito | G_cyt del leucocito es el de MDA-MB-231 | reclasificado: medidos 0.712, transferidos 0.349 |
| leucocito de a = 8.5 µm (HL-60) en microfiltro y Fig. 5 | `paso123_analisis` (poro), `diseno_optimizado.csv` | compara la tumoral con un blasto, no con un leucocito primario | **pendiente de re-ejecutar** con radios primarios |
| poros etiquetados por radio nominal | `paso123_analisis` (poro) | apertura real 10.9–16.3 µm, no 7.5–14 | **pendiente** |
| `figuras_bmmb.py` no ejecutable | `panel_letter` sin definir | — | definida |
| falta el script de `diseno_optimizado.csv` | — | Fig. 5 no reproducible | **pendiente** (no está en el depósito) |

Ejecución:
```
cd 07_correccion_reduccion
python reejecutar_fem_mcf7.py        # opcional: regenera datos/fem_mcf7_referencia_reejecucion.jsonl (~4 min)
python corregir_tabla2.py            # pruebas T1–T4 + factores de corrección
python tabla2_directa.py             # Tabla 2 corregida desde FEM directo por línea
```
Para el barrido completo: en `paso123_analisis.py`, sustituir `splines()` por
`reduccion.spline_por_cuello(...)` y `L.presion_critica(...)` por `reduccion.presion_critica(...)`;
en la rama del microfiltro evaluar `spline_por_cuello` en R_p y usar radios primarios.
