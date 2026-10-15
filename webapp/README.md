# Web app: paso de una célula por un poro

`index.html` es una página autocontenida (sin servidor, sin dependencias) que calcula al instante la
presión crítica de entrada de un fenotipo celular en un poro o constricción de 5–15 µm y anima su paso.

## Cómo calcula
* La energía de confinamiento colapsa como `W = Gs a² g(ρ/a)` (verificado en el artículo entre líneas, 1–6 %).
* `g` está tabulada en `curvas_fem.json` desde el FEM de `ctcfem` (MCF-7, hasta a/R = 2.4) para cada dispositivo.
* Para otro fenotipo, `Gs = E_app · a / 2.608` (calibración del córtex).
* Canal: reducción axial indexada por el cuello. Membrana: trabajo sobre el volumen que entra en el poro.
* Cota independiente del FEM: gota líquida `2γ(1/R − 1/a)`.
* La animación es un esquema de conservación de volumen, **no** la forma calculada por el FEM.

## Archivos
| archivo | función |
|---|---|
| `plantilla.html` | interfaz y dibujo |
| `nucleo.js` | cálculo (mismo algoritmo que `src/ctcfem/reduccion.py`) |
| `curvas_fem.json` | curvas `g(ρ/a)` del FEM |
| `exportar_web.py` | genera `index.html` desde `data/poblaciones.yaml` + curvas |
| `generar_curvas.py` | regenera `curvas_fem.json` con el FEM (~40 s) |
| `prueba.js` | `node webapp/prueba.js`: el JS reproduce a Python (< 0.1 % dentro de los datos) |

## Añadir un tipo celular
1. Añadir la población en `data/poblaciones.yaml` (D, an_a, E_app, gamma, eta, tau), cada parámetro con `estado` y `fuente`.
2. `python webapp/exportar_web.py`.

## Límites
* El colapso se verificó con la misma razón núcleo/célula; para núcleos mayores (linfocitos, blastos) es aproximación declarada.
* a/R > 1.636: extrapolación; a/R > 2.4: no se calcula.
* Solo dos dispositivos tabulados (membrana PC, canal PDMS). Para PET, PDMS-membrana o PMMA, añadirlos en `generar_curvas.py`.
* Sin fluido resuelto, adhesión ni fricción.
