# ctcfem — célula multicompartimento en constricciones y membranas

Modelo axisimétrico de elementos finitos (córtex Skalak, citoplasma y nucleoplasma neo-Hookeanos casi
incompresibles, envoltura nuclear) de una célula **muestreada de una población** dentro de un
dispositivo de captura (membrana de PC, PET o PDMS; canal de PDMS o PMMA), con aperturas de **5 a
15 µm** y atención especial a **8, 10 y 12 µm**.

## Qué cambia respecto al código del manuscrito

| antes | ahora |
|---|---|
| radio de garganta *nominal*: la apertura real era mayor (g₀ = 0.168) | perfil suave con **apertura exacta**: el diámetro declarado es el que ve la célula |
| quiebre de tanh(\|z\|) en el cuello | perfil par analítico, derivada nula en el cuello |
| W indexado por R_t nominal y evaluado en R_w(z) | W indexado por el radio real del cuello (exacto en z = 0) |
| una corrida FEM por apertura | **una** continuación por célula da W(ρ) para todas las aperturas |
| leucocito = HL-60 (a = 8.5 µm) | neutrófilo, linfocito, monocito, eosinófilo y basófilo **primarios**, muestreados |
| pared rígida | rigidez del material (PC, PET, PDMS, PMMA) en serie con la penalización |
| parámetros fijos | muestreo por hipercubo latino con semilla y procedencia de cada parámetro |

## Instalación (una vez)

```bash
conda env create -f environment.yml      # versiones exactas verificadas
conda activate ctcfem
make pruebas                             # 10 pruebas; deben pasar todas
make humo                                # ~1 min en 1 núcleo: cadena completa
```

## Corrida en el servidor (40 núcleos)

```bash
bash scripts/correr_servidor.sh configs/convergencia_malla_h10.yaml 20   # 1) convergencia de malla
bash scripts/correr_servidor.sh configs/convergencia_malla_h14.yaml 20   #    (pareadas: misma muestra)
bash scripts/correr_servidor.sh configs/completo.yaml 40                 # 2) producción
tail -f resultados/completo/logs/progreso.log                            #    progreso y ETA
python -m ctcfem.cli analizar --config configs/completo.yaml             # 3) análisis
```

* Corre con `nohup setsid`: sobrevive al cierre de la sesión SSH.
* **Reanuda**: si se corta, se vuelve a lanzar el mismo comando y salta lo ya hecho.
* Un hilo BLAS por proceso (40 procesos = 40 núcleos). Memoria por proceso < 1 GB.
* `python -m ctcfem.cli planificar --config configs/completo.yaml` da el número de tareas y el coste.

## Reproducibilidad

* `muestras.csv` se escribe una vez y el FEM la lee; los ids de tarea son hashes de entradas.
* `manifest.json` registra hashes de configuración, datos y cada archivo de código, versiones de
  librerías, máquina y semilla. Dos corridas son la misma corrida si y solo si coinciden sus hashes.
* El solver (`src/ctcfem/solver_axisimetrico.py`) es el verificado en el Apéndice A, sin cambios
  (sha256 `c5261e18…53dfe`).

## Salidas (`resultados/<nombre>/analisis/`)

`tabla_larga.csv` · `resumen_aperturas.csv` · `ventanas.csv` · `informe_8_10_12.md` · `fig_<dispositivo>.png`

## Límites que hay que leer antes de interpretar

1. **Envolvente validada: confinamiento efectivo a/R ≤ 1.636.** Una MCF-7 (17.5 µm) a 12 µm está en
   1.46 (dentro); a 10 µm en 1.75 (fuera); a 8 µm en 2.19 (muy fuera). Fuera de la envolvente el
   resultado es extrapolación y se marca (`dentro_envolvente = False`). Si la continuación no llega
   a la apertura, el valor es **cota inferior** (`completo = False`).
2. A 8 µm una célula tumoral tendría que alargarse a ~55 µm para pasar, más que el espesor de la
   membrana: el modelo simétrico no representa esa extrusión. Lo que sí resuelve bien a 8 µm es el
   lado leucocitario.
3. Sin fluido resuelto, sin adhesión ni fricción (relevante para PDMS, hidrófobo).
4. En membranas con borde agudo, `dP_local` escala como 1/radio de borde; el estimador recomendado
   es `dP_trabajo` (declara `l_factor`).
5. La procedencia de cada parámetro (medido / transferido / declarado) está en
   `procedencia_parametros.csv`. La tensión cortical depende mucho de la temperatura (ver
   `data/poblaciones.yaml: notas`).

Detalle del método y fuentes: `docs/METODO.md`.

## Robustez de la ejecución paralela (`src/ctcfem/ejecucion.py`)

Soluciones trasladadas de T2D_landscape (`[FIX-threads]`, `[FIX-worker]`):

| problema | solución | prueba |
|---|---|---|
| un `OMP_NUM_THREADS` heredado (p. ej. de `.bashrc`) hacía que N workers × M hilos sobresuscribieran la máquina | 1 hilo **forzado** en la CLI y **dentro de cada worker** con `threadpoolctl` | `test_un_hilo_por_worker_aunque_el_entorno_herede_otro_valor` |
| si el sistema mata un worker (OOM killer, fallo nativo) un `Pool` se cuelga para siempre | `ProcessPoolExecutor`: se detecta, se recrea el pool y las tareas en vuelo se reintentan **aisladas**; solo la culpable se registra como error | `test_worker_muerto_y_tarea_colgada_no_cuelgan_el_lote` |
| una tarea desbocada | límite de tiempo (`timeout_tarea_s`) y de memoria por worker (`limite_memoria_gb`, RLIMIT_DATA) | ídem |
| fallo del entorno que mata a todos los workers | se **aborta** sin marcar tareas (`FalloSistematico`) | `test_fallo_sistematico_aborta_sin_marcar_tareas` |
| poca memoria libre | se pausa el envío por debajo de `memoria_minima_libre_gb`; workers limitados por CPU y memoria | — |

Cada crudo guarda `diag.hilos_os` y `diag.rss_pico_MB`. Para comprobar en el servidor que no hubo sobresuscripción:

```bash
python -c "import json,glob;d=[json.load(open(f)).get('diag',{}) for f in glob.glob('resultados/completo/crudos/*.json')];print('hilos',sorted({x.get('hilos_os') for x in d}),'RSS max MB',max(x.get('rss_pico_MB',0) for x in d))"
```

Reintentar las tareas que fallaron: `python -m ctcfem.cli correr --config configs/completo.yaml --reintentar-errores`.
