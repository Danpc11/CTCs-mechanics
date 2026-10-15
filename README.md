# CTCs-mechanics

Code and data accompanying the manuscript

> **Cortical resistance and entry geometry govern the entry pressure of circulating tumour cells**
> Pérez-Calixto D, Vazquez-Victorio G, Escutia-Guadarrama L (submitted to the *Journal of the Mechanical Behavior of Biomedical Materials*).

A computational proof of concept, parameterised with published measurements of breast cancer cell
lines and leucocytes. It is not a validated device.

## Contents

| folder | what it is |
|---|---|
| `manuscript_analysis/` | The analysis behind every number of the manuscript: verified solver (`01_solver`), parameter table with provenance (`02_parametros`), constriction/DLD/rate analysis (`03_modalidades`), per-population Sobol (`04_correcciones`), figures (`05_figuras`), chip hydraulics (`06_hidraulica_chip`), **corrected reduction to critical pressure and re-run FEM data** (`07_correccion_reduccion`), size selection with **patient-derived CTC diameters** (`08_ctc_paciente`), all tables (`datos/`) and audit reports (`informes/`). |
| `src/ctcfem/` | Reproducible package for population and device runs (primary leucocyte subpopulations, membranes and channels, apertures 5–15 µm), with resumable parallel execution and a run manifest. |
| `data/` | Population distributions and device/material database, every entry with status (measured / transferred / declared) and source. |
| `webapp/` | Self-contained web app: entry pressure and a schematic of passage for any phenotype (`webapp/index.html`). |
| `figures/` | Manuscript figures regenerated from `manuscript_analysis/datos/` by `05_figuras/generar_figuras.py` (600 dpi). |
| `tests/` | Automated tests (solver patch test, exact-aperture geometry, sampling reproducibility, reduction, fault-tolerant parallel execution). |

## Quick start

```bash
conda env create -f environment.yml && conda activate ctcfem   # pinned versions
make pruebas                                                   # 14 tests
cd manuscript_analysis/07_correccion_reduccion && python corregir_tabla2.py   # corrected Table 2
node webapp/prueba.js                                          # web app reproduces Python
```

Full population run on a workstation (e.g. 20 physical cores):
`python -m ctcfem.cli correr --config configs/completo.yaml --workers 20` (in byobu/tmux) or
`bash scripts/correr_servidor.sh configs/completo.yaml 20` (detached). One BLAS thread per worker is enforced;
a dead or runaway worker is recovered without hanging the run (see `docs/USAGE_ctcfem_es.md`).

All manuscript figures (600 dpi): `cd manuscript_analysis/05_figuras && python generar_figuras.py` → `figures/`.

## Notes for reuse

* The Skalak parameter follows `Ks = (1 + 2C) Gs`, as implemented in `solver_axisimetrico.py`.
* Absolute critical pressures depend on the neck shape (corner vs smooth neck, 1.6–2.1×); ratios do not (< 5 %).
* Results above an effective confinement of 1.636 are extrapolation.
* Mechanical parameters are those of cultured cell lines. Patient-derived breast CTCs are smaller (median 12.4 µm on CellSearch, Mendelaar et al. 2021); `data/poblaciones.yaml` includes them as `CTC_mama_paciente`, with a transferred, unmeasured apparent modulus.

## Citation

See `CITATION.cff`. A licence has not yet been chosen; until one is added, all rights are reserved by the authors.
