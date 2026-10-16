# CTCs-mechanics

Code and data for the manuscript

> **Cortical resistance and entry geometry govern the entry pressure of circulating tumour cells**
> Pérez-Calixto D, Vazquez-Victorio G, Escutia-Guadarrama L. Submitted to the *Journal of the Mechanical Behavior of Biomedical Materials*.

A computational proof of concept parameterised with published measurements of breast cancer cell lines and leucocytes.

**Web app:** <https://danpc11.github.io/CTCs-mechanics/webapp/> — computed cell shape and critical entry pressure for any phenotype, device and aperture.

## Contents

| folder | contents |
|---|---|
| `manuscript_analysis/` | Analysis behind the manuscript: solver and its verification (`01_solver`), parameter table with provenance (`02_parametros`), constriction, DLD and rate analysis (`03_modalidades`), per-population sensitivity analysis (`04_correcciones`), figures (`05_figuras`), chip hydraulics (`06_hidraulica_chip`), reduction to critical pressure and re-run finite-element data (`07_correccion_reduccion`), patient-derived CTC size analysis (`08_ctc_paciente`), operating windows as flow rates in µL/min (`09_flow_units`) and all tables (`datos/`). |
| `src/ctcfem/` | Package for population and device runs (primary leucocyte subpopulations, membranes and channels, apertures of 5–15 µm) with resumable, fault-tolerant parallel execution and a run manifest. |
| `data/` | Population distributions and device/material data; every entry has a status (measured, transferred, declared, derived) and a source. |
| `configs/` | Run configurations (smoke test, mesh convergence, full run). |
| `webapp/` | Self-contained web app and the scripts that build it. |
| `figures/` | Manuscript figures (Fig. 1–7, Fig. S1) and graphical abstract at 600 dpi, regenerated from `manuscript_analysis/datos/` by `05_figuras/generar_figuras.py`. |
| `tests/` | Automated tests: solver patch test, geometry, sampling, reduction, volume conservation, fault-tolerant execution. |

## Quick start

```bash
conda env create -f environment.yml && conda activate ctcfem
make test          # Python tests and the web-app core
make smoke         # end-to-end run on one core (~1 min)
make figures       # manuscript figures -> figures/
```

Full run (e.g. 20 cores, in byobu or tmux): `python -m ctcfem.cli run --config configs/full.yaml --workers 20`,
then `python -m ctcfem.cli analyze --config configs/full.yaml`. See `docs/USAGE.md`.

## Re-running the manuscript analysis

`bash scripts/rerun_manuscript.sh 20 1e5` re-runs every finite-element sweep and analysis of `manuscript_analysis/` with the
volume-conserving bulk modulus (kappa/mu = 1e5) in a separate working copy (`rerun_k1e5/`), regenerates the figures and writes
`rerun_k1e5.tar.gz` with a summary (`SUMMARY.txt`), the tables, the figures and the logs. The deposited results are not modified.

## Model notes

* Skalak cortex with `Ks = (1 + 2C) Gs`; near-incompressible neo-Hookean cytoplasm and nucleoplasm (bulk/shear ratio 1e5).
* Results above an effective confinement of 1.636 are extrapolation.
* Absolute critical pressures depend on the neck shape (1.6–2.1×); tumour/leucocyte ratios do not (< 5 %).
* Mechanical parameters are those of cultured cell lines. Patient-derived breast CTCs are smaller (Mendelaar et al. 2021) and their apparent modulus has not been measured.

## Citation and licence

See `CITATION.cff`. A licence has not yet been chosen; until one is added, all rights are reserved by the authors.
