# Changelog

## v0.1.1 — 2026-10-15
* **Patient-derived breast CTCs**: population `CTC_mama_paciente` in `data/poblaciones.yaml` (median from
  Mendelaar et al. 2021, nucleus-to-cell ratio from Sebastian et al. 2021, apparent modulus transferred and flagged
  as not measured) and `manuscript_analysis/08_ctc_paciente/` (size selection and contact census; Table S2).
* **Fault-tolerant parallel execution** (`src/ctcfem/ejecucion.py`): one BLAS thread forced inside each worker,
  recovery from dead workers without hanging, per-task time and memory limits, memory-aware throttling,
  abort on systematic failure, `--reintentar-errores`. Default 20 workers; `scripts/medir_workers.sh`.
* **All figures reproducible**: `manuscript_analysis/05_figuras/generar_figuras.py` → `figures/` at 600 dpi.
* Skalak parameter documented consistently as `Ks = (1 + 2C) Gs` in code comments and the parameter table.
* Web app: patient CTC phenotype.

## v0.1.0 — 2026-10-15
* First release: verified solver, corrected reduction to critical pressure, sensitivity analysis, `ctcfem`
  population/device pipeline, web app.
