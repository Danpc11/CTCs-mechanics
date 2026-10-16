# Manuscript analysis

Scripts and tables behind every number of the manuscript. Run each script from its own folder.

| folder | contents |
|---|---|
| `01_solver/` | Axisymmetric multi-compartment finite-element solver and its verification against exact solutions (`verificacion.py`). |
| `02_parametros/` | Parameter table with source, technique and transfer flag (`build_params.py` → `datos/parametros_compartimento.csv`). |
| `03_modalidades/` | Constriction sweep, compartment split, DLD and rate-dependent analysis. |
| `04_correcciones/` | Design sweep and global sensitivity analysis with independent per-population sampling. |
| `05_figuras/` | `generar_figuras.py` regenerates all figures into `../figures/` at 600 dpi. |
| `06_hidraulica_chip/` | Hydraulic network of the chip with a track-etched membrane. |
| `07_correccion_reduccion/` | Reduction of the energy landscape to critical pressure indexed by the neck radius; tests T1–T4 (`corregir_tabla2.py`) and the re-run finite-element data for Table 2. |
| `08_ctc_paciente/` | Size selection and contact census with patient-derived breast CTC diameters (Table S2). |
| `datos/` | All input and result tables. |
