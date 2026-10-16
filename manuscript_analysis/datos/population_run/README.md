# Population run (configs/full.yaml)

720 cells (11 populations) x 5 device families, apertures 5-15 µm, kappa/mu = 1e5: 3600 finite-element continuations,
0 errors (`progress.log`). `manifest.json` records the hashes of the code, inputs and configuration, library versions
and machine. `long_table.csv.gz` has one row per cell, device and aperture; `aperture_summary.csv` and `windows.csv`
(operating points, Table S4) are produced from it by `python -m ctcfem.cli analyze --config configs/full.yaml`.
The 3600 raw JSON results are archived with the Zenodo release.
