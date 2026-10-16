# Web app

`index.html` is self-contained (no server). It shows the computed finite-element shape of a cell held at the
neck of a pore or constriction (5–15 µm) and its critical entry pressure.

| file | role |
|---|---|
| `template.html` | interface and drawing |
| `core.js` | critical-pressure calculation (same algorithm as `src/ctcfem/reduction.py`) |
| `energy_curves.json` | normalised confinement energy g(ρ/a) from the finite-element model (`build_energy_curves.py`) |
| `fem_shapes.json` | deformed meshes and element strain from the finite-element model (`build_shapes.py`) |
| `build_app.py` | builds `index.html` from `data/populations.yaml` and the two JSON files |
| `test_core.js` | `node webapp/test_core.js`: the JavaScript core reproduces the Python reduction |

To add a cell type: add it to `data/populations.yaml` (each parameter with status and source) and run
`python webapp/build_app.py`.
