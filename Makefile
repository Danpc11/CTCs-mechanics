# make env | make test | make smoke | make convergence | make full | make analysis | make figures | make app
PY := PYTHONPATH=src python
env:         ; conda env create -f environment.yml
test:        ; $(PY) -m pytest -q && node webapp/test_core.js
smoke:       ; $(PY) -m ctcfem.cli run --config configs/smoke_test.yaml && $(PY) -m ctcfem.cli analyze --config configs/smoke_test.yaml
convergence: ; $(PY) -m ctcfem.cli run --config configs/mesh_convergence_h10.yaml && $(PY) -m ctcfem.cli run --config configs/mesh_convergence_h14.yaml
full:        ; $(PY) -m ctcfem.cli run --config configs/full.yaml
analysis:    ; $(PY) -m ctcfem.cli analyze --config configs/full.yaml
figures:     ; cd manuscript_analysis/05_figuras && python generar_figuras.py
app:         ; python webapp/build_app.py
