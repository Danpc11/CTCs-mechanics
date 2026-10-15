# make entorno | make pruebas | make humo | make convergencia | make completo | make analisis
PY := PYTHONPATH=src python
entorno:      ; conda env create -f environment.yml
pruebas:      ; $(PY) -m pytest -q
humo:         ; $(PY) -m ctcfem.cli correr --config configs/prueba_rapida.yaml && $(PY) -m ctcfem.cli analizar --config configs/prueba_rapida.yaml
convergencia: ; bash scripts/correr_servidor.sh configs/convergencia_malla_h10.yaml 20 && bash scripts/correr_servidor.sh configs/convergencia_malla_h14.yaml 20
completo:     ; bash scripts/correr_servidor.sh configs/completo.yaml 40
analisis:     ; $(PY) -m ctcfem.cli analizar --config configs/completo.yaml
