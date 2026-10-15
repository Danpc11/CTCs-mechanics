#!/usr/bin/env bash
# Corrida en el servidor del INMEGEN (local, 40 nucleos, 250 GB).
# Uso:  bash scripts/correr_servidor.sh configs/completo.yaml [workers]
# Sobrevive al cierre de la sesion SSH (nohup + setsid). Reanuda si se vuelve a lanzar.
set -euo pipefail
CFG="${1:-configs/completo.yaml}"
W="${2:-40}"
cd "$(dirname "$0")/.."
export PYTHONPATH="$PWD/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0
SAL=$(python -c "import yaml,sys;print(yaml.safe_load(open('$CFG'))['salida'])")
mkdir -p "$SAL/logs"
echo "== pruebas =="          ; python -m pytest -q
echo "== plan =="             ; python -m ctcfem.cli planificar --config "$CFG" --workers "$W"
echo "== lanzando con $W procesos; log: $SAL/logs/progreso.log =="
nohup setsid python -m ctcfem.cli correr --config "$CFG" --workers "$W" \
      > "$SAL/logs/stdout.log" 2>&1 < /dev/null &
echo $! > "$SAL/logs/pid"
echo "PID $(cat "$SAL/logs/pid").  Seguir:  tail -f $SAL/logs/progreso.log"
echo "Al terminar:  python -m ctcfem.cli analizar --config $CFG"
