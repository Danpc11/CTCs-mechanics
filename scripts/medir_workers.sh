#!/usr/bin/env bash
# Mide el rendimiento con distinto numero de procesos sobre las MISMAS tareas (las mas pesadas).
# Uso: bash scripts/medir_workers.sh [n_tareas] [lista_de_workers]
#   ej.: bash scripts/medir_workers.sh 60 "20 30 40"
# El solver lineal (SuperLU) es de un hilo: el paralelismo es por procesos, 1 hilo BLAS cada uno.
set -euo pipefail
N="${1:-60}"; LISTA="${2:-20 30 40}"
cd "$(dirname "$0")/.."
export PYTHONPATH="$PWD/src" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
echo "CPU: $(nproc) hilos logicos; $(lscpu 2>/dev/null | awk -F: '/^Core\(s\) per socket|^Socket\(s\)|^Thread\(s\) per core/{gsub(/ /,"",$2); printf "%s=%s ", $1, $2}')"
printf "%8s %10s %14s\n" workers segundos tareas_por_hora
for W in $LISTA; do
  CFG="configs/_bench_${W}.yaml"
  sed "s#^salida: .*#salida: resultados/_bench_${W}#" configs/completo.yaml > "$CFG"
  rm -rf "resultados/_bench_${W}"
  T0=$(date +%s)
  python -m ctcfem.cli correr --config "$CFG" --workers "$W" --limite "$N" > /dev/null
  T=$(( $(date +%s) - T0 ))
  printf "%8s %10s %14s\n" "$W" "$T" "$(python3 -c "print(round($N*3600/max($T,1)))")"
  rm -f "$CFG"
done
echo "Usa el numero de workers con mas tareas_por_hora. Borra resultados/_bench_* al terminar."
