#!/usr/bin/env bash
# Throughput with different numbers of workers on the SAME (heaviest) tasks.
# Usage: bash scripts/benchmark_workers.sh [n_tasks] ["20 30 40"]
set -euo pipefail
N="${1:-60}"; LIST="${2:-20 30 40}"
cd "$(dirname "$0")/.."; export PYTHONPATH="$PWD/src"
echo "logical CPUs: $(nproc)"; printf "%8s %10s %14s\n" workers seconds tasks_per_hour
for W in $LIST; do
  CFG="configs/_bench_${W}.yaml"; sed "s#^output: .*#output: results/_bench_${W}#" configs/full.yaml > "$CFG"
  rm -rf "results/_bench_${W}"; T0=$(date +%s)
  python -m ctcfem.cli run --config "$CFG" --workers "$W" --limit "$N" > /dev/null
  T=$(( $(date +%s) - T0 )); printf "%8s %10s %14s\n" "$W" "$T" "$(python3 -c "print(round($N*3600/max($T,1)))")"
  rm -f "$CFG"
done
echo "Use the number of workers with the most tasks_per_hour; then delete results/_bench_*."
