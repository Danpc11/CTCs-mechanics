#!/usr/bin/env bash
# Run detached from the terminal (survives closing SSH). In byobu/tmux, run the CLI directly instead.
# Usage: bash scripts/run_detached.sh configs/full.yaml [workers]
set -euo pipefail
CFG="${1:-configs/full.yaml}"; W="${2:-20}"
cd "$(dirname "$0")/.."
export PYTHONPATH="$PWD/src" PYTHONHASHSEED=0
OUT=$(python -c "import yaml;print(yaml.safe_load(open('$CFG'))['output'])")
mkdir -p "$OUT/logs"
python -m pytest -q
python -m ctcfem.cli plan --config "$CFG" --workers "$W"
nohup setsid python -m ctcfem.cli run --config "$CFG" --workers "$W" > "$OUT/logs/stdout.log" 2>&1 < /dev/null &
echo $! > "$OUT/logs/pid"
echo "PID $(cat "$OUT/logs/pid"). Follow: tail -f $OUT/logs/progress.log. Stop: kill -- -\$(cat $OUT/logs/pid)"
