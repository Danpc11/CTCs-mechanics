#!/usr/bin/env bash
# Re-run every number of the manuscript_analysis/ pipeline with a given bulk-to-shear ratio, in a separate working
# copy (rerun_k<kappa>/) so that the deposited results are never overwritten.
# Usage (in byobu/tmux): bash scripts/rerun_manuscript.sh [workers=20] [kappa=1e5]
set -euo pipefail
W="${1:-20}"; K="${2:-1e5}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; WORK="$ROOT/rerun_k${K}"
export CTC_KAPPA_MU="$K" CTC_WORKERS="$W" PYTHONHASHSEED=0 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
echo "== working copy $WORK  (kappa/mu = $K, $W workers)"
rm -rf "$WORK"; mkdir -p "$WORK/webapp" "$WORK/logs"
cp -a "$ROOT/manuscript_analysis" "$WORK/"; ln -s "$ROOT/src" "$WORK/src"; ln -s "$ROOT/data" "$WORK/data"
cp "$ROOT/webapp/energy_curves.json" "$WORK/webapp/"
M="$WORK/manuscript_analysis"; L="$WORK/logs"; T0=$(date +%s)
step() { echo "== [$(( ($(date +%s) - T0) / 60 )) min] $*"; }
shards() {  # run a sharded script with W independent processes
  local script="$1" tag="$2"; seq 0 $((W - 1)) | xargs -P "$W" -I{} sh -c "python $script --shard {} --de $W > $L/${tag}_{}.log 2>&1"; }

step "1/6 constriction sweep (paso1_fem, 300 FEM points)"; cd "$M/03_modalidades"; rm -f paso1_fem_s*.json; shards paso1_fem.py paso1
step "2/6 analysis of the sweep (paso123) and rate/purity (paso2b)"
python paso123_analisis.py > "$L/paso123.log" 2>&1; python paso2b_tumoral_leucocito.py > "$L/paso2b.log" 2>&1
step "3/6 design hypercube (paso5, 640 FEM points) and Sobol"; rm -f paso5b_diseno_s*.json; shards paso5_diseno.py paso5
python paso5_sobol.py > "$L/paso5_sobol.log" 2>&1
cd "$M/04_correcciones"; python diseno_f6.py > "$L/diseno_f6.log" 2>&1; python sobol_f6.py > "$L/sobol_f6.log" 2>&1
step "4/6 Table 2 and entry geometries (64 FEM points)"; cd "$M/07_correccion_reduccion"
python rerun_fem_points.py > "$L/fem_points.log" 2>&1; python tabla2_directa.py > "$L/tabla2.log" 2>&1; python geometrias_170.py > "$L/geometrias.log" 2>&1
step "5/6 copy results into datos/ and regenerate figures"
for f in "$M"/03_modalidades/*.csv "$M"/03_modalidades/paso123.json "$M"/03_modalidades/paso2b.json "$M"/04_correcciones/*.csv "$M"/07_correccion_reduccion/*.csv; do
  [ -f "$f" ] && cp "$f" "$M/datos/"; done
cd "$M/05_figuras"; python compute_fem_shapes.py > "$L/shapes.log" 2>&1; CTC_DIRECT=1 python generar_figuras.py > "$L/figures.log" 2>&1
step "6/6 summary"; cd "$M"; python summarize_rerun.py | tee "$WORK/SUMMARY.txt"
tar czf "$ROOT/rerun_k${K}.tar.gz" -C "$ROOT" "rerun_k${K}/SUMMARY.txt" "rerun_k${K}/logs" "rerun_k${K}/figures" \
    $(cd "$ROOT" && ls -d rerun_k${K}/manuscript_analysis/datos/*.csv rerun_k${K}/manuscript_analysis/datos/*.json rerun_k${K}/manuscript_analysis/datos/*.jsonl 2>/dev/null)
step "done -> $ROOT/rerun_k${K}.tar.gz"
