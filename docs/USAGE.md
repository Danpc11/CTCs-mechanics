# Usage

## Commands
```bash
python -m ctcfem.cli sample  --config configs/full.yaml     # samples.csv, parameter_provenance.csv
python -m ctcfem.cli plan    --config configs/full.yaml     # number of tasks and cost estimate
python -m ctcfem.cli run     --config configs/full.yaml --workers 20 [--limit N] [--retry-errors]
python -m ctcfem.cli analyze --config configs/full.yaml     # results/full/analysis/
```
`run` resumes where it stopped: finished tasks are skipped. To run detached from the terminal instead of in
byobu/tmux: `bash scripts/run_detached.sh configs/full.yaml 20`. Follow progress with
`tail -f results/full/logs/progress.log`.

## Choosing the number of workers
The linear solver (SuperLU) is single-threaded, so parallelism is by processes with one thread each.
`bash scripts/benchmark_workers.sh 60 "20 30 40"` measures throughput on the heaviest tasks.

## Robust execution
* One BLAS/OpenMP thread is enforced in the command line and inside every worker.
* A worker killed by the system does not hang the run: the pool is rebuilt and only the task that kills its
  worker is recorded as failed. Per-task limits: `task_timeout_s`, `memory_limit_gb`.
* If workers die without any task finishing, the run stops with a message instead of marking tasks.
* Every raw result stores `diag.os_threads` and `diag.peak_rss_MB`.

## Outputs (`results/<name>/`)
`samples.csv`, `parameter_provenance.csv`, `manifest.json` (hashes of code, inputs and configuration),
`raw/*.json`, and in `analysis/`: `long_table.csv`, `aperture_summary.csv`, `windows.csv`,
`report_8_10_12.md`, `fig_<device>.png`.

## Reproducibility
The sample is written once and read back, task ids are hashes of their inputs, writes are atomic, and the
manifest records versions and hashes. The number of workers changes speed, not results.
