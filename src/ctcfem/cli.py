"""Command line:  python -m ctcfem.cli <command> --config configs/<file>.yaml

commands:
  sample    write samples.csv and parameter_provenance.csv
  plan      count tasks and estimate cost
  run       run (resumes if interrupted)   [--workers N] [--limit N] [--retry-errors]
  analyze   raw results -> tables, pressure windows, report for 8/10/12 um and figures
"""
from .execution import force_single_thread
force_single_thread()        # forced (not setdefault): an inherited OMP_NUM_THREADS would oversubscribe the machine
import argparse, json, pathlib
from . import runner as R
from . import analysis as A


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ctcfem")
    ap.add_argument("command", choices=["sample", "plan", "run", "analyze"])
    ap.add_argument("--config", required=True)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--retry-errors", action="store_true", help="delete raw results with status error and recompute them")
    a = ap.parse_args(argv)
    cfg = R.load_config(a.config)
    if a.command == "sample":
        print(R.sample(cfg).groupby("population").size().to_string())
    elif a.command == "plan":
        ts = R.tasks(cfg, R.sample(cfg)); w = a.workers or cfg["execution"]["workers"]
        print(f"tasks: {len(ts)}")
        for s in (30, 120, 400):
            print(f"  if each task takes {s:4d} s -> {len(ts) * s / w / 3600:6.1f} h with {w} workers")
    elif a.command == "run":
        if a.retry_errors:
            raw = R.output_dir(cfg) / "raw"; nb = 0
            for f in raw.glob("*.json"):
                if json.loads(f.read_text(encoding="utf-8")).get("status") == "error":
                    f.unlink(); nb += 1
            print(f"{nb} failed tasks will be recomputed")
        n, err = R.run(cfg, workers=a.workers, limit=a.limit); print(f"done {n}, errors {err}")
    elif a.command == "analyze":
        T, _, _ = A.analyze(cfg, R.ROOT); print(f"rows {len(T)}; see {cfg['output']}/analysis/")


if __name__ == "__main__":
    main()
