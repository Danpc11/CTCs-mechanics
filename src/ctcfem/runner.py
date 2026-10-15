"""Reproducible pipeline: sampling -> tasks -> parallel execution -> raw JSON results.

* The cell sample is written ONCE (samples.csv) and the FEM reads it back, so a restart never resamples.
* Each task id is sha256(relevant configuration + sample row + device); a task whose JSON exists is skipped,
  which makes every run resumable.
* Atomic writes (tmp + rename). manifest.json records hashes of configuration, data and code, library
  versions and machine.
"""
from __future__ import annotations
import hashlib, json, os, platform, socket, sys, time, traceback
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import populations as P
from . import geometry as G
from . import reduction as RD
from . import execution as E
from .cell_model import close_aperture, CONF_EFF_ASSESSED

ROOT = Path(__file__).resolve().parents[2]


def _sha(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


def _sha_file(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_config(path):
    cfg = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    cfg["_path"] = str(Path(path).resolve())
    return cfg


def apertures(cfg):
    a = cfg["apertures_um"]
    v = np.round(np.arange(a["start"], a["stop"] + 1e-9, a["step"]), 4).tolist()
    return sorted(set(v) | set(float(x) for x in a.get("extra", [])))


def output_dir(cfg) -> Path:
    d = ROOT / cfg["output"]
    (d / "raw").mkdir(parents=True, exist_ok=True); (d / "logs").mkdir(parents=True, exist_ok=True)
    return d


def sample(cfg) -> pd.DataFrame:
    d = output_dir(cfg); f = d / "samples.csv"
    if not f.exists():
        data = P.load(ROOT / cfg["data"]["populations"])
        df = pd.concat([P.sample(data, nm, int(n), int(cfg["seed"])) for nm, n in cfg["populations"].items()], ignore_index=True)
        df.to_csv(f, index=False, float_format="%.17g")
        P.provenance_table(data).to_csv(d / "parameter_provenance.csv", index=False)
    return pd.read_csv(f)          # always from disk: identical task ids on first run and on restart


def tasks(cfg, samples: pd.DataFrame):
    mats = yaml.safe_load((ROOT / cfg["data"]["materials"]).read_text(encoding="utf-8"))
    radii = [0.5 * D for D in apertures(cfg)]
    rel = {k: cfg[k] for k in ("model", "apertures_um", "reduction")}
    out = []
    for nd in cfg["devices"]:
        dev = mats["devices"][nd]; mat = mats["materials"][dev["material"]]
        for _, row in samples.iterrows():
            c = row.to_dict()
            tid = _sha([rel, nd, dev, mat, {k: c[k] for k in sorted(c)}])[:20]
            out.append(dict(id=tid, device=nd, dev=dev, mat=mat, cell=c, radii=radii, model=cfg["model"], reduction=cfg["reduction"]))
    return out


def run_task(t, raw_dir):
    target = Path(raw_dir) / f"{t['id']}.json"
    if target.exists():
        return t["id"], "skipped"
    c, dev, mat, mo = t["cell"], t["dev"], t["mat"], t["model"]
    prof = G.profile_for_device(dev)
    rec = dict(id=t["id"], device=t["device"], material=dev["material"], type=dev["type"], population=c["population"],
               group=c["group"], sample=int(c["sample"]), cell=c, profile=dict(zt=prof.zt, L=prof.L, R_ext=prof.R_ext))
    try:
        E_kPa = float(mat["E_GPa"]["nominal"]) * 1e6 if mo.get("compliant_wall", True) else None
        l_w = dev["thickness_um"]["nominal"] if dev["type"] == "membrane" else None
        r = close_aperture(c, prof, t["radii"], conf_max_attempt=mo["conf_max_attempt"], hrel=mo["hrel"], kn_f=mo["kn_f"],
                           E_wall_kPa=E_kPa, nu_wall=float(mat["nu"]), l_wall=l_w, use_gamma=bool(mo.get("use_gamma", False)))
        rec["fem"] = r
        sp = RD.spline_W(r["curve"], c["a"]) if len(r["curve"]) > 1 else None
        rows = []
        for Rn in t["radii"]:
            a = c["a"]
            row = dict(D_aperture=2 * Rn, conf_eff=a / Rn, contact=bool(Rn < a), inside_assessed_range=bool(a / Rn <= CONF_EFF_ASSESSED))
            if Rn >= a:
                row.update(dP_local=0.0, dP_work=0.0, complete=True)
            elif sp is None:
                row.update(dP_local=np.nan, dP_work=np.nan, complete=False)
            else:
                pl, ok1 = RD.dP_local(prof, Rn, a, sp, r["rho_min_reached"])
                pw, ok2 = RD.dP_work(Rn, a, sp, r["rho_min_reached"], t["reduction"]["l_factor"])
                row.update(dP_local=pl, dP_work=pw, complete=bool(ok1 and ok2))
            row["dP_drop"] = RD.dP_drop(Rn, a, c["gamma"])
            row["L_load_um"] = prof.contact_half_length(a, Rn)
            row["entry_angle_deg"] = prof.entry_angle_deg(Rn)
            rows.append(row)
        rec["apertures"] = rows; rec["status"] = "ok"
        import resource
        rec["diag"] = dict(os_threads=E.process_threads(), peak_rss_MB=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    except Exception as exc:                       # recorded, never hidden
        rec["status"] = "error"; rec["error"] = repr(exc)[:500]; rec["trace"] = traceback.format_exc()[-3000:]
    tmp = target.with_suffix(".tmp"); tmp.write_text(json.dumps(rec, default=float), encoding="utf-8"); os.replace(tmp, target)
    return t["id"], rec["status"]


def _worker(args):
    return run_task(*args)


def _manifest(cfg, d, n_tasks):
    import scipy, skfem, gmsh, matplotlib
    m = dict(started=time.strftime("%Y-%m-%dT%H:%M:%S%z"), host=socket.gethostname(), platform=platform.platform(),
             cpu=os.cpu_count(), python=sys.version,
             versions=dict(numpy=np.__version__, scipy=scipy.__version__, scikit_fem=skfem.__version__, gmsh=gmsh.__version__,
                           pandas=pd.__version__, matplotlib=matplotlib.__version__, pyyaml=yaml.__version__),
             config=_sha_file(cfg["_path"]), data={k: _sha_file(ROOT / v) for k, v in cfg["data"].items()},
             code={p.name: _sha_file(p) for p in sorted((ROOT / "src" / "ctcfem").glob("*.py"))},
             n_tasks=n_tasks, seed=cfg["seed"])
    (d / "manifest.json").write_text(json.dumps(m, indent=1), encoding="utf-8")
    (d / "config_used.yaml").write_text(Path(cfg["_path"]).read_text(encoding="utf-8"), encoding="utf-8")


def run(cfg, workers=None, limit=None):
    d = output_dir(cfg)
    ts = tasks(cfg, sample(cfg))
    ts.sort(key=lambda t: -t["cell"]["a"])          # long tasks first: better load balance
    if limit is not None:
        ts = ts[:int(limit)]
    _manifest(cfg, d, len(ts))
    done = {p.stem for p in (d / "raw").glob("*.json")}
    pend = [t for t in ts if t["id"] not in done]
    ex = cfg["execution"]
    log = open(d / "logs" / "progress.log", "a", encoding="utf-8")

    def L(msg):
        log.write(f"{time.ctime()}  {msg}\n"); log.flush(); print(msg, flush=True)

    w = E.safe_workers(int(workers or ex["workers"]), float(ex.get("memory_per_task_gb", 1.0)), float(ex.get("memory_reserve_gb", 8.0)), log=L)
    L(f"tasks={len(ts)} pending={len(pend)} workers={w} threads_per_worker=1 "
      f"timeout={ex.get('task_timeout_s', 3600)}s memory_limit={ex.get('memory_limit_gb', 4.0)}GB")
    raw = str(d / "raw")

    def on_death(tid, arg, reason):
        t = arg[0]
        rec = dict(id=tid, device=t["device"], material=t["dev"]["material"], type=t["dev"]["type"], population=t["cell"]["population"],
                   group=t["cell"]["group"], sample=int(t["cell"]["sample"]), cell=t["cell"], status="error", error=reason)
        p = Path(raw) / f"{tid}.json"; tmp = p.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, default=float), encoding="utf-8"); os.replace(tmp, p)

    t0, n, err = time.time(), 0, 0
    items = ((t["id"], (t, raw)) for t in pend)
    gen = (((k, _worker(a)) for k, a in items) if w <= 1 else
           E.run_batch(items, _worker, w, timeout_s=float(ex.get("task_timeout_s", 3600)),
                       memory_limit_gb=float(ex.get("memory_limit_gb", 4.0)), nice=int(ex.get("nice", 0)),
                       min_free_gb=float(ex.get("min_free_memory_gb", 4.0)), tasks_per_child=int(ex.get("tasks_per_child", 25)),
                       on_death=on_death, log=L))
    for tid, res in gen:
        n += 1
        status = res[1] if isinstance(res, tuple) and res and res[0] == tid else (res[0] if isinstance(res, tuple) else res)
        err += status in ("error", "died")
        if n % 10 == 0 or n == len(pend):
            el = time.time() - t0
            L(f"{n}/{len(pend)}  errors={err}  {el/60:.1f} min  ETA {el/n*(len(pend)-n)/60:.1f} min  free memory {E.available_memory_gb():.0f} GB")
    L(f"finished  done={n} errors={err}")
    log.close()
    return n, err
