"""Aggregation of raw results into tables and figures.

Recommended estimator by device type (see reduction.py): channel -> dP_local, membrane -> dP_work.

Outputs in <output>/analysis/:
  long_table.csv        one row per (cell, device, aperture)
  aperture_summary.csv  quantiles by population, device and aperture
  windows.csv           leucocyte/tumour pressure window by device and aperture
  report_8_10_12.md     direct reading of the apertures of interest
  fig_<device>.png      critical pressure against aperture (median and 5th-95th percentile) by population
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

ESTIMATOR = {"channel": "dP_local", "membrane": "dP_work"}


def long_table(d: Path) -> pd.DataFrame:
    rows = []
    for f in sorted((d / "raw").glob("*.json")):
        r = json.loads(f.read_text(encoding="utf-8"))
        base = dict(id=r["id"], device=r["device"], material=r["material"], type=r["type"], population=r["population"],
                    group=r["group"], sample=r["sample"], status=r["status"], D_cell=r["cell"]["D"],
                    E_app=r["cell"]["E_app"], tau=r["cell"]["tau"])
        if r["status"] != "ok":
            rows.append(base); continue
        fem = r["fem"]
        base.update(stop_reason=fem.get("stop_reason"), rho_min=fem.get("rho_min_reached"), seconds=fem.get("seconds", 0.0))
        for a in r["apertures"]:
            row = dict(base, **a)
            row["dP"] = a[ESTIMATOR[r["type"]]]
            row["beta_c"] = row["dP"] * base["tau"] / a["L_load_um"] * 1e3 if a["L_load_um"] > 0 else 0.0   # kPa per mm/s
            rows.append(row)
    return pd.DataFrame(rows)


def _q(x, p):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.quantile(x, p)) if len(x) else np.nan


def summary(T: pd.DataFrame) -> pd.DataFrame:
    g = T[T.status == "ok"].groupby(["device", "population", "group", "D_aperture"])
    return g.apply(lambda s: pd.Series(dict(
        n=len(s), frac_contact=s.contact.mean(), frac_complete=s.complete.mean(),
        frac_inside_range=s[s.contact].inside_assessed_range.mean() if s.contact.any() else 1.0,
        dP_p5=_q(s.dP, .05), dP_p50=_q(s.dP, .5), dP_p95=_q(s.dP, .95), dP_drop_p50=_q(s.dP_drop, .5))),
        include_groups=False).reset_index()


def windows(T: pd.DataFrame, differential: dict) -> pd.DataFrame:
    """Window [P99 of the differential-weighted leucocyte mixture, P1 of each tumour population]."""
    out = []
    for (nd, D), s in T[T.status == "ok"].groupby(["device", "D_aperture"]):
        leu = s[s.group == "leucocyte"]
        if leu.empty:
            continue
        w = leu.population.map(differential).fillna(0.0) / leu.groupby("population").population.transform("count")
        x, ww = leu.dP.to_numpy(float), w.to_numpy(float)
        o = np.argsort(x); cw = np.cumsum(ww[o]) / ww.sum()
        p99_leu = float(x[o][np.searchsorted(cw, 0.99)])
        for tum, st in s[s.group == "tumour"].groupby("population"):
            p1 = _q(st.dP, 0.01)
            P_op = np.sqrt(max(p99_leu, 1e-6) * max(p1, 1e-6)) if p1 > p99_leu else np.nan
            out.append(dict(device=nd, D_aperture=D, tumour=tum, P99_leucocyte=p99_leu, P1_tumour=p1,
                            window_exists=bool(p1 > p99_leu), window_factor=p1 / p99_leu if p99_leu > 0 else np.inf,
                            P_operating=P_op,
                            tumour_yield=float((st.dP > P_op).mean()) if np.isfinite(P_op) else np.nan,
                            leucocyte_carryover_per_million=float(ww[x > P_op].sum() / ww.sum() * 1e6) if np.isfinite(P_op) else np.nan,
                            frac_tumour_complete=float(st.complete.mean()), frac_tumour_inside_range=float(st.inside_assessed_range.mean())))
    return pd.DataFrame(out)


def report(R: pd.DataFrame, V: pd.DataFrame, focus=(8.0, 10.0, 12.0)) -> str:
    L = ["# Apertures of interest", "",
         "Estimator: dP_local (channels) and dP_work (membranes). 'complete' = the continuation reached that aperture;",
         "otherwise the value is a LOWER BOUND. 'inside range' = effective confinement <= 1.636 (assessed range);",
         "outside it the result is an extrapolation.", ""]
    for D in focus:
        L += [f"## Aperture {D:g} µm", ""]
        r = R[np.isclose(R.D_aperture, D)]
        if r.empty:
            L += ["(no data)", ""]; continue
        L += [r[["device", "population", "frac_contact", "frac_complete", "frac_inside_range", "dP_p5", "dP_p50", "dP_p95"]]
              .round(4).to_markdown(index=False), ""]
        v = V[np.isclose(V.D_aperture, D)] if not V.empty else V
        if not v.empty:
            L += ["Pressure windows:", "", v.drop(columns="D_aperture").round(4).to_markdown(index=False), ""]
    return "\n".join(L)


def figures(R: pd.DataFrame, d: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for nd, s in R.groupby("device"):
        fig, ax = plt.subplots(figsize=(6.85, 3.6))
        for pop, q in s.groupby("population"):
            q = q.sort_values("D_aperture")
            ls = {"tumour": "-", "leucocyte": "--"}.get(q.group.iloc[0], ":")
            ax.plot(q.D_aperture, q.dP_p50.clip(lower=1e-5), ls, label=pop)
            ax.fill_between(q.D_aperture, q.dP_p5.clip(lower=1e-5), q.dP_p95.clip(lower=1e-5), alpha=0.15)
        for D in (8, 10, 12):
            ax.axvline(D, color="0.6", lw=0.6)
        ax.set_yscale("log"); ax.set_xlabel("aperture diameter (µm)"); ax.set_ylabel("critical entry pressure (kPa)")
        ax.legend(fontsize=6, ncol=2, frameon=False); fig.tight_layout()
        fig.savefig(d / f"fig_{nd}.png", dpi=300); plt.close(fig)


def analyze(cfg, root: Path):
    d = root / cfg["output"]; out = d / "analysis"; out.mkdir(exist_ok=True)
    T = long_table(d); T.to_csv(out / "long_table.csv", index=False)
    if T.empty or "D_aperture" not in T:
        return T, None, None
    R = summary(T); R.to_csv(out / "aperture_summary.csv", index=False)
    diff = yaml.safe_load((root / cfg["data"]["populations"]).read_text(encoding="utf-8"))["differential"]
    V = windows(T, diff); V.to_csv(out / "windows.csv", index=False)
    (out / "report_8_10_12.md").write_text(report(R, V), encoding="utf-8")
    figures(R, out)
    return T, R, V
