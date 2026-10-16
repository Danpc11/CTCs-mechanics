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
    """Contact-based operating points.

    Size decides which cells touch the wall; mechanics acts only on those that do. For each device, aperture and
    tumour population, the operating pressure P_op is set to retain a target fraction of ALL tumour cells (cells
    that do not touch the wall pass whatever the pressure, so the attainable yield is capped by the contact
    fraction). Quantiles use only complete states: continuations that stopped before the aperture give lower
    bounds and are excluded and counted. Leucocyte retention is reported as the observed differential-weighted
    fraction and the observed number of retained sampled leucocytes; with a stratified, weighted sample no
    binomial bound is attached to a zero count.
    """
    out = []
    for (nd, D), s in T[T.status == "ok"].groupby(["device", "D_aperture"]):
        leu = s[s.group == "leucocyte"].copy()
        if leu.empty:
            continue
        leu["w"] = leu.population.map(differential).fillna(0.0) / leu.groupby("population").population.transform("count")
        for tum, st in s[s.group == "tumour"].groupby("population"):
            cap = float(st.contact.mean())
            usable = st[(~st.contact) | st.complete]
            for target in (0.25, 0.50, 0.70, 0.90, 0.99):
                if target > cap or usable.empty:
                    continue
                P_op = float(np.quantile(usable.dP, 1.0 - target))
                if P_op <= 0:
                    continue
                kept = usable[usable.dP > P_op]
                out.append(dict(device=nd, D_aperture=D, tumour=tum, tumour_contact_fraction=cap,
                                tumour_contact_inside_range=float((st.contact & st.inside_assessed_range).mean()),
                                tumour_contact_above_range=float((st.contact & ~st.inside_assessed_range).mean()),
                                tumour_lower_bounds_excluded=int(len(st) - len(usable)),
                                target_yield=target, P_operating=P_op,
                                leucocytes_retained_weighted_fraction=float(leu.w[leu.dP > P_op].sum() / leu.w.sum()),
                                leucocytes_retained_count=int((leu.dP > P_op).sum()), leucocytes_sampled=int(len(leu)),
                                retained_tumour_inside_range=float(kept.inside_assessed_range.mean()) if len(kept) else np.nan))
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
            L += ["Operating points (target fraction of ALL tumour cells retained; yield is capped by contact):", "",
                  v.drop(columns="D_aperture").round(4).to_markdown(index=False), ""]
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
