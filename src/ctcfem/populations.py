"""Reproducible sampling of cell populations from data/populations.yaml.

* Latin hypercube with a seed derived from (global seed, population name): the sample of one
  population does not depend on which other populations are requested or on execution order.
* Truncated marginals by inverse transform (normal or lognormal).
* Parameters are sampled independently (declared: no joint measurements of size and stiffness
  per subpopulation were found).
"""
from __future__ import annotations
import hashlib
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from scipy.stats import qmc

PARAMS = ["D", "an_a", "E_app", "gamma", "eta", "tau"]


def load(path):
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _seed(global_seed: int, name: str) -> int:
    return int(hashlib.sha256(f"{global_seed}:{name}".encode()).hexdigest()[:15], 16)


def _ppf(spec: dict, u: np.ndarray) -> np.ndarray:
    d = spec["dist"]
    if d == "fixed":
        return np.full_like(u, float(spec["value"]))
    m = float(spec["mean"])
    if d == "normal":
        dist = stats.norm(m, float(spec["sd"]) if "sd" in spec else float(spec["cv"]) * m)
    elif d == "lognormal":
        cv = float(spec["cv"]); s = np.sqrt(np.log(1.0 + cv ** 2))
        dist = stats.lognorm(s, scale=m / np.sqrt(1.0 + cv ** 2))      # arithmetic mean = m
    else:
        raise ValueError(d)
    lo, hi = spec.get("trunc", [-np.inf, np.inf])
    Flo, Fhi = dist.cdf(lo), dist.cdf(hi)
    return dist.ppf(Flo + u * (Fhi - Flo))


def sample(data: dict, name: str, n: int, global_seed: int) -> pd.DataFrame:
    pop, mod = data["populations"][name], data["model"]
    U = qmc.LatinHypercube(d=len(PARAMS), seed=_seed(global_seed, name)).random(n)
    df = pd.DataFrame({p: _ppf(pop[p], U[:, j]) for j, p in enumerate(PARAMS)})
    df.insert(0, "population", name); df.insert(1, "group", pop["group"]); df.insert(2, "sample", np.arange(n))
    # derived quantities used by the FEM (units: um, mN/m, kPa)
    df["a"] = 0.5 * df.D
    df["an"] = df.an_a * df.a
    df["Gs"] = df.E_app * df.a / mod["K_E"]
    df["G_cyt"] = mod["Pi"] * df.Gs / df.a
    df["C"], df["r_nuc"], df["r_env"] = mod["C_skalak"], mod["r_nuc"], mod["r_env"]
    return df


def provenance_table(data: dict) -> pd.DataFrame:
    rows = [dict(population=nm, parameter=p, dist=pop[p]["dist"], status=pop[p].get("status", ""),
                 source=pop[p].get("source", "")) for nm, pop in data["populations"].items() for p in PARAMS]
    return pd.DataFrame(rows)
