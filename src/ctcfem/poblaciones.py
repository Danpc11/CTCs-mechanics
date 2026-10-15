"""Modelo poblacional: muestreo reproducible de celulas a partir de data/poblaciones.yaml.

* Hipercubo latino con semilla derivada de (semilla_global, nombre de poblacion)
  -> la muestra de una poblacion NO depende de que otras poblaciones se pidan ni
  del orden de ejecucion.
* Marginales truncadas por transformada inversa (normal o lognormal).
* Parametros independientes entre si (DECLARADO: no hay medidas conjuntas de
  tamano y rigidez por subpoblacion).
"""
from __future__ import annotations
import hashlib
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from scipy.stats import qmc

PARAMS = ["D", "an_a", "E_app", "gamma", "eta", "tau"]


def cargar(path):
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _semilla(semilla_global: int, nombre: str) -> int:
    h = hashlib.sha256(f"{semilla_global}:{nombre}".encode()).hexdigest()
    return int(h[:15], 16)


def _ppf(spec: dict, u: np.ndarray) -> np.ndarray:
    """Cuantil truncado de la marginal declarada."""
    d = spec["dist"]
    if d == "fijo":
        return np.full_like(u, float(spec["valor"]))
    m = float(spec["media"])
    if d == "normal":
        sd = float(spec["sd"]) if "sd" in spec else float(spec["cv"]) * m
        dist = stats.norm(m, sd)
    elif d == "lognormal":
        cv = float(spec["cv"])
        s = np.sqrt(np.log(1.0 + cv ** 2))
        dist = stats.lognorm(s, scale=m / np.sqrt(1.0 + cv ** 2))   # media aritmetica = m
    else:
        raise ValueError(d)
    lo, hi = spec.get("trunc", [-np.inf, np.inf])
    Flo, Fhi = dist.cdf(lo), dist.cdf(hi)
    return dist.ppf(Flo + u * (Fhi - Flo))


def muestrear(datos: dict, nombre: str, n: int, semilla_global: int) -> pd.DataFrame:
    pob = datos["poblaciones"][nombre]
    mod = datos["modelo"]
    lhs = qmc.LatinHypercube(d=len(PARAMS), seed=_semilla(semilla_global, nombre))
    U = lhs.random(n)
    df = pd.DataFrame({p: _ppf(pob[p], U[:, j]) for j, p in enumerate(PARAMS)})
    df.insert(0, "poblacion", nombre)
    df.insert(1, "grupo", pob["grupo"])
    df.insert(2, "muestra", np.arange(n))
    # magnitudes derivadas que usa el FEM (unidades: um, mN/m, kPa)
    df["a"] = 0.5 * df.D
    df["an"] = df.an_a * df.a
    df["Gs"] = df.E_app * df.a / mod["K_E"]          # mN/m  (E_app kPa * um / K_E)
    df["G_cyt"] = mod["Pi"] * df.Gs / df.a          # kPa
    df["C"] = mod["C_skalak"]
    df["r_nuc"] = mod["r_nuc"]
    df["r_env"] = mod["r_env"]
    return df


def tabla_de_estados(datos: dict) -> pd.DataFrame:
    """Tabla de procedencia: poblacion x parametro -> estado y fuente."""
    filas = []
    for nm, pob in datos["poblaciones"].items():
        for p in PARAMS:
            s = pob[p]
            filas.append(dict(poblacion=nm, parametro=p, dist=s["dist"],
                              estado=s.get("estado", ""), fuente=s.get("fuente", "")))
    return pd.DataFrame(filas)
