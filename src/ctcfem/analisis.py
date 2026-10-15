"""Agregacion de crudos -> tablas y figuras.

Estimador recomendado por tipo de dispositivo (ver reduccion.py):
    canal    -> dP_local      membrana -> dP_trabajo
Las otras estimaciones se conservan en la tabla larga.

Salidas en <salida>/analisis/:
  tabla_larga.csv           una fila por (celula, dispositivo, apertura)
  resumen_aperturas.csv     cuantiles por poblacion, dispositivo y apertura
  ventanas.csv              ventana de presion leucocito/tumoral por dispositivo y apertura
  informe_8_10_12.md        lectura directa de las aperturas de interes
  fig_<dispositivo>.png     dP frente a apertura (mediana y p5-p95) por poblacion
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

ESTIMADOR = {"canal": "dP_local", "membrana": "dP_trabajo"}


def tabla_larga(d: Path) -> pd.DataFrame:
    filas = []
    for f in sorted((d / "crudos").glob("*.json")):
        r = json.loads(f.read_text(encoding="utf-8"))
        base = dict(id=r["id"], dispositivo=r["dispositivo"], material=r["material"], tipo=r["tipo"],
                    poblacion=r["poblacion"], grupo=r["grupo"], muestra=r["muestra"], estado=r["estado"],
                    D_celula=r["celula"]["D"], E_app=r["celula"]["E_app"], tau=r["celula"]["tau"])
        if r["estado"] != "ok":
            filas.append(base); continue
        fem = r["fem"]
        base.update(motivo_parada=fem.get("motivo_parada"), rho_min=fem.get("rho_min_alcanzado"),
                    segundos=fem.get("segundos", 0.0))
        for a in r["aperturas"]:
            fila = dict(base, **a)
            fila["dP"] = a[ESTIMADOR[r["tipo"]]]
            # umbral de presion de dispositivo (Kelvin-Voigt), kPa por mm/s
            fila["beta_c"] = (fila["dP"] * base["tau"] / a["L_carga_um"] * 1e3
                              if a["L_carga_um"] > 0 else 0.0)
            filas.append(fila)
    return pd.DataFrame(filas)


def _q(x, p):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.quantile(x, p)) if len(x) else np.nan


def resumen(T: pd.DataFrame) -> pd.DataFrame:
    g = T[T.estado == "ok"].groupby(["dispositivo", "poblacion", "grupo", "D_apertura"])
    return g.apply(lambda s: pd.Series(dict(
        n=len(s), frac_contacto=s.contacto.mean(), frac_completo=s.completo.mean(),
        frac_dentro_envolvente=s[s.contacto].dentro_envolvente.mean() if s.contacto.any() else 1.0,
        dP_p5=_q(s.dP, .05), dP_p50=_q(s.dP, .5), dP_p95=_q(s.dP, .95),
        dP_gota_p50=_q(s.dP_gota, .5))), include_groups=False).reset_index()


def ventanas(T: pd.DataFrame, diferencial: dict, prevalencia=1e-6) -> pd.DataFrame:
    """Ventana [P99 de la mezcla leucocitaria ponderada, P1 de cada linea tumoral]."""
    ok = T[(T.estado == "ok")]
    out = []
    for (nd, D), s in ok.groupby(["dispositivo", "D_apertura"]):
        leu = s[s.grupo == "leucocito"]
        if leu.empty:
            continue
        # cuantil ponderado de la mezcla: remuestreo por pesos del diferencial
        w = leu.poblacion.map(diferencial).fillna(0.0)
        w = w / leu.groupby("poblacion").poblacion.transform("count")
        x, ww = leu.dP.to_numpy(float), w.to_numpy(float)
        o = np.argsort(x); cw = np.cumsum(ww[o]) / ww.sum()
        p99_leu = float(x[o][np.searchsorted(cw, 0.99)])
        for tum, st in s[s.grupo == "tumoral"].groupby("poblacion"):
            p1_tum = _q(st.dP, 0.01)
            P_op = np.sqrt(max(p99_leu, 1e-6) * max(p1_tum, 1e-6)) if p1_tum > p99_leu else np.nan
            ret_tum = float((st.dP > P_op).mean()) if np.isfinite(P_op) else np.nan
            ret_leu = float((ww[x > P_op].sum() / ww.sum())) if np.isfinite(P_op) else np.nan
            out.append(dict(dispositivo=nd, D_apertura=D, tumoral=tum, P99_leucocito=p99_leu,
                            P1_tumoral=p1_tum, ventana_existe=bool(p1_tum > p99_leu),
                            ventana_factor=p1_tum / p99_leu if p99_leu > 0 else np.inf,
                            P_operacion=P_op, rendimiento_tumoral=ret_tum,
                            arrastre_leu_por_millon=ret_leu * 1e6 if np.isfinite(ret_leu) else np.nan,
                            frac_tumoral_completo=float(st.completo.mean()),
                            frac_tumoral_dentro_envolvente=float(st.dentro_envolvente.mean())))
    return pd.DataFrame(out)


def informe(R: pd.DataFrame, V: pd.DataFrame, foco=(8.0, 10.0, 12.0)) -> str:
    L = ["# Lectura directa de las aperturas de interes", "",
         "Estimador: dP_local (canales) y dP_trabajo (membranas). 'completo' = la continuacion",
         "llego a esa apertura; si no, el valor es COTA INFERIOR. 'dentro_envolvente' = confinamiento",
         "efectivo <= 1.636 (validado). Fuera de envolvente el resultado es extrapolacion.", ""]
    for D in foco:
        L += [f"## Apertura {D:g} um", ""]
        r = R[np.isclose(R.D_apertura, D)]
        if r.empty:
            L += ["(sin datos)", ""]; continue
        L += [r[["dispositivo", "poblacion", "frac_contacto", "frac_completo", "frac_dentro_envolvente",
                 "dP_p5", "dP_p50", "dP_p95"]].round(4).to_markdown(index=False), ""]
        v = V[np.isclose(V.D_apertura, D)]
        if not v.empty:
            L += ["Ventanas:", "", v.drop(columns="D_apertura").round(4).to_markdown(index=False), ""]
    return "\n".join(L)


def figuras(R: pd.DataFrame, d: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for nd, s in R.groupby("dispositivo"):
        fig, ax = plt.subplots(figsize=(6.85, 3.6))
        for pob, q in s.groupby("poblacion"):
            q = q.sort_values("D_apertura")
            ls = "-" if q.grupo.iloc[0] == "tumoral" else ("--" if q.grupo.iloc[0] == "leucocito" else ":")
            ax.plot(q.D_apertura, q.dP_p50.clip(lower=1e-5), ls, label=pob)
            ax.fill_between(q.D_apertura, q.dP_p5.clip(lower=1e-5), q.dP_p95.clip(lower=1e-5), alpha=0.15)
        for D in (8, 10, 12):
            ax.axvline(D, color="0.6", lw=0.6)
        ax.set_yscale("log"); ax.set_xlabel("diametro de apertura (um)"); ax.set_ylabel("presion critica (kPa)")
        ax.legend(fontsize=6, ncol=2, frameon=False); fig.tight_layout()
        fig.savefig(d / f"fig_{nd}.png", dpi=300); plt.close(fig)


def analizar(cfg, raiz: Path):
    d = raiz / cfg["salida"]
    out = d / "analisis"; out.mkdir(exist_ok=True)
    T = tabla_larga(d); T.to_csv(out / "tabla_larga.csv", index=False)
    if T.empty or "D_apertura" not in T:
        return T, None, None
    R = resumen(T); R.to_csv(out / "resumen_aperturas.csv", index=False)
    dif = yaml.safe_load((raiz / cfg["datos"]["poblaciones"]).read_text(encoding="utf-8"))["diferencial"]
    V = ventanas(T, dif); V.to_csv(out / "ventanas.csv", index=False)
    (out / "informe_8_10_12.md").write_text(informe(R, V), encoding="utf-8")
    figuras(R, out)
    return T, R, V
