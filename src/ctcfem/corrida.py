"""Orquestacion reproducible: muestreo -> tareas -> ejecucion paralela -> crudos JSON.

Garantias:
* La muestra de celulas se escribe UNA vez (muestras.csv) y el FEM la lee: re-ejecutar no
  re-muestrea.
* Cada tarea tiene un id determinista = sha256(config relevante + fila de la muestra +
  dispositivo).  Si el JSON existe, la tarea se salta (reanudacion tras corte de luz o kill).
* Escritura atomica (tmp + rename): nunca queda un JSON a medias.
* manifest.json guarda hashes de config, datos y codigo, versiones de librerias y maquina.
"""
from __future__ import annotations
import hashlib, json, os, platform, socket, sys, time, traceback
from pathlib import Path
import multiprocessing as mp

import numpy as np
import pandas as pd
import yaml

from . import poblaciones as P
from . import geometria as G
from . import reduccion as RD
from .modelo_celula import cerrar_apertura, CONF_EF_VALIDADO

RAIZ = Path(__file__).resolve().parents[2]


def _sha(obj) -> str:
    s = obj if isinstance(obj, (bytes, bytearray)) else json.dumps(obj, sort_keys=True, default=str).encode()
    return hashlib.sha256(s).hexdigest()


def _sha_archivo(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def cargar_config(path):
    cfg = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    cfg["_ruta"] = str(Path(path).resolve())
    return cfg


def aperturas(cfg):
    a = cfg["aperturas_um"]
    v = np.round(np.arange(a["inicio"], a["fin"] + 1e-9, a["paso"]), 4).tolist()
    v = sorted(set(v) | set(float(x) for x in a.get("extra", [])))
    return v


def dir_salida(cfg) -> Path:
    d = RAIZ / cfg["salida"]
    (d / "crudos").mkdir(parents=True, exist_ok=True)
    (d / "logs").mkdir(parents=True, exist_ok=True)
    return d


# ------------------------------------------------------------------ muestreo
def muestrear(cfg) -> pd.DataFrame:
    d = dir_salida(cfg)
    f = d / "muestras.csv"
    if f.exists():
        return pd.read_csv(f)
    datos = P.cargar(RAIZ / cfg["datos"]["poblaciones"])
    df = pd.concat([P.muestrear(datos, nm, int(n), int(cfg["semilla"]))
                    for nm, n in cfg["poblaciones"].items()], ignore_index=True)
    df.to_csv(f, index=False, float_format="%.17g")
    P.tabla_de_estados(datos).to_csv(d / "procedencia_parametros.csv", index=False)
    return pd.read_csv(f)          # SIEMPRE desde disco: ids identicos en la primera corrida y al reanudar


def tareas(cfg, muestras: pd.DataFrame):
    mats = yaml.safe_load((RAIZ / cfg["datos"]["materiales"]).read_text(encoding="utf-8"))
    radios = [0.5 * D for D in aperturas(cfg)]
    rel = {k: cfg[k] for k in ("modelo", "aperturas_um", "reduccion")}
    out = []
    for nd in cfg["dispositivos"]:
        disp = mats["dispositivos"][nd]
        mat = mats["materiales"][disp["material"]]
        for _, fila in muestras.iterrows():
            c = fila.to_dict()
            tid = _sha([rel, nd, disp, mat, {k: c[k] for k in sorted(c)}])[:20]
            out.append(dict(id=tid, dispositivo=nd, disp=disp, mat=mat, celula=c, radios=radios,
                            modelo=cfg["modelo"], reduccion=cfg["reduccion"]))
    return out


# ------------------------------------------------------------------ una tarea
def ejecutar_tarea(t, dir_crudos):
    destino = Path(dir_crudos) / f"{t['id']}.json"
    if destino.exists():
        return t["id"], "saltada"
    c, disp, mat, mo = t["celula"], t["disp"], t["mat"], t["modelo"]
    perfil = G.perfil_de_dispositivo(disp)
    rec = dict(id=t["id"], dispositivo=t["dispositivo"], material=disp["material"], tipo=disp["tipo"],
               poblacion=c["poblacion"], grupo=c["grupo"], muestra=int(c["muestra"]),
               celula=c, perfil=dict(zt=perfil.zt, L=perfil.L, R_ext=perfil.R_ext))
    try:
        E_kPa = float(mat["E_GPa"]["nominal"]) * 1e6 if mo.get("pared_compliante", True) else None
        l_w = disp["espesor_um"]["nominal"] if disp["tipo"] == "membrana" else None
        r = cerrar_apertura(c, perfil, t["radios"], conf_max_intento=mo["conf_max_intento"],
                            hrel=mo["hrel"], kn_f=mo["kn_f"], E_pared_kPa=E_kPa,
                            nu_pared=float(mat["nu"]), l_pared=l_w,
                            usar_gamma=bool(mo.get("usar_gamma", False)))
        rec["fem"] = r
        sp = RD.spline_W(r["curva"], c["a"]) if len(r["curva"]) > 1 else None
        filas = []
        for Rn in t["radios"]:
            a = c["a"]
            fila = dict(D_apertura=2 * Rn, conf_ef=a / Rn, contacto=bool(Rn < a),
                        dentro_envolvente=bool(a / Rn <= CONF_EF_VALIDADO))
            if Rn >= a:
                fila.update(dP_local=0.0, dP_trabajo=0.0, completo=True)
            elif sp is None:
                fila.update(dP_local=np.nan, dP_trabajo=np.nan, completo=False)
            else:
                pl, ok1 = RD.dP_local(perfil, Rn, a, sp, r["rho_min_alcanzado"])
                pw, ok2 = RD.dP_trabajo(Rn, a, sp, r["rho_min_alcanzado"], t["reduccion"]["l_factor"])
                fila.update(dP_local=pl, dP_trabajo=pw, completo=bool(ok1 and ok2))
            fila["dP_gota"] = RD.dP_gota(Rn, c["a"], c["gamma"])
            fila["L_carga_um"] = perfil.longitud_contacto(c["a"], Rn)
            fila["angulo_entrada_deg"] = perfil.angulo_entrada_deg(Rn)
            filas.append(fila)
        rec["aperturas"] = filas
        rec["estado"] = "ok"
    except Exception as exc:                                   # se registra, no se oculta
        rec["estado"] = "error"
        rec["error"] = repr(exc)[:500]
        rec["traza"] = traceback.format_exc()[-3000:]
    tmp = destino.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec, default=float), encoding="utf-8")
    os.replace(tmp, destino)
    return t["id"], rec["estado"]


def _worker(args):
    return ejecutar_tarea(*args)


def _manifest(cfg, d, n_tareas):
    import scipy, skfem, gmsh, matplotlib
    codigo = {p.name: _sha_archivo(p) for p in sorted((RAIZ / "src" / "ctcfem").glob("*.py"))}
    m = dict(fecha_inicio=time.strftime("%Y-%m-%dT%H:%M:%S%z"), maquina=socket.gethostname(),
             plataforma=platform.platform(), cpu=os.cpu_count(), python=sys.version,
             versiones=dict(numpy=np.__version__, scipy=scipy.__version__, scikit_fem=skfem.__version__,
                            gmsh=gmsh.__version__, pandas=pd.__version__, matplotlib=matplotlib.__version__,
                            pyyaml=yaml.__version__),
             config=_sha_archivo(cfg["_ruta"]),
             datos={k: _sha_archivo(RAIZ / v) for k, v in cfg["datos"].items()},
             codigo=codigo, n_tareas=n_tareas, semilla=cfg["semilla"])
    (d / "manifest.json").write_text(json.dumps(m, indent=1), encoding="utf-8")
    (d / "config_usada.yaml").write_text(Path(cfg["_ruta"]).read_text(encoding="utf-8"), encoding="utf-8")
    return m


def correr(cfg, workers=None, limite=None):
    d = dir_salida(cfg)
    mu = muestrear(cfg)
    ts = tareas(cfg, mu)
    # las celulas sin contacto en ninguna apertura no necesitan FEM, pero se registran igual
    ts.sort(key=lambda t: -t["celula"]["a"])          # las largas primero: mejor balance de carga
    if limite:
        ts = ts[:int(limite)]
    _manifest(cfg, d, len(ts))
    hechas = {p.stem for p in (d / "crudos").glob("*.json")}
    pend = [t for t in ts if t["id"] not in hechas]
    w = int(workers or cfg["ejecucion"]["workers"])
    log = open(d / "logs" / "progreso.log", "a", encoding="utf-8")
    log.write(f"{time.ctime()}  tareas={len(ts)} pendientes={len(pend)} workers={w}\n"); log.flush()
    t0, n, err = time.time(), 0, 0
    args = [(t, str(d / "crudos")) for t in pend]
    if w <= 1:
        it = map(_worker, args)
        pool = None
    else:
        pool = mp.get_context("fork").Pool(w, maxtasksperchild=int(cfg["ejecucion"].get("tareas_por_proceso", 25)))
        it = pool.imap_unordered(_worker, args, chunksize=1)
    for tid, est in it:
        n += 1; err += est == "error"
        if n % 10 == 0 or n == len(pend):
            el = time.time() - t0
            log.write(f"{time.ctime()}  {n}/{len(pend)}  errores={err}  {el/60:.1f} min  "
                      f"ETA {el/n*(len(pend)-n)/60:.1f} min\n"); log.flush()
    if pool:
        pool.close(); pool.join()
    log.write(f"{time.ctime()}  fin  hechas={n} errores={err}\n"); log.close()
    return n, err
