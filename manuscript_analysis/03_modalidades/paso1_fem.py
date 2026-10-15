"""Paso 1-A: paisaje energetico W(R) por estirpe y por geometria de
constriccion.  Es el insumo FEM de todo lo demas (pasos 1, 2, 3 y 5)."""
import sys, os, json, time, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fase2_lib as L

# Gs calibrado contra el modulo aparente medido (paso0_calibracion.py).
# MCF-10A no tiene E_app medido: se TRANSFIERE el Gs de MCF-7 y se declara.
GS_CAL = {nm: (L.Gs_de_Eapp(e.E_app, e.a) if e.E_app is not None else None)
          for nm, e in L.ESTIRPES.items()}
GS_CAL["MCF-10A"] = GS_CAL["MCF-7"]      # E_app no medido: transferido
GS_TRANSFERIDO = {"MCF-10A": "Gs transferido de MCF-7 (E_app no medido)"}

# confinamientos: 1.05 y 1.12 por DEBAJO de la envolvente (1.2), 2.65 por
# ENCIMA del techo validado (2.43).  Se calculan y se excluyen, no se omiten.
CONFS = [1.05, 1.12, 1.20, 1.30, 1.42, 1.55, 1.70, 1.90, 2.10, 2.30, 2.43, 2.65]

# (etiqueta, semilongitud de garganta zt, longitud de transicion Ltap) en um
GEOMS = [("referencia", 4.0, 5.0),
         ("garganta_corta", 2.0, 2.5),
         ("garganta_larga", 8.0, 10.0),
         ("entrada_abrupta", 4.0, 2.0),
         ("entrada_suave", 4.0, 12.0)]


def tarea(arg):
    nm, geo, zt, Ltap, c = arg
    e = L.ESTIRPES[nm]
    Rt = e.a / c
    t0 = time.time()
    try:
        gs = GS_CAL[nm]
        rec = L.resolver_constriccion(
            e.a, e.an, gs, e.C_skalak, L.PI_REG * gs / e.a, e.r_nuc, e.r_env,
            Rt, zt=zt, Ltap=Ltap)
    except Exception as exc:
        rec = dict(Rt=Rt, conf=c, converge=False, error=repr(exc)[:200])
    gs = GS_CAL[nm]
    rec.update(estirpe=nm, geometria=geo, zt=zt, Ltap=Ltap, Gs=gs,
               a=e.a, an=e.an, G_cyt=L.PI_REG * gs / e.a, Pi=L.PI_REG,
               G_cyt_medido=e.G_cyt, Pi_medido=e.G_cyt * e.a / gs,
               conf_nominal=c, Rwide=1.6 * e.a,
               conf_efectivo=L.conf_efectivo(e.a, Rt, 1.6 * e.a, zt, Ltap),
               theta_entrada=L.angulo_de_Ltap(Ltap, Rt, 1.6 * e.a, zt),
               L_carga_um=L.longitud_de_carga(e.a, Rt, 1.6 * e.a, zt, Ltap),
               dentro_envolvente=bool(L.dentro_envolvente(c)),
               segundos=round(time.time() - t0, 2))
    if rec.get("converge"):
        rec["f_adim"] = rec["W_total"] / (gs * e.a**2)
    return rec


if __name__ == "__main__":
    # el sandbox no permite semaforos POSIX: se reparte por procesos
    # independientes lanzados desde el shell (--shard i --de n)
    shard, nsh = 0, 1
    if "--shard" in sys.argv:
        shard = int(sys.argv[sys.argv.index("--shard") + 1])
        nsh = int(sys.argv[sys.argv.index("--de") + 1])
    args = [(nm, geo, zt, Lt, c)
            for nm in L.ESTIRPES
            for (geo, zt, Lt) in GEOMS
            for c in CONFS]
    mios = args[shard::nsh]
    t0 = time.time()
    res = []
    for i, a in enumerate(mios):
        res.append(tarea(a))
        if (i + 1) % 10 == 0:
            print("  shard %d: %3d/%d  (%.0f s)"
                  % (shard, i + 1, len(mios), time.time() - t0), flush=True)
    json.dump({"GS_CAL": GS_CAL, "CONFS": CONFS, "GEOMS": GEOMS,
               "registros": res, "segundos": round(time.time() - t0, 1)},
              open(os.path.join(HERE, "paso1_fem_s%d.json" % shard), "w"), indent=1)
    print("shard %d listo: %d registros, %.0f s"
          % (shard, len(res), time.time() - t0), flush=True)
