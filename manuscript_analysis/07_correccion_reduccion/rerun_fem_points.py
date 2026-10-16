"""Finite-element points behind Table 2 and the entry-geometry comparison, in parallel (one thread per worker).
Writes the three *_reejecucion.jsonl files of ../datos/ (kappa/mu from CTC_KAPPA_MU)."""
import sys, os, json, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(HERE, "..", "01_solver"), os.path.join(HERE, "..", "03_modalidades"), HERE]
REF = [1.12, 1.20, 1.30, 1.42, 1.55, 1.70, 1.90, 2.10, 2.43]
GEO = [1.05, 1.12, 1.20, 1.30, 1.42, 1.55, 1.70]
JOBS = ([("fem_mcf7_referencia_reejecucion.jsonl", ("MCF-7", "referencia", 4.0, 5.0, c)) for c in REF] +
        [("fem_lineas_referencia_reejecucion.jsonl", (nm, "referencia", 4.0, 5.0, c)) for nm in ("MDA-MB-231", "leucocito", "SK-BR-3") for c in REF] +
        [("fem_mcf7_geometrias_reejecucion.jsonl", ("MCF-7", g, zt, Lt, c)) for g, zt, Lt in
         (("garganta_larga", 8.0, 10.0), ("garganta_corta", 2.0, 2.5), ("entrada_abrupta", 4.0, 2.0), ("entrada_suave", 4.0, 12.0)) for c in GEO])


def run(job):
    import paso1_fem as P
    f, args = job; t = time.time(); r = P.tarea(args)
    keep = {k: v for k, v in r.items() if isinstance(v, (int, float, str, bool)) or v is None}
    keep["seg_total"] = time.time() - t
    return f, keep


if __name__ == "__main__":
    D = os.path.join(HERE, "..", "datos")
    for f in {j[0] for j in JOBS}:
        open(os.path.join(D, f), "w").close()
    with Pool(int(os.environ.get("CTC_WORKERS", "8"))) as po:
        for i, (f, rec) in enumerate(po.imap_unordered(run, JOBS), 1):
            open(os.path.join(D, f), "a").write(json.dumps(rec) + "\n")
            print(f"  {i}/{len(JOBS)}  {rec.get('estirpe')} {rec.get('geometria')} {rec.get('conf_nominal')}  converge={rec.get('converge')}", flush=True)
