"""Measured (0.68) against transferred (0.50) nucleus-to-cell ratio: share of stored energy per compartment and critical
entry pressure of an MCF-7-like cell in the PDMS channel constriction, kappa/mu = 1e5, nucleus 2x or 4x stiffer than
the cytoplasm. Output: datos/nucleus_ratio.csv (about 1 min)."""
import sys, json, time, numpy as np, yaml
import os; ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')); os.chdir(ROOT); sys.path.insert(0, 'src')
from ctcfem.execution import force_single_thread; force_single_thread()
from ctcfem import geometry as G, reduction as RD
from ctcfem.cell_model import close_aperture
prof = G.profile_for_device(yaml.safe_load(open('data/materials.yaml'))['devices']['channel_PDMS'])
a = 8.74; Gs = 2.65 * a / 2.608
CONFS = [1.20, 1.35, 1.45, 1.55, 1.636, 1.80, 2.00]
out = []
for an_a in (0.50, 0.68):
    for r_nuc in (2.0, 4.0):
        cell = dict(a=a, an=an_a * a, Gs=Gs, G_cyt=0.1 * Gs / a, C=10.0, r_nuc=r_nuc, r_env=2.0, gamma=0.0)
        t = time.time()
        r = close_aperture(cell, prof, [a / c for c in CONFS], conf_max_attempt=2.01)
        sp = RD.spline_W(r['curve'], a) if len(r['curve']) > 1 else None
        for cp in r['checkpoints']:
            W = cp['W_total']; Rn = cp['rho']
            dP = RD.dP_local(prof, Rn, a, sp, r['rho_min_reached'])[0] if sp is not None else float('nan')
            out.append(dict(an_a=an_a, r_nuc=r_nuc, conf=round(cp['conf_eff'], 3), dP_kPa=dP,
                            cortex=cp['W_cortex'] / W, cytoplasm=cp['W_cytoplasm'] / W,
                            nucleus_plus_envelope=(cp['W_nucleus'] + cp['W_envelope']) / W))
        print(f'an/a={an_a} r_nuc={r_nuc}: reached a/R={a / r["rho_min_reached"]:.2f} ({r["stop_reason"]}) in {time.time() - t:.0f} s', flush=True)
import pandas as pd; pd.DataFrame(out).to_csv('manuscript_analysis/datos/nucleus_ratio.csv', index=False)
