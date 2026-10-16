"""Recompute webapp/energy_curves.json with the finite-element model (~40 s): normalised confinement energy
g(rho/a) = W / (Gs a^2) of an MCF-7 cell from first contact to a/R = 2.4, for each tabulated device."""
import json, os, sys, yaml
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.dirname(H)
sys.path.insert(0, os.path.join(R, "src"))
from ctcfem.execution import force_single_thread; force_single_thread()
from ctcfem import geometry as G, populations as P
from ctcfem.cell_model import close_aperture
c = P.sample(P.load(os.path.join(R, "data", "populations.yaml")), "MCF-7", 1, 20261015).iloc[0].to_dict()
mats = yaml.safe_load(open(os.path.join(R, "data", "materials.yaml"))); out = {}
for nd in ("membrane_PC", "channel_PDMS"):
    dev = mats["devices"][nd]; mat = mats["materials"][dev["material"]]; prof = G.profile_for_device(dev)
    r = close_aperture(c, prof, [c["a"] / x for x in (1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.4)], conf_max_attempt=2.4,
                       E_wall_kPa=mat["E_GPa"]["nominal"] * 1e6, nu_wall=mat["nu"],
                       l_wall=dev["thickness_um"]["nominal"] if dev["type"] == "membrane" else None)
    cur = sorted(r["curve"], key=lambda p: p["rho"])
    out[nd] = dict(u=[p["rho"] / c["a"] for p in cur], g=[p["W_total"] / (c["Gs"] * c["a"] ** 2) for p in cur],
                   profile=dict(zt=prof.zt, L=prof.L, R_ext=prof.R_ext))
    print(nd, "reached a/R =", round(c["a"] / r["rho_min_reached"], 2))
json.dump(out, open(os.path.join(H, "energy_curves.json"), "w"), indent=0)
