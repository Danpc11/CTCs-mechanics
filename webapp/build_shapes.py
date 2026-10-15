"""Compute the deformed FEM shapes shown by the web app -> webapp/fem_shapes.json (~1 min).

For a tumour-like (a = 8.74 um, a_n/a = 0.50) and a leucocyte-like (a = 4.20 um, a_n/a = 0.625) reference cell,
in each tabulated device, the aperture is closed quasi-statically from first contact to a/R = 2.4 and the
deformed mesh is stored at every 0.05 of confinement, with the equivalent Green-Lagrange strain per triangle.
Coordinates are normalised by the cell radius, so the app rescales them to any phenotype at the same
confinement (the shape depends on confinement and on the dimensionless groups, not on stiffness).
"""
import json, os, sys
import numpy as np, yaml
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
from ctcfem.execution import force_single_thread; force_single_thread()
from ctcfem import geometry as G
from ctcfem.cell_model import close_aperture
from ctcfem.solver_axisimetrico import axi_F, det3

mats = yaml.safe_load(open(os.path.join(ROOT, "data", "materials.yaml")))
CELLS = {"tumour": dict(a=8.74, an_a=0.50, E=2.65), "leucocyte": dict(a=4.20, an_a=0.625, E=0.43)}
CONFS = np.round(np.arange(1.05, 2.4001, 0.05), 3)


def fe_strain(mdl):
    """Equivalent Green-Lagrange strain from the finite-element field itself (quadrature points), averaged per
    element, and the element-mean volume ratio J as a check (near-incompressible: J ~ 1)."""
    b = mdl.basis; R = b.global_coordinates()[0]; uq = b.interpolate(mdl.u)
    F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], R)
    E = 0.5 * (np.einsum("kiep,kjep->ijep", F, F) - np.eye(3)[:, :, None, None])
    eq = np.sqrt(2.0 / 3.0 * np.einsum("ijep,ijep->ep", E, E))
    return eq.mean(axis=1), det3(F).mean(axis=1)


JCHECK = []
out = dict(confinements=[1.0] + CONFS.tolist(), cells={})
for cn, cc in CELLS.items():
    a = cc["a"]; Gs = cc["E"] * a / 2.608
    cell = dict(a=a, an=cc["an_a"] * a, Gs=Gs, G_cyt=0.10 * Gs / a, C=10.0, r_nuc=2.0, r_env=2.0, gamma=0.0)
    entry = dict(an_a=cc["an_a"], devices={})
    for dev in ("membrane_PC", "channel_PDMS"):
        prof = G.profile_for_device(mats["devices"][dev]); frames = []
        def grab(m, mdl, cp, frames=frames):
            u = mdl.x[mdl.basis.nodal_dofs]; X = m.p
            eq, J = fe_strain(mdl); JCHECK.append((cp["conf_eff"], float(J.min()), float(J.max())))
            frames.append((cp["conf_eff"], u.copy(), eq))
            if "mesh" not in entry:
                entry["mesh"] = dict(p=np.round(X / a, 4).tolist(), t=m.t.tolist(),
                                     nucleus=np.isin(np.arange(m.t.shape[1]), mdl.subdomains["nucleus"]).astype(int).tolist())
        r = close_aperture(cell, prof, [a / c for c in CONFS], conf_max_attempt=2.401, on_checkpoint=grab)
        frames.sort(key=lambda f: f[0])
        entry["devices"][dev] = dict(conf=[1.0] + [round(f[0], 3) for f in frames],
                                     u=[[[0.0] * len(frames[0][1][0])] * 2] + [np.round(f[1] / a, 4).tolist() for f in frames],
                                     strain=[[0.0] * len(frames[0][2])] + [np.round(f[2], 4).tolist() for f in frames])
        print(f"{cn:9s} {dev:12s} frames {len(frames)+1}  reached a/R = {a / r['rho_min_reached']:.2f}  ({r['stop_reason']})", flush=True)
    out["cells"][cn] = entry
print("finite-element volume ratio J per element: min %.5f  max %.5f" % (min(j[1] for j in JCHECK), max(j[2] for j in JCHECK)))
json.dump(out, open(os.path.join(HERE, "fem_shapes.json"), "w"), separators=(",", ":"))
print("fem_shapes.json", os.path.getsize(os.path.join(HERE, "fem_shapes.json")) // 1024, "kB")
