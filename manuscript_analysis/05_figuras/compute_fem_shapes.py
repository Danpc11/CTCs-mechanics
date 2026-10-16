"""Deformed finite-element states used in Figs. 1 and 3 -> datos/fem_shapes_channel.json (~1 min).
Tumour-like (a = 8.74 um, a_n/a = 0.50) and leucocyte-like (a = 4.20 um, a_n/a = 0.625) cells in the PDMS channel
constriction at effective confinements 1.20, 1.45 and 1.636, with element strain and the share of stored energy per
compartment. kappa/mu = 1e3, the value used for the numbers reported in the manuscript."""
import sys, os, json, numpy as np, yaml
R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")); sys.path.insert(0, os.path.join(R, "src"))
from ctcfem.execution import force_single_thread; force_single_thread()
from ctcfem import geometry as G, cell_model as CM
from ctcfem.solver_axisimetrico import axi_F
CM.KAPPA_MU = float(os.environ.get("CTC_KAPPA_MU", "1e3"))
prof = G.profile_for_device(yaml.safe_load(open(os.path.join(R, "data", "materials.yaml")))["devices"]["channel_PDMS"])
out = {}
for name, a, an_a, E in [("tumour", 8.74, 0.50, 2.65), ("leucocyte", 4.20, 0.625, 0.43)]:
    Gs = E * a / 2.608; cell = dict(a=a, an=an_a * a, Gs=Gs, G_cyt=0.1 * Gs / a, C=10., r_nuc=2., r_env=2., gamma=0.)
    fr, mesh = [], {}
    def grab(m, mdl, cp):
        b = mdl.basis; Rq = b.global_coordinates()[0]; uq = b.interpolate(mdl.u); F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], Rq)
        Ee = 0.5 * (np.einsum("kiep,kjep->ijep", F, F) - np.eye(3)[:, :, None, None])
        fr.append(dict(conf=cp["conf_eff"], u=(mdl.x[mdl.basis.nodal_dofs] / a).tolist(),
                       strain=np.sqrt(2 / 3 * np.einsum("ijep,ijep->ep", Ee, Ee)).mean(1).tolist(), W=cp["W_total"],
                       Wc={k: cp[k] for k in ("W_cortex", "W_cytoplasm", "W_nucleus", "W_envelope")}))
        mesh.update(p=(m.p / a).tolist(), t=m.t.tolist(), nuc=np.isin(np.arange(m.t.shape[1]), mdl.subdomains["nucleus"]).astype(int).tolist())
    CM.close_aperture(cell, prof, [a / c for c in (1.2, 1.45, 1.636)], conf_max_attempt=1.7, on_checkpoint=grab)
    out[name] = dict(a=a, an_a=an_a, mesh=mesh, frames=fr)
D = os.path.join(R, "manuscript_analysis", "datos")
json.dump(out, open(os.path.join(D, "fem_shapes_channel.json"), "w"))
json.dump(dict(profile=dict(zt=prof.zt, L=prof.L, R_ext=prof.R_ext)), open(os.path.join(D, "fem_shapes_channel_profile.json"), "w"))
