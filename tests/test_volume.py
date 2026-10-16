"""The confined cell must keep its volume: with a soft cytoplasm, a bulk modulus tied to it must be large."""
import numpy as np
from ctcfem.geometry import Profile
from ctcfem.cell_model import close_aperture
from ctcfem.solver_axisimetrico import axi_F, det3


def test_volume_conserved_at_the_top_of_the_assessed_range():
    a = 8.74; Gs = 2.65 * a / 2.608
    cell = dict(a=a, an=0.5 * a, Gs=Gs, G_cyt=0.1 * Gs / a, C=10.0, r_nuc=2.0, r_env=2.0, gamma=0.0)
    J = []
    def grab(m, mdl, cp):
        b = mdl.basis; R = b.global_coordinates()[0]; uq = b.interpolate(mdl.u)
        F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], R); w = 2 * np.pi * R * b.dx
        J.append(float((det3(F) * w).sum() / w.sum()))
    r = close_aperture(cell, Profile(zt=5.0, L=5.0, R_ext=15.0), [a / 1.636], conf_max_attempt=1.7, on_checkpoint=grab)
    assert r["complete"] and abs(J[-1] - 1.0) < 0.005
