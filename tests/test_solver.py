"""Patch test del manuscrito (V1): extension uniaxial homogenea, error < 1e-10."""
import numpy as np
from scipy.optimize import brentq
from skfem import MeshTri
from ctcfem.solver_axisimetrico import NeoHookean, AxiModel, piola

def test_patch_extension_uniaxial():
    R0, H = 1.0, 2.0
    m = MeshTri.init_tensor(np.linspace(0, R0, 5), np.linspace(0, H, 7))
    mdl = AxiModel(m, {"default": NeoHookean(C10=0.5, kappa=1e4)})
    mdl.add_dirichlet(0, lambda x: x[0] < 1e-12, 0.0)
    mdl.add_dirichlet(1, lambda x: x[1] < 1e-12, 0.0)
    lt = 1.3
    mdl.add_dirichlet(1, lambda x: np.abs(x[1] - H) < 1e-12, (lt - 1.0) * H)
    mdl.solve([0.25, 0.5, 0.75, 1.0])
    def Prr(lr):
        F = np.zeros((3, 3, 1, 1)); F[0, 0] = lr; F[1, 1] = lt; F[2, 2] = lr
        return piola(F, np.full((1, 1), 1e4 * (lr * lr * lt - 1.0)), np.full((1, 1), 0.5))[0, 0, 0, 0]
    lr = brentq(Prr, 0.3, 1.5)
    ud = mdl.x[mdl.basis.nodal_dofs]
    lr_fe = float(np.mean(1.0 + ud[0, np.abs(m.p[0] - R0) < 1e-12] / R0))
    assert abs(lr_fe - lr) / abs(lr - 1) < 1e-10
