import numpy as np
from ctcfem.geometry import Profile
from ctcfem import reduction as RD


def _curve(a, rmin):
    return [dict(rho=r, W_total=(a - r) ** 2) for r in np.linspace(a, rmin, 40)]


def test_zero_without_contact_and_positive_with_contact():
    a, p = 4.2, Profile(zt=5.0, L=1.0, R_ext=25.0); sp = RD.spline_W(_curve(a, 3.0), a)
    assert RD.dP_local(p, 5.0, a, sp, 3.0)[0] == 0.0
    v, ok = RD.dP_local(p, 4.0, a, sp, 3.0)
    assert v > 0 and ok                                  # slight contact: finite barrier, not zero


def test_lower_bound_flag_when_not_reached():
    a, p = 8.74, Profile(zt=5.0, L=5.0, R_ext=15.0)
    assert not RD.dP_local(p, 4.0, a, RD.spline_W(_curve(a, 6.0), a), 6.0)[1]


def test_liquid_drop_laplace():
    assert abs(RD.dP_drop(4.0, 4.2, 0.037) - 2 * 0.037 * (1 / 4.0 - 1 / 4.2)) < 1e-15
