import numpy as np
from ctcfem.geometry import Profile


def test_exact_aperture_and_no_corner():
    for zt, L in [(5.0, 1.0), (10.0, 2.0), (5.0, 10.0)]:
        p = Profile(zt=zt, L=L, R_ext=25.0)
        for Rn in (2.5, 4.0, 6.0, 7.5):
            assert abs(p.Rw(0.0, Rn) - Rn) < 1e-12           # the declared aperture is the real one
            assert abs(p.dRw(0.0, Rn)) < 1e-12               # zero slope at the neck
            assert np.all(p.Rw(np.linspace(-40, 40, 8001), Rn) >= Rn - 1e-12)
            assert abs(p.Rw(1e3, Rn) - 25.0) < 1e-9


def test_analytic_derivative():
    p = Profile(zt=5.0, L=2.0, R_ext=20.0); z = np.linspace(0.1, 30, 500); e = 1e-6
    assert np.max(np.abs((p.Rw(z + e, 4.0) - p.Rw(z - e, 4.0)) / (2 * e) - p.dRw(z, 4.0))) < 1e-6
