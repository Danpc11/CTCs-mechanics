import numpy as np
from ctcfem.geometria import Perfil

def test_apertura_exacta_y_sin_quiebre():
    for zt, L in [(5.0, 1.0), (10.0, 2.0), (5.0, 10.0)]:
        p = Perfil(zt=zt, L=L, R_ext=25.0)
        for Rn in (2.5, 4.0, 6.0, 7.5):
            assert abs(p.Rw(0.0, Rn) - Rn) < 1e-12           # el diametro declarado es el real
            assert abs(p.dRw(0.0, Rn)) < 1e-12               # derivada nula en el cuello
            z = np.linspace(-40, 40, 8001)
            assert np.all(p.Rw(z, Rn) >= Rn - 1e-12)         # minimo global en el cuello
            assert abs(p.Rw(1e3, Rn) - 25.0) < 1e-9

def test_derivada_analitica():
    p = Perfil(zt=5.0, L=2.0, R_ext=20.0)
    z = np.linspace(0.1, 30, 500); e = 1e-6
    num = (p.Rw(z + e, 4.0) - p.Rw(z - e, 4.0)) / (2 * e)
    assert np.max(np.abs(num - p.dRw(z, 4.0))) < 1e-6
