import numpy as np
from ctcfem.geometria import Perfil
from ctcfem import reduccion as RD

def _curva(a, rmin):
    rho = np.linspace(a, rmin, 40)
    return [dict(rho=r, W_total=(a - r) ** 2) for r in rho]

def test_sin_contacto_es_cero_y_con_contacto_positivo():
    a, p = 4.2, Perfil(zt=5.0, L=1.0, R_ext=25.0)
    sp = RD.spline_W(_curva(a, 3.0), a)
    assert RD.dP_local(p, 5.0, a, sp, 3.0)[0] == 0.0          # apertura mayor que la celula
    v, ok = RD.dP_local(p, 4.0, a, sp, 3.0)
    assert v > 0 and ok                                       # contacto leve: barrera finita, no cero

def test_cota_inferior_si_no_se_alcanza():
    a, p = 8.74, Perfil(zt=5.0, L=5.0, R_ext=15.0)
    sp = RD.spline_W(_curva(a, 6.0), a)
    _, ok = RD.dP_local(p, 4.0, a, sp, 6.0)
    assert not ok                                             # se marca como incompleto

def test_gota_liquida_laplace():
    # neutrofilo (a = 4.2 um, gamma = 0.037 mN/m) en poro de 8 um: dP = 2 gamma (1/4 - 1/4.2)
    assert abs(RD.dP_gota(4.0, 4.2, 0.037) - 2 * 0.037 * (1 / 4.0 - 1 / 4.2)) < 1e-15
