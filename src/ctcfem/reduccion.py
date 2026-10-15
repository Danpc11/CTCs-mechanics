"""Reduccion de W(rho) a presion critica de paso (kPa).

Tres estimaciones, reportadas por separado:

1. dP_local  : max_z |dW/drho (R_w(z))| |dR_w/dz| / (pi R_w(z)^2), con W indexado por el
               radio REAL del cuello (correccion del manuscrito, Sect. 2.6).  Apropiada para
               canales con transicion larga frente a la celula.  En membranas con borde
               agudo crece como 1/radio_borde: se reporta, pero no se recomienda.
2. dP_trabajo: W(R_n) / (pi R_n^2 l),  l = l_factor * a  (manuscrito, microfiltro).  Trabajo
               de la presion sobre el volumen que entra en el poro.  Recomendada para membranas.
3. dP_gota   : 2 gamma (1/R_n - 1/a)  (gota liquida con tension cortical; aspiracion por
               micropipeta).  Contraste analitico independiente del FEM para leucocitos.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator


def spline_W(curva, a):
    rho = np.array([p["rho"] for p in curva] + [a], float)
    W = np.array([p["W_total"] for p in curva] + [0.0], float)
    rho, idx = np.unique(rho, return_index=True)
    W = np.maximum.accumulate(np.maximum(W[idx][::-1], 0.0))[::-1]   # monotona decreciente en rho
    return PchipInterpolator(rho, W, extrapolate=False)


def dP_local(perfil, Rn, a, sp, rho_min_alcanzado, nz=40001):
    if Rn >= a:
        return 0.0, True
    z = np.linspace(0.0, perfil.zt + 8.0 * perfil.L, nz)[1:]
    Rw = perfil.Rw(z, Rn)
    ok = (Rw <= a) & (Rw >= max(Rn, rho_min_alcanzado) - 1e-9)
    q = np.abs(sp.derivative()(Rw)) * np.abs(perfil.dRw(z, Rn)) / (np.pi * Rw ** 2)
    q = np.where(ok & np.isfinite(q), q, 0.0)
    completo = rho_min_alcanzado <= Rn + 1e-9
    return float(q.max()), completo


def dP_trabajo(Rn, a, sp, rho_min_alcanzado, l_factor=1.0):
    if Rn >= a:
        return 0.0, True
    r = max(Rn, rho_min_alcanzado)
    W = float(sp(r)) if np.isfinite(sp(r)) else 0.0
    return W / (np.pi * Rn ** 2 * l_factor * a), rho_min_alcanzado <= Rn + 1e-9


def dP_gota(Rn, a, gamma_mN_m):
    if Rn >= a or gamma_mN_m <= 0:
        return 0.0
    return 2.0 * gamma_mN_m * (1.0 / Rn - 1.0 / a)     # mN/m / um = kPa
