"""Reduction of the energy landscape W(rho) to a critical entry pressure (kPa).

* dP_local: max_z |dW/drho(R_w(z))| |dR_w/dz| / (pi R_w(z)^2), with W indexed by the radius actually
  reached at the neck. Suited to channels whose transition is long compared with the cell. In membranes
  with a sharp edge it grows as 1/edge radius; it is reported but not recommended there.
* dP_work: W(R_n) / (pi R_n^2 l), l = l_factor * a: work of the pressure on the volume that enters the
  pore. Recommended for membranes.
* dP_drop: 2 gamma (1/R_n - 1/a): liquid drop with cortical tension (micropipette aspiration); an
  analytical check independent of the finite-element model, for leucocytes.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator


def spline_W(curve, a):
    rho = np.array([p["rho"] for p in curve] + [a], float)
    W = np.array([p["W_total"] for p in curve] + [0.0], float)
    rho, idx = np.unique(rho, return_index=True)
    W = np.maximum.accumulate(np.maximum(W[idx][::-1], 0.0))[::-1]       # monotone decreasing in rho
    return PchipInterpolator(rho, W, extrapolate=False)


def dP_local(profile, Rn, a, sp, rho_min_reached, nz=40001):
    if Rn >= a:
        return 0.0, True
    z = np.linspace(0.0, profile.zt + 8.0 * profile.L, nz)[1:]
    Rw = profile.Rw(z, Rn)
    ok = (Rw <= a) & (Rw >= max(Rn, rho_min_reached) - 1e-9)
    q = np.abs(sp.derivative()(Rw)) * np.abs(profile.dRw(z, Rn)) / (np.pi * Rw ** 2)
    q = np.where(ok & np.isfinite(q), q, 0.0)
    return float(q.max()), bool(rho_min_reached <= Rn + 1e-9)


def dP_work(Rn, a, sp, rho_min_reached, l_factor=1.0):
    if Rn >= a:
        return 0.0, True
    r = max(Rn, rho_min_reached)
    W = float(sp(r)) if np.isfinite(sp(r)) else 0.0
    return W / (np.pi * Rn ** 2 * l_factor * a), bool(rho_min_reached <= Rn + 1e-9)


def dP_drop(Rn, a, gamma_mN_m):
    if Rn >= a or gamma_mN_m <= 0:
        return 0.0
    return 2.0 * gamma_mN_m * (1.0 / Rn - 1.0 / a)      # mN/m / um = kPa
