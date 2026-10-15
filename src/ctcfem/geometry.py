"""Wall profiles with an EXACT aperture and no corner.

    R_w(z) = R_n + (R_ext - R_n) h(z),   h(z) = (g(z) - g(0)) / (1 - g(0)),
    g(z)   = 1 + [tanh((z - z_t)/L) - tanh((z + z_t)/L)] / 2

h is even and analytic (zero slope at z = 0), h(0) = 0 and h(inf) = 1, so the neck radius is exactly
R_n: the declared aperture diameter is the one the cell meets.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


def _g(z, zt, L):
    return 1.0 + 0.5 * (np.tanh((z - zt) / L) - np.tanh((z + zt) / L))


@dataclass(frozen=True)
class Profile:
    zt: float      # throat half-length (channel) or half-thickness (membrane), um
    L: float       # transition length (channel) or edge radius (membrane), um
    R_ext: float   # radius far from the throat, um

    def h(self, z):
        g0 = _g(0.0, self.zt, self.L)
        return (_g(z, self.zt, self.L) - g0) / (1.0 - g0)

    def dh(self, z):
        g0 = _g(0.0, self.zt, self.L)
        s1 = 1.0 / np.cosh((z - self.zt) / self.L) ** 2
        s2 = 1.0 / np.cosh((z + self.zt) / self.L) ** 2
        return 0.5 * (s1 - s2) / self.L / (1.0 - g0)

    def Rw(self, z, Rn):
        return Rn + (self.R_ext - Rn) * self.h(z)

    def dRw(self, z, Rn):
        return (self.R_ext - Rn) * self.dh(z)

    def entry_angle_deg(self, Rn):
        z = np.linspace(0.0, self.zt + 8.0 * self.L, 20001)
        return float(np.degrees(np.arctan(np.abs(self.dRw(z, Rn)).max())))

    def contact_half_length(self, a, Rn):
        """Axial half-length over which R_w < a (potential contact zone)."""
        z = np.linspace(0.0, self.zt + 8.0 * self.L, 20001)
        i = np.where(self.Rw(z, Rn) < a)[0]
        return float(z[i.max()]) if len(i) else 0.0


def profile_for_device(dev: dict) -> Profile:
    if dev["type"] == "membrane":
        return Profile(zt=0.5 * dev["thickness_um"]["nominal"], L=dev["edge_radius_um"]["nominal"], R_ext=dev["R_ext_um"])
    if dev["type"] == "channel":
        return Profile(zt=dev["throat_half_length_um"]["nominal"], L=dev["transition_um"]["nominal"], R_ext=dev["R_ext_um"])
    raise ValueError(dev["type"])
