"""Perfiles de pared con apertura EXACTA y sin quiebre.

    R_w(z) = R_n + (R_ext - R_n) h(z)
    h(z)   = (g(z) - g(0)) / (1 - g(0))
    g(z)   = 1 + [tanh((z - z_t)/L) - tanh((z + z_t)/L)] / 2

h es par, analitica (derivada nula en z = 0), h(0) = 0 y h(inf) = 1, de modo que
el radio en el cuello es EXACTAMENTE R_n: el diametro de apertura que se declara
es el que ve la celula.  Esto elimina la distincion nominal/efectivo del
manuscrito original y el quiebre del perfil tanh(|z|).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


def _g(z, zt, L):
    return 1.0 + 0.5 * (np.tanh((z - zt) / L) - np.tanh((z + zt) / L))


@dataclass(frozen=True)
class Perfil:
    zt: float          # semilongitud de garganta (canal) o semiespesor (membrana), um
    L: float           # longitud de transicion (canal) o radio de borde (membrana), um
    R_ext: float       # radio lejos de la garganta, um

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

    def angulo_entrada_deg(self, Rn):
        z = np.linspace(0.0, self.zt + 8.0 * self.L, 20001)
        return float(np.degrees(np.arctan(np.abs(self.dRw(z, Rn)).max())))

    def longitud_contacto(self, a, Rn):
        """Semilongitud axial con R_w < a (zona de contacto potencial)."""
        z = np.linspace(0.0, self.zt + 8.0 * self.L, 20001)
        i = np.where(self.Rw(z, Rn) < a)[0]
        return float(z[i.max()]) if len(i) else 0.0


def perfil_de_dispositivo(disp: dict) -> Perfil:
    """Construye el perfil desde una entrada de materiales.yaml:dispositivos."""
    if disp["tipo"] == "membrana":
        return Perfil(zt=0.5 * disp["espesor_um"]["nominal"],
                      L=disp["radio_borde_um"]["nominal"], R_ext=disp["R_ext_um"])
    if disp["tipo"] == "canal":
        return Perfil(zt=disp["semilongitud_garganta_um"]["nominal"],
                      L=disp["transicion_um"]["nominal"], R_ext=disp["R_ext_um"])
    raise ValueError(disp["tipo"])
