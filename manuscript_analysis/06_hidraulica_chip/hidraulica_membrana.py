"""
Red hidraulica del chip de membrana de dos canales separados por membrana de policarbonato.

Resuelve la presion transmembrana en cada cruce del canal superior con el serpentin inferior,
tratando el canal superior como una escalera 1D con derivaciones a traves de la membrana y el
serpentin como presion uniforme (su resistencia es 12 veces menor que la de la membrana).

Geometria y supuestos DECLARADOS (sustituir por los del CAD real):
  - canal superior ("rojo"): w = 150 um (estimado del render), h = 200 um, L = 10 mm
  - serpentin inferior ("naranja"): w = 700 um (estimado), h = 600 um, L = 55 mm
  - 12 cruces entre ambos canales
  - membrana de policarbonato track-etched: espesor 10 um, poro 8 um, 1e5 poros/cm2
  - viscosidades: tampon 1,0 / sangre 1:10 1,3 / sangre entera 3,5 mPa.s

Resultado principal: la dispersion de presion transmembrana entre los 12 cruces consume casi
entera la ventana de selectividad del poro de 8 um (2,67x), dejando un margen de caudal de solo
1,02x con el canal de 150 um. Ensanchar a 400-500 um lo sube a 2,19-2,30x.

Verificacion incluida: conservacion de caudal exacta (Q_in = Q_membrana + Q_salida).
"""
import numpy as np
from scipy.optimize import brentq

ETA = {"tampon": 1.0e-3, "sangre_1a10": 1.3e-3, "sangre_entera": 3.5e-3}
W_RED, H_RED, L_RED = 150e-6, 200e-6, 10e-3
W_OR, H_OR, L_OR = 700e-6, 600e-6, 55e-3
N_CRUCES = 12
T_MEM, R_PORO = 10e-6, 4e-6
DENS_POROS_CM2 = 1e5
VENTANA_8UM = (0.406, 1.084)          # kPa, de ventana_microfiltro.csv

A_OVERLAP = N_CRUCES * W_RED * W_OR
N_POROS = DENS_POROS_CM2 * (A_OVERLAP * 1e4)


def R_rect(w, h, L, eta):
    """Resistencia de conducto rectangular (Bruus); a = dimension menor."""
    a, b = min(w, h), max(w, h)
    return 12 * eta * L / (b * a**3 * (1 - 0.63 * a / b))


def R_poro(r, t, eta):
    """Poro cilindrico: Hagen-Poiseuille + entrada/salida de Sampson."""
    return 8 * eta * t / (np.pi * r**4) + 3 * eta / r**3


def red_restr(Q_in, eta, R_restr, N=N_CRUCES, w=W_RED, h=H_RED, p_serp=0.0):
    """Escalera 1D con restrictor a la salida. Nodo 0 = fuente de caudal."""
    Rseg = R_rect(w, h, L_RED / (N + 1), eta)
    Rm1 = R_poro(R_PORO, T_MEM, eta) / (N_POROS / N)
    A = np.zeros((N, N)); b = np.zeros(N)
    for i in range(N):
        A[i, i] += 1.0 / Rm1; b[i] += p_serp / Rm1
        if i > 0:   A[i, i] += 1.0 / Rseg; A[i, i - 1] -= 1.0 / Rseg
        if i < N-1: A[i, i] += 1.0 / Rseg; A[i, i + 1] -= 1.0 / Rseg
    b[0] += Q_in
    A[N-1, N-1] += 1.0 / (Rseg + R_restr)        # tramo + restrictor, en serie a atmosfera
    p = np.linalg.solve(A, b)
    dP = p - p_serp
    return p, dP, dP / Rm1, p[-1] / (Rseg + R_restr)


def ventana_caudal(eta_nm="tampon", mult=30, w=W_RED):
    """Caudal (uL/min) que mete los 12 cruces dentro de la ventana de retencion."""
    eta = ETA[eta_nm]
    Rm_tot = R_poro(R_PORO, T_MEM, eta) / N_POROS
    def f(Q, cual):
        _, dP, _, _ = red_restr(Q * 1e-9 / 60, eta, mult * Rm_tot, w=w)
        d = dP / 1e3
        return (d.min() - VENTANA_8UM[0]) if cual == "lo" else (d.max() - VENTANA_8UM[1])
    return brentq(lambda Q: f(Q, "lo"), 20, 3000), brentq(lambda Q: f(Q, "hi"), 20, 3000)


if __name__ == "__main__":
    print(f"area funcional de membrana: {A_OVERLAP*1e6:.3f} mm2 | poros: {N_POROS:.0f}")
    print(f"fraccion de area abierta  : {N_POROS*np.pi*R_PORO**2/A_OVERLAP*100:.2f} %\n")
    for nm, eta in ETA.items():
        Rr, Ro = R_rect(W_RED, H_RED, L_RED, eta), R_rect(W_OR, H_OR, L_OR, eta)
        Rm = R_poro(R_PORO, T_MEM, eta) / N_POROS
        print(f"{nm:<14} rojo/membrana {Rr/Rm:.2f}x | membrana/serpentin {Rm/Ro:.2f}x "
              f"| rojo/serpentin {Rr/Ro:.2f}x")
    print()
    # conservacion de caudal
    _, _, Qm, Qo = red_restr(10e-9/60, ETA["tampon"], 0.0)
    print(f"conservacion: Q_in 10,000 = Q_mem {Qm.sum()*60e9:.3f} + Q_out {Qo*60e9:.3f} uL/min\n")
    for nm in ETA:
        lo, hi = ventana_caudal(nm)
        print(f"{nm:<14} ventana de caudal: {lo:.0f}-{hi:.0f} uL/min (margen {hi/lo:.2f}x)")
    print()
    for w_um in [150, 250, 400, 500, 700, 1000]:
        lo, hi = ventana_caudal("tampon", w=w_um*1e-6)
        print(f"w_rojo {w_um:>5} um -> margen {hi/lo:.2f}x  ({lo:.0f}-{hi:.0f} uL/min)")
