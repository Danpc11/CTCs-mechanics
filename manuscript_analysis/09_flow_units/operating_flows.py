"""Operating windows of the population run expressed as flow rates (uL/min).

Membranes: pores are hydraulically in parallel, so a cell blocking one pore feels the pressure drop of the whole
membrane; flow per open pore Q = dP / R_pore, with R_pore = 8 mu L/(pi r^4) + 3 mu / r^3 (Poiseuille + Sampson
entrance correction for short pores). Flow per cm^2 = pore density x Q (all pores open; clogging raises dP at fixed flow).
Channel constriction: R = integral of 8 mu / (pi R_w(z)^4) dz over the wall profile (axisymmetric throat).
Window: lower edge = differential-weighted 99th percentile of leucocyte critical pressure (leucocytes pass above it);
upper edge = pressure that still retains the target fraction of tumour cells. Quasi-static: at higher velocity the
viscoelastic term raises the thresholds (Section 3.5 of the manuscript).
Input: results/full/analysis/long_table.csv of the population run. Output: ../datos/operating_flows.csv
"""
import os, sys, numpy as np, pandas as pd, yaml
H = os.path.dirname(os.path.abspath(__file__)); R0 = os.path.abspath(os.path.join(H, "..", ".."))
sys.path.insert(0, os.path.join(R0, "src"))
from ctcfem.geometry import profile_for_device
MU = 1.0e-3                                  # Pa s, buffer or diluted blood (scale flows by 1/mu)
LT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(R0, "results", "full", "analysis", "long_table.csv")
mats = yaml.safe_load(open(os.path.join(R0, "data", "materials.yaml")))
dif = yaml.safe_load(open(os.path.join(R0, "data", "populations.yaml")))["differential"]
T = pd.read_csv(LT, usecols=["device", "population", "group", "D_aperture", "contact", "dP"])
UL_MIN = 1e9 * 60                            # m^3/s -> uL/min


def resistance(dev, D):
    d = mats["devices"][dev]; r = 0.5 * D * 1e-6
    if d["type"] == "membrane":
        L = d["thickness_um"]["nominal"] * 1e-6
        n = d.get("pore_density_cm2", {}).get("nominal"); return 8 * MU * L / (np.pi * r ** 4) + 3 * MU / r ** 3, (float(n) if n is not None else None)
    p = profile_for_device(d); z = np.linspace(-(p.zt + 6 * p.L), p.zt + 6 * p.L, 20001)
    Rw = p.Rw(np.abs(z), 0.5 * D) * 1e-6
    return float(np.trapezoid(8 * MU / (np.pi * Rw ** 4), z * 1e-6)), None


rows = []
for (dev, D), s in T[T.D_aperture.isin([8.0, 10.0, 12.0])].groupby(["device", "D_aperture"]):
    leu = s[s.group == "leucocyte"].copy()
    leu["w"] = leu.population.map(dif) / leu.groupby("population").population.transform("count")
    o = np.argsort(leu.dP.values); cw = np.cumsum(leu.w.values[o]) / leu.w.sum()
    P_low = float(leu.dP.values[o][np.searchsorted(cw, 0.99)])
    Rh, n = resistance(dev, D)
    for pop, target in [("MDA-MB-231", 0.99), ("MCF-7", 0.99), ("patient_breast_CTC", 0.50)]:
        st = s[s.population == pop]
        if target > st.contact.mean():
            continue
        P_high = float(np.quantile(st.dP, 1 - target))
        q = lambda P: P * 1e3 / Rh * UL_MIN            # kPa -> Pa
        rows.append(dict(device=dev, D_um=D, tumour=pop, target_yield=target, P_low_kPa=P_low, P_high_kPa=P_high,
                         window_exists=P_high > P_low, Q_low_per_unit_uL_min=q(P_low), Q_high_per_unit_uL_min=q(P_high),
                         unit="pore" if n else "constriction",
                         Q_low_per_cm2_uL_min=q(P_low) * n if n else np.nan, Q_high_per_cm2_uL_min=q(P_high) * n if n else np.nan,
                         pore_density_cm2=n))
out = pd.DataFrame(rows); out.to_csv(os.path.join(H, "..", "datos", "operating_flows.csv"), index=False)
pd.set_option("display.width", 220)
print(out.round(5).to_string(index=False))
