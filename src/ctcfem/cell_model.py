"""A multi-compartment cell held at the neck of a device: one continuation closes the aperture from first
contact to the smallest requested aperture and records the energy at EVERY accepted step, so a single run
gives W(rho) for all apertures.

Internal units: um, mN/m (Gs, gamma), kPa (G_cyt, pressures); W in mN/m um^2.
"""
from __future__ import annotations
import time
import numpy as np

from .solver_axisimetrico import (NeoHookean, SurfaceLayer, ProfileWall, AxiModel, mesh_cell_quarter,
                                  boundary_facets_where, interface_facets_where, axi_F, det3)
from .geometry import Profile

KAPPA_MU = 1.0e5            # bulk/shear ratio of cytoplasm and nucleoplasm; 1e3 lets the cortex compress the cell by ~6% at a/R = 1.636 (tests/test_volume.py)
CONF_EFF_ASSESSED = 1.636   # upper limit of the assessed range (effective confinement a/R)


def _materials(G_cyt, r_nuc):
    cyt = NeoHookean(name="cytoplasm", C10=0.5 * G_cyt, kappa=KAPPA_MU * G_cyt)
    gn = r_nuc * G_cyt
    return cyt, NeoHookean(name="nucleus", C10=0.5 * gn, kappa=KAPPA_MU * gn)


def energies(mdl, layers):
    b, pb = mdl.basis, mdl.pbasis
    R = b.global_coordinates()[0]
    uq = b.interpolate(mdl.u); pq = np.asarray(pb.interpolate(mdl.p))
    F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], R); J = det3(F)
    I1 = np.einsum("ij...,ij...->...", F, F)
    Psi = mdl.C10 * (J ** (-2.0 / 3.0) * I1 - 3.0) + pq * (J - 1.0) - pq ** 2 / (2.0 * mdl.kappa)
    wq = 2.0 * np.pi * R * b.dx
    out = {f"W_{nm}": float(np.sum(Psi[els] * wq[els])) for nm, els in mdl.subdomains.items()}
    for layer, (_, fb) in zip(layers, mdl.surfaces):
        X = np.asarray(fb.global_coordinates()); Rref = X[0]
        n = np.asarray(fb.normals); t = np.stack([-n[1], n[0]])
        uqf = fb.interpolate(mdl.u); gu = np.asarray(uqf.grad)
        a0 = t[0] + gu[0, 0] * t[0] + gu[0, 1] * t[1]; a1 = t[1] + gu[1, 0] * t[0] + gu[1, 1] * t[1]
        ls = np.sqrt(a0 ** 2 + a1 ** 2); lt = (Rref + np.asarray(uqf)[0]) / Rref
        out[f"W_{layer.name}"] = float(np.sum(layer.w(ls, lt) * 2.0 * np.pi * Rref * fb.dx))
        out[f"ls_max_{layer.name}"] = float(ls.max())
    out["W_total"] = sum(v for k, v in out.items() if k.startswith("W_"))
    return out


def relative_area(mdl, outer):
    fb = mdl.facet_basis(outer)
    X = np.asarray(fb.global_coordinates()); n = np.asarray(fb.normals); t = np.stack([-n[1], n[0]])
    uq = fb.interpolate(mdl.u); gu = np.asarray(uq.grad)
    a0 = t[0] + gu[0, 0] * t[0] + gu[0, 1] * t[1]; a1 = t[1] + gu[1, 0] * t[0] + gu[1, 1] * t[1]
    ls = np.sqrt(a0 ** 2 + a1 ** 2); lt = (X[0] + np.asarray(uq)[0]) / X[0]
    A0 = float(np.sum(2.0 * np.pi * X[0] * fb.dx))
    return float(np.sum(ls * lt * 2.0 * np.pi * X[0] * fb.dx)) / A0 - 1.0


def wall_stiffness(Gs, a, h, kn_f, E_wall_kPa=None, nu_wall=0.4, l_wall=None):
    """Numerical penalty in series with a Winkler wall stiffness E/((1 - nu^2) l), kPa/um."""
    k_num = kn_f * (Gs / a) / h
    if not E_wall_kPa:
        return k_num, k_num, np.inf
    k_w = E_wall_kPa / ((1.0 - nu_wall ** 2) * (l_wall if l_wall else a))
    return 1.0 / (1.0 / k_num + 1.0 / k_w), k_num, k_w


def build(c: dict, profile: Profile, rho_start, rho_end, hrel=10.0, kn_f=200.0, E_wall_kPa=None,
          nu_wall=0.4, l_wall=None, use_gamma=False):
    a, an = c["a"], c["an"]
    h, hn = a / hrel, 0.7 * a / hrel
    cyt, nuc = _materials(c["G_cyt"], c["r_nuc"])
    m, sub = mesh_cell_quarter(Rc=a, Rn=an, h=h, hn=hn)
    sub = {("cytoplasm" if k == "citoplasma" else "nucleus" if k == "nucleo" else k): v for k, v in sub.items()}
    mdl = AxiModel(m, {"cytoplasm": cyt, "nucleus": nuc}, subdomains=sub)
    mdl.add_dirichlet(0, lambda x: np.abs(x[0]) < 1e-9, 0.0)
    mdl.add_dirichlet(1, lambda x: np.abs(x[1]) < 1e-9, 0.0)
    outer = boundary_facets_where(m, lambda x: np.hypot(x[0], x[1]) > a - 1e-7)
    iface = interface_facets_where(m, lambda x: np.abs(np.hypot(x[0], x[1]) - an) < 1e-7)
    cortex = SurfaceLayer(name="cortex", Gs=c["Gs"], C=c["C"], law="skalak", gamma=(c["gamma"] if use_gamma else 0.0))
    envelope = SurfaceLayer(name="envelope", Gs=c["r_env"] * c["Gs"], C=c["C"], law="skalak")
    mdl.add_surface(cortex, outer); mdl.add_surface(envelope, iface)
    k_eff, k_num, k_w = wall_stiffness(c["Gs"], a, h, kn_f, E_wall_kPa, nu_wall, l_wall)
    mdl.add_wall(ProfileWall(Rw=lambda z, lam: profile.Rw(z, rho_start + lam * (rho_end - rho_start)), kn=k_eff), outer)
    return m, mdl, outer, (cortex, envelope), dict(k_eff=k_eff, k_num=k_num, k_wall=k_w, h=h)


def close_aperture(c: dict, profile: Profile, target_radii, conf_max_attempt=2.4, dlam0=0.05, dlam_min=2e-4,
                   dlam_max=0.08, maxit=30, max_attempts=400, on_checkpoint=None, **kw):
    """Continuation from rho = 1.02 a to the smallest target radius.

    Returns 'curve' (rho and energies at every accepted step), 'checkpoints' (full metrics at every target
    radius reached) and diagnostics. on_checkpoint(mesh, model, checkpoint) is called at each target.
    """
    t0 = time.time(); a = c["a"]
    targets = sorted([r for r in target_radii if r < a and a / r <= conf_max_attempt], reverse=True)
    res = dict(a=a, targets=targets, curve=[], checkpoints=[], complete=False, rho_min_reached=a, stop_reason="no_contact")
    if not targets:
        res["complete"] = True
        return res
    rho0, rho_end = 1.02 * a, targets[-1]
    m, mdl, outer, layers, kinfo = build(c, profile, rho0, rho_end, **kw)
    res.update(ne=int(m.t.shape[1]), ndof=int(mdl.N), **kinfo)
    lam_cp = [(rho0 - r) / (rho0 - rho_end) for r in targets]
    lam, dlam, x_ok, k, attempts = 0.0, dlam0, mdl.x.copy(), 0, 0
    res["curve"].append(dict(rho=rho0, W_total=0.0))
    while k < len(lam_cp):
        attempts += 1
        if attempts > max_attempts:
            res["stop_reason"] = "max_attempts"; break
        trial = min(lam + dlam, lam_cp[k])
        mdl.x = x_ok.copy()
        rep = mdl.newton(trial, maxit=maxit)
        if rep.converged and np.all(np.isfinite(mdl.x)):
            lam, x_ok = trial, mdl.x.copy(); mdl.x_prev = x_ok.copy()
            rho = rho0 + lam * (rho_end - rho0)
            e = energies(mdl, layers)
            res["curve"].append(dict(rho=rho, **{q: e[q] for q in e if q.startswith("W_")}))
            res["rho_min_reached"] = rho
            if rep.iters <= 6:
                dlam = min(dlam * 1.5, dlam_max)
            if abs(lam - lam_cp[k]) < 1e-12:
                cp = dict(rho=rho, D_aperture=2 * rho, conf_eff=a / rho, inside_assessed_range=bool(a / rho <= CONF_EFF_ASSESSED),
                          dA_A0=relative_area(mdl, outer), Fz_contact=float(mdl.contact_force_z()), iters=int(rep.iters), **e)
                res["checkpoints"].append(cp)
                if on_checkpoint:
                    on_checkpoint(m, mdl, cp)
                k += 1
        else:
            dlam *= 0.4
            if dlam < dlam_min:
                res["stop_reason"] = "not_converged"; break
    else:
        res["complete"] = True; res["stop_reason"] = "ok"
    res["seconds"] = round(time.time() - t0, 2)
    return res
