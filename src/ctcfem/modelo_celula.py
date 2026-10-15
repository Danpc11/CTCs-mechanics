"""Celula multicompartimento en un dispositivo: una continuacion que cierra la
apertura desde el primer contacto hasta la apertura mas pequena pedida y
registra la energia en CADA paso aceptado.  Una sola corrida da W(rho) para
todas las aperturas.

Unidades internas: um, mN/m (Gs, gamma), kPa (G_cyt, presiones); W en mN/m*um^2.
"""
from __future__ import annotations
import time
import numpy as np

from .solver_axisimetrico import (NeoHookean, SurfaceLayer, ProfileWall, AxiModel,
                                  mesh_cell_quarter, boundary_facets_where,
                                  interface_facets_where, axi_F, det3)
from .geometria import Perfil

KAPPA_MU = 1.0e3          # nu = 0.4995 (manuscrito; J medio = 1.0000000 verificado)
CONF_EF_VALIDADO = 1.636  # techo validado (confinamiento efectivo)


def _materiales(G_cyt, r_nuc):
    cyt = NeoHookean(name="citoplasma", C10=0.5 * G_cyt, kappa=KAPPA_MU * G_cyt)
    gn = r_nuc * G_cyt
    nuc = NeoHookean(name="nucleo", C10=0.5 * gn, kappa=KAPPA_MU * gn)
    return cyt, nuc


def energias(mdl, capas):
    b, pb = mdl.basis, mdl.pbasis
    R = b.global_coordinates()[0]
    uq = b.interpolate(mdl.u)
    pq = np.asarray(pb.interpolate(mdl.p))
    F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], R)
    J = det3(F)
    I1 = np.einsum("ij...,ij...->...", F, F)
    Psi = (mdl.C10 * (J ** (-2.0 / 3.0) * I1 - 3.0) + pq * (J - 1.0) - pq ** 2 / (2.0 * mdl.kappa))
    wq = 2.0 * np.pi * R * b.dx
    out = {f"W_{nm}": float(np.sum(Psi[els] * wq[els])) for nm, els in mdl.subdomains.items()}
    for capa, (_, fb) in zip(capas, mdl.surfaces):
        X = np.asarray(fb.global_coordinates()); Rref = X[0]
        n = np.asarray(fb.normals); t = np.stack([-n[1], n[0]])
        uqf = fb.interpolate(mdl.u); gu = np.asarray(uqf.grad)
        a0 = t[0] + gu[0, 0] * t[0] + gu[0, 1] * t[1]
        a1 = t[1] + gu[1, 0] * t[0] + gu[1, 1] * t[1]
        ls = np.sqrt(a0 ** 2 + a1 ** 2); lt = (Rref + np.asarray(uqf)[0]) / Rref
        out[f"W_{capa.name}"] = float(np.sum(capa.w(ls, lt) * 2.0 * np.pi * Rref * fb.dx))
        out[f"ls_max_{capa.name}"] = float(ls.max())
    out["W_total"] = sum(v for k, v in out.items() if k.startswith("W_"))
    return out


def area_relativa(mdl, outer):
    fb = mdl.facet_basis(outer)
    X = np.asarray(fb.global_coordinates()); n = np.asarray(fb.normals)
    t = np.stack([-n[1], n[0]]); uq = fb.interpolate(mdl.u); gu = np.asarray(uq.grad)
    a0 = t[0] + gu[0, 0] * t[0] + gu[0, 1] * t[1]
    a1 = t[1] + gu[1, 0] * t[0] + gu[1, 1] * t[1]
    ls = np.sqrt(a0 ** 2 + a1 ** 2); lt = (X[0] + np.asarray(uq)[0]) / X[0]
    A0 = float(np.sum(2.0 * np.pi * X[0] * fb.dx))
    return float(np.sum(ls * lt * 2.0 * np.pi * X[0] * fb.dx)) / A0 - 1.0


def rigidez_pared(Gs, a, h, kn_f, E_pared_kPa=None, nu_pared=0.4, l_pared=None):
    """Penalizacion numerica en serie con la rigidez de Winkler de la pared."""
    k_num = kn_f * (Gs / a) / h
    if not E_pared_kPa:
        return k_num, k_num, np.inf
    l_w = l_pared if l_pared else a
    k_w = E_pared_kPa / ((1.0 - nu_pared ** 2) * l_w)
    return 1.0 / (1.0 / k_num + 1.0 / k_w), k_num, k_w


def construir(c: dict, perfil: Perfil, rho_ini, rho_fin, hrel=10.0, kn_f=200.0,
              E_pared_kPa=None, nu_pared=0.4, l_pared=None, usar_gamma=False):
    a, an = c["a"], c["an"]
    h, hn = a / hrel, 0.7 * a / hrel
    cyt, nuc = _materiales(c["G_cyt"], c["r_nuc"])
    m, sub = mesh_cell_quarter(Rc=a, Rn=an, h=h, hn=hn)
    mdl = AxiModel(m, {"citoplasma": cyt, "nucleo": nuc}, subdomains=sub)
    mdl.add_dirichlet(0, lambda x: np.abs(x[0]) < 1e-9, 0.0)
    mdl.add_dirichlet(1, lambda x: np.abs(x[1]) < 1e-9, 0.0)
    outer = boundary_facets_where(m, lambda x: np.hypot(x[0], x[1]) > a - 1e-7)
    iface = interface_facets_where(m, lambda x: np.abs(np.hypot(x[0], x[1]) - an) < 1e-7)
    cortex = SurfaceLayer(name="cortex", Gs=c["Gs"], C=c["C"], law="skalak",
                          gamma=(c["gamma"] if usar_gamma else 0.0))
    envol = SurfaceLayer(name="envoltura", Gs=c["r_env"] * c["Gs"], C=c["C"], law="skalak")
    mdl.add_surface(cortex, outer)
    mdl.add_surface(envol, iface)
    k_eff, k_num, k_w = rigidez_pared(c["Gs"], a, h, kn_f, E_pared_kPa, nu_pared, l_pared)

    def Rw(z, lam):
        return perfil.Rw(z, rho_ini + lam * (rho_fin - rho_ini))

    mdl.add_wall(ProfileWall(Rw=Rw, kn=k_eff), outer)
    return m, mdl, outer, (cortex, envol), dict(k_eff=k_eff, k_num=k_num, k_pared=k_w, h=h)


def cerrar_apertura(c: dict, perfil: Perfil, radios_objetivo, conf_max_intento=2.4,
                    dlam0=0.05, dlam_min=2e-4, dlam_max=0.08, maxit=30, max_intentos=400,
                    **kw):
    """Continuacion desde rho = 1.02 a hasta el radio objetivo mas pequeno.

    Devuelve dict con 'curva' (rho, energias en cada paso aceptado), 'checkpoints'
    (metricas completas en cada radio objetivo alcanzado) y diagnostico.
    """
    t0 = time.time()
    a = c["a"]
    objetivos = sorted([r for r in radios_objetivo if r < a and a / r <= conf_max_intento],
                       reverse=True)
    res = dict(a=a, objetivos=objetivos, curva=[], checkpoints=[], completo=False,
               rho_min_alcanzado=a, motivo_parada="sin_contacto")
    if not objetivos:
        res["completo"] = True
        return res
    rho0, rho_end = 1.02 * a, objetivos[-1]
    m, mdl, outer, capas, kinfo = construir(c, perfil, rho0, rho_end, **kw)
    res.update(ne=int(m.t.shape[1]), ndof=int(mdl.N), **kinfo)
    lam_de = lambda r: (rho0 - r) / (rho0 - rho_end)
    lam_cp = [lam_de(r) for r in objetivos]
    lam, dlam, x_ok, k, intentos = 0.0, dlam0, mdl.x.copy(), 0, 0
    res["curva"].append(dict(rho=rho0, W_total=0.0))
    while k < len(lam_cp):
        intentos += 1
        if intentos > max_intentos:
            res["motivo_parada"] = "max_intentos"; break
        trial = min(lam + dlam, lam_cp[k])
        mdl.x = x_ok.copy()
        rep = mdl.newton(trial, maxit=maxit)
        if rep.converged and np.all(np.isfinite(mdl.x)):
            lam, x_ok = trial, mdl.x.copy()
            mdl.x_prev = x_ok.copy()
            rho = rho0 + lam * (rho_end - rho0)
            e = energias(mdl, capas)
            res["curva"].append(dict(rho=rho, **{q: e[q] for q in e if q.startswith("W_")}))
            res["rho_min_alcanzado"] = rho
            if rep.iters <= 6:
                dlam = min(dlam * 1.5, dlam_max)
            if abs(lam - lam_cp[k]) < 1e-12:
                cp = dict(rho=rho, D_apertura=2 * rho, conf_ef=a / rho,
                          dentro_envolvente=bool(a / rho <= CONF_EF_VALIDADO),
                          dA_A0=area_relativa(mdl, outer),
                          Fz_contacto=float(mdl.contact_force_z()), iters=int(rep.iters), **e)
                res["checkpoints"].append(cp)
                k += 1
        else:
            dlam *= 0.4
            if dlam < dlam_min:
                res["motivo_parada"] = "no_converge"; break
    else:
        res["completo"] = True
        res["motivo_parada"] = "ok"
    res["segundos"] = round(time.time() - t0, 2)
    return res
