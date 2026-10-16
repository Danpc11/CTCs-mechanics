import os, sys; sys.path.append(os.getcwd())
import numpy as np, json
from solver_axisimetrico import *
# Unidades: longitud um, fuerza nN  ->  tension kPa, modulo de superficie mN/m

a, an = 8.5, 4.25                 # radios celula y nucleo (um)
Gs, Cs = 0.39, 10.0               # membrana HL-60 ajustada por Wang 2023
Gsn = 2.0*Gs                      # envoltura nuclear = 2x membrana
mu_cyt = 0.0046                   # kPa  (= 0.1 * Gs/a, interior casi fluido)
cyt = NeoHookean(C10=0.5*mu_cyt, kappa=1e4*mu_cyt)
nuc = NeoHookean(C10=0.5*mu_cyt, kappa=1e4*mu_cyt)

def build(h, hn, Rt, Rwide=1.6*8.5, zt=4.0, Ltap=5.0, kn=None):
    m, sub = mesh_cell_quarter(Rc=a, Rn=an, h=h, hn=hn)
    mdl = AxiModel(m, {"citoplasma": cyt, "nucleo": nuc}, subdomains=sub)
    mdl.add_dirichlet(0, lambda x: np.abs(x[0]) < 1e-9, 0.0)
    mdl.add_dirichlet(1, lambda x: np.abs(x[1]) < 1e-9, 0.0)
    outer = boundary_facets_where(m, lambda x: np.hypot(x[0],x[1]) > a-1e-7)
    iface = interface_facets_where(m, lambda x: np.abs(np.hypot(x[0],x[1])-an) < 1e-7)
    mdl.add_surface(SurfaceLayer(name="cortex", Gs=Gs, C=Cs, law="skalak"), outer)
    mdl.add_surface(SurfaceLayer(name="envoltura", Gs=Gsn, C=Cs, law="skalak"), iface)
    def Rw(z, lam):
        Rth = Rwide + lam*(Rt - Rwide)
        return Rth + (Rwide-Rth)*0.5*(1.0+np.tanh((np.abs(z)-zt)/Ltap))
    kn = kn if kn is not None else 50.0*Gs/h**2
    mdl.add_wall(ProfileWall(Rw=Rw, kn=kn), outer)
    return m, mdl, outer, iface

def metrics(m, mdl, outer):
    """DI = 1 - 2 sqrt(pi A)/P sobre el contorno meridiano del cuarto de
    celula reflejado a la celula completa; y deformacion de area de membrana."""
    nod = mdl.basis.nodal_dofs
    P = m.p + mdl.x[nod]
    v = m.facets[:, outer]
    verts = np.unique(v)
    rr, zz = P[0, verts], P[1, verts]
    ang = np.arctan2(zz, rr); o = np.argsort(ang)
    rr, zz = rr[o], zz[o]
    # SECCION MERIDIANA COMPLETA: reflejar en z=0 y despues en r=0, igual
    # que el plano de simetria y=0 de Wang et al. (seccion llena, no media)
    Rh = np.concatenate([rr[::-1], rr[1:]])        # de (0,-a) a (0,+a), r>=0
    Zh = np.concatenate([-zz[::-1], zz[1:]])
    Rfull = np.concatenate([Rh, -Rh[::-1][1:-1]])  # cerrar por r<0
    Zfull = np.concatenate([Zh,  Zh[::-1][1:-1]])
    Aproj = 0.5*np.abs(np.sum(Rfull*np.roll(Zfull,-1) - np.roll(Rfull,-1)*Zfull))
    seg = np.hypot(np.diff(np.append(Rfull,Rfull[0])),
                   np.diff(np.append(Zfull,Zfull[0]))).sum()
    DI = 1.0 - 2.0*np.sqrt(np.pi*Aproj)/seg
    # area de membrana (superficie de revolucion) deformada / referencia
    fb = mdl.facet_basis(outer); X = np.asarray(fb.global_coordinates())
    n = np.asarray(fb.normals); t = np.stack([-n[1], n[0]])
    uq = fb.interpolate(mdl.u); gu = np.asarray(uq.grad)
    a0 = t[0] + gu[0,0]*t[0] + gu[0,1]*t[1]; a1 = t[1] + gu[1,0]*t[0] + gu[1,1]*t[1]
    ls = np.sqrt(a0**2+a1**2); lt = (X[0]+np.asarray(uq)[0])/X[0]
    from skfem import Functional
    @Functional
    def Adef(w): return 2*np.pi*X[0]*ls*lt
    @Functional
    def A0(w): return 2*np.pi*X[0]
    return dict(DI=float(DI), dA=float(Adef.assemble(fb)/A0.assemble(fb)-1.0),
                ls_max=float(ls.max()), lt_max=float(lt.max()),
                rmax=float(P[0].max()), zmax=float(P[1].max()))

if __name__ == "__main__":
    res = {"convergencia_malla": [], "barrido_confinamiento": []}
    Rt = 6.0   # radio de garganta um -> confinamiento 2a/D = 8.5/6.0 = 1.417
    for h, hn in [(1.7,1.4),(1.2,1.0),(0.85,0.7),(0.6,0.5)]:
        m, mdl, outer, iface = build(h, hn, Rt)
        reps = mdl.solve(np.linspace(0.1,1.0,12))
        ok = all(r.converged for r in reps)
        mt = metrics(m, mdl, outer) if ok else {}
        rec = dict(h=h, ne=int(m.t.shape[1]), ndof=int(mdl.N), conv=ok,
                   Fc=float(mdl.contact_force_z()), res=float(reps[-1].residuals[-1]),
                   its=[r.iters for r in reps], **mt)
        res["convergencia_malla"].append(rec)
        print(f"MALLA h={h} ne={rec['ne']:5d} ndof={rec['ndof']:6d} conv={ok} "
              f"DI={mt.get('DI',float('nan')):.5f} dA={mt.get('dA',float('nan')):+.4f} "
              f"rmax={mt.get('rmax',float('nan')):.4f} res={rec['res']:.2e}", flush=True)
    for Rt in [8.6, 8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0, 4.6, 4.25, 3.8, 3.5]:
        m, mdl, outer, iface = build(0.85, 0.7, Rt)
        reps = mdl.solve(np.linspace(0.08,1.0,14))
        ok = all(r.converged for r in reps)
        mt = metrics(m, mdl, outer) if ok else {}
        rec = dict(Rt=Rt, conf=a/Rt, conv=ok, Fc=float(mdl.contact_force_z()), **mt)
        res["barrido_confinamiento"].append(rec)
        print(f"CONF Rt={Rt} 2a/D={a/Rt:.3f} conv={ok} DI={mt.get('DI',float('nan')):.5f} "
              f"dA={mt.get('dA',float('nan')):+.4f} ls={mt.get('ls_max',float('nan')):.3f}", flush=True)
    json.dump(res, open("constriccion.json","w"), indent=1)
