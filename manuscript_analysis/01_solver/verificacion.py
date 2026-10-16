import os, sys; sys.path.append(os.getcwd())
import numpy as np, json
from scipy.optimize import brentq
from scipy.integrate import quad
from solver_axisimetrico import *
from constriccion import build, metrics, a as a_cel
OUT = {}

# ---------- V1 patch test ----------
from skfem import MeshTri
R0,H = 1.0,2.0
m = MeshTri.init_tensor(np.linspace(0,R0,5), np.linspace(0,H,7))
mat = NeoHookean(C10=0.5, kappa=1e4)
mdl = AxiModel(m, {"default": mat})
mdl.add_dirichlet(0, lambda x: x[0]<1e-12, 0.0)
mdl.add_dirichlet(1, lambda x: x[1]<1e-12, 0.0)
lam_t=1.3
mdl.add_dirichlet(1, lambda x: np.abs(x[1]-H)<1e-12, (lam_t-1.0)*H)
reps = mdl.solve([0.25,0.5,0.75,1.0])
def Prr(lr,lz):
    F=np.zeros((3,3,1,1)); F[0,0]=lr; F[1,1]=lz; F[2,2]=lr
    p=1e4*(lr*lr*lz-1.0)
    return piola(F,np.full((1,1),p),np.full((1,1),0.5))
lr=brentq(lambda l: Prr(l,lam_t)[0,0,0,0],0.3,1.5)
nod=mdl.basis.nodal_dofs; ud=mdl.x[nod]
lr_fe=float(np.mean(1.0+ud[0,np.abs(m.p[0]-R0)<1e-12]/R0))
OUT["V1_patch"]=dict(lam_axial=lam_t, lam_r_exacto=lr, lam_r_fe=lr_fe,
    err_rel=abs(lr_fe-lr)/abs(lr-1), J_fe=mdl.mean_J(), J_ex=lr*lr*lam_t,
    p_fe=mdl.mean_pressure(), p_ex=1e4*(lr*lr*lam_t-1), 
    res=[r.residuals for r in reps], its=[r.iters for r in reps])
print("V1", OUT["V1_patch"]["err_rel"], flush=True)

# ---------- V2 Hertz ----------
R,E,nu = 1.0,1.0,0.3; Es=E/(1-nu**2); mat=NeoHookean.from_E_nu(E,nu)
hz=[]
for delta in [0.0025,0.005,0.01,0.02,0.04]:
    a_h=np.sqrt(R*delta); F_h=4/3*Es*np.sqrt(R)*delta**1.5; hp=max(a_h/10,0.004)
    mm=mesh_quarter_disk(R=R,h=0.08,h_pole=hp,pole_radius=max(3*a_h,0.12))
    md=AxiModel(mm,{"default":mat})
    md.add_dirichlet(0,lambda x: np.abs(x[0])<1e-10,0.0)
    md.add_dirichlet(1,lambda x: np.abs(x[1])<1e-10,0.0)
    top=boundary_facets_where(mm,lambda x: np.hypot(x[0],x[1])>R-1e-8)
    md.add_wall(PlaneWall(z0=R,v=-delta,kn=100*E/hp),top)
    rp=md.solve(np.linspace(0.25,1.0,4)); F=md.contact_force_z()
    hz.append(dict(delta=delta,aR=a_h/R,ne=int(mm.t.shape[1]),conv=bool(all(r.converged for r in rp)),
                   F=F,F_hertz=F_h,err=(F-F_h)/F_h,a_fe=md.contact_radius(),a_hertz=a_h,
                   res=float(rp[-1].residuals[-1])))
    print("V2",delta,hz[-1]["err"],flush=True)
c=np.polyfit([d["aR"] for d in hz],[d["err"] for d in hz],1)
OUT["V2_hertz"]=dict(casos=hz, ajuste_err_vs_aR=list(c), err_extrapolado_aR0=float(c[1]))

# sensibilidad a la rigidez de penalizacion (delta/R=0.01)
delta=0.01; a_h=np.sqrt(R*delta); hp=a_h/10; F_h=4/3*Es*np.sqrt(R)*delta**1.5; pen=[]
for fac in [10,30,100,300,1000]:
    mm=mesh_quarter_disk(R=R,h=0.08,h_pole=hp,pole_radius=0.3)
    md=AxiModel(mm,{"default":mat})
    md.add_dirichlet(0,lambda x: np.abs(x[0])<1e-10,0.0)
    md.add_dirichlet(1,lambda x: np.abs(x[1])<1e-10,0.0)
    top=boundary_facets_where(mm,lambda x: np.hypot(x[0],x[1])>R-1e-8)
    md.add_wall(PlaneWall(z0=R,v=-delta,kn=fac*E/hp),top)
    rp=md.solve(np.linspace(0.25,1.0,4)); F=md.contact_force_z()
    pen.append(dict(fac=fac,kn=fac*E/hp,F=F,err=(F-F_h)/F_h,conv=bool(all(r.converged for r in rp))))
    print("PEN",fac,pen[-1]["err"],flush=True)
OUT["V2_penalizacion"]=pen

# ---------- V3 esfera hueca ----------
A,B,C1=0.6,1.0,0.5
def P_an(aa):
    c=aa**3-A**3; bb=(B**3+c)**(1/3)
    f=lambda r: 4*C1*((r/(r**3-c)**(1/3))**2-(r/(r**3-c)**(1/3))**-4)/r
    return quad(f,aa,bb,limit=200)[0], bb
a_t=0.9; P,b_t=P_an(a_t); hs=[]
def ann(nr,nt):
    mm=MeshTri.init_tensor(np.linspace(A,B,nr),np.linspace(0,np.pi/2,nt))
    rho,th=mm.p[0].copy(),mm.p[1].copy()
    return MeshTri(np.ascontiguousarray(np.vstack([rho*np.sin(th),rho*np.cos(th)])),mm.t)
for nr,nt in [(4,8),(6,12),(9,18),(13,26),(19,38)]:
    mm=ann(nr,nt); md=AxiModel(mm,{"default":NeoHookean(C10=C1,kappa=1e4*2*C1)})
    md.add_dirichlet(0,lambda x: np.abs(x[0])<1e-10,0.0)
    md.add_dirichlet(1,lambda x: np.abs(x[1])<1e-10,0.0)
    md.add_pressure(P, mm.facets_satisfying(lambda x: np.hypot(x[0],x[1])<A+0.02*(B-A)))
    rp=md.solve(np.linspace(0.1,1.0,10))
    nod=md.basis.nodal_dofs; rr=np.hypot(mm.p[0],mm.p[1]); s=np.abs(rr-A)<1e-8
    a_fe=float(np.mean(np.hypot(mm.p[0][s]+md.x[nod[0]][s], mm.p[1][s]+md.x[nod[1]][s])))
    hs.append(dict(n=[nr,nt],ndof=int(md.N),h=1.0/nt,conv=bool(all(r.converged for r in rp)),
                   a_fe=a_fe,a_ex=a_t,err=(a_fe-a_t)/a_t,err_desp=(a_fe-a_t)/(a_t-A),J=md.mean_J(),
                   res=float(rp[-1].residuals[-1])))
    print("V3",nr,hs[-1]["err"],flush=True)
OUT["V3_esfera_hueca"]=dict(P_analitico=P, a_objetivo=a_t, b=b_t, casos=hs)

# ---------- V4 Laplace ----------
Rs,C10L,kapL,gam=1.0,0.2,1.0,0.02
lam_ex=brentq(lambda l: kapL*(l**3-1.0)+2*gam/(l*Rs),0.5,1.2); p_ex=kapL*(lam_ex**3-1.0)
lp=[]
for h in [0.20,0.12,0.07,0.045]:
    mm=mesh_quarter_disk(R=Rs,h=h); md=AxiModel(mm,{"default":NeoHookean(C10=C10L,kappa=kapL)})
    md.add_dirichlet(0,lambda x: np.abs(x[0])<1e-10,0.0)
    md.add_dirichlet(1,lambda x: np.abs(x[1])<1e-10,0.0)
    md.add_surface(SurfaceLayer(gamma=gam), boundary_facets_where(mm,lambda x: np.hypot(x[0],x[1])>Rs-1e-8))
    rp=md.solve([0.3,0.6,1.0])
    nod=md.basis.nodal_dofs; s=np.hypot(mm.p[0],mm.p[1])>Rs-1e-10
    a_fe=float(np.mean(np.hypot(mm.p[0][s]+md.x[nod[0]][s],mm.p[1][s]+md.x[nod[1]][s])))
    lp.append(dict(h=h,ne=int(mm.t.shape[1]),conv=bool(all(r.converged for r in rp)),lam=a_fe/Rs,
                   err_lam=(a_fe/Rs-lam_ex)/(1-lam_ex),p=md.mean_pressure(),err_p=(md.mean_pressure()-p_ex)/abs(p_ex)))
    print("V4",h,lp[-1]["err_p"],flush=True)
OUT["V4_laplace"]=dict(lam_exacto=lam_ex,p_exacto=p_ex,casos=lp)

# ---------- V5 esfera casi incompresible ----------
OUT["V5"]=json.load(open("v5.json"))
# ---------- constriccion ----------
OUT["constriccion"]=json.load(open("constriccion.json"))

# ---------- sensibilidad al modulo de cizalla sustituto del citoplasma ----
import constriccion as CC
sens=[]
for fac in [0.1,0.3,1.0,3.0,10.0]:
    mu=0.0046*fac
    CC.cyt=NeoHookean(C10=0.5*mu,kappa=1e4*mu); CC.nuc=NeoHookean(C10=0.5*mu,kappa=1e4*mu)
    mm,md,outer,_=CC.build(0.85,0.7,5.0)
    rp=md.solve(np.linspace(0.08,1.0,14)); mt=CC.metrics(mm,md,outer)
    sens.append(dict(factor=fac,mu_cyt=mu,conv=bool(all(r.converged for r in rp)),**mt))
    print("SENS",fac,mt["DI"],mt["dA"],flush=True)
OUT["sensibilidad_mu_cyt"]=sens
json.dump(OUT, open("validacion.json","w"), indent=1)
print("LISTO")
