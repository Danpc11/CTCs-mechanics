"""
fase2_lib.py -- Fisica de las modalidades (Fase 2) sobre el solver
axisimetrico validado en la Fase 1.

UNIDADES (identicas a constriccion.py de la Fase 1)
    longitud  um        fuerza  nN        ->  tension kPa
    modulo de superficie  mN/m (= nN/um)
    energia   nN*um (= fJ)                    presion kPa
    tiempo    s         viscosidad de fluido nN*s/um^2 (= kPa*s)

ENVOLVENTE DE VALIDEZ (Fase 1, validacion_solver.md Sec. 5.4)
    confinamiento  c = 2a/D = a/R_garganta  debe cumplir  1.2 <= c <= 2.43.
    Por debajo de 1.2 no hay contacto y el modelo cuasi-estatico predice
    deformacion identicamente nula (subestimacion conocida).
    Por encima de 2.43 no se ha explorado.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from skfem import Functional

from solver_axisimetrico import (
    NeoHookean, SurfaceLayer, ProfileWall, PlaneWall, AxiModel,
    mesh_cell_quarter, boundary_facets_where, interface_facets_where,
    axi_F, det3,
)

# ---------------------------------------------------------------------
#  Envolvente de validez
# ---------------------------------------------------------------------
CONF_MIN = 1.2
CONF_MAX = 2.43


def dentro_envolvente(conf):
    c = np.asarray(conf, dtype=float)
    return (c >= CONF_MIN) & (c <= CONF_MAX)


# ---------------------------------------------------------------------
#  Estirpes: SOLO lo que esta en parametros_compartimento.csv
# ---------------------------------------------------------------------

@dataclass
class Estirpe:
    nombre: str
    a: float                      # radio celular en suspension (um)
    a_src: str                    # 'medido' | 'transferido'
    E_app: float | None           # modulo de Young aparente en suspension (kPa)
    E_app_src: str
    G_cyt: float                  # modulo elastico del citoplasma (kPa)
    G_cyt_src: str
    tau: float | None             # tiempo de relajacion citoplasmatico (s)
    tau_src: str
    eta_cyt: float | None         # viscosidad estatica del citoplasma (Pa s)
    # rangos publicados (min, max) para la propagacion de dispersion
    a_rng: tuple = (None, None)
    E_app_rng: tuple = (None, None)
    G_cyt_rng: tuple = (None, None)
    tau_rng: tuple = (None, None)
    # transferidos identicos para las cinco estirpes (auditoria Sec. 2)
    an_a: float = 0.50            # relacion de radios nucleo/celula
    C_skalak: float = 10.0        # Ks/Gs = 1 + 2C
    r_nuc: float = 1.4            # G_nucleoplasma / G_citoplasma
    r_env: float = 2.0            # Gs_envoltura / Gs_cortex
    Gs_lit: float | None = None   # Gs transferido de la literatura (mN/m)
    Gs_lit_src: str = ""
    Gs: float | None = None       # Gs CALIBRADO contra E_app (mN/m)

    @property
    def an(self):
        return self.an_a * self.a


# tau y G_cyt: microrreologia rotacional 10.1039/d4na00003j
# E_app: pinzas acusticas 10.1038/s41378-020-0150-6
# Gs_lit tumoral: PC-3 por ajuste IB-LBM 10.1007/s10237-024-01920-9
# Gs_lit leucocito: HL-60 por ajuste IB-LBM 10.1007/s10237-023-01705-6
ESTIRPES = {
    "MCF-10A": Estirpe(
        "MCF-10A", a=8.74, a_src="transferido (MCF-7, mismo panel)",
        E_app=None, E_app_src="NO MEDIDO",
        G_cyt=0.0793, G_cyt_src="medido", G_cyt_rng=(0.0720, 0.0866),
        tau=0.77, tau_src="medido", tau_rng=(0.40, 1.14), eta_cyt=41.6,
        a_rng=(8.51, 10.0), Gs_lit=13.25, Gs_lit_src="PC-3 transferido"),
    "MCF-7": Estirpe(
        "MCF-7", a=8.74, a_src="medido", a_rng=(8.59, 8.89),
        E_app=2.650, E_app_src="medido", E_app_rng=(2.650, 2.650),
        G_cyt=0.0329, G_cyt_src="medido", G_cyt_rng=(0.0269, 0.0389),
        tau=3.23, tau_src="medido", tau_rng=(2.70, 3.76), eta_cyt=56.4,
        Gs_lit=13.25, Gs_lit_src="PC-3 transferido"),
    "SK-BR-3": Estirpe(
        "SK-BR-3", a=10.0, a_src="transferido (media del panel)",
        a_rng=(8.6, 10.0),
        E_app=2.772, E_app_src="medido", E_app_rng=(2.772, 2.772),
        G_cyt=0.0793, G_cyt_src="transferido (MCF-10A)",
        G_cyt_rng=(0.0329, 0.0793),
        tau=None, tau_src="NO MEDIDO", tau_rng=(None, None), eta_cyt=None,
        Gs_lit=13.25, Gs_lit_src="PC-3 transferido"),
    "MDA-MB-231": Estirpe(
        "MDA-MB-231", a=8.62, a_src="medido", a_rng=(8.51, 8.73),
        E_app=1.527, E_app_src="medido", E_app_rng=(1.527, 1.527),
        G_cyt=0.0386, G_cyt_src="medido", G_cyt_rng=(0.0328, 0.0444),
        tau=0.57, tau_src="medido", tau_rng=(0.37, 0.77), eta_cyt=10.7,
        Gs_lit=13.25, Gs_lit_src="PC-3 transferido"),
    "leucocito": Estirpe(
        "leucocito", a=8.5, a_src="transferido (HL-60)", a_rng=(8.5, 8.5),
        E_app=0.620, E_app_src="medido (PBMC)", E_app_rng=(0.400, 0.660),
        G_cyt=0.0386, G_cyt_src="NO MEDIDO - sustituto declarado",
        G_cyt_rng=(0.0329, 0.0793),
        tau=None, tau_src="NO MEDIDO", tau_rng=(None, None), eta_cyt=None,
        Gs_lit=0.39, Gs_lit_src="HL-60 ajuste IB-LBM"),
}

TAU_DISPONIBLE = [k for k, v in ESTIRPES.items() if v.tau is not None]

# REGULARIZACION DEL CITOPLASMA (heredada de la Fase 1).
#   Wang et al. modelan citoplasma y nucleoplasma como fluidos viscosos SIN
#   elasticidad; un solver cuasi-estatico de solido necesita rigidez no nula.
#   La Fase 1 fijo  mu_cyt = 0.1 * Gs / a,  es decir  Pi = G_cyt a / Gs = 0.1,
#   y verifico que DI cambia solo +4.5% y dA/A0 +22.5% sobre 0.3x-10x en mu_cyt
#   (Pi de 0.03 a 1.0).  Se mantiene ese valor.
#   Los G_cyt MEDIDOS por microrreologia implican Pi = 0.024-0.203 (ver
#   resultados_fase2.md): son un modulo LOCAL dependiente de la escala de sonda
#   (el propio trabajo publica 15-160 Pa s al variar el hilo de 2 a 6 um), no el
#   modulo continuo de la gota.  Lo que SI es robusto a la escala es tau = eta/G,
#   cociente de dos medidas de la misma sonda, y es tau lo que fija el eje
#   temporal de este trabajo.
PI_REG = 0.10
PI_REG_RANGO = (0.03, 1.00)
K_EAPP = 2.608      # E_app = K_EAPP * Gs / a  (exacto a Pi fijo; ver paso 0)


def Gs_de_Eapp(E_app, a):
    return E_app * a / K_EAPP


# ---------------------------------------------------------------------
#  Continuacion adaptativa en la carga
# ---------------------------------------------------------------------

def continuacion(mdl, dlam0=0.12, dlam_min=4e-4, maxit=30, dlam_max=0.2,
                 max_intentos=160):
    """Continuacion con bisección del incremento de carga.  El contacto por
    penalizacion encadena iteraciones cuando la zona activa crece de golpe;
    reducir el incremento lo resuelve sin tocar el solver validado."""
    lam, dlam, x_ok = 0.0, dlam0, mdl.x.copy()
    pasos, fallos, its, intentos = [], 0, [], 0
    while lam < 1.0 - 1e-12:
        intentos += 1
        if intentos > max_intentos:
            mdl.x = x_ok
            mdl._last_lam = lam
            return False, lam, len(pasos), fallos, max(its or [0])
        trial = min(lam + dlam, 1.0)
        mdl.x = x_ok.copy()
        rep = mdl.newton(trial, maxit=maxit)
        if rep.converged and np.all(np.isfinite(mdl.x)):
            lam, x_ok = trial, mdl.x.copy()
            mdl.x_prev = x_ok.copy()
            pasos.append(trial)
            its.append(rep.iters)
            if rep.iters <= 6:
                dlam = min(dlam * 1.5, dlam_max)
        else:
            fallos += 1
            dlam *= 0.4
            if dlam < dlam_min:
                mdl.x = x_ok
                mdl._last_lam = lam
                return False, lam, len(pasos), fallos, max(its or [0])
    mdl.x = x_ok
    mdl._last_lam = 1.0
    return True, 1.0, len(pasos), fallos, max(its or [0])


# ---------------------------------------------------------------------
#  Montaje del modelo
# ---------------------------------------------------------------------

KAPPA_MU = 1e3      # nu = 0.4995 ; la Fase 1 valido 1e4 (nu = 0.49995) y
                    # el J medio se mantuvo en 1.0000000, asi que bajar a 1e3
                    # no cambia la fisica y mejora el condicionamiento del
                    # sistema cortex-rigido / citoplasma-blando


def _materiales(G_cyt, r_nuc, kappa_mu=None):
    km = kappa_mu or KAPPA_MU
    cyt = NeoHookean(name="citoplasma", C10=0.5 * G_cyt, kappa=km * G_cyt)
    gn = r_nuc * G_cyt
    nuc = NeoHookean(name="nucleo", C10=0.5 * gn, kappa=km * gn)
    return cyt, nuc


def build_constriccion(a, an, Gs, C, G_cyt, r_nuc, r_env, Rt,
                       zt=4.0, Ltap=5.0, Rwide=None, h=None, hn=None,
                       kn_f=200.0, hrel=10.0):
    """Cuarto de celula de tres compartimentos dentro de una constriccion
    axisimetrica de perfil tanh.  Identico al montaje validado en Fase 1,
    con los parametros abiertos."""
    Rwide = Rwide if Rwide is not None else 1.6 * a
    h = h if h is not None else a / hrel
    hn = hn if hn is not None else 0.7 * a / hrel
    cyt, nuc = _materiales(G_cyt, r_nuc)
    m, sub = mesh_cell_quarter(Rc=a, Rn=an, h=h, hn=hn)
    mdl = AxiModel(m, {"citoplasma": cyt, "nucleo": nuc}, subdomains=sub)
    mdl.add_dirichlet(0, lambda x: np.abs(x[0]) < 1e-9, 0.0)
    mdl.add_dirichlet(1, lambda x: np.abs(x[1]) < 1e-9, 0.0)
    outer = boundary_facets_where(m, lambda x: np.hypot(x[0], x[1]) > a - 1e-7)
    iface = interface_facets_where(
        m, lambda x: np.abs(np.hypot(x[0], x[1]) - an) < 1e-7)
    cortex = SurfaceLayer(name="cortex", Gs=Gs, C=C, law="skalak")
    envol = SurfaceLayer(name="envoltura", Gs=r_env * Gs, C=C, law="skalak")
    mdl.add_surface(cortex, outer)
    mdl.add_surface(envol, iface)

    # perfil EXACTO de la Fase 1 (tanh sin renormalizar): es la geometria en
    # la que se validaron la convergencia de malla y la calibracion 2D-3D.
    g0 = 0.5 * (1.0 + np.tanh(-zt / Ltap))
    # el camino de carga arranca en el PRIMER CONTACTO (R_min = 1.02 a) para
    # repartir los pasos de continuacion sobre el aplastamiento y no sobre el
    # tramo sin contacto; el estado final en lam = 1 es el mismo.
    Rth0 = min(Rwide, (1.02 * a - g0 * Rwide) / (1.0 - g0))

    def Rw(z, lam):
        Rth = Rth0 + lam * (Rt - Rth0)
        g = 0.5 * (1.0 + np.tanh((np.abs(z) - zt) / Ltap))
        return Rth + (Rwide - Rth) * g

    mdl.add_wall(ProfileWall(Rw=Rw, kn=kn_f * (Gs / a) / h), outer)
    return m, mdl, outer, iface, (cortex, envol)


def build_placas(a, an, Gs, C, G_cyt, r_nuc, r_env, delta1,
                 h=None, hn=None, kn_f=200.0, hrel=20.0):
    """Compresion entre placas rigidas (indentacion delta1 por cara) para
    extraer el modulo de Young aparente.  Reutiliza la geometria de los
    casos V2/V5 de la validacion."""
    h = h if h is not None else a / hrel
    hn = hn if hn is not None else 0.7 * a / hrel
    cyt, nuc = _materiales(G_cyt, r_nuc)
    m, sub = mesh_cell_quarter(Rc=a, Rn=an, h=h, hn=hn)
    mdl = AxiModel(m, {"citoplasma": cyt, "nucleo": nuc}, subdomains=sub)
    mdl.add_dirichlet(0, lambda x: np.abs(x[0]) < 1e-9, 0.0)
    mdl.add_dirichlet(1, lambda x: np.abs(x[1]) < 1e-9, 0.0)
    outer = boundary_facets_where(m, lambda x: np.hypot(x[0], x[1]) > a - 1e-7)
    iface = interface_facets_where(
        m, lambda x: np.abs(np.hypot(x[0], x[1]) - an) < 1e-7)
    mdl.add_surface(SurfaceLayer(name="cortex", Gs=Gs, C=C, law="skalak"), outer)
    mdl.add_surface(SurfaceLayer(name="envoltura", Gs=r_env * Gs, C=C,
                                 law="skalak"), iface)
    mdl.add_wall(PlaneWall(z0=a, v=-delta1, kn=kn_f * (Gs / a) / h, side=+1), outer)
    return m, mdl, outer


# ---------------------------------------------------------------------
#  Energia almacenada, por compartimento
# ---------------------------------------------------------------------

def build_P2(a, an, Gs, C, G_cyt, r_nuc, r_env, p0, hrel=14.0):
    """Celula libre (SIN contacto) cargada con una traccion normal seguidora
    de patron cuadrupolar  sigma_rr = p0 * P2(cos theta).  Es el modo de
    deformacion que impone una pinza acustica de haz unico y, a diferencia
    de la indentacion de Hertz, no tiene conjunto activo de contacto, por lo
    que Newton converge limpio tambien con cortex rigido sobre nucleo
    blando."""
    h, hn = a / hrel, 0.7 * a / hrel
    cyt, nuc = _materiales(G_cyt, r_nuc)
    m, sub = mesh_cell_quarter(Rc=a, Rn=an, h=h, hn=hn)
    mdl = AxiModel(m, {"citoplasma": cyt, "nucleo": nuc}, subdomains=sub)
    mdl.add_dirichlet(0, lambda x: np.abs(x[0]) < 1e-9, 0.0)
    mdl.add_dirichlet(1, lambda x: np.abs(x[1]) < 1e-9, 0.0)
    outer = boundary_facets_where(m, lambda x: np.hypot(x[0], x[1]) > a - 1e-7)
    iface = interface_facets_where(
        m, lambda x: np.abs(np.hypot(x[0], x[1]) - an) < 1e-7)
    if Gs > 0:
        mdl.add_surface(SurfaceLayer(name="cortex", Gs=Gs, C=C, law="skalak"),
                        outer)
        mdl.add_surface(SurfaceLayer(name="envoltura", Gs=r_env * Gs, C=C,
                                     law="skalak"), iface)
    fb = mdl.facet_basis(outer)
    X = np.asarray(fb.global_coordinates())
    ct = X[1] / np.hypot(X[0], X[1])
    patron = p0 * 0.5 * (3.0 * ct**2 - 1.0)
    mdl.add_pressure(lambda lam, arr=patron: lam * arr, outer)
    return m, mdl, outer


def _amplitud_P2(m, mdl, a):
    """(u_r ecuador - u_z polo)/a  del estado deformado."""
    nod = mdl.basis.nodal_dofs
    u = mdl.x[nod]
    d_pol = np.hypot(m.p[0], m.p[1] - a)
    d_eq = np.hypot(m.p[0] - a, m.p[1])
    ip, ie = int(np.argmin(d_pol)), int(np.argmin(d_eq))
    return float((u[0, ie] - u[1, ip]) / a)


K_P2 = None     # constante adimensional  S*E  de la sonda P2 (se calibra)


def calibrar_sonda_P2(a=8.74, an=4.37, E_ref=2.0, hrel=14.0, p0_rel=0.01):
    """Determina K_P2 = (amplitud/p0)*E sobre una esfera HOMOGENEA casi
    incompresible de modulo conocido.  Es una constante adimensional (por
    invariancia de escala) que define operativamente E_app."""
    G = E_ref / 3.0                      # incompresible: mu = E/3
    p0 = p0_rel * E_ref
    m, mdl, outer = build_P2(a, an, 0.0, 10.0, G, 1.0, 0.0, p0, hrel=hrel)
    if not _respuesta_lineal(mdl):
        return None
    return _amplitud_P2(m, mdl, a) / p0 * E_ref


def energias(mdl, capas):
    """Energia de deformacion almacenada (nN*um) desglosada por
    compartimento.  Volumen:  int Psi 2 pi R dR dZ ;
    superficies:  int w(ls,lt) 2 pi R_ref dS_ref."""
    b, pb = mdl.basis, mdl.pbasis
    R = b.global_coordinates()[0]
    uq = b.interpolate(mdl.u)
    pq = np.asarray(pb.interpolate(mdl.p))
    F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], R)
    J = det3(F)
    I1 = np.einsum('ij...,ij...->...', F, F)
    Psi = (mdl.C10 * (J ** (-2.0 / 3.0) * I1 - 3.0)
           + pq * (J - 1.0) - pq**2 / (2.0 * mdl.kappa))
    wq = 2.0 * np.pi * R * b.dx
    out = {}
    for nm, els in mdl.subdomains.items():
        out[f"W_{nm}"] = float(np.sum(Psi[els] * wq[els]))

    for capa, (_, fb) in zip(capas, mdl.surfaces):
        X = np.asarray(fb.global_coordinates())
        Rref = X[0]
        n = np.asarray(fb.normals)
        t = np.stack([-n[1], n[0]])
        uqf = fb.interpolate(mdl.u)
        gu = np.asarray(uqf.grad)
        a0 = t[0] + gu[0, 0] * t[0] + gu[0, 1] * t[1]
        a1 = t[1] + gu[1, 0] * t[0] + gu[1, 1] * t[1]
        ls = np.sqrt(a0**2 + a1**2)
        lt = (Rref + np.asarray(uqf)[0]) / Rref
        ws = capa.w(ls, lt)
        out[f"W_{capa.name}"] = float(np.sum(ws * 2.0 * np.pi * Rref * fb.dx))
        out[f"ls_max_{capa.name}"] = float(ls.max())
        out[f"lt_max_{capa.name}"] = float(lt.max())
    out["W_total"] = sum(v for k, v in out.items() if k.startswith("W_"))
    return out


def area_membrana(mdl, outer):
    """Deformacion relativa de area de la membrana externa."""
    fb = mdl.facet_basis(outer)
    X = np.asarray(fb.global_coordinates())
    n = np.asarray(fb.normals)
    t = np.stack([-n[1], n[0]])
    uq = fb.interpolate(mdl.u)
    gu = np.asarray(uq.grad)
    a0 = t[0] + gu[0, 0] * t[0] + gu[0, 1] * t[1]
    a1 = t[1] + gu[1, 0] * t[0] + gu[1, 1] * t[1]
    ls = np.sqrt(a0**2 + a1**2)
    lt = (X[0] + np.asarray(uq)[0]) / X[0]
    A0 = float(np.sum(2.0 * np.pi * X[0] * fb.dx))
    Ad = float(np.sum(ls * lt * 2.0 * np.pi * X[0] * fb.dx))
    return Ad / A0 - 1.0


# ---------------------------------------------------------------------
#  Resolucion de un punto del barrido
# ---------------------------------------------------------------------

def resolver_constriccion(a, an, Gs, C, G_cyt, r_nuc, r_env, Rt,
                          zt=4.0, Ltap=5.0, cont_kw=None, **kw):
    """Devuelve el dict de energias + metricas, o None si no converge."""
    m, mdl, outer, iface, capas = build_constriccion(
        a, an, Gs, C, G_cyt, r_nuc, r_env, Rt, zt=zt, Ltap=Ltap, **kw)
    ok, lam, npas, nfal, itmax = continuacion(mdl, **(cont_kw or {}))
    rec = dict(Rt=Rt, conf=a / Rt, converge=bool(ok), lam_alcanzado=float(lam),
               ne=int(m.t.shape[1]), ndof=int(mdl.N),
               pasos_carga=int(npas), reducciones=int(nfal),
               iters_max=int(itmax))
    if not ok:
        return rec
    rec.update(energias(mdl, capas))
    rec["dA_A0"] = area_membrana(mdl, outer)
    rec["Fz_contacto"] = float(mdl.contact_force_z())
    return rec


def _respuesta_lineal(mdl):
    """Respuesta LINEALIZADA exacta: el primer paso de Newton desde el estado
    sin deformar resuelve K(0) dx = f, que es la solucion del problema lineal.
    Evita la bifurcacion de arrugamiento de la membrana (que es un fenomeno de
    carga finita) y basta un unico sistema disperso."""
    import scipy.sparse.linalg as spl
    mdl.x = np.zeros(mdl.N)
    mdl._last_lam = 1.0
    D, Dval = mdl.constrained_dofs(1.0)
    free = np.setdiff1d(np.arange(mdl.N), D)
    R, K = mdl.residual_and_tangent(mdl.x, 1.0)
    dx = spl.spsolve(K[free][:, free].tocsc(), -R[free])
    if not np.all(np.isfinite(dx)):
        return False
    mdl.x[free] = dx
    return True


def modulo_aparente(a, an, Gs, C, G_cyt, r_nuc, r_env, p0_rel=0.01,
                    K=1.18559, hrel=20.0):
    """Modulo de Young aparente en suspension (kPa), definido operativamente
    como el modulo de la esfera homogenea casi incompresible que mostraria la
    MISMA amplitud de deformacion cuadrupolar bajo la MISMA traccion P2.

    Se usa la respuesta LINEALIZADA (_respuesta_lineal): es exacta en el regimen
    de pequena deformacion en que se mide E_app y evita la bifurcacion de
    arrugamiento de la membrana sin rigidez de flexion, que impide converger el
    problema no lineal con cortex rigido.  K = 1.18559 es el valor por defecto,
    obtenido con calibrar_sonda_P2(hrel=20) sobre esfera homogenea de modulo
    conocido; es una constante adimensional (invariante en a, E y p0) y difiere
    en 0.13% del valor no lineal convergido 1.18687 y en 0.02% del obtenido a
    hrel=14.  Si se cambia hrel hay que recalibrar K con la misma malla."""
    p0 = p0_rel * 2.0
    m, mdl, outer = build_P2(a, an, Gs, C, G_cyt, r_nuc, r_env, p0, hrel=hrel)
    if not _respuesta_lineal(mdl):
        return None
    S = _amplitud_P2(m, mdl, a) / p0
    return float(K / S) if S > 0 else None


def calibrar_Gs(est, E_obj=None, K=1.18508, hrel=20.0, bracket=(0.05, 200.0)):
    """Resuelve Gs tal que el modulo aparente del modelo iguale el medido.
    E_app(Gs) es monotona y suave (funcion racional: la rigidez es afin en Gs,
    la compliancia no), asi que se invierte por biseccion de Brent."""
    from scipy.optimize import brentq
    E_obj = E_obj if E_obj is not None else est.E_app
    if E_obj is None:
        return None, None, 0
    n = [0]

    def g(x):
        n[0] += 1
        e = modulo_aparente(est.a, est.an, x, est.C_skalak, est.G_cyt,
                            est.r_nuc, est.r_env, K=K, hrel=hrel)
        return (e if e is not None else np.nan) - E_obj

    lo, hi = bracket
    if g(lo) > 0:
        return float(lo), float(E_obj + g(lo)), n[0]
    gs = brentq(g, lo, hi, xtol=1e-4, rtol=1e-8)
    e = modulo_aparente(est.a, est.an, gs, est.C_skalak, est.G_cyt,
                        est.r_nuc, est.r_env, K=K, hrel=hrel)
    return float(gs), float(e), n[0]


# =====================================================================
#  REDUCCION A UN GRADO DE LIBERTAD
#  ------------------------------------------------------------------
#  El FEM entrega el paisaje energetico W(R) de confinamiento radial
#  (celula centrada en una garganta de radio R).  La entrada se reduce a
#  una coordenada axial s = posicion del ecuador celular, con la
#  aproximacion declarada
#        W(s) ~= W_FEM(R_w(s))
#  (la celula a la posicion s tiene el estado elastico del equilibrio
#  simetrico al radio local de pared).  La fuerza axial resistente es
#        dW/ds = |dW/dR| * |dR_w/ds|
#  y el avance exige   dP * pi R_w(s)^2 >= dW/ds.
# =====================================================================

def perfil_Rw(z, Rt, Rwide, zt, Ltap):
    """Perfil de pared de la Fase 1 (tanh sin renormalizar).

    OJO con la convencion heredada: en z = 0 vale
        R_min = Rt + g0 (Rwide - Rt),   g0 = 0.5 (1 + tanh(-zt/Ltap)),
    o sea que el radio nominal Rt NO se alcanza en el cuello (con zt/Ltap=0.8,
    g0 = 0.168).  El 'confinamiento 2a/D = a/Rt' de la Fase 1 es por tanto
    NOMINAL; el efectivo es a/R_min y siempre es menor.  Se reportan ambos."""
    g = 0.5 * (1.0 + np.tanh((np.abs(z) - zt) / Ltap))
    return Rt + (Rwide - Rt) * g


def R_minimo(Rt, Rwide, zt, Ltap):
    return float(perfil_Rw(0.0, Rt, Rwide, zt, Ltap))


def conf_efectivo(a, Rt, Rwide, zt, Ltap):
    return float(a / R_minimo(Rt, Rwide, zt, Ltap))


def dRw_dz(z, Rt, Rwide, zt, Ltap):
    ch = np.cosh((np.abs(z) - zt) / Ltap)
    return (Rwide - Rt) * 0.5 / (Ltap * ch**2) * np.sign(z)


def angulo_de_Ltap(Ltap, Rt, Rwide, zt=4.0):
    """Semiangulo de entrada efectivo: tan(theta) = max|dR_w/dz|."""
    z = np.linspace(0.0, zt + 8.0 * Ltap, 4001)
    return float(np.degrees(np.arctan(
        np.abs(dRw_dz(z, Rt, Rwide, zt, Ltap)).max())))


def longitud_de_carga(a, Rt, Rwide, zt, Ltap, nz=20001):
    """Distancia axial entre el PRIMER CONTACTO (R_w = a) y el confinamiento
    pleno (R_w = Rt).  Es el recorrido sobre el que se desarrolla la
    deformacion, y por tanto el que fija el tiempo de carga impuesto."""
    z = np.linspace(0.0, zt + 10.0 * Ltap, nz)
    Rw = perfil_Rw(z, Rt, Rwide, zt, Ltap)
    i = np.where(Rw <= a)[0]
    return float(z[i.max()]) if len(i) else 0.0


def presion_critica(a, Rt, Rwide, zt, Ltap, W_of_R, dW_dR, nz=6001):
    """Presion critica de paso (kPa) y geometria del camino de entrada.

    El maximo se alcanza en el propio cuello (z -> 0+): el perfil tanh(|z|)
    tiene un QUIEBRE ahi, de modo que |dR_w/dz| salta de 0 a su valor finito y
    el radio es minimo.  Por eso hace falta rejilla fina: de nz=301 a nz=12001
    el valor sube de 0.1368 a 0.1609 kPa (caso MCF-7, c=1.70); con nz=6001 el
    residuo frente al valor converguido es del orden del 1%.

    Solo se maximiza DENTRO de la envolvente de validez: el tramo del
    camino con R_w fuera de [a/CONF_MAX, a/CONF_MIN] se excluye.  Como en
    el tramo 2a/D < 1.2 el modelo predice energia nula (no hay contacto),
    la presion critica que sale es una COTA INFERIOR.
    """
    Rmin_val = R_minimo(a / CONF_MAX, Rwide, zt, Ltap)
    Rmax_val = R_minimo(a / CONF_MIN, Rwide, zt, Ltap)
    z = np.linspace(0.0, zt + 6.0 * Ltap, nz)
    Rw = perfil_Rw(z, Rt, Rwide, zt, Ltap)
    dR = np.abs(dRw_dz(z, Rt, Rwide, zt, Ltap))
    val = (Rw >= Rmin_val) & (Rw <= Rmax_val)
    if not np.any(val):
        return dict(dP_crit=np.nan, L_entrada=np.nan, z_crit=np.nan,
                    R_crit=np.nan, n_camino_valido=0)
    dWds = np.abs(dW_dR(Rw)) * dR
    q = dWds / (np.pi * Rw**2)
    q = np.where(val, q, -np.inf)
    i = int(np.argmax(q))
    zv = z[val]
    return dict(dP_crit=float(q[i]),
                L_banda_valida=float(zv.max() - zv.min()),
                L_entrada=longitud_de_carga(a, Rt, Rwide, zt, Ltap),
                z_crit=float(z[i]), R_crit=float(Rw[i]),
                conf_crit=float(a / Rw[i]), n_camino_valido=int(val.sum()))


def tiempo_transito(a, Rt, Rwide, zt, Ltap, W_of_R, dW_dR, dP, tau, nz=2001):
    """Tiempo de transito de la ETAPA DE ENTRADA (s) bajo presion impuesta
    dP, con respuesta de Kelvin-Voigt:  sigma = G (eps + tau d(eps)/dt).
    Como la energia es lineal en los modulos a campo de deformacion fijo,
    la fuerza generalizada se multiplica por (1 + tau * eps_punto/eps):

        dP pi R_w^2 = (dW/ds) [1 + tau (deps/ds) s_punto / eps]

    Devuelve None si dP no basta en algun punto del camino valido."""
    Rmin_val = R_minimo(a / CONF_MAX, Rwide, zt, Ltap)
    Rmax_val = R_minimo(a / CONF_MIN, Rwide, zt, Ltap)
    z = np.linspace(0.0, zt + 6.0 * Ltap, nz)
    Rw = perfil_Rw(z, Rt, Rwide, zt, Ltap)
    val = (Rw >= Rmin_val) & (Rw <= Rmax_val)
    if val.sum() < 5:
        return None
    z, Rw = z[val], Rw[val]
    o = np.argsort(-z)                    # avance: de z grande hacia el cuello
    z, Rw = z[o], Rw[o]
    s = z[0] - z
    dRds = np.abs(np.gradient(Rw, s))
    dWds = np.abs(dW_dR(Rw)) * dRds
    eps = 1.0 - Rw / a
    deps_ds = dRds / a
    with np.errstate(divide='ignore', invalid='ignore'):
        m = dP * np.pi * Rw**2 / dWds
    # dW/ds -> 0 en el cuello (dR_w/dz = 0 por simetria): ahi no hay
    # resistencia elastica, m = inf y el avance es instantaneo; eso NO es un
    # fallo, aporta cero al tiempo.  Solo hay bloqueo si m <= 1 donde dW/ds > 0.
    finito = np.isfinite(m)
    if not finito.any() or np.nanmin(m[finito]) <= 1.0:
        return None
    with np.errstate(divide='ignore', invalid='ignore'):
        sdot = (m - 1.0) * eps / (tau * deps_ds)
        inv = np.where(np.isfinite(sdot) & (sdot > 0), 1.0 / sdot, 0.0)
    if np.count_nonzero(inv) < 5:
        return None
    return float(np.trapezoid(inv, s))


# ---------------------------------------------------------------------
#  Operacion a CAUDAL impuesto  (nucleo del articulo)
# ---------------------------------------------------------------------

def deborah(tau, U, L_entrada):
    """De = tau / t_carga  con  t_carga = L_entrada / U.
    U en um/s, L_entrada en um, tau en s."""
    return tau * U / L_entrada


def margen_caudal(U, dP_el, tau, L_entrada, beta):
    """Margen de paso a velocidad de garganta U.

    Bajo rampa de duracion t_carga = L_entrada/U un solido de Kelvin-Voigt
    responde con rigidez aparente G(1 + tau/t_carga) = G(1 + De), luego
        dP_crit(U) = dP_el (1 + De).
    La presion motriz en operacion a caudal constante crece linealmente con
    el caudal,  dP_motriz = beta * U  (beta = resistencia hidraulica del
    dispositivo x area de garganta; parametro de diseno declarado).
        m(U) = beta U / dP_crit(U) = m_inf * De/(1+De),
        m_inf = beta L_entrada / (dP_el tau).
    La celula PASA si m > 1; queda CAPTURADA si m <= 1.
    Consecuencia: si m_inf <= 1 la celula no pasa a NINGUN caudal."""
    De = deborah(tau, U, L_entrada)
    m_inf = beta * L_entrada / (dP_el * tau)
    return m_inf * De / (1.0 + De), De, m_inf


def U_critica(dP_el, tau, L_entrada, beta):
    """Velocidad de garganta a la que m(U) = 1.  np.inf si m_inf <= 1."""
    m_inf = beta * L_entrada / (dP_el * tau)
    if m_inf <= 1.0:
        return np.inf
    De_c = 1.0 / (m_inf - 1.0)
    return De_c * L_entrada / tau


# ---------------------------------------------------------------------
#  Microfiltro: barrera de energia de entrada en poro de canto vivo
# ---------------------------------------------------------------------

def presion_poro(W_Rp, Rp, a, ell_factor=1.0):
    """Umbral de presion de paso por un poro de canto vivo de radio Rp.
    El canto vivo no define una pendiente dR_w/ds, asi que la barrera se
    estima como la energia total del estado confinado repartida sobre la
    longitud de deglucion  ell = ell_factor * a  (la celula debe avanzar
    del orden de un radio para quedar engullida).  ell_factor se declara y
    entra en el analisis de sensibilidad."""
    return W_Rp / (ell_factor * a * np.pi * Rp**2)
