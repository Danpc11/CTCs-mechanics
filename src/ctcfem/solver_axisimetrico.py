"""
solver_axisimetrico.py
=======================
Solver cuasi-estatico AXISIMETRICO de gran deformacion para celulas
multicompartimento (membrana/cortex + citoplasma + envoltura nuclear +
nucleoplasma) en canales microfluidicos.

Formulacion
-----------
* Coordenadas de referencia  X = (R, Z);  desplazamiento  u = (u_r, u_z).
* Gradiente de deformacion axisimetrico (3x3, componente circunferencial
  explicita):

        F = [[1+du_r/dR,   du_r/dZ ,      0      ],
             [  du_z/dR , 1+du_z/dZ,      0      ],
             [     0    ,     0    , 1 + u_r/R  ]]

* Volumen: dV = 2*pi*R dR dZ (configuracion de referencia).
* Material: Neo-Hookeano casi incompresible en formulacion MIXTA
  desplazamiento-presion (Lagrangiano perturbado / Simo-Taylor-Pister):

        Pi = int [ C10*(I1bar - 3) + p*(J-1) - p^2/(2*kappa) ] dV

  con I1bar = J^(-2/3) tr(F^T F).  Primer tensor de Piola-Kirchhoff:

        P = 2*C10*J^(-2/3)*(F - (I1/3)*F^-T) + p*J*F^-T

  Elementos: desplazamiento P2+burbuja (ElementTriP2B) vectorial,
  presion P1 discontinua (ElementTriP1DG)  -> par de Crouzeix-Raviart,
  inf-sup estable y con salto de presion permitido en las interfases
  entre compartimentos (nucleo/citoplasma).

* Superficies (cortex y envoltura nuclear): capa delgada con energia por
  unidad de area de REFERENCIA funcion de los estiramientos principales
  de superficie (meridional lambda_s y circunferencial lambda_theta):

        w_s = gamma*lambda_s*lambda_theta                       (tension)
            + (Gs/2)*(ls^2 + lt^2 + 1/(ls*lt)^2 - 3)            (cizalla 2D)
            + (Ks/2)*(ls*lt - 1)^2                              (dilatacion)

  integrada como int w_s * 2*pi*R dS.  Opcionalmente un termino viscoso
  de Kelvin-Voigt  (eta_s/dt)*((ls-ls_n)^2+(lt-lt_n)^2)/2  para los
  transitorios.

* Contacto con la pared del canal por PENALIZACION:
        w_c = (k_n/2) * <g>_+^2 ,  g = violacion de la restriccion,
  con la pared definida por un perfil r = Rw(z) o un plano z = cte.

* Carga por presion seguidora (follower) con la relacion de Nanson
  exacta,  t dA = p*J*F^-T N dA_0.

* Newton-Raphson con continuacion en la carga; tangente consistente.
  Las derivadas del material se obtienen por PASO COMPLEJO
  (complex-step, error de truncamiento nulo) y las de los terminos de
  superficie/contacto por paso complejo + diferencia central en el
  espacio (lambda_s, lambda_theta), que esta bien condicionado.

Autor: generado para el proyecto CTC-mama.  Unidades SI salvo aviso
(longitudes en m, modulos en Pa, fuerzas en N).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

from skfem import (Basis, FacetBasis, BilinearForm, LinearForm, Functional,
                   ElementVector, ElementTriP2B, ElementTriP1DG,
                   ElementTriP2, ElementTriP1, MeshTri)

__all__ = ["piola", "piola_tangent", "dpiola_dp", "axi_F", "axi_G", "det3", "inv3", "NeoHookean", "SurfaceLayer", "Wall", "PlaneWall", "ProfileWall",
           "AxiModel", "NewtonReport", "mesh_quarter_disk", "mesh_cell_with_nucleus", "constriction_profile", "boundary_facets_where", "interface_facets_where", "mesh_cell_quarter"]

# ---------------------------------------------------------------------
#  Algebra 3x3 compatible con numeros complejos (necesario para el
#  complex-step; numpy.linalg.det/inv no vectorizan como queremos)
# ---------------------------------------------------------------------


def det3(A):
    return (A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
            - A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
            + A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]))


def inv3(A):
    d = det3(A)
    B = np.empty_like(A)
    B[0, 0] = (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) / d
    B[0, 1] = (A[0, 2] * A[2, 1] - A[0, 1] * A[2, 2]) / d
    B[0, 2] = (A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]) / d
    B[1, 0] = (A[1, 2] * A[2, 0] - A[1, 0] * A[2, 2]) / d
    B[1, 1] = (A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0]) / d
    B[1, 2] = (A[0, 2] * A[1, 0] - A[0, 0] * A[1, 2]) / d
    B[2, 0] = (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]) / d
    B[2, 1] = (A[0, 1] * A[2, 0] - A[0, 0] * A[2, 1]) / d
    B[2, 2] = (A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]) / d
    return B


def transpose3(A):
    return np.einsum('ij...->ji...', A)


# indices de F con contribucion no nula en axisimetria
AXI_IDX = ((0, 0), (0, 1), (1, 0), (1, 1), (2, 2))


def axi_F(gradu, ur, R):
    """Gradiente de deformacion 3x3 a partir de grad(u) (2x2) y u_r/R."""
    sh = gradu.shape[2:]
    F = np.zeros((3, 3) + sh, dtype=gradu.dtype)
    F[0, 0] = 1.0 + gradu[0, 0]
    F[0, 1] = gradu[0, 1]
    F[1, 0] = gradu[1, 0]
    F[1, 1] = 1.0 + gradu[1, 1]
    F[2, 2] = 1.0 + ur / R
    return F


def axi_G(v, R):
    """'Gradiente' 3x3 de una funcion de forma vectorial (operador B)."""
    gv = v.grad
    sh = gv.shape[2:]
    G = np.zeros((3, 3) + sh)
    G[0, 0] = gv[0, 0]
    G[0, 1] = gv[0, 1]
    G[1, 0] = gv[1, 0]
    G[1, 1] = gv[1, 1]
    G[2, 2] = v[0] / R
    return G


# ---------------------------------------------------------------------
#  Material
# ---------------------------------------------------------------------

@dataclass
class NeoHookean:
    """Neo-Hookeano casi incompresible (mixto u-p).

    C10 = mu/2 con mu el modulo de cizalla en pequena deformacion.
    kappa = modulo volumetrico.
    eta  = viscosidad de Kelvin-Voigt (Pa s), opcional, solo transitorios.
    """
    name: str = "material"
    C10: float = 1.0
    kappa: float = 1.0e4
    eta: float = 0.0

    @property
    def mu(self):
        return 2.0 * self.C10

    @property
    def E(self):
        mu, k = self.mu, self.kappa
        return 9.0 * k * mu / (3.0 * k + mu)

    @property
    def nu(self):
        mu, k = self.mu, self.kappa
        return (3.0 * k - 2.0 * mu) / (2.0 * (3.0 * k + mu))

    @staticmethod
    def from_E_nu(E, nu, name="material", eta=0.0):
        mu = E / (2.0 * (1.0 + nu))
        kappa = E / (3.0 * (1.0 - 2.0 * nu))
        return NeoHookean(name=name, C10=0.5 * mu, kappa=kappa, eta=eta)


def piola(F, p, C10):
    """Primer Piola-Kirchhoff del Neo-Hookeano desviatorico + parte p.

    Valido para F complejo (complex-step)."""
    J = det3(F)
    Fit = transpose3(inv3(F))
    I1 = np.einsum('ij...,ij...->...', F, F)
    Jm23 = J ** (-2.0 / 3.0)
    return 2.0 * C10 * Jm23 * (F - (I1 / 3.0) * Fit) + p * J * Fit


def piola_tangent(F, p, C10, h=1e-20):
    """A_{iJkL} = dP_iJ/dF_kL por paso complejo (exacto a precision maquina)."""
    sh = F.shape[2:]
    A = np.zeros((3, 3, 3, 3) + sh)
    Fc = F.astype(complex)
    for (k, L) in AXI_IDX:
        Fp = Fc.copy()
        Fp[k, L] = Fp[k, L] + 1j * h
        A[:, :, k, L] = piola(Fp, p, C10).imag / h
    return A


def dpiola_dp(F):
    """dP/dp = J F^-T."""
    return det3(F) * transpose3(inv3(F))


# ---------------------------------------------------------------------
#  Capas superficiales (cortex, envoltura nuclear)
# ---------------------------------------------------------------------

@dataclass
class SurfaceLayer:
    """Energia de superficie por unidad de area de referencia."""
    name: str = "cortex"
    gamma: float = 0.0      # tension cortical [N/m]
    Gs: float = 0.0         # modulo de cizalla de superficie [N/m]
    Ks: float = 0.0         # modulo de dilatacion de area [N/m]
    eta_s: float = 0.0      # viscosidad de superficie [N s/m]
    law: str = "neohookean"  # 'neohookean' | 'skalak'
    C: float = None          # parametro de dureza de Skalak; Ks=(1+2C)Gs

    def __post_init__(self):
        if self.law == "skalak":
            if self.C is None and self.Gs:
                self.C = (self.Ks / self.Gs - 1.0) / 2.0
            if self.Ks == 0.0 and self.C is not None:
                self.Ks = (1.0 + 2.0 * self.C) * self.Gs

    def w(self, ls, lt, ls_n=None, lt_n=None, dt=None):
        """Densidad de energia; admite argumentos complejos."""
        out = self.gamma * ls * lt
        if self.law == "skalak" and self.Gs:
            # Skalak et al. (1973):
            #   W = Gs/4 * (I1^2 + 2 I1 - 2 I2) + C*Gs/4 * I2^2
            I1 = ls**2 + lt**2 - 2.0
            I2 = (ls * lt)**2 - 1.0
            out = out + 0.25 * self.Gs * (I1**2 + 2.0 * I1 - 2.0 * I2) \
                + 0.25 * self.C * self.Gs * I2**2
        elif self.Gs:
            out = out + 0.5 * self.Gs * (ls**2 + lt**2
                                         + 1.0 / (ls * lt)**2 - 3.0)
        if self.Ks and self.law != "skalak":
            out = out + 0.5 * self.Ks * (ls * lt - 1.0)**2
        if self.eta_s and dt:
            out = out + 0.5 * (self.eta_s / dt) * ((ls - ls_n)**2
                                                  + (lt - lt_n)**2)
        return out


def _grad_hess_2(fun, x, y, d=1e-5, h=1e-20):
    """Gradiente (paso complejo, exacto) y Hessiano (diferencia central del
    gradiente complejo, error O(d^2)) de fun(x, y) evaluada en arrays."""
    def gx(xx, yy):
        return fun(xx + 1j * h, yy).imag / h

    def gy(xx, yy):
        return fun(xx, yy + 1j * h).imag / h

    g = np.stack([gx(x, y), gy(x, y)])
    dx = d * np.maximum(np.abs(x), 1.0)
    dy = d * np.maximum(np.abs(y), 1.0)
    H = np.empty((2, 2) + x.shape)
    H[0, 0] = (gx(x + dx, y) - gx(x - dx, y)) / (2 * dx)
    H[1, 1] = (gy(x, y + dy) - gy(x, y - dy)) / (2 * dy)
    H01 = (gx(x, y + dy) - gx(x, y - dy)) / (2 * dy)
    H10 = (gy(x + dx, y) - gy(x - dx, y)) / (2 * dx)
    H[0, 1] = H[1, 0] = 0.5 * (H01 + H10)
    return g, H


# ---------------------------------------------------------------------
#  Paredes rigidas (contacto por penalizacion)
# ---------------------------------------------------------------------

class Wall:
    """Interfaz: gap(r, z, t) > 0 significa penetracion."""
    kn: float = 0.0

    def gap(self, r, z, lam):
        raise NotImplementedError


@dataclass
class PlaneWall(Wall):
    """Plano rigido z = z0 + lam*v.  side=+1: el solido debe quedar z<=zw."""
    z0: float = 0.0
    v: float = 0.0
    kn: float = 1.0e6
    side: int = +1

    def gap(self, r, z, lam):
        zw = self.z0 + lam * self.v
        return self.side * (z - zw)


@dataclass
class ProfileWall(Wall):
    """Pared axisimetrica r = Rw(z).  El solido debe quedar r <= Rw(z)."""
    Rw: Callable = None
    kn: float = 1.0e6

    def gap(self, r, z, lam):
        return r - self.Rw(z, lam)


def constriction_profile(R_wide, R_throat, z_throat, L_taper, sharp=1.0):
    """Perfil suave de constriccion (tanh) centrado en z=0, garganta de
    semiancho R_throat y longitud 2*z_throat."""
    def Rw(z, lam=1.0):
        s = (np.abs(z) - z_throat) / L_taper
        f = 0.5 * (1.0 + np.tanh(sharp * s))
        return R_throat + (R_wide - R_throat) * f
    return Rw


# ---------------------------------------------------------------------
#  Informe de Newton
# ---------------------------------------------------------------------

@dataclass
class NewtonReport:
    converged: bool = False
    iters: int = 0
    residuals: list = field(default_factory=list)
    message: str = ""


# ---------------------------------------------------------------------
#  Modelo
# ---------------------------------------------------------------------

class AxiModel:
    """Modelo axisimetrico multicompartimento.

    Parameters
    ----------
    mesh : skfem.MeshTri en coordenadas (R, Z), R >= 0.
    materials : dict nombre -> NeoHookean.
    subdomains : dict nombre -> array de indices de elementos.  Si es None
        todo el dominio usa materials['default'].
    surfaces : lista de (SurfaceLayer, facets)
    walls : lista de Wall con facets asociadas: (Wall, facets)
    pressure_loads : lista de (valor, facets) presion seguidora.
    dirichlet : lista de (component, facet_selector_or_facets, value_fn)
        component en {0,1} (r, z); value_fn(lam, x) -> valor del dof.
    element : 'P2B-P1DG' (por defecto) o 'P2-P1' (Taylor-Hood).
    """

    def __init__(self, mesh, materials, subdomains=None, element="P2B-P1DG",
                 intorder=None):
        self.mesh = mesh
        self.materials = materials
        if element == "P2B-P1DG":
            eu, ep, order = ElementTriP2B(), ElementTriP1DG(), 4
        elif element == "P2-P1":
            eu, ep, order = ElementTriP2(), ElementTriP1(), 4
        else:
            raise ValueError(element)
        self.intorder = intorder or order
        self.eu = ElementVector(eu)
        self.basis = Basis(mesh, self.eu, intorder=self.intorder)
        self.pbasis = Basis(mesh, ep, intorder=self.intorder)
        self.nu_dofs = self.basis.N
        self.np_dofs = self.pbasis.N
        self.N = self.nu_dofs + self.np_dofs

        ne = mesh.t.shape[1]
        nq = self.basis.X.shape[1]
        self.C10 = np.zeros((ne, nq))
        self.kappa = np.zeros((ne, nq))
        self.eta = np.zeros((ne, nq))
        if subdomains is None:
            subdomains = {"default": np.arange(ne)}
        self.subdomains = subdomains
        for nm, els in subdomains.items():
            mat = materials[nm]
            self.C10[els, :] = mat.C10
            self.kappa[els, :] = mat.kappa
            self.eta[els, :] = mat.eta

        self.surfaces = []      # (SurfaceLayer, FacetBasis)
        self.walls = []         # (Wall, FacetBasis)
        self.pressure_loads = []
        self.dirichlet = []
        self.x = np.zeros(self.N)
        self.x_prev = np.zeros(self.N)
        self.dt = None

    # ---------------- construccion -----------------
    def facet_basis(self, facets):
        return FacetBasis(self.mesh, self.eu, facets=facets,
                          intorder=self.intorder)

    def add_surface(self, layer, facets):
        self.surfaces.append((layer, self.facet_basis(facets)))

    def add_wall(self, wall, facets):
        self.walls.append((wall, self.facet_basis(facets)))

    def add_pressure(self, value, facets):
        """Presion seguidora.  value > 0 = presion que COMPRIME el solido
        a traves de esa superficie (sentido -n).  Para presurizar una
        cavidad interna basta con pasar el valor positivo sobre las
        facetas de la cavidad."""
        self.pressure_loads.append((value, self.facet_basis(facets)))

    def add_dirichlet(self, component, facets, value=0.0):
        key = 'u^1' if component == 0 else 'u^2'
        dofs = self.basis.get_dofs(facets).all(key)
        self.dirichlet.append((np.asarray(dofs), value))

    # ---------------- ensamblado -----------------
    def _bulk(self, u, p):
        """Residuo y tangentes del volumen."""
        b, pb = self.basis, self.pbasis
        R = b.global_coordinates()[0]          # (ne, nq) coordenada radial ref
        uq = b.interpolate(u)
        pq = pb.interpolate(p)
        F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], R)
        pv = np.asarray(pq)
        P = piola(F, pv, self.C10)
        A = piola_tangent(F, pv, self.C10)
        dPdp = dpiola_dp(F)
        J = det3(F)

        if self.dt is not None:
            # Kelvin-Voigt (aproximacion de rotacion pequena)
            uq_n = b.interpolate(self.x_prev[:self.nu_dofs])
            Fn = axi_F(np.asarray(uq_n.grad), np.asarray(uq_n)[0], R)
            P = P + self.eta * (F - Fn) / self.dt
            for (k, L) in AXI_IDX:
                A[k, L, k, L] = A[k, L, k, L] + self.eta / self.dt

        w2pi = 2.0 * np.pi * R

        @LinearForm
        def res_u(v, w):
            Gv = axi_G(v, R)
            return w2pi * np.einsum('ij...,ij...->...', P, Gv)

        @LinearForm
        def res_p(q, w):
            return w2pi * (J - 1.0 - pv / self.kappa) * np.asarray(q)

        @BilinearForm
        def k_uu(uu, vv, w):
            Gu = axi_G(uu, R)
            Gv = axi_G(vv, R)
            return w2pi * np.einsum('ijkl...,ij...,kl...->...', A, Gv, Gu)

        @BilinearForm
        def k_up(qq, vv, w):
            Gv = axi_G(vv, R)
            return w2pi * np.asarray(qq) * np.einsum('ij...,ij...->...',
                                                     dPdp, Gv)

        @BilinearForm
        def k_pp(qq, rr, w):
            return -w2pi * np.asarray(qq) * np.asarray(rr) / self.kappa

        Ru = res_u.assemble(b)
        Rp = res_p.assemble(pb)
        Kuu = k_uu.assemble(b)
        Kup = k_up.assemble(pb, b)
        Kpp = k_pp.assemble(pb)
        return Ru, Rp, Kuu, Kup, Kpp

    def _surface_terms(self, u, lam):
        """Residuo y tangente de las capas superficiales."""
        Ru = np.zeros(self.nu_dofs)
        rows, cols, vals = [], [], []
        Kglob = sp.csr_matrix((self.nu_dofs, self.nu_dofs))
        for layer, fb in self.surfaces:
            x = fb.global_coordinates()
            Rref = np.asarray(x)[0]
            n = np.asarray(fb.normals)
            t = np.stack([-n[1], n[0]])        # tangente unitaria (r,z)
            uq = fb.interpolate(u)
            gu = np.asarray(uq.grad)
            # a = (I + grad u) . t
            a0 = t[0] + gu[0, 0] * t[0] + gu[0, 1] * t[1]
            a1 = t[1] + gu[1, 0] * t[0] + gu[1, 1] * t[1]
            ls = np.sqrt(a0**2 + a1**2)
            ur = np.asarray(uq)[0]
            lt = (Rref + ur) / Rref

            if self.dt is not None and layer.eta_s:
                uq_n = fb.interpolate(self.x_prev[:self.nu_dofs])
                gun = np.asarray(uq_n.grad)
                b0 = t[0] + gun[0, 0] * t[0] + gun[0, 1] * t[1]
                b1 = t[1] + gun[1, 0] * t[0] + gun[1, 1] * t[1]
                ls_n = np.sqrt(b0**2 + b1**2)
                lt_n = (Rref + np.asarray(uq_n)[0]) / Rref
                fun = lambda s, th: layer.w(s, th, ls_n, lt_n, self.dt)
            else:
                fun = lambda s, th: layer.w(s, th)

            g2, H2 = _grad_hess_2(fun, ls, lt)
            w_s, w_t = g2[0], g2[1]

            # cadena: dls/da_i = a_i/ls ;  dlt/du_r = 1/Rref
            inv_ls = 1.0 / ls
            na0, na1 = a0 * inv_ls, a1 * inv_ls
            S0 = w_s * na0
            S1 = w_s * na1
            Q = w_t / Rref
            meas = 2.0 * np.pi * Rref

            # Hessiano en (a0, a1, ur)
            Haa00 = (H2[0, 0] * na0 * na0
                     + w_s * (1.0 - na0 * na0) * inv_ls)
            Haa11 = (H2[0, 0] * na1 * na1
                     + w_s * (1.0 - na1 * na1) * inv_ls)
            Haa01 = (H2[0, 0] * na0 * na1
                     + w_s * (-na0 * na1) * inv_ls)
            Hau0 = H2[0, 1] * na0 / Rref
            Hau1 = H2[0, 1] * na1 / Rref
            Huu = H2[1, 1] / Rref**2

            tt = t

            @LinearForm
            def rs(v, w):
                gv = v.grad
                dv0 = gv[0, 0] * tt[0] + gv[0, 1] * tt[1]
                dv1 = gv[1, 0] * tt[0] + gv[1, 1] * tt[1]
                return meas * (S0 * dv0 + S1 * dv1 + Q * v[0])

            @BilinearForm
            def ks(uu, vv, w):
                gu_ = uu.grad
                gv_ = vv.grad
                du0 = gu_[0, 0] * tt[0] + gu_[0, 1] * tt[1]
                du1 = gu_[1, 0] * tt[0] + gu_[1, 1] * tt[1]
                dv0 = gv_[0, 0] * tt[0] + gv_[0, 1] * tt[1]
                dv1 = gv_[1, 0] * tt[0] + gv_[1, 1] * tt[1]
                return meas * (Haa00 * dv0 * du0 + Haa11 * dv1 * du1
                               + Haa01 * (dv0 * du1 + dv1 * du0)
                               + Hau0 * (dv0 * uu[0] + vv[0] * du0)
                               + Hau1 * (dv1 * uu[0] + vv[0] * du1)
                               + Huu * vv[0] * uu[0])

            Ru += rs.assemble(fb)
            Kglob = Kglob + ks.assemble(fb)
        return Ru, Kglob

    def _contact_terms(self, u, lam):
        Ru = np.zeros(self.nu_dofs)
        Kglob = sp.csr_matrix((self.nu_dofs, self.nu_dofs))
        self._contact_force = 0.0
        self._contact_active = 0
        for wall, fb in self.walls:
            x = np.asarray(fb.global_coordinates())
            uq = fb.interpolate(u)
            ur = np.asarray(uq)[0]
            uz = np.asarray(uq)[1]
            r = x[0] + ur
            z = x[1] + uz
            kn = wall.kn

            # Penalizacion  E = 0.5*kn*<g>_+^2  con derivadas ANALITICAS en
            # g: dE/du = kn*<g>_+ * dg/du ;
            #    d2E/du2 = kn*H(g)*(dg/du (x) dg/du) + kn*<g>_+ * d2g/du2.
            # dg/du se obtiene por paso complejo (exacto) y d2g/du2 por
            # diferencia central sobre la GEOMETRIA (bien escalada), nunca
            # sobre la penetracion, que puede ser ~1e-7 veces menor.
            h = 1e-20
            gr = wall.gap(r + 1j * h, z + 0j, lam).imag / h
            gz = wall.gap(r + 0j, z + 1j * h, lam).imag / h
            scale = float(np.max(np.abs(x))) or 1.0
            d = 1e-5 * scale
            grr = (wall.gap(r + d + 1j * h, z + 0j, lam).imag
                   - wall.gap(r - d + 1j * h, z + 0j, lam).imag) / h / (2 * d)
            gzz = (wall.gap(r + 0j, z + d + 1j * h, lam).imag
                   - wall.gap(r + 0j, z - d + 1j * h, lam).imag) / h / (2 * d)
            grz = (wall.gap(r + 0j, z + d + 1j * h, lam).imag
                   - wall.gap(r + 0j, z - d + 1j * h, lam).imag) * 0.0
            grz = (wall.gap(r + d + 0j, z + 1j * h, lam).imag
                   - wall.gap(r - d + 0j, z + 1j * h, lam).imag) / h / (2 * d)

            g = wall.gap(r, z, lam)
            active = g > 0.0
            gp = np.where(active, g, 0.0)
            act = active.astype(float)
            Fr = kn * gp * gr
            Fz = kn * gp * gz
            meas = 2.0 * np.pi * x[0]
            self._contact_active += int(np.sum(active))

            Hrr = kn * (act * gr * gr + gp * grr)
            Hzz = kn * (act * gz * gz + gp * gzz)
            Hrz = kn * (act * gr * gz + gp * grz)

            @LinearForm
            def rc(v, w):
                return meas * (Fr * v[0] + Fz * v[1])

            @BilinearForm
            def kc(uu, vv, w):
                return meas * (Hrr * vv[0] * uu[0]
                               + Hrz * (vv[0] * uu[1] + vv[1] * uu[0])
                               + Hzz * vv[1] * uu[1])

            Ru += rc.assemble(fb)
            Kglob = Kglob + kc.assemble(fb)
        return Ru, Kglob

    def _pressure_terms(self, u, lam):
        Ru = np.zeros(self.nu_dofs)
        Kglob = sp.csr_matrix((self.nu_dofs, self.nu_dofs))
        for value, fb in self.pressure_loads:
            pv = value * lam if np.isscalar(value) else value(lam)
            x = np.asarray(fb.global_coordinates())
            Rref = x[0]
            n2 = np.asarray(fb.normals)
            uq = fb.interpolate(u)
            gu = np.asarray(uq.grad)
            ur = np.asarray(uq)[0]
            F = axi_F(gu, ur, Rref)
            JFit = dpiola_dp(F)
            N = np.zeros((3,) + Rref.shape)
            N[0], N[1] = n2[0], n2[1]
            tN = np.einsum('ij...,j...->i...', JFit, N)
            meas = 2.0 * np.pi * Rref

            # tangente: derivada de JFit.N respecto de F (paso complejo)
            dtN = np.zeros((3, 3, 3) + Rref.shape)
            Fc = F.astype(complex)
            for (k, L) in AXI_IDX:
                Fp = Fc.copy()
                Fp[k, L] = Fp[k, L] + 1j * 1e-20
                JF = dpiola_dp(Fp)
                dtN[:, k, L] = np.einsum(
                    'ij...,j...->i...', JF, N.astype(complex)).imag / 1e-20

            @LinearForm
            def rp(v, w):
                vv3 = np.stack([v[0], v[1], np.zeros_like(v[0])])
                return meas * pv * np.einsum('i...,i...->...', tN, vv3)

            @BilinearForm
            def kp(uu, vv, w):
                Gu = axi_G(uu, Rref)
                vv3 = np.stack([vv[0], vv[1], np.zeros_like(vv[0])])
                return meas * pv * np.einsum(
                    'ikl...,i...,kl...->...', dtN, vv3, Gu)

            Ru += rp.assemble(fb)
            Kglob = Kglob + kp.assemble(fb)
        return Ru, Kglob

    # ---------------- Newton -----------------
    def residual_and_tangent(self, x, lam):
        u = x[:self.nu_dofs]
        p = x[self.nu_dofs:]
        Ru, Rp, Kuu, Kup, Kpp = self._bulk(u, p)
        for fn in (self._surface_terms, self._contact_terms,
                   self._pressure_terms):
            ru, kuu = fn(u, lam)
            Ru = Ru + ru
            Kuu = Kuu + kuu
        R = np.concatenate([Ru, Rp])
        K = sp.bmat([[Kuu, Kup], [Kup.T, Kpp]], format='csr')
        return R, K

    def constrained_dofs(self, lam):
        idx, vals = [], []
        for dofs, value in self.dirichlet:
            idx.append(dofs)
            v = value(lam) if callable(value) else value * lam
            vals.append(np.full(len(dofs), v))
        if idx:
            return np.concatenate(idx), np.concatenate(vals)
        return np.array([], dtype=int), np.array([])

    def newton(self, lam, tol=1e-9, maxit=30, verbose=False,
               rel_tol=1e-8, armijo=True):
        self._last_lam = lam
        D, Dval = self.constrained_dofs(lam)
        x = self.x.copy()
        if len(D):
            x[D] = Dval
        free = np.setdiff1d(np.arange(self.N), D)
        rep = NewtonReport()
        r0 = None
        for it in range(maxit):
            R, K = self.residual_and_tangent(x, lam)
            rn = np.linalg.norm(R[free])
            rep.residuals.append(float(rn))
            if r0 is None:
                r0 = max(rn, 1e-300)
            if verbose:
                print(f"    it {it:2d}  |R| = {rn:.6e}")
            if rn < tol or rn / r0 < rel_tol:
                rep.converged = True
                rep.iters = it
                self.x = x
                return rep
            Kff = K[free][:, free].tocsc()
            try:
                dx = spl.spsolve(Kff, -R[free])
            except Exception as exc:                      # pragma: no cover
                rep.message = f"fallo del solver lineal: {exc}"
                return rep
            if not np.all(np.isfinite(dx)):
                rep.message = "incremento no finito (matriz singular)"
                return rep
            step = 1.0
            if armijo:
                for _ in range(12):
                    xt = x.copy()
                    xt[free] += step * dx
                    Rt, _ = self.residual_and_tangent(xt, lam)
                    if np.linalg.norm(Rt[free]) < (1.0 - 1e-4 * step) * rn:
                        break
                    step *= 0.5
                else:
                    step = 1.0
            x[free] += step * dx
        rep.iters = maxit
        rep.message = "no converge en maxit"
        self.x = x
        return rep

    def solve(self, lams, verbose=False, **kw):
        """Continuacion en el parametro de carga."""
        reports = []
        for lam in lams:
            if verbose:
                print(f"  lam = {lam:.4f}")
            rep = self.newton(lam, verbose=verbose, **kw)
            reports.append(rep)
            if not rep.converged:
                break
            self.x_prev = self.x.copy()
        return reports

    # ---------------- postproceso -----------------
    @property
    def u(self):
        return self.x[:self.nu_dofs]

    @property
    def p(self):
        return self.x[self.nu_dofs:]

    def deformed_nodes(self):
        ub = self.basis.split_bases() if hasattr(self.basis, 'split_bases') \
            else None
        nodal = self.basis.nodal_dofs           # (2, nvert)
        return self.mesh.p + self.x[nodal]

    def volume(self, deformed=True):
        b = self.basis
        R = b.global_coordinates()[0]
        uq = b.interpolate(self.u)
        F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], R)
        J = det3(F) if deformed else np.ones_like(R)

        @Functional
        def vol(w):
            return 2.0 * np.pi * R * J
        return float(vol.assemble(b))

    def mean_J(self):
        b = self.basis
        R = b.global_coordinates()[0]
        uq = b.interpolate(self.u)
        F = axi_F(np.asarray(uq.grad), np.asarray(uq)[0], R)
        J = det3(F)

        @Functional
        def num(w):
            return 2.0 * np.pi * R * J

        @Functional
        def den(w):
            return 2.0 * np.pi * R
        return float(num.assemble(b) / den.assemble(b))

    def mean_pressure(self):
        pb = self.pbasis
        R = pb.global_coordinates()[0]
        pq = np.asarray(pb.interpolate(self.p))

        @Functional
        def num(w):
            return 2.0 * np.pi * R * pq

        @Functional
        def den(w):
            return 2.0 * np.pi * R
        return float(num.assemble(pb) / den.assemble(pb))

    def reaction(self, dofs_component, facets=None):
        """Fuerza de reaccion total (N) sobre un conjunto de dofs."""
        R, _ = self.residual_and_tangent(self.x, self._last_lam)
        return float(np.sum(R[dofs_component]))

    def contact_force_z(self, wall_index=0):
        """Fuerza axial total transmitida por la pared `wall_index`."""
        wall, fb = self.walls[wall_index]
        x = np.asarray(fb.global_coordinates())
        uq = fb.interpolate(self.u)
        ur = np.asarray(uq)[0]
        uz = np.asarray(uq)[1]
        lam = getattr(self, "_last_lam", 1.0)

        r = x[0] + ur
        z = x[1] + uz
        h = 1e-20
        gz = wall.gap(r + 0j, z + 1j * h, lam).imag / h
        g = wall.gap(r, z, lam)
        Fz = wall.kn * np.where(g > 0.0, g, 0.0) * gz

        @Functional
        def f(w):
            return 2.0 * np.pi * w.x[0] * Fz
        return float(f.assemble(fb))

    def contact_radius(self, wall_index=0):
        wall, fb = self.walls[wall_index]
        x = np.asarray(fb.global_coordinates())
        uq = fb.interpolate(self.u)
        r = x[0] + np.asarray(uq)[0]
        z = x[1] + np.asarray(uq)[1]
        lam = getattr(self, "_last_lam", 1.0)
        g = wall.gap(r, z, lam)
        act = g > 0.0
        return float(np.max(r[act])) if np.any(act) else 0.0


# ---------------------------------------------------------------------
#  Generadores de malla (gmsh)  --  dominios en el semiplano (R, Z)
# ---------------------------------------------------------------------

def boundary_facets_where(m, pred, facets=None):
    """Facetas de frontera cuyos DOS vertices cumplen pred(x) -> bool.
    Evita el problema de evaluar el predicado en el punto medio de la
    cuerda, que no esta sobre la frontera curva."""
    bf = m.boundary_facets() if facets is None else np.asarray(facets)
    v = m.facets[:, bf]
    return bf[pred(m.p[:, v[0]]) & pred(m.p[:, v[1]])]


def interface_facets_where(m, pred):
    """Igual pero sobre TODAS las facetas (incluye interiores)."""
    v = m.facets
    sel = pred(m.p[:, v[0]]) & pred(m.p[:, v[1]])
    return np.nonzero(sel)[0]


def _gmsh_to_skfem(dim_tags_surface_names=None):
    import gmsh
    import numpy as _np
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)
    ntags, ncoord, _ = gmsh.model.mesh.getNodes()
    order = _np.argsort(ntags)
    ntags = ntags[order]
    pts = ncoord.reshape(-1, 3)[order, :2].T
    remap = {int(t): i for i, t in enumerate(ntags)}
    tris = []
    subdom = {}
    for (dim, tag) in gmsh.model.getEntities(2):
        et, etags, enodes = gmsh.model.mesh.getElements(dim, tag)
        for e, en in zip(et, enodes):
            if e != 2:
                continue
            tt = _np.array([remap[int(k)] for k in en]).reshape(-1, 3).T
            i0 = len(tris[0][0]) if tris else 0
            start = sum(t.shape[1] for t in tris)
            tris.append(tt)
            subdom[tag] = _np.arange(start, start + tt.shape[1])
    t = _np.ascontiguousarray(_np.hstack(tris))
    pts = _np.ascontiguousarray(pts)
    return MeshTri(pts, t), subdom


def mesh_quarter_disk(R=1.0, h=0.1, h_pole=None, pole_radius=0.2):
    """Cuarto de disco en (r,z): r>=0, z>=0, de radio R.
    Refinado opcional cerca del polo (0, R) donde ocurre el contacto."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("qd")
    occ = gmsh.model.occ
    p0 = occ.addPoint(0, 0, 0, h)
    p1 = occ.addPoint(R, 0, 0, h)
    p2 = occ.addPoint(0, R, 0, h)
    l1 = occ.addLine(p0, p1)
    arc = occ.addCircleArc(p1, p0, p2)
    l2 = occ.addLine(p2, p0)
    cl = occ.addCurveLoop([l1, arc, l2])
    s = occ.addPlaneSurface([cl])
    occ.synchronize()
    if h_pole is not None:
        f1 = gmsh.model.mesh.field.add("Ball")
        gmsh.model.mesh.field.setNumber(f1, "Radius", pole_radius)
        gmsh.model.mesh.field.setNumber(f1, "XCenter", 0.0)
        gmsh.model.mesh.field.setNumber(f1, "YCenter", R)
        gmsh.model.mesh.field.setNumber(f1, "ZCenter", 0.0)
        gmsh.model.mesh.field.setNumber(f1, "VIn", h_pole)
        gmsh.model.mesh.field.setNumber(f1, "VOut", h)
        gmsh.model.mesh.field.setNumber(f1, "Thickness", pole_radius)
        gmsh.model.mesh.field.setAsBackgroundMesh(f1)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    m, _ = _gmsh_to_skfem()
    gmsh.finalize()
    return m


def mesh_cell_quarter(Rc=1.0, Rn=0.5, h=0.08, hn=0.05):
    """Cuarto de celula axisimetrica (r>=0, z>=0): citoplasma de radio Rc
    con nucleo CONCENTRICO de radio Rn.  Para problemas con plano de
    simetria z=0 (constriccion simetrica).  Devuelve (malla, subdominios)."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cellq")
    occ = gmsh.model.occ

    def quarter(R, hh):
        o = occ.addPoint(0, 0, 0, hh)
        b = occ.addPoint(R, 0, 0, hh)
        c = occ.addPoint(0, R, 0, hh)
        arc = occ.addCircleArc(b, o, c)
        l1 = occ.addLine(o, b)
        l2 = occ.addLine(c, o)
        return occ.addPlaneSurface([occ.addCurveLoop([l1, arc, l2])])

    sc = quarter(Rc, h)
    sn = quarter(Rn, hn)
    out, _ = occ.fragment([(2, sc)], [(2, sn)])
    occ.synchronize()
    tags = [t for (d, t) in out if d == 2]
    areas = {t: occ.getMass(2, t) for t in tags}
    nuc = min(areas, key=lambda t: areas[t])
    cyt = [t for t in tags if t != nuc]
    m, sub = _gmsh_to_skfem()
    gmsh.finalize()
    return m, {"nucleo": sub[nuc],
               "citoplasma": np.concatenate([sub[t] for t in cyt])}


def mesh_cell_with_nucleus(Rc=1.0, Rn=0.5, zn=0.0, h=0.08, hn=0.05,
                           half=True):
    """Celula axisimetrica: semidisco (r>=0) de radio Rc con un nucleo
    esferico de radio Rn centrado en (0, zn).  Devuelve (malla, subdominios)
    con claves 'citoplasma' y 'nucleo'."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cell")
    occ = gmsh.model.occ
    zmin = -Rc if half else 0.0

    def half_disk(R, zc, hh):
        a = occ.addPoint(0, zc - R, 0, hh)
        b = occ.addPoint(R, zc, 0, hh)
        c = occ.addPoint(0, zc + R, 0, hh)
        o = occ.addPoint(0, zc, 0, hh)
        a1 = occ.addCircleArc(a, o, b)
        a2 = occ.addCircleArc(b, o, c)
        ln = occ.addLine(c, a)
        return occ.addPlaneSurface([occ.addCurveLoop([a1, a2, ln])])

    sc = half_disk(Rc, 0.0, h)
    sn = half_disk(Rn, zn, hn)
    out, _ = occ.fragment([(2, sc)], [(2, sn)])
    occ.synchronize()
    tags = [t for (d, t) in out if d == 2]
    coms = {t: occ.getCenterOfMass(2, t) for t in tags}
    areas = {t: gmsh.model.occ.getMass(2, t) for t in tags}
    nuc = min(areas, key=lambda t: areas[t])
    cyt = [t for t in tags if t != nuc]
    m, sub = _gmsh_to_skfem()
    gmsh.finalize()
    subdomains = {"nucleo": sub[nuc],
                  "citoplasma": np.concatenate([sub[t] for t in cyt])}
    return m, subdomains
