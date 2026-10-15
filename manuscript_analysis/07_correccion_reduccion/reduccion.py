"""
reduccion.py -- reduccion CORREGIDA del paisaje FEM W a la presion critica de paso.

Corrige tres defectos de fase2_lib.presion_critica / paso123_analisis.splines:

1. MAPEO (el importante).  El FEM devuelve W como funcion del radio NOMINAL de
   garganta R_t, pero presion_critica lo evalua en el radio GEOMETRICO local
   R_w(z).  En el cuello R_w(0) = R_min = R_t + g0 (R_wide - R_t) > R_t, asi que
   la aproximacion declarada W(s) ~ W_FEM(R_w(s)) no reproduce W_FEM(R_t) ni
   siquiera en z = 0, el unico punto donde es exacta por construccion.
   Aqui W se indexa por el radio REAL del cuello del estado simulado, R_min(R_t),
   y se ancla en W = 0 en el contacto (R_min = a).

2. CERO ARTIFICIAL a 2a/D = 1.20.  La banda valida se definia en radio nominal;
   a 2a/D = 1.20 solo el punto z = 0 caia dentro y alli np.sign(0) = 0 anula la
   pendiente: dP = 0 para CUALQUIER W, aunque el confinamiento efectivo 1.039
   implica contacto.  Aqui z > 0 estricto y la banda se define en confinamiento
   EFECTIVO: 1 <= a/R_w <= techo validado.

3. QUIEBRE DEL PERFIL.  tanh(|z|) tiene pendiente discontinua en el cuello y el
   maximo cae siempre en el primer punto de la rejilla (de ahi la "convergencia
   en nz").  Se ofrece un perfil par y suave con el mismo radio de cuello, para
   medir la sensibilidad a la forma del cuello.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

CONF_EF_MAX = 1.636          # techo validado, confinamiento efectivo (2a/D = 2.43 nominal)


def g0(zt, Ltap):
    return 0.5 * (1.0 + np.tanh(-zt / Ltap))


def perfil_quiebre(z, Rt, Rwide, zt, Ltap):
    """Perfil de la Fase 1 (el que simula el FEM): quiebre en z = 0."""
    return Rt + (Rwide - Rt) * 0.5 * (1.0 + np.tanh((np.abs(z) - zt) / Ltap))


def perfil_suave(z, Rt, Rwide, zt, Ltap):
    """Perfil par con derivada nula en z = 0 y la misma forma lejana."""
    return Rt + (Rwide - Rt) * (1.0 + 0.5 * (np.tanh((z - zt) / Ltap) - np.tanh((z + zt) / Ltap)))


def R_cuello(Rt, Rwide, zt, Ltap):
    return Rt + g0(zt, Ltap) * (Rwide - Rt)


def spline_por_cuello(Rt, W, a, Rwide, zt, Ltap):
    """W(rho) con rho = radio REAL del cuello del estado FEM; W(a) = 0 (contacto)."""
    rho = np.append(R_cuello(np.asarray(Rt, float), Rwide, zt, Ltap), a)
    Wv = np.append(np.asarray(W, float), 0.0)
    o = np.argsort(rho)
    return PchipInterpolator(rho[o], np.maximum(Wv[o], 0.0), extrapolate=False)


def spline_original(Rt, W, a):
    """Reproduce paso123_analisis.splines (indexado por R_t nominal, ancla 1.02a)."""
    R = np.append(np.asarray(Rt, float), 1.02 * a)
    Wv = np.append(np.asarray(W, float), 0.0)
    o = np.argsort(R)
    return PchipInterpolator(R[o], np.maximum(Wv[o], 0.0), extrapolate=False)


def presion_critica(a, Rt, Rwide, zt, Ltap, spW, perfil="quiebre", nz=200001,
                    conf_ef_max=CONF_EF_MAX):
    """Presion critica (kPa si W en mN/m*um^2 y longitudes en um).

    perfil = "quiebre": el dispositivo simulado (R_t es el nominal).
    perfil = "suave":   cuello suave con el MISMO radio de cuello que el de quiebre.
    """
    f = perfil_quiebre
    Rt_use = Rt
    if perfil == "suave":
        rmin = R_cuello(Rt, Rwide, zt, Ltap)
        Rt_use = brentq(lambda x: perfil_suave(0.0, x, Rwide, zt, Ltap) - rmin, 1e-3, rmin)
        f = perfil_suave
    z = np.linspace(0.0, zt + 6.0 * Ltap, nz)[1:]            # z > 0 estricto
    Rw = f(z, Rt_use, Rwide, zt, Ltap)
    dR = np.abs(np.gradient(Rw, z))
    ok = (Rw <= a) & (a / Rw <= conf_ef_max + 1e-9)          # contacto y dentro del techo
    q = np.abs(spW.derivative()(Rw)) * dR / (np.pi * Rw ** 2)
    q = np.where(ok & np.isfinite(q), q, -np.inf)
    if not np.isfinite(q.max()):
        return dict(dP_crit=0.0, z_crit=np.nan, conf_crit=np.nan)
    i = int(np.argmax(q))
    return dict(dP_crit=float(q[i]), z_crit=float(z[i]), conf_crit=float(a / Rw[i]),
                R_crit=float(Rw[i]))
