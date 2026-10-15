import numpy as np, pathlib
from ctcfem import poblaciones as P
DATOS = P.cargar(pathlib.Path(__file__).resolve().parents[1] / "data" / "poblaciones.yaml")

def test_reproducible_e_independiente_del_orden():
    a = P.muestrear(DATOS, "neutrofilo", 32, 7)
    b = P.muestrear(DATOS, "neutrofilo", 32, 7)
    assert np.array_equal(a.D.values, b.D.values)
    c = P.muestrear(DATOS, "neutrofilo", 32, 8)
    assert not np.array_equal(a.D.values, c.D.values)

def test_medias_y_truncamiento():
    for nm, pob in DATOS["poblaciones"].items():
        df = P.muestrear(DATOS, nm, 4000, 1)
        lo, hi = pob["D"]["trunc"]
        assert df.D.between(lo, hi).all()
        assert abs(df.D.mean() - pob["D"]["media"]) / pob["D"]["media"] < 0.03

def test_calibracion_cortex_reproduce_manuscrito():
    df = P.muestrear(DATOS, "MCF-7", 1, 1)
    # con E_app = 2.65 kPa y a = 8.74 um el manuscrito da Gs = 8.881 mN/m
    assert abs(2.65 * 8.74 / DATOS["modelo"]["K_E"] - 8.881) < 1e-3

def test_toda_fila_tiene_estado_y_fuente():
    t = P.tabla_de_estados(DATOS)
    assert (t.estado.str.len() > 0).all() and (t.fuente.str.len() > 0).all()
