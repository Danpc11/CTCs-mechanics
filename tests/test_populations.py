import numpy as np, pathlib
from ctcfem import populations as P
DATA = P.load(pathlib.Path(__file__).resolve().parents[1] / "data" / "populations.yaml")


def test_reproducible_and_seed_dependent():
    a, b = P.sample(DATA, "neutrophil", 32, 7), P.sample(DATA, "neutrophil", 32, 7)
    assert np.array_equal(a.D.values, b.D.values)
    assert not np.array_equal(a.D.values, P.sample(DATA, "neutrophil", 32, 8).D.values)


def test_means_and_truncation():
    for nm, pop in DATA["populations"].items():
        df = P.sample(DATA, nm, 4000, 1); lo, hi = pop["D"]["trunc"]
        assert df.D.between(lo, hi).all()
        assert abs(df.D.mean() - pop["D"]["mean"]) / pop["D"]["mean"] < 0.03


def test_cortex_calibration_matches_manuscript():
    assert abs(2.65 * 8.74 / DATA["model"]["K_E"] - 8.881) < 1e-3     # Gs of MCF-7, mN/m


def test_every_parameter_has_status_and_source():
    t = P.provenance_table(DATA)
    assert (t.status.str.len() > 0).all() and (t.source.str.len() > 0).all()
