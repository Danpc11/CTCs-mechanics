"""Seleccion por tamano y censo de contacto con CTC de mama de PACIENTE (no lineas celulares).

Entradas: data/poblaciones.yaml (poblacion CTC_mama_paciente y leucocitos primarios).
Salidas:  ../datos/tamano_ctc_paciente.csv  y  ../datos/censo_ctc_paciente.csv
"""
import os, sys
import numpy as np, pandas as pd
R = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(R, 'src'))
from ctcfem import poblaciones as P
d = P.cargar(os.path.join(R, 'data', 'poblaciones.yaml')); dif = d['diferencial']
OUT = os.path.join(R, 'manuscript_analysis', 'datos')
N = 400_000
leu = np.concatenate([P.muestrear(d, k, max(int(N * w), 1), 3).D.values for k, w in dif.items()])
ENV = 1.636


def ctc(Dm, cv, n=200_000, seed=1):
    x = np.random.default_rng(seed).normal(Dm, cv * Dm, n)
    return x[(x > 6) & (x < 24)]


def metricas(c):
    th = np.linspace(6, 16, 2001)
    sens = np.array([(c > t).mean() for t in th]); spec = np.array([(leu <= t).mean() for t in th])
    J = sens + spec - 1; i = int(J.argmax())
    t99 = np.quantile(c, 0.01)
    return dict(youden=J[i], D_critico_um=th[i], arrastre_99_por_millon=(leu > t99).mean() * 1e6,
                enriquecimiento_99=0.99 / max((leu > t99).mean(), 1e-12))


filas = []
for nombre, Dm in [('MCF-7 (linea, articulo)', 17.48), ('CTC paciente, via lineas', 11.8),
                   ('CTC paciente, punto medio', 11.5), ('CTC paciente, via leucocitos', 11.1)]:
    for cv in (0.10, 0.15, 0.20, 0.25):
        filas.append(dict(poblacion=nombre, D_mediana_um=Dm, cv=cv, **metricas(ctc(Dm, cv))))
T = pd.DataFrame(filas); T.to_csv(os.path.join(OUT, 'tamano_ctc_paciente.csv'), index=False)
pd.set_option('display.width', 200)
print(T.round(3).to_string(index=False))

c = P.muestrear(d, 'CTC_mama_paciente', N, 1).D.values
cen = []
for ap in (8.0, 10.0, 12.0):
    cen.append(dict(apertura_um=ap, sin_contacto=(c <= ap).mean(),
                    contacto_dentro_rango=((c > ap) & (c / ap <= ENV)).mean(), fuera_rango=(c / ap > ENV).mean(),
                    confinamiento_mediano=np.median(c) / ap,
                    leucocitos_en_contacto_por_millon=sum(dif[k] * (P.muestrear(d, k, 200_000, 1).D.values > ap).mean() for k in dif) * 1e6))
C = pd.DataFrame(cen); C.to_csv(os.path.join(OUT, 'censo_ctc_paciente.csv'), index=False)
print('\nCenso CTC de paciente (D = 11.5 um, cv 20 %):'); print(C.round(4).to_string(index=False))
