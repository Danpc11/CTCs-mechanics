# Auditoría de la Fase 1 — dos hallazgos que obligan a decidir antes de seguir

## 1. El solver pasó, y con holgura

Cuatro soluciones analíticas reproducidas, error máximo 3,7e-3; tangente consistente verificado
contra diferencias finitas a 2,1e-10; convergencia de Newton cuadrática. El punto de corte del plan
queda superado sin reservas. Detalle en `validacion_solver.md`.

Pero la calibración contra el 3D de Wang et al. **no se pudo hacer como estaba planeada**, y la
razón importa: su caso tiene confinamiento 2a/l = 0,567, es decir la célula **no toca la pared** y la
deformación se la impone el campo de esfuerzos del fluido. Un solver cuasi-estático de sólido sin
acoplamiento fluido no puede generar eso. Se calibró en su lugar la relación hiperelástica entre
invariantes de forma, que es independiente del modo de carga, con una discrepancia de +16% —
provisional, porque el valor 3D de referencia se extrajo gráficamente de una figura (rango real
-1% a +42%).

**Consecuencia de diseño:** el modelo solo es válido para 2a/D >= 1,2, o sea en régimen de contacto
con la pared. Constricción y microfiltro operan ahí; **DLD no**. La vía DLD del plan, tal como estaba
escrita, no es ejecutable con este solver.

## 2. Los dos ejes del mapa de fases son degenerados

Este es el hallazgo grave. De las 86 filas de la tabla de parámetros, 53 son transferidas, y la
transferencia se concentra exactamente donde el artículo prometía discriminar:

- **Relación núcleo/célula (a_n/a): 0,50 para las cinco estirpes**, incluido el leucocito, todas
  transferidas del ajuste de Wang et al. a HL-60/PC-3/K-562. Cero fuentes publicadas para líneas de
  mama, ni en PubMed ni en OpenAlex, ni en el barrido dedicado posterior (el dato de tamaño nuclear
  existe en la literatura de imagen pero vive en figuras, no en texto extraíble).
- **Todas las propiedades nucleares: idénticas para las cinco estirpes**, todas transferidas de las
  mismas fuentes.
- **Membrana/córtex de cualquier línea de mama: cero fuentes.** Todo transferido de PC-3.
- **SK-BR-3 (HER2+): cero fuentes de citoplasma o núcleo** pese a búsqueda dedicada sobre 123
  registros.

Dicho sin rodeos: en los ejes planeados (rigidez nuclear × N/C) las cinco estirpes caen en **un solo
punto**, con un único valor distinto en cada eje. Un mapa de fases sobre esos ejes no mediría nada;
sería una figura que aparenta información y no la tiene.

## 3. Lo que sí está medido y sí discrimina

| Propiedad | MCF-10A | MCF-7 | SK-BR-3 | MDA-MB-231 | leucocito |
|---|---|---|---|---|---|
| Módulo aparente en suspensión (kPa) | — | 2,65 | 2,77 | 1,53 | 0,40–0,62 |
| Tiempo de relajación citoplasmático τ (s) | 0,77 | 3,23 | — | 0,57 | — |
| Viscosidad estática citoplasmática (Pa·s) | 41,6 | 56,4 | — | 10,7 | — |
| Módulo elástico citoplasmático (Pa) | 79,3 | 32,9 | — | 38,6 | — |

Dos cosas destacan. Primero, **las células tumorales en suspensión son 2,5–7× MÁS RÍGIDAS que los
leucocitos** (pinzas acústicas, 10.1038/s41378-020-0150-6), lo contrario del relato habitual de
"la célula cancerosa es más blanda" — que procede de medidas en célula adherida. Segundo, y es el
dato aprovechable: **τ difiere 5,7× entre luminal (MCF-7, 3,23 s) y basal/mesenquimal
(MDA-MB-231, 0,57 s)** (microrreología rotacional, 10.1039/d4na00003j).

## 4. Resolución del conflicto de viscosidad citoplasmática

El barrido marcó como incompatibles 32,5 mPa·s (Wang et al.) y 10–60 Pa·s (microrreología). No son
la misma magnitud: el primero es una viscosidad de gota continua dentro de una formulación fluida
IB-LBM, el segundo es microrreología de la red citoesquelética. Para un modelo de sólido
viscoelástico los valores correctos son los de microrreología, y la comprobación de consistencia
interna lo respalda por dos vías (`consistencia_viscoelastica.csv`):

| línea | η (Pa·s) | E (Pa) | η/E (s) | τ medido (s) | (η/E)/τ |
|---|---|---|---|---|---|
| MDA-MB-231 | 10,7 | 38,6 | 0,28 | 0,57 | 0,49 |
| MCF-10A | 41,6 | 79,3 | 0,52 | 0,77 | 0,68 |
| MCF-7 | 56,4 | 32,9 | 1,71 | 3,23 | 0,53 |

Primero, **η/E reproduce exactamente el orden de los τ medidos** (MDA-MB-231 < MCF-10A < MCF-7), y el
cociente (η/E)/τ se mantiene en 0,49–0,68 con media 0,57 en las tres líneas: un factor constante,
que es lo que cabe esperar de una diferencia de convención entre la definición de Maxwell y la de
Kelvin-Voigt, no un desacuerdo entre magnitudes. Segundo, la viscosidad de Wang con esos mismos
módulos daría τ entre 4,1e-4 y 9,9e-4 s, unas **2.000 veces por debajo** de los τ medidos: queda
descartada como viscosidad de red para este modelo. El valor de Wang se usará solo al reproducir su
caso dentro de su propia formulación. Hueco cerrado.

## 5. La decisión que hay que tomar

La afirmación central aprobada — "el fenotipo mecánico por compartimento determina qué modalidad
captura, y la óptima cambia con el subtipo", sobre ejes de rigidez nuclear y N/C — **no es
sostenible con los datos que existen**. Hay que reformularla o hay que generar el dato que falta.
