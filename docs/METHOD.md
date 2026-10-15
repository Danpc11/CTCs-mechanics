# Method of the population and device runs

**Populations** (`data/populations.yaml`). Diameter, nucleus-to-cell ratio, apparent modulus, cortical tension,
viscosity and relaxation time per subpopulation, sampled by Latin hypercube with a seed derived from the
population name. Gs = E_app a / 2.608 and G_cyt = 0.1 Gs / a.

**Devices** (`data/materials.yaml`). Track-etched PC and PET membranes, a moulded PDMS membrane and PDMS and
PMMA channel constrictions, with an even analytic wall profile that reaches the declared aperture exactly at
the neck. The wall material enters through the geometry and a Winkler stiffness in series with the contact
penalty.

**Finite elements** (`src/ctcfem/cell_model.py`). One continuation per cell and device closes the aperture
from first contact and records the energy at every step, giving W for all apertures (5–15 µm).

**Critical pressure** (`src/ctcfem/reduction.py`). Channels: axial reduction with W indexed by the neck radius.
Membranes: work on the volume entering the pore. Leucocytes: liquid-drop bound 2γ(1/R − 1/a) as an independent check.

**Analysis** (`src/ctcfem/analysis.py`). Pressure window between the 99th percentile of the differential-weighted
leucocyte mixture and the 1st percentile of each tumour population; operating pressure, tumour yield and
leucocyte carry-over per million.
