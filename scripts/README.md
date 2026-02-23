# scripts/

Standalone one-off plotting scripts for visualising simulation output.

## `plot_tree.py`

Reads a FLUKA mgdraw output file (e.g. `nEXO_OD001_fort.72`) and renders an
interactive 3-D event-tree plot using Plotly.

**Input:** a FLUKA fort output file (text format written by the mgdraw routine).

**Usage:**
```bash
python scripts/plot_tree.py
```

Edit the filename near the top of the script to point at your output file.

---

## `plot_with_TPC.py`

Plots the 3-D positions of neutron hits scored in the TPC alongside the TPC
cylinder geometry.

**Input:** `od_master.h5` — a merged HDF5 file produced by
`flukatools.analysis.merge_hdf5_files()`.

**Usage:**
```bash
python scripts/plot_with_TPC.py
```

---

## `plot_with_OC.py`

Plots the 3-D positions of neutron hits alongside the outer cryostat sphere
geometry.

**Input:** `od_master.h5` — a merged HDF5 file produced by
`flukatools.analysis.merge_hdf5_files()`.

**Usage:**
```bash
python scripts/plot_with_OC.py
```
