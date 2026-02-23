# flukaSims

#### Purpose and Scope
This is a repository for the code that interfaces with [FLUKA](https://fluka.cern)— particularly, for [nEXO](https://nexo.llnl.gov)'s simulations with FLUKA. The simulations were created for the purpose of propagating muons through the nEXO configuration and counting the resultant activation products. To that end, this project has become larger and this directory contains a lot of more general tools that can be used in future simulations of nEXO using FLUKA, or maybe something else entirely. That, or this is never touched again and this documentation is entirely a matter of due diligence and will never be read. Regardless, I'll try to make it sufficiently thorough.

---

## Installation / Setup

### Python environment

```bash
pip install -r requirements.txt
```

Python 3.8+ is recommended.

### FLUKA binaries

FLUKA must be obtained separately from [fluka.cern](https://fluka.cern) (requires registration). The simulation is designed to run inside a Singularity container that bundles the FLUKA binaries — see the [Docker/Singularity documentation](./Docker_tools/documentation.md).

### Singularity / Docker

Build or pull the Singularity image using the files in `Docker_tools/`. The resulting `.sif` file path must be set in `simconfig.yaml` as `SIFPath`.

---

## Configuration

1. Copy the template configuration file:

   ```bash
   cp simconfig.yaml.template simconfig.yaml
   ```

2. Edit `simconfig.yaml`:
   - Set `SIFPath` to the **absolute path** of your flukaSims directory on the cluster (i.e. where the `.sif` container file lives).
   - Adjust `Muons` to the number of primary muons to simulate.
   - Choose `MGDrawFile` — `mgdraw_resnuc.f` scores residual nuclei; `mgdraw_neutron_count.f` counts neutrons.

---

## Running the Simulation

The simulation is submitted as a SLURM job array on the S3DF cluster:

```bash
sbatch job_script.slurm
```

Each array task runs `runsim.py`, which performs the following steps in order:

1. Reads `simconfig.yaml`.
2. Copies input files to a temporary working directory.
3. Creates the output directory.
4. Generates the muon phase-space file (Monte Carlo sampling from the Mei & Hime parameterization).
5. Copies the muon file to the output directory.
6. Updates the path to the phase-space file inside the FLUKA source routine.
7. Updates the number of muons in the `.inp` file.
8. Links and compiles the Fortran user routines into a FLUKA executable.
9. Runs the FLUKA simulation.
10. Organises the output files.

---

## Analysis Workflow

After the simulation finishes:

1. **Merge HDF5 outputs** from all job-array tasks into a single master file:

   ```python
   from flukatools import analysis as at
   import os

   data_dir = os.environ.get('FLUKA_DATA_PATH', './data')
   h5files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.endswith('.h5')]
   at.merge_hdf5_files(h5files, 'master.h5')
   ```

2. **Open the analysis notebooks** in `analysis/`:
   - `main_analysis.ipynb` — general overview plots and summary statistics.
   - `od_neutrons.ipynb` — neutron analysis for the Outer Detector.
   - `tpc_neutrons.ipynb` — neutron analysis for the TPC.

3. **Standalone plotting scripts** are available in `scripts/` — see [`scripts/README.md`](./scripts/README.md).

---

## What's here?

##### Tools for FLUKA I/O and Execution `(./flukatools)` 
When you're going to use the same tools over and over, it's best to pack 'em up into a module and iteratively improve them. The `flukatools` module contains basically everything required to run and parse through the outputs of the simulations run in this work. Check that out [here](./flukatools/documentation.md).

##### Containerization `(./Docker_tools)`
Code deployed on the [S3DF](https://s3df.slac.stanford.edu) cluster, as this was, must be inside a container. You can find the documentation for this [here](./Docker_tools/documentation.md). 

##### FLUKA Simulation Input
FLUKA is written mostly in FORTRAN 77— that is 1977. She's an old monster and takes some funky input files. These files include the file that is called the *input file* which ends in `.inp` and some other FORTRAN routines that request specific quantities to be scored in the simulation. These other files end in `.f` because they are FORTRAN files that are compiled with the rest of the simulation. You can find these files and their descriptions [here](./input/documentation.md).