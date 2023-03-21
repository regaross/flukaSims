# flukaSims
A FLUKA-based simulation toolkit for nEXO's Outer Detector (OD)
## Purpose
This toolkit is for simulating cosmogenic muons and the neutrons they produce in proximity to nEXO's OD. Input options may be altered to score different regions, adjust geometry and more.

## What's included?

Everything is here to run the simulation remotely on a cluster (SDF) using the singularity container built from the included `.def` file or locally (with Mac OSX).

Included are 
- A FLUKA input `.inp` file with a simple nEXO OD geometry
- The definition `.def` file for a singularity image with all dependencies 
- Run scripts to permit laziness and run simulations with few operations
- YAML files to configure the run scripts
- A directory including analysis tools for analyzing hdf5 and other output
- `.flair` files for modifying input geometry \& scoring in the Flair graphical user interface


## Getting Started with the Singularity Container

The Singularity container has within it the FLUKA binaries necessary to compile the simulation (this is done every time a user routine is changed i.e. `mgdraw.f` files or the source file: `muon_from_file.f`. It also contains some python packages in order to perform (simple) analysis in *one* job as opposed to several.

In order to build the singularity container, one needs either a Linux machine or a VM. Using MacOSX, Sylabs recommends a VM created with Vagrant and Virtualbox which is what the author used
Within a Linux machine (or a VM), get the [FLUKA binary executables](https://flukafiles.web.cern.ch/flukafiles/fluka-4-3.1/fluka_4-3.1.x86-Linux-gfor7_amd64.deb). NOTE: you will need a CERN licence to access the file.
Next, you should be able to build the singularity container using <br>
``singularity build fluka_nEXO.sif fluka_definition.def``


## The `src` directory

The `src` directory contains the `muon_functions.py` python module which is called from the `runsim` script to generate a new file of muons each time the simulation is called. There are other tools within this module the user may find useful to explore; see [this repo](https://github.com/regaross/muons/) for more about that. The generated `muon_file.txt` is also stored in this directory for later analysis. 

## The `misc_tools` directory

This directory currently contains only a python and a fortran script for interpreting the ASCII resnucle(i) scoring output from a simulation and spitting out human-readable output as a table of isotopes. More (individual) tools will be added here, but these tools will largely be incorporated into the general code. In essence, they may be used on their own, but will likely also be used in the general sim.

## The `old` directory 

The `old` directory contains several pieces that were useful in the construction of the code along the way. These are *not* used in the simulation and are deprecated. However, a future user may find something useful within them, therefore, they will stay. For instance, early versions of event storing scripts, the FLUKA original `mgdraw.f` file and a jupyter notebook used for making muon files are found here. The user may review these files to better understand each of these individual steps or to construct their own analysis code or scripts.


