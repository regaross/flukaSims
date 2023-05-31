# flukaSims

A FLUKA-based simulation toolkit for nEXO. This repo contains all the accessory code that with the FLUKA binaries, composes a near-complete toolkit.

## Purpose

This toolkit is for simulating cosmogenic muons and the neutrons they produce in proximity to nEXO's OD. Input options may be altered to score different regions, adjust geometry and more. Any FLUKA features may be enabled by subsequent users to deploy the full force of FLUKA with respect to nEXO.

## What's included?

Almost everything is here to run the simulation remotely on a cluster (SDF) using the singularity container built from the included `.def` file.

Included are 
- A FLUKA input `.inp` file with a simple nEXO OD geometry
- The definition `.def` file for a singularity image with all dependencies 
- Run scripts to permit laziness and run simulations with few operations using SLURM
- YAML files to configure the simulation run
- A directory including analysis tools for analyzing hdf5 and other output

## What's not included?

The FLUKA source code or FLUKA binaries. These require particular permissions from CERN. Users without CERN accounts can easily obtain individual accounts and thus permissions to access the FLUKA binaries (but not the source code). The source is not needed, only the tar ball with the binaries are necessary. Once licensed, users may obtain these files by making an account on the [CERN website](https://auth.cern.ch/auth/realms/cern/protocol/openid-connect/auth?client_id=webframeworks-drupal-fluka2019&response_type=code&scope=openid%20email&redirect_uri=https%3A//fluka.cern/openid-connect/generic&state=TDYggvaxKQy-JGzvCjrm7N__iwLy9pJUqkhzHcoGi0s).


## Getting Started with the Singularity Container

The Singularity container will have within it the FLUKA binaries necessary to compile the simulation (this is done every time a user routine is changed i.e. `mgdraw.f` files or the source file: `muon_from_file.f`. It also contains some python packages in order to perform (simple) data management after the simulation has been run. 

In order to build the singularity container, one needs either a Linux machine or a VM. Using MacOSX, Sylabs recommends a VM created with Vagrant and Virtualbox which is what the author used.

Within a Linux machine (or a VM), get the [FLUKA binary executables](https://flukafiles.web.cern.ch/flukafiles/fluka-4-3.1/fluka_4-3.1.x86-Linux-gfor7_amd64.deb). NOTE: as mentioned, you will need a CERN licence to access the file.

Next, you should be able to build the singularity container using:
<br>
``singularity build fluka_nEXO.sif simfiles/fluka_definition.def``


## The `simfiles` directory

### `muon_functions.py`

The `simfiles` directory contains the `muon_functions.py` python module which is called from the `run_simulation.py` script to generate a new file of muons each time the simulation is called. There are other tools within this module the user may find useful to explore; see [this repo](https://github.com/regaross/muons/) for more about that. This module is used to produce a `.txt` phase space file of muons for each simulation.

### `mgdraw.f`

The `mgdraw.f` file is a routine called by FLUKA for each event in the stack. In its current state, it is used to print out tabulated data about neutrons, their parent muons, and their attributes when they're counted in the OD and or the TPC. This file is used for advanced customizable scoring for specific types of simulations. It allows users to access all data from the simulation (if they know how).

### `muon_from_file.f`

The `muon_from_file.f` file is the file used to provide FLUKA a custom particle source. While the `.inp` file can specify a mono-energetic beam in a single direction, the muon spectra at SNOLAB is much more complicated. This particular file calls one function that reads in the muons from the file produced by `muon_functions.py`.

### `nEXO_OD.inp`

This file contains the entirety of the physical configuration of the simulation as well as a few builtin scoring cards currently used to score isotopic activation in various regions. It is written in a very strict format with character limits, specific spacing, and ought to be edited with the [flair](https://flair.web.cern.ch/flair/). Flair is a free and open source input editor for FLUKA. It is recommended that if you must edit the input cards, you do so via Flair which takes care of the input formats and contains information for each type of FLUKA card.

### `nEXO_FLUKA.def`

This is the definition file used to build the singularity container. It contains tells singularity what pieces to put in the container so that the simulation will run properly. Note that you will need the FLUKA binaries, `fluka_4-3.2.x86-Linux-gfor9_amd64.deb`, and the ENDF cross section database configured for FLUKA in the build directory, [fluka-pw-endf-viii0_1-0_all.deb](https://flukafiles.web.cern.ch/flukafiles/neutron/fluka-pw-endf-viii0_1-0_all.deb).

## The `simtools` directory

The simtools directory contains tools for running and analyzing the simulation data (including muon_functions.py).
