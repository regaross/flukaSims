# flukaSims
A FLUKA-based simulation toolkit for nEXO's Outer Detector (OD)
## Purpose
This toolkit is for simulating cosmogenic muons and the neutrons they produce in proximity to nEXO's OD. Input options may be altered to score different regions, adjust geometry and more.

## What's included?

Everything is here to run the simulation locally (with Mac OSX) or remotely on a cluster using the singularity container built from the included `.def` file.

Included are 
- A FLUKA input `.inp` file with a simple nEXO OD geometry
- The definition `.def` file for a singularity image with all dependencies 
- Run scripts to permit laziness and run simulations with few operations
- YAML files to configure the run scripts


## Getting Started

In order to build the singularity container, one needs either a Linux machine or a VM. Using MacOSX, Sylabs recommends a VM created with Vagrant and Virtualbox which is what the author used
Within a Linux machine (or a VM), get the [FLUKA binary executables](https://flukafiles.web.cern.ch/flukafiles/fluka-4-3.1/fluka_4-3.1.x86-Linux-gfor7_amd64.deb). NOTE: you will need a CERN licence to access the file.
Next, you should be able to build the singularity container using <br>
``singularity build fluka_nEXO.sif fluka_definition.def``


# Soon to be updated
