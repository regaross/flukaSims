# flukaSims

#### Purpose and Scope
This is a repository for the code that interfaces with [FLUKA](https://fluka.cern)— particularly, for [nEXO](https://nexo.llnl.gov)'s simulations with FLUKA. The simulations were created for the purpose of propagating muons through the nEXO configuration and counting the resultant activation products. To that end, this project has become larger and this directory contains a lot of more general tools that can be used in future simulations of nEXO using FLUKA, or maybe something else entirely. That, or this is never touched again and this documentation is entirely a matter of due diligence and will never be read. Regardless, I'll try to make it sufficiently thorough.

#### What's here?
There are a few things here. Some of which are analysis tools, some of which are absolutely necessary to compile this simulation and there are other things too.

##### Tools for FLUKA I/O and Execution `(./flukatools)` 
When you're going to use the same tools over and over, it's best to pack 'em up into a module and iteratively improve them. The `flukatools` module contains basically everything required to run and parse through the outputs of the simulations run in this work. Check that out [here](./flukatools/documentation.md).

##### Containerization `(./Docker_tools)`
Code deployed on the [S3DF](https://s3df.slac.stanford.edu) cluster, as this was, must be inside a container. You can find the documentation for this [here](./Docker_tools/documentation.md). 

##### FLUKA Simulation Input
FLUKA is written mostly in FORTRAN 77— that is 1977. She's an old monster and takes some funky input files. These files include the file that is called the *input file* which ends in `.inp` and some other FORTRAN routines that request specific quantities to be scored in the simulation. These other files end in `.f` because they are FORTRAN files that are compiled with the rest of the simulation. You can find these files and their descriptions [here](./input/documentation.md).