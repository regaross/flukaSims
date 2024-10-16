# input
These are the files that are input to the simulation. The `.f` ones are Fortran files that must be linked and compiled with the other FLUKA files before run time. There's another file which is the actual *input file* to the simulation. It has a suffix `.inp`. 

### `nEXO_2024.inp`
This is the input file to the simulation. It contains the nEXO configuration, physics options, material properties, and more. If a future user wishes to simulated something else with respect to nEXO, this is probably the first file to look at. Configuration changes will be made here, builtin scoring options are enabled/disabled herein.

### `muon_from_file.f`
This is the source file. This file points to a "phase space file" which is basically a 2D array in a text file. It then sequentially selects muons from it (row by row) to feed the simulation. These phase space files are created by [muons.py](../flukatools/muons.py).

### `mgdraw_neutron_count.f`
This is a custom scoring file. It uses a few routines to print to an output file attributes of neutrons scored in particular regions. It's not currently in use. Neutrons have the FLUKA variable `JTRACK = 8`. Regions are numbered sequentially in the input file. So there are routines in this file that score neutrons in the TPC, and anywhere inside the nEXO volume (region number is smaller than the number of the volume outside the OD). 

### `mgdraw_resnuc.f`
This is a user routine that scores activation products or residual nuclei (but not really). It works in conjunction with the next file `usrrnc.f` to print out specific kinds of events— particularly inelastic scattering events. Those that produce "residual nuclei" are handled with...

### `usrrnc.f`
This is the FLUKA routine that handles so-called "residual nuclei". This is how we are able to score activation and relate activated nuclei to their parent muon. This routine is what is called by the residual nuclei card in the input files as well. The data are almost the same; this custom routine can give information on an event by event basis, whereas the builtin residual nuclei scoring bins things per run.