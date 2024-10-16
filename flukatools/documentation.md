# flukatools
The flukatools directory is a python module. The tools herein can be imported into other python scripts— and they are. This is basically a module to house what is needed to run the simulation without having to think too hard about it.

### `__init__.py`
This is the "init" file. It's basically the place where global variables can be assigned, or declared. I learned this the *hard* way: if you want modifications of a variable in this file from another file to be accessible via yet another file, that variable cannot be a primitive in python. It must be an object. Let me rephrase. If you have `__init__.py`, `a.py`, and `b.py` and some globally accessible variable `v` in `__init__.py`, and you want your change of that variable in `a.py` to be accessible in `b.py`, `v` must be an (pointer to an) object. ANYWAY this file contains pieces that are read in at the beginning of the simulation— system name, seed for the simulation, and so forth. These variables are used elsewhere. Basically anything that needs to be accessed by several of the files in the directory should be put here.

### `runtools.py`
This file contains commands that interface with the OS for convenient running of the simulation. Basically, it has functions that link and compile the user routines with the fluka binaries, and a command to run the simulation.

### `filemanip.py`
This file is used for file inputs and outputs. Not all file inputs and outputs are performed in this file, but most are. It generally makes sense that you'd put the I/O in one file, as opposed to doing it all over in many files. 

### `analysis.py`
This one is likely to become the largest of the files. When you create a plot that's handy, or you manipulate the data in a useful way, you should save that code somewhere other than a jupyter notebook. This file is meant to house the code that is used for analysis (obviously). 

### `muons.py`
The simulations for which this repo was made were originally (probably still are) for propagating *cosmic* muons through the nEXO configuration to count the isotopes produced as a consequence. This file is called in order to make a [phase space file](https://indico.cern.ch/event/1200922/contributions/5411811/attachments/2659968/4607569/04_Source_routine_2023_Advanced_ANL.pdf#page=19) from which FLUKA will sample the "muons".

This file uses Monte Carlo techniques to sample from various probability density functions to determine muon zenith angles and their energies.