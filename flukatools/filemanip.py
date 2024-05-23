#!/usr/bin/python

__version__ = 0.0
__author__ = 'Regan Ross'
## Last Edited May 23, 2024

'''
filemanip.py

This part of the module is for interfacing with any of the simulation files. Copying them, renaming them, changing particular lines as required for the given simulation run.

'''

# Import all the constants and dictionaries from the __init__ file
from . import *

# Others (only the necessary ones)
from os import path, makedirs, system
from shutil import copy

def copy_input_to_workdir() -> None:
    '''This function makes copies of the input files that are required to run each 
    simulation and puts them in a directory from which the simulation will run.
    See WORKPATH & WORKDIR in the __init__.py file'''

    # Create the working directory
    if not path.exists(WORKPATH):
        # If it doesn't exist, create it
        makedirs(WORKPATH)

    for entry in FLUKA_FILES:
        filename = FLUKA_FILES[entry]

        # If the file exists, copy it with a new name- the ones that don't yet exist are created elsewhere (compilation, etc...)
        if path.exists(INPUT_PATH + filename):
            # Copy the file to the working directory
            copy(INPUT_PATH + filename, WORKPATH + filename)

            # Rename file with the SEED number created at first module call
            ## Separate the filename into the base and the extension
            base, ext = path.splitext(filename)
            ## Add the SEED string between the two
            filename = base + str(SEED) + ext

            # Change the global variable for this filename
            FLUKA_FILES[entry] = WORKPATH + filename

def change_number_of_muons() -> None:
    '''Changes the number of muons in the input file for a given simulation'''

    num_muons = YAML_PARAMS['num_muons']
    input_filename = FLUKA_FILES['input']

    space_string = '               '
    num_spaces = len(space_string) - len(str(num_muons))
    num_string = 'START' + space_string[:num_spaces] + str(num_muons)
    system('sed -i \'s/^START.*/' + num_string + '/\' ' + input_filename)

def change_seed() -> None:
    '''Changes the seed to the simulation in a given input file to the SEED defined by the first call to the module'''

    input_filename = FLUKA_FILES['input']
    seed = SEED
    system('echo RR: Changing seed in input file to: ' + seed)
    space_string = '                      '
    num_spaces = len(space_string) - len(str(seed))
    num_string = 'RANDOMIZ' + space_string[:num_spaces] + seed
    system('sed -i \'s/^RANDOMIZ.*/' + num_string + '/\' ' + input_filename)

def remove_leftovers() -> None:
    '''FLUKA simulations leave a bunch of useless files at the end that from compilation and running the simulation. These ones can't be read, and are effectively useless to the user. We'll remove them.'''

    file_list = ['*mod', '*.o', '*.exe', 'ran*']

    for ext in file_list:
        system('rm ' +  WORKPATH + ext )
        system('rm ' + ext )
