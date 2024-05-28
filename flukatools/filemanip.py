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
from yaml import safe_load

def read_in_config_yaml(yaml_filename : str) -> None :
    '''Reads in a particular yaml file with the appropriate parameters for the simulation'''


    # Must use the global keyword as the dictionaries are being modified.

    global PATHS
    global YAML_PARAMS
    
    with open(yaml_filename) as yaml_file:

        input_yaml = safe_load(yaml_file)

        Simulation = input_yaml.get('Simulation')
        Input = input_yaml.get('Input')

        PATHS['input'] = Input.get('InputPath')
        PATHS['workpath'] = PATHS['input'] + '.temp/'
        PATHS['workdir'] = '.temp/'

        # Simulation Parameters
        YAML_PARAMS['num_muons']     =  Simulation.get('Muons')
        YAML_PARAMS['intersecting']  =  Simulation.get('Intersecting')
        YAML_PARAMS['make_new']      =  Simulation.get('MakeNewFile')
        YAML_PARAMS['roi_radius']    =  Simulation.get('ROI_Radius')
        YAML_PARAMS['roi_height']    =  Simulation.get('ROI_Height')

        # Input Parameters
        YAML_PARAMS['input_file']        =  Input.get('InputFile')
        YAML_PARAMS['source_routine']    =  Input.get('SourceFile')
        YAML_PARAMS['mgdraw_file']       =  Input.get('MGDrawFile')

        # Source Parameters
        YAML_PARAMS['source_path'] =  input_yaml.get('Source').get('FlukaPath')

def copy_input_to_workdir() -> None:
    '''This function makes copies of the input files that are required to run each 
    simulation and puts them in a directory from which the simulation will run.
    See WORKPATH & WORKDIR in the __init__.py file'''

    # Only this dictionary is being modified.
    global FLUKA_JOB_FILES

    WORKPATH = PATHS['workpath']
    INPUT_PATH = PATHS['input']

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
            FLUKA_JOB_FILES[entry] = WORKPATH + filename

def change_number_of_muons() -> None:
    '''Changes the number of muons in the input file for a given simulation'''

    num_muons = YAML_PARAMS['num_muons']
    input_filename = FLUKA_JOB_FILES['input']

    space_string = '               '
    num_spaces = len(space_string) - len(str(num_muons))
    num_string = 'START' + space_string[:num_spaces] + str(num_muons)
    system('sed -i \'s/^START.*/' + num_string + '/\' ' + input_filename)

def change_seed() -> None:
    '''Changes the seed to the simulation in a given input file to the SEED defined by the first call to the module'''

    input_filename = FLUKA_JOB_FILES['input']
    seed = SEED
    system('echo RR: Changing seed in input file to: ' + seed)
    space_string = '                      '
    num_spaces = len(space_string) - len(str(seed))
    num_string = 'RANDOMIZ' + space_string[:num_spaces] + seed
    system('sed -i \'s/^RANDOMIZ.*/' + num_string + '/\' ' + input_filename)

def remove_leftovers() -> None:
    '''FLUKA simulations leave a bunch of useless files at the end that from compilation and running the simulation. These ones can't be read, and are effectively useless to the user. We'll remove them.'''
    WORKPATH = PATHS['workpath']

    file_list = ['*mod', '*.o', '*.exe', 'ran*']

    for ext in file_list:
        system('rm ' +  WORKPATH + ext )
        system('rm ' + ext )

def change_muon_filepath() -> None:
    '''Changes the path to the muon_file in the provided fluka source file. This function is very likely to throw an error in the future if the line number that reads the muon phase space file is ever changed!!!'''

    source_routine = FLUKA_JOB_FILES['source_routine']
    muon_file = FLUKA_JOB_FILES['muon_file']

    # Open the copied source routine file in read mode from FLUKA_JOB_FILES. Read it line by line.
    with open(source_routine, 'r') as source:
        lines = source.readlines()
        
        replace_string = '      call read_phase_space_file(\"'+ muon_file + '\", \'GeV\', \'m\', phase_space_entry, .true. , nomore )'
        lines[527] = replace_string

    # Open the copied source routine file in write mode, and replace the text within to include the appropriate path to the muon phase space file
    with open(source_routine, 'w') as source:
        source.writelines(lines)

def manage_output_files() -> None:
    '''This function is meant to move the relevant and useful output files to a particular directory.
    It should also remove the files that are no longer relevant: the copies of input files, the compiled binaries'''

    # Create an output directory

    # Probably already exists
    system('mkdir ./data/')
    
    output_dir = './data/' + TODAY + '/'
    # This one probably does not
    system('mkdir ' + output_dir)

    # Remove the temporary job specific files
    remove = ['*.o', '*ran*', '*.exe', ]
    for entry in remove:
        system('rm ' + entry)

        
    for entry in FLUKA_JOB_FILES:
        system('rm ' + FLUKA_JOB_FILES[entry])

    # The name of the copy of the input file, followed by the seed number
    input_prefix = FLUKA_JOB_FILES['input'][:-4]

    for entry in FLUKA_OUTPUT_CHANNELS:
        # Rename the fort.## files to something sensible for later parsing
        filename = input_prefix + '_001.' + entry
        new_filename = output_dir + FLUKA_OUTPUT_CHANNELS[entry] + '.asc'

    extras = ['.log', '.err', '.out']