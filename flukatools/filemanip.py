#!/usr/bin/python

__version__ = 0.1
__author__ = 'Regan Ross'
## Last Edited May 23, 2024

'''
filemanip.py

This part of the module is for interfacing with any of the simulation files. Copying them, renaming them, changing particular lines as required for the given simulation run. It is also used for parsing through the datafiles and output and putting them into manageable formats for analysis.

'''
################################################################################
#                                                                              #
#                                   IMPORTS                                    #
#                                                                              #
################################################################################

# Import all the constants and dictionaries from the __init__ file
from . import *

# Others (only the necessary ones)
from os import path, makedirs, system, rename
from shutil import copy
from yaml import safe_load

################################################################################
#                                                                              #
#                         SIMULATION FILE HANDLING                             #
#                                                                              #
################################################################################


def read_in_config_yaml(yaml_filename : str) -> None :
    '''Reads in a particular yaml file with the appropriate parameters for the simulation'''
    # Must use the global keyword as the dictionaries are being modified globally— not merely read.

    global PATHS
    global YAML_PARAMS
    
    with open(yaml_filename) as yaml_file:

        input_yaml = safe_load(yaml_file)

        Simulation = input_yaml.get('Simulation')
        Input = input_yaml.get('Input')

        PATHS['SIF'] = Input.get('SIFPath')
        PATHS['input'] = Input.get('InputPath')

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


    FLUKA_FILES['input'] = YAML_PARAMS['input_file']
    FLUKA_FILES['source_routine'] = YAML_PARAMS['source_routine']
    FLUKA_FILES['mgdraw'] = YAML_PARAMS['mgdraw_file']

def copy_input_to_workdir() -> None:
    '''This function makes copies of the input files that are required to run each 
    simulation and puts them in a directory from which the simulation will run.
    See WORKPATH & WORKDIR in the __init__.py file'''

    # Only this dictionary is being modified.
    global FLUKA_JOB_FILES

    WORKPATH = PATHS['SIF'] + PATHS['workdir']
    INPUT_PATH = PATHS['SIF'] + PATHS['input']

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
            new_filename = base + str(SEED) + ext

            rename(WORKPATH + filename, WORKPATH + new_filename)

            # Change the global variable for this filename
            FLUKA_JOB_FILES[entry] = WORKPATH + new_filename

def change_number_of_muons() -> None:
    '''Changes the number of muons in the input file for a given simulation'''
    
    # Read the number of muons as specified in the yaml file
    num_muons = YAML_PARAMS['num_muons']
    input_filename = FLUKA_JOB_FILES['input']

    # THESE CANNOT be changed. The FLUKA input cards require a very specific format...
    space_string = '               '
    num_spaces = len(space_string) - len(str(num_muons))
    num_string = 'START' + space_string[:num_spaces] + str(num_muons)

    # Parse the file with sed inline and change the required line.
    system('sed -i \'s/^START.*/' + num_string + '/\' ' + input_filename)

def change_seed() -> None:
    '''Changes the seed to the simulation in a given input file to the SEED defined by the first call to the module'''

    input_filename = FLUKA_JOB_FILES['input']
    system('echo RR: Changing seed in input file to: ' + str(SEED))
    space_string = '                      '
    num_spaces = len(space_string) - len(str(SEED))
    num_string = 'RANDOMIZ' + space_string[:num_spaces] + str(SEED)
    system('sed -i \'s/^RANDOMIZ.*/' + num_string + '/\' ' + input_filename)

def remove_leftovers() -> None:
    '''FLUKA simulations leave a bunch of useless files at the end that from compilation and running the simulation. These ones can't be read, and are effectively useless to the user. We'll remove them.'''
    WORKPATH = PATHS['SIF'] + PATHS['workdir']

    file_list = ['*mod', '*.o', '*.exe', 'ran*']

    for ext in file_list:
        system('rm ' +  WORKPATH + ext )
        system('rm ' + ext )

def change_muon_filepath() -> None:
    '''Changes the path to the muon_file in the provided fluka source file. This function is very likely to throw an error in the future if the line number that reads the muon phase space file is ever changed!!!'''

    source_routine = FLUKA_JOB_FILES['source_routine']
    muon_file = FLUKA_JOB_FILES['muons']

    # Open the copied source routine file in read mode from FLUKA_JOB_FILES. Read it line by line.
    with open(source_routine, 'r') as source:
        lines = source.readlines()
        
        # Again, this string CANNOT be modified. FORTRAN is very touchy about spaces
        replace_string = '      call read_phase_space_file(\"'+ muon_file + '\", \'GeV\', \'m\', phase_space_entry, .true. , nomore )'
        lines[527] = replace_string

    # Open the copied source routine file in write mode, and replace the text within to include the appropriate path to the muon phase space file
    with open(source_routine, 'w') as source:
        source.writelines(lines)

def make_output_directory() -> None:
    '''Produces a subdirectory in ./data/ that is labelled with the date of the simulation. Makes copies of the input file to put there.'''

    output_dir = PATHS['output']
    system('mkdir data')
    system('mkdir ' + output_dir)
    system('cp ' + FLUKA_FILES['input']  + ' ' + output_dir + FLUKA_FILES['input'])

def manage_output_files() -> None:
    '''This function is meant to move the relevant and useful output files to a particular directory.
    It should also remove the files that are no longer relevant: the copies of input files, the compiled binaries
    In particular it should:
        1. put the appropriate FLUKA output files in PATHS['output'] (those that are named fort##)
        2. copy the muon.txt file from .temp/ into PATHS['output']
        3. produce a concatenated std error & std output & logfile for meta data in PATHS['outout']
        4. remove files that are not needed (don't directly delete .temp/ as it may be in use)
        '''
    
    output = PATHS['output']
    input_prefix = FLUKA_FILES['input'][:-4] + str(SEED)

    # Step ONE 
    # Move the output data files to the proper directory
    for entry in FLUKA_OUTPUT_CHANNELS:
        # Rename the fort.## files to something sensible for later parsing
        filename = input_prefix + '001_fort.' + str(entry)
        new_filename = output + FLUKA_OUTPUT_CHANNELS[entry] + str(SEED) + '.asc'
        system('mv ' + PATHS['SIF'] + filename + ' ' + new_filename)
    
    # Step TWO
    # Copy the muon file from .temp 
    system('mv ' + FLUKA_JOB_FILES['muons'] + ' ' + output + 'muons' + str(SEED) + '.txt')

    # Step THREE
    # Produce a concatenated std error, std output and log files
    slurm_output_filename = 'slurm_outlogerr' + str(SEED) + '.txt'
    fluka_output_filename = 'fluka_outlogerr' + str(SEED) + '.txt'

    # Copy the SLURM simulation output files 
    system('cat ' + SLURM['SLURM_PREFIX'] + '.out >> ' + output + slurm_output_filename + ' && rm ' + PATHS['SIF'] + SLURM['SLURM_PREFIX'] + '.out')
    system('cat ' + SLURM['SLURM_PREFIX'] + '.log >> ' + output + slurm_output_filename + ' && rm ' + PATHS['SIF'] + SLURM['SLURM_PREFIX'] + '.log')
    system('cat ' + SLURM['SLURM_PREFIX'] + '.err >> ' + output + slurm_output_filename + ' && rm ' + PATHS['SIF'] + SLURM['SLURM_PREFIX'] + '.err')

    # Copy the FLUKA simulation output files
    fluka_output_filename = 'foutlogerr' + str(SEED)
    system('cat ' + input_prefix + '001.out >> ' + output + fluka_output_filename + ' && rm ' + PATHS['SIF'] + input_prefix + '001.out')
    system('cat ' + input_prefix + '001.log >> ' + output + fluka_output_filename + ' && rm ' + PATHS['SIF'] + input_prefix + '001.log')
    system('cat ' + input_prefix + '001.err >> ' + output + fluka_output_filename + ' && rm ' + PATHS['SIF'] + input_prefix + '001.err')

    # Step FOUR: remove what remains
    system('rm .temp/*' + str(SEED) + '* *' + str(SEED) + '*')

################################################################################
#                                                                              #
#                      DATA FILE HANDLING FOR ANALYSIS                         #
#                                                                              #
################################################################################

'''Native FLUKA scoring functions produce an ASCII file with a header of 14 lines
followed by an array of dimension n x 10; that is, there are 10 columns and n rows.'''

def read_resnuclei_file(filepath) -> dict:
    ''' One FLUKA scoring card is called RESNUCLEI for scoring residual nuclei.
    As with all native FLUKA scoring outputs, there is a header for the file
    that can be parsed for meta data, and a subsequent 2D array of numbers 
    representing the number of nuclei counted in a given region. The dimensions
    of this array are not straightforward and obvious. Edit this function with 
    caution.
    '''
    raw =  np.loadtxt(filepath, skiprows = 14) 

def read_muon_phase_space_file(filepath)-> dict:
    ''' The customized muon source for these FLUKA simulations requires producing 
    a 'phase space' file with information containing the muon energies, direction
    cosines and so forth. As the simulation currently works, the phase space file
    is saved along with the simulation output as a '.txt' file. This way, the muon
    file isn't just saved in volatile memory and can be stored and, if necessary, 
    re-used. This function will read in that file and return the data.'''
    pass

def read_neutron_file(filepath)-> dict:
    ''' The customized mgdraw output deployed to score neutrons in various regions
    produces a file that is effectively a 2D array of numbers providing information
    about the neutrons as they are scored by FLUKA. The parameters and their order 
    are defined in the mgdraw file. This function will read in that _fort.## file and
    return the relevant data.'''
    pass