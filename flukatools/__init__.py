#!/usr/bin/python
from os import environ, uname, getenv
from datetime import datetime
import numpy as np
from .constants import (
    OD_RADIUS, OD_HEIGHT, OD_CENTER, OC_RADIUS, OC_POSITION,
    TPC_RADIUS, TPC_HEIGHT,
    ROI_RADIUS, ROI_HEIGHT, GEN_OFFSET, GEN_RADIUS,
    SNOLAB_MU_FLUX, SNOLAB_MU_E_AVG, SNOLAB_DEPTH,
    MIN_ENERGY, MAX_ENERGY,
    ICODE_DICTIONARY, JTRACK_DICTIONARY, JTRACK_REST_ENERGIES,
    JTRACK_LABELS, NEXO_FLUKA_REGIONS,
    icode_dictionary, jtrack_dictionary, jtrack_rest_energies,
    jtrack_labels, fluka_nEXO_regions,
)


################################################################################
#                                                                              #
#               CONSTANTS AND DICTIONARIES- FILES, PATHS AND ENV               #
#                                                                              #
################################################################################
# Note: Any variable declared or assigned here will be available through each sub 
# file of this module. HOWEVER, once changed, a change to any of these variables 
# is only global if the variable is of a mutable type AND the change is not a 
# reassignment,but rather a modification like appending to an array, adding a new 
# key : val in a dictionary and so on. STRINGS ARE IMMUTABLE. Therefore, any STRING 
# must not be altered in subordinate files with the change assumed global. Wrap 
# things in dictionaries, and life will be okay.

# This is a nicely formatted date-time string that is used for labelling output directories
TODAY = datetime.now().strftime("%B%d-%Hh%M").lower()





# This is used to seed the random number generator and also to name files. Each simulation
# run will get its own seed value as defined on the system that will be appended to file
# names. It is unlikely that there will be an overlap... but it's possible

# If we're loading the module on the server— in particular S3DF, these are important
# environment variables for subsequent file manipulation.
if uname().nodename[:3] == 'sdf':

    SEED = int(getenv('FLUKA_RANDOM_SEED'))
    np.random.seed(SEED)

    # SLURM environment variables for running job arrays
    SLURM = {
    'SLURM_JOB_ID' : int(environ["SLURM_ARRAY_JOB_ID"]),
    'SLURM_TASK_ID' : int(environ["SLURM_ARRAY_TASK_ID"]),
    'SLURM_PREFIX' : 'simrun-' + str(environ["SLURM_ARRAY_JOB_ID"]) + '-' + str(environ["SLURM_ARRAY_TASK_ID"])
    }


# Declare a global dictionary to house the YAML configuration file parameters
# These are set when the YAML file is read in from within the filemanip script; 
# these are declared here to be overwritten later. Dictionaries are mutable!
YAML_PARAMS = {}

# These are set and forgotten in the filemanip functions
PATHS = {
    # System dependent (absolute path)— where is the SIF file?
    'SIF'      :   '',
    # Relative path to the simulation files
    'input' :   '',
    # Relative path to the copies of simulation files and compiled FLUKA code
    'workdir'  :   '.temp/',
    'output'   :   './data/' + TODAY + '/'
    }       

FLUKA_OUTPUT_CHANNELS = {
    # The FLUKA data emerge in fort.## files. The <##> tells us what the output files are. This will have to change to something more general.
    21  : 'resnucTPC',
    22  : 'resnucTPCCu',
    23  : 'resnucCryo',
    70  : 'neutronsOD',
    72  : 'neutronsTPC',
    96  : 'new_resnuclei_tpc'
}

FLUKA_FILES = {
    # These files will be found in the relative path PATHS['simfiles'] but may have different names...
    'input'             :   '',
    'source_routine'    :   '',
    'mgdraw'            :   '',
}

FLUKA_JOB_FILES = { # This dictionary will be populated with files that are only required for the particular run; the copies with unique names.
    'muons'             : '',
    'executable'        : '',
    'input'             : '',
    'source_routine'    : '',
    'mgdraw'            : '',
    'resnuc'            : '',
    
}

# Configuring Input and Output
HDF5_STRUCTURE = {
    ###---->  Meta data about the respective simulation
    'meta' : {   
        'seed':             {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'year':             {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'month':            {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'day':              {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'hour':             {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'minute':           {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'second':           {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'muons_simulated':  {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'roi_radius':       {'shape' : (0,), 'dtype' : float, 'maxshape': (None,)},
        'roi_height':       {'shape' : (0,), 'dtype' : float, 'maxshape': (None,)},
        },

    ###----> Data points for each neutron counted in the TPC
    'tpc_data' : {   
        'muon_energy':          {'shape' :   (0,),   'dtype' : float,    'maxshape': (None,)},
        'muon_impact':          {'shape' :   (0,),   'dtype' : float,    'maxshape': (None,)},
        'muon_initial':         {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'muon_direction':       {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'muon_pn':              {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},

        'neutron_icode':        {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'neutron_generation':   {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'neutron_region':       {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'neutron_energy':       {'shape' :   (0,),   'dtype' : float,    'maxshape': (None,)},
        'neutron_xyz':          {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'neutron_direction':    {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'neutron_parent':       {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'neutron_birth_icode':  {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        },

    ###----> Total data points for after each run tabulated in the TPC
    'tpc_totals' : {   
        'neutrons_counted':     {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'muons_simulated':      {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'muon_parents':         {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        },

    ###----> Data points for each neutron counted in the OD or somewhere within
    'od_data' : {   
        'muon_energy':          {'shape' :   (0,),   'dtype' : float,    'maxshape': (None,)},
        'muon_impact':          {'shape' :   (0,),   'dtype' : float,    'maxshape': (None,)},
        'muon_initial':         {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'muon_direction':       {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'muon_pn':              {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},

        'neutron_icode':        {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'neutron_generation':   {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'neutron_region':       {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'neutron_energy':       {'shape' :   (0,),   'dtype' : float,    'maxshape': (None,)},
        'neutron_xyz':          {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'neutron_direction':    {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'neutron_parent':       {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'neutron_birth_icode':  {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        },
    ###----> Total data points for the entire OD set of events
    'od_totals' : {   
        'neutrons_counted':     {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'muons_simulated':      {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'muon_parents':         {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        },

    'resnuclei': {
        'resnuclei':        {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},
        'resnuclei_cu':     {'shape' :   (0,3),  'dtype' : float,    'maxshape': (None,3)},}
                     
                    }

# For various plotting tools (not super important)
PARTICLE_COLOUR_DICTIONARY = {
    None:   None,
    -6:     'darkorange',   #Alpha
    -5:     'deeppink',     #Helium-3
    -4:     'maroon',       #Triton
    -3:     'darkviolet',   #Deuteron
    -2:     'blue',         #Heavy Ion
    -1:     'ivory',        #Optical P
    0:      'white',        #Ray pseudo straight-line particle
    1:      'red',          #Proton
    2:      'white',        #Antiproton
    3:      'darkcyan',     #Electron
    4:      'magenta',      #Positron
    5:      'lightgrey',    #Electron neutrino
    6:      'slategrey',    #Electron antineutrino
    7:      'gold',         #Photon
    8:      'green',        #Neutron
    9:      'white',        #Antineutron
    10:     'turquoise',    #Positive muon
    11:     'turquoise',    #Negative muon
    12:     'mistyrose',    #Kaon-zero long
    13:     'sienna',       #Positive Pion
    14:     'sienna',       #Negative Pion
    15:     'olive',        #Positive Kaon
    16:     'olive',        #Negative Kaon
    23:     'white',        #Pion Zero
    208:    'aliceblue',    #Heavy recoil
    211:    'linen',        #EM Energy deposition
    308:    'honeydew',     #Low energy neutron kerma
}

PARTICLE_COLOURS_CAPITALIZED = {
    -6:     'DarkOrange',   #Alpha
    -5:     'DeepPink',     #Helium-3
    -4:     'Maroon',       #Triton
    -3:     'DarkViolet',   #Deuteron
    -2:     'Blue',         #Heavy Ion
    -1:     'Ivory',        #Optical P
    0:      'White',        #Ray pseudo straight-line particle
    1:      'Red',          #Proton
    2:      'White',        #Antiproton
    3:      'Cyan',         #Electron
    4:      'Magenta',      #Positron
    5:      'LightGrey',    #Electron neutrino
    6:      'SlateGrey',    #Electron antineutrino
    7:      'Yellow',       #Photon
    8:      'Green',        #Neutron
    9:      'White',        #Antineutron
    10:     'Turquoise',    #Positive muon
    11:     'Turquoise',    #Negative muon
    12:     'MistyRose',    #Kaon-zero long
    13:     'Sienna',       #Positive Pion
    14:     'Sienna',       #Negative Pion
    15:     'Olive',        #Positive Kaon
    16:     'Olive',        #Negative Kaon
    23:     'White',        #Pion Zero
    208:    'AliceBlue',    #Heavy recoil
    211:    'Linen',        #EM Energy deposition
    308:    'Honeydew',     #Low energy neutron kerma
}
