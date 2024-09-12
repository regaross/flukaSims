#!/usr/bin/python
from os import environ, uname, getenv
from datetime import datetime
import numpy as np
from scipy.constants import physical_constants

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

#### nEXO Specific Constants ####

# # nEXO OUTER DETECTOR PARAMETERS
OD_RADIUS = 6.1722      # m
OD_HEIGHT = 12.800      # m
OD_CENTER = (0,0,0)     # m         defines coordinate system with respect to literal centre of OD
OC_RADIUS = 2.270       # m
ROI_RADIUS = OD_RADIUS + 2 # m
ROI_HEIGHT = OD_HEIGHT + 4 # m
GEN_OFFSET = OD_HEIGHT
GEN_RADIUS = np.tan(1)*(OD_HEIGHT + GEN_OFFSET) + OD_RADIUS

# OC_POSITION = (0,0,0.40) # m         positions OC with respect to OD centre
# TPC_RADIUS = 0.575      # m         from the pre-conceptual design report
# TPC_HEIGHT = 0.625      # m

# MUON FLUX PARAMETERS AT SNOLAB
SNOLAB_MU_FLUX = 3.31e-10       # \pm (0.01 (stat) \pm 0.09 (sys))e-10 mu/cm^2/s # arXiv:0902.2776v1
SNOLAB_MU_E_AVG = 363.0         # \pm 1.2 GeV # arXiv:1909.11728v1
SNOLAB_DEPTH = 5.890    #km.w.e        # \pm 94 km.w.e.  # arXiv:1909.11728v1

MIN_ENERGY = 0      #GeV
MAX_ENERGY = 25e3   #GeV





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

################################################################################
#                                                                              #
#                    DICTIONARIES- ANALYSIS & OTHER TOOLS                      #
#                                                                              #
################################################################################

# Here the FLUKA variable dictionaries are taken from entries in the FLUKA manual

# FLUKA variable information
ICODE_DICTIONARY = {
    None: None,
    -1: 'event not completed',
    0: 'normal event termination',
    4: 'stack overflow', 

    # Icode = 1x: call from Kaskad
    10: 'elastic interaction recoil', 
    11: 'inelastic interaction recoil',
    12: 'stopping particle',
    13: 'pseudo-neutron deposition',
    14: 'escape', 
    15: 'time kill',

    # Icode = 2x: call from Emfsco  
    20: 'local energy deposition (i.e. photoelectric)',
    21: 'below threshold, iarg=1',
    22: 'below threshold, iarg=2',
    23: 'escape',
    24: 'time kill',

    # Icode = 3x: call from Kasneu
    30: 'target recoil',
    31: 'below threshold',
    32: 'escape',
    33: 'time kill',

    # Icode = 4x: call from Kashea
    40: 'escape',
    41: 'time kill',
    42: 'delta ray stack overflow',

    # Icode = 5x: call from Kasoph
    50: 'optical photon absorption',
    51: 'escape',
    52: 'time kill',

    # Icode = 10x: call from Kaskad 
    100: 'elastic interaction secondaries',
    101: 'inelastic interaction secondaries',
    102: 'particle decay  secondaries',
    103: 'delta ray  generation secondaries',
    104: 'pair production secondaries',
    105: 'bremsstrahlung  secondaries',
    110: 'decay products',

    # Icode = 20x: call from Emfsco
    208: 'bremsstrahlung secondaries',
    210: 'Moller secondaries',
    212: 'Bhabha secondaries',
    214: 'in-flight annihilation secondaries',
    215: 'annihilation at rest secondaries',
    217: 'pair production secondaries',
    219: 'Compton scattering secondaries',
    221: 'photoelectric secondaries',
    225: 'Rayleigh scattering secondaries',

    # Icode = 30x: call from Kasneu
    300: 'interaction secondaries',

    # Icode = 40x: call from Kashea
    400: 'delta ray generation secondaries',
}

#### Need to add scipy into the Singularity container to use this!

jtrack_rest_energies = { # Expressed in MeV if known, 0 otherwise (for heavy ions for instance)
    None:   None,
    -6:     physical_constants['alpha particle mass energy equivalent in MeV'][0],
    -5:     physical_constants['helion mass energy equivalent in MeV'][0],
    -4:     physical_constants['triton mass energy equivalent in MeV'][0],
    -3:     physical_constants['deuteron mass energy equivalent in MeV'][0],
    -2:     0,  # Generic heavy ion (see FLUKA Manual)
    -1:     0,  # Optical photon
    0:      0,  # Pseudo-particle (see FLUKA Manual)
    1:      physical_constants['proton mass energy equivalent in MeV'][0],
    2:      physical_constants['proton mass energy equivalent in MeV'][0],
    3:      physical_constants['electron mass energy equivalent in MeV'][0],
    4:      physical_constants['electron mass energy equivalent in MeV'][0],
    5:      0,  # Electron neutrino
    6:      0,  # Electron anti-neutrino
    7:      0,  # Photon
    8:      physical_constants['neutron mass energy equivalent in MeV'][0],
    9:      physical_constants['neutron mass energy equivalent in MeV'][0],
    10:     physical_constants['muon mass energy equivalent in MeV'][0],
    11:     physical_constants['muon mass energy equivalent in MeV'][0],
    12:     497.611, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    13:     139.57039, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    14:     139.57039, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    15:     493.677, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    16:     493.677, # https://pdg.lbl.gov/2022/listings/contents_listings.html

    23:     134.9768, # https://pdg.lbl.gov/2022/listings/contents_listings.html

    208:    0,  # Energy for dose scoring (see FLUKA Manual)
    211:    0,  
    308:    0,

    # We don't care about many of these masses; for instance strange mesons, or anything else we won't be counting in the analysis.

}

JTRACK_DICTIONARY = {
    None:   None,
    -6:     'alpha, 4-HELIUM',
    -5:     'helium-3, 3-HELIUM',
    -4:     'triton, TRITON',
    -3:     'deuteron, DEUTERON',
    -2:     'generic heavy ion with Z > 2, HEAVYION',
    -1:     'optical photon, OPTIPHOT',
    0:      'ray: pseudo straight-line particle, RAY',
    1:      'proton, PROTON',
    2:      'antiproton, APROTON',
    3:      'electron, ELECTRON',
    4:      'positron, POSITRON',
    5:      'electron neutrino, NEUTRIE',
    6:      'electron antineutrino, ANEUTRIE',
    7:      'photon, PHOTON',
    8:      'neutron, NEUTRON',
    9:      'antineutron, ANEUTRON',
    10:     'positive muon, MUON+',
    11:     'negative muon, MUON-',
    12:     'Kaon-zero long, KAONLONG',
    13:     'Positive Pion, PION+',
    14:     'Negative Pion, PION-',
    15:     'Positive Kaon, KAON+',
    16:     'Negative Kaon, KAON-',

    23:     'Pion-zero, PIZERO',

    208:    'Energy Deposition',

    211:    'EM Energy Deposition',

    308:    'Low Energy Neutron Kerma'

}

JTRACK_LABELS = {
    None:   None,
    -6:     'Alpha',
    -5:     'Helium-3',
    -4:     'Triton',
    -3:     'Deuteron',
    -2:     'Generic heavy ion with Z > 2',
    -1:     'Optical Photon',
    0:      'Ray: pseudo straight-line particle',
    1:      'Proton',
    2:      'Antiproton',
    3:      'Electron',
    4:      'Positron',
    5:      'Electron Neutrino',
    6:      'Electron Antineutrino',
    7:      'Photon',
    8:      'Neutron',
    9:      'Antineutron',
    10:     'Positive Muon',
    11:     'Negative Muon',
    12:     'Kaon-zero Long',
    13:     'Positive Pion',
    14:     'Negative Pion',
    15:     'Positive Kaon',
    16:     'Negative Kaon',
    17:     'Lambda',
    18:     'Antilambda',
    19:     'Kaon zero short',
    20:     'Negative Sigma',
    21:     'Positive Sigma',
    22:     'Sigma-zero',
    23:     'Pion-zero',
    24:     'Kaon-zero',
    25:     'Antikaon-zero',
    27:     'Muon neutrino',
    28:     'Muon antineutrino',
    31:     'Antisigma-minus',
    32:     'Antisigma-zero',
    33:     'Antisigma-plus',
    34:     'Xi-zero',
    35:     'Antixi-zero',
    36:     'Negative Xi',
    37:     'Positive Xi',
    38:     'Omega-minus',
    39:     'Antiomega',
    41:     'Positive Tau',
    


    208:    'Energy Deposition',

    211:    'EM Energy Deposition',

    308:    'Low Energy Neutron Kerma'

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

NEXO_FLUKA_REGIONS = {
                1:     'tpc_in',
                2:     'tpc',
                3:     'tpc_con',
                4:     'icryo_in',
                5:     'icryo',
                6:     'ic_con',
                7:     'ocryo_in',
                8:     'ocryo',
                9:     'oc_con',
                10:    'oc_sup',
                11:    'od_in',
                12:    'cov_gas',
                13:    'od',
                14:    'cpit',
                15:    'rock',
                16:    'blkhole'
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
