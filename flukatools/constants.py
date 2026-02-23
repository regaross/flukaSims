#!/usr/bin/python

'''
constants.py

Authoritative source for physical/geometry constants and FLUKA dictionaries
shared across the flukatools module and analysis scripts.
'''

import numpy as np
from scipy.constants import physical_constants

#### nEXO Specific Constants ####

# nEXO OUTER DETECTOR PARAMETERS
OD_RADIUS = 6.1722      # m
OD_HEIGHT = 12.800      # m
OD_CENTER = (0, 0, 0)   # m  defines coordinate system w.r.t. centre of OD
OC_RADIUS = 2.270       # m
OC_POSITION = (0, 0, 0.40)  # m  positions OC with respect to OD centre
TPC_RADIUS = 0.575      # m  from the pre-conceptual design report
TPC_HEIGHT = 0.625      # m

ROI_RADIUS = OD_RADIUS + 2   # m
ROI_HEIGHT = OD_HEIGHT + 4   # m
GEN_OFFSET = 10 + 4.6 + 0.43  # How far above the top of the targeting region?
GEN_RADIUS = np.tan(1) * (ROI_HEIGHT + GEN_OFFSET) + ROI_RADIUS  # Approximately 60 m

# MUON FLUX PARAMETERS AT SNOLAB
SNOLAB_MU_FLUX = 3.31e-10    # \pm (0.01 (stat) \pm 0.09 (sys))e-10 mu/cm^2/s # arXiv:0902.2776v1
SNOLAB_MU_E_AVG = 363.0      # \pm 1.2 GeV # arXiv:1909.11728v1
SNOLAB_DEPTH = 5.890         # km.w.e  \pm 94 km.w.e.  # arXiv:1909.11728v1

MIN_ENERGY = 0       # GeV
MAX_ENERGY = 25e3    # GeV

# FLUKA variable information — from the FLUKA manual
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

NEXO_FLUKA_REGIONS = {
    1:      'tpc_in',
    2:      'tpc',
    3:      'tpc_con',
    4:      'icryo_in',
    5:      'icryo',
    6:      'ic_con',
    7:      'ocryo_in',
    8:      'ocryo',
    9:      'oc_con',
    10:     'oc_sup',
    11:     'od_in',
    12:     'cov_gas',
    13:     'od',
    14:     'cpit',
    15:     'rock',
    16:     'blkhole',
}

JTRACK_REST_ENERGIES = {  # Expressed in MeV if known, 0 otherwise (for heavy ions for instance)
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
    12:     497.611,    # https://pdg.lbl.gov/2022/listings/contents_listings.html
    13:     139.57039,  # https://pdg.lbl.gov/2022/listings/contents_listings.html
    14:     139.57039,  # https://pdg.lbl.gov/2022/listings/contents_listings.html
    15:     493.677,    # https://pdg.lbl.gov/2022/listings/contents_listings.html
    16:     493.677,    # https://pdg.lbl.gov/2022/listings/contents_listings.html
    23:     134.9768,   # https://pdg.lbl.gov/2022/listings/contents_listings.html
    208:    0,  # Energy for dose scoring (see FLUKA Manual)
    211:    0,
    308:    0,
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
    308:    'Low Energy Neutron Kerma',
}

# Lowercase aliases for backward compatibility with analysis scripts
icode_dictionary = ICODE_DICTIONARY
jtrack_dictionary = JTRACK_DICTIONARY
jtrack_rest_energies = JTRACK_REST_ENERGIES
jtrack_labels = JTRACK_LABELS
fluka_nEXO_regions = NEXO_FLUKA_REGIONS
