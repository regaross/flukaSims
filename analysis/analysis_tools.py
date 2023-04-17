#!/usr/bin/python

__version__ = 2.2
__author__ = 'Regan Ross'
## Last Edited March 8, 2023

'''
Contact:
Regan Ross
rross@laurentian.ca
'''
#################################################
#                    IMPORTS                    }
#                                               }
#################################################

import numpy as np
from scipy import constants as sc
import h5py as h5
import numpy as np
import os
import matplotlib.pyplot as plt
params = {'text.usetex' : True,
        'font.family' : 'lmodern'}
plt.rcParams.update(params)
import re

#################################################
#                   CONSTANTS                    }
#               (And Dictionaries)               }
#################################################

icode_dictionary = {
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

jtrack_dictionary = {
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

jtrack_rest_energies = { # Expressed in MeV if known, 0 otherwise (for heavy ions for instance)
    None:   None,
    -6:     sc.physical_constants['alpha particle mass energy equivalent in MeV'][0],
    -5:     sc.physical_constants['helion mass energy equivalent in MeV'][0],
    -4:     sc.physical_constants['triton mass energy equivalent in MeV'][0],
    -3:     sc.physical_constants['deuteron mass energy equivalent in MeV'][0],
    -2:     0,
    -1:     0,
    0:      0,
    1:      sc.physical_constants['proton mass energy equivalent in MeV'][0],
    2:      sc.physical_constants['proton mass energy equivalent in MeV'][0],
    3:      sc.physical_constants['electron mass energy equivalent in MeV'][0],
    4:      sc.physical_constants['electron mass energy equivalent in MeV'][0],
    5:      0,
    6:      0,
    7:      0,
    8:      sc.physical_constants['neutron mass energy equivalent in MeV'][0],
    9:      sc.physical_constants['neutron mass energy equivalent in MeV'][0],
    10:     sc.physical_constants['muon mass energy equivalent in MeV'][0],
    11:     sc.physical_constants['muon mass energy equivalent in MeV'][0],
    12:     497.611, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    13:     139.57039, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    14:     139.57039, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    15:     493.677, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    16:     493.677, # https://pdg.lbl.gov/2022/listings/contents_listings.html

    23:     134.9768, # https://pdg.lbl.gov/2022/listings/contents_listings.html

    208:    0,

    211:    0,

    308:    0,

}

jtrack_labels = {
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

    23:     'Pion-zero',

    208:    'Energy Deposition',

    211:    'EM Energy Deposition',

    308:    'Low Energy Neutron Kerma'

}

particle_colour_dictionary = {
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

particle_colors_capitalized = {
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


#################################################
#             Analyses of H5 Files              }
#                                               }
#################################################


def plot_e_vs_e(h5_file):
    ''' A simple scatter plot of muon vs neutron energies'''

    h5_file = h5.File(h5_file)

    data = h5_file['data']
    meta = h5_file['meta']

    muon_energies = data['muon_energy']
    neutron_energies = data['neutron_energy']
    parents = meta['muon_parents'][0]

    plt.scatter(muon_energies, neutron_energies)
    plt.xscale('log'); plt.yscale('log')
    plt.title('Neutron vs. Muon Energies')
    plt.xlabel('Parent Muon Energy [GeV]')
    plt.ylabel('Neutron Energy [GeV]')
    #sum_string = str(parents) + ' unique parent muons\nproducing ' + str(len(neutron_energies)) +  ' neutrons'
    #plt.text(3000,10, sum_string)

def plot_neutron_energy_histogram(h5_file, bins = 100):
    ''' A simple histogram of the neutron energies from the file'''
    h5_file = h5.File(h5_file)
    data = h5_file['data']
    neutron_energies = data['neutron_energy']
    logbins = np.logspace(-1, 1, bins)
    plt.hist(neutron_energies, bins = logbins)
    plt.yscale('log'); plt.xscale('log')
    plt.title('Neutron Energy Spectrum for N = ' + str(len(neutron_energies)))
    plt.xlabel('Energy [GeV]'); plt.ylabel('Count')
    plt.show()

def plot_impact_hist(h5_file, bins = 20, label = ''):

    h5_file = h5.File(h5_file)
    data = h5_file['data']
    plt.hist(np.unique(data['muon_impact']), bins=bins, histtype='step', label = label)
    plt.title('Muon Impact Parameters')
    plt.xlabel('Impact Parameter [cm]'); plt.ylabel('Count')

def plot_coz_neutrons(h5_file, bins=50):
    h5_file = h5.File(h5_file)
    data = h5_file['data']
    plt.hist(data['neutron_direction'][:,2], bins=bins)
    plt.title('Neutron Zenith Angles')
    plt.xlabel(r'$\cos \theta$'); plt.ylabel('Count')

def initialize_h5_file(h5_filename):
    '''Creates an h5 file with the correct datasets and group names and layouts for storing neutron data from a FLUKA run'''

# If the file does not exist
    if not os.path.isfile(h5_filename): # The file must be created.
        file = h5.File(h5_filename,'a')

        # Create groups for the data
        tpc_data = file.create_group('tpc_data')
        od_data = file.create_group('od_data')
        tpc_totals = file.create_group('tpc_totals')
        od_totals = file.create_group('od_totals')
        # Create a group for the meta data (to dilineate different simulation data sets)
        meta = file.create_group('meta')

        # TPC Data

        # So too must the datasets be created and instantiated with zero size.
        tpc_data.create_dataset("muon_energy", (0,), dtype=float, maxshape=(None,))
        tpc_data.create_dataset("muon_impact", (0,), dtype=float, maxshape=(None,))
        tpc_data.create_dataset("muon_initial", (0,3), dtype=float, maxshape=(None,3))
        tpc_data.create_dataset("muon_direction", (0,3), dtype=float, maxshape=(None,3))
        tpc_data.create_dataset("muon_pn", (0,), dtype=int, maxshape=(None,))

        # ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK
        tpc_data.create_dataset("neutron_icode", (0,), dtype=int, maxshape=(None,))
        tpc_data.create_dataset("neutron_region", (0,), dtype=int, maxshape=(None,))
        tpc_data.create_dataset("neutron_generation", (0,), dtype=int, maxshape=(None,))
        tpc_data.create_dataset("neutron_energy", (0,), dtype=float, maxshape=(None,))
        tpc_data.create_dataset("neutron_xyz", (0,3), dtype=float, maxshape=(None,3))
        tpc_data.create_dataset("neutron_direction", (0,3), dtype=float, maxshape=(None,3))
        tpc_data.create_dataset("parent", (0,), dtype=int, maxshape=(None,))
        tpc_data.create_dataset("birth_icode", (0,), dtype=int, maxshape=(None,))

        # TPC Totals

        tpc_totals.create_dataset("neutrons_counted",(0,), dtype=int, maxshape=(None,))
        # number of muons simulated
        tpc_totals.create_dataset("muons_simulated",(0,), dtype=int, maxshape=(None,))
        # number of muons responsible for creating neutrons
        tpc_totals.create_dataset("muon_parents",(0,), dtype=int, maxshape=(None,))
        # a dataset for the resnuclei output
        tpc_totals.create_dataset("resnuclei",(0,3), dtype=float, maxshape=(None,3))
        # a dataset for the resnuclei output scoring the TPC copper
        tpc_totals.create_dataset("resnuclei_cu",(0,3), dtype=float, maxshape=(None,3))


        ## Meta data

        # The integer seed used in the simulation
        meta.create_dataset("seed",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("year",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("month",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("day",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("hour",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("minute",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("second",(0,), dtype=int, maxshape=(None,))

        # OD Data

        od_data.create_dataset("muon_energy", (0,), dtype=float, maxshape=(None,))
        od_data.create_dataset("muon_impact", (0,), dtype=float, maxshape=(None,))
        od_data.create_dataset("muon_initial", (0,3), dtype=float, maxshape=(None,3))
        od_data.create_dataset("muon_direction", (0,3), dtype=float, maxshape=(None,3))
        od_data.create_dataset("muon_pn", (0,), dtype=int, maxshape=(None,))

        # ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK
        od_data.create_dataset("neutron_icode", (0,), dtype=int, maxshape=(None,))
        od_data.create_dataset("neutron_region", (0,), dtype=int, maxshape=(None,))
        od_data.create_dataset("neutron_generation", (0,), dtype=int, maxshape=(None,))
        od_data.create_dataset("neutron_energy", (0,), dtype=float, maxshape=(None,))
        od_data.create_dataset("neutron_xyz", (0,3), dtype=float, maxshape=(None,3))
        od_data.create_dataset("neutron_direction", (0,3), dtype=float, maxshape=(None,3))
        od_data.create_dataset("parent", (0,), dtype=int, maxshape=(None,))
        od_data.create_dataset("birth_icode", (0,), dtype=int, maxshape=(None,))

        od_totals.create_dataset("neutrons_counted",(0,), dtype=int, maxshape=(None,))
        # number of muons simulated
        od_totals.create_dataset("muons_simulated",(0,), dtype=int, maxshape=(None,))
        # number of muons responsible for creating neutrons
        od_totals.create_dataset("muon_parents",(0,), dtype=int, maxshape=(None,))

        file.close()
    else:
        # send the data to an hdf5 file with a different name, and merge them with the merge function
        initialize_h5_file('retry_' + h5_filename)
        ### NEED TO ADD THE MERGE HERE LATER

def merge_hdf5_files(h5_output, *args):

    '''A function to combine multiple h5 neutron files into a larger file (which may already exist).
    Checks for no duplicates by ensuring seeds are different'''

    # Load in the h5 files in read-only mode
    file_list = [h5.File(arg, 'r') for arg in args]

    # Must first confirm these files can be merged. Namely, we need to know that they are not identical.
    meta_list = [file['meta'] for file in file_list]

    # Checking seeds
    seeds = [int(m['seed'][0]) for m in meta_list]

    # Check if the output file already exists and append its attributes to the lists for uniqueness checks

    if not os.path.isfile(h5_output): # The file must be created.
        initialize_h5_file(h5_output)

    else:
        file = h5.File(h5_output, 'a')
        for seed in file['meta']['seed']:
            seeds.append(seed)
    
        file.close()

    # Make sure the input files are unique simulations with the same region number
    if len(np.unique(seeds)) < len(seeds):
        print('INCOMPATIBLE SCORING REGIONS or SAME SEED!')
        return None

    
    # We should be safe to add the attributes of all the other files to the output h5 file now.

    # How much longer does 'data' have to be?
    tpc_data_length = 0
    od_data_length = 0
    tpc_totals_length = 0
    od_totals_length = 0

    for file in file_list:
        tpc_data_length = tpc_data_length + len(file['tpc_data']['neutron_energy'])
        od_data_length = od_data_length + len(file['od_data']['neutron_energy'])
        tpc_totals_length = tpc_totals_length + len(file['tpc_totals']['neutrons_counted'])
        od_totals_length = od_totals_length + len(file['od_totals']['neutrons_counted'])

    # How much longer does 'meta' have to be?
    meta_length = 0
    for meta in meta_list:
        meta_length = meta_length + len(meta['seed'])

    # Now we can re-open the output file and resize it as appropriate.

    file = h5.File(h5_output,'a')

    # grabbing the file by the data groups
    tpc_data = file['tpc_data']
    od_data = file['od_data']
    tpc_totals = file['tpc_totals']
    od_totals = file['od_totals']

    meta = file['meta']

    # Current length of the data in the file
    current_od_data_size = len(od_data['neutron_energy'])
    current_tpc_data_size = len(tpc_data['neutron_energy'])
    current_od_tot_data_size = len(od_totals['neutrons_counted'])
    current_tpc_tot_data_size = len(tpc_totals['neutrons_counted'])

    current_meta_size = len(meta['year'])

    # Resizing the datasets in the file by the number of elements we need to append therein.
    for dset in tpc_data:
        tpc_data[dset].resize(current_tpc_data_size + tpc_data_length, axis = 0)
    for dset in od_data:
        od_data[dset].resize(current_od_data_size + od_data_length, axis = 0)
    for dset in od_totals:
        od_totals[dset].resize(current_od_tot_data_size + od_totals_length, axis = 0)
    for dset in tpc_totals:
        tpc_totals[dset].resize(current_tpc_tot_data_size + tpc_totals_length, axis = 0)

    for dset in meta:
        meta[dset].resize(current_meta_size + meta_length, axis = 0)

    meta_index = 0

    for temp_file in file_list:
        temp_tpc_data = temp_file['tpc_data']
        temp_od_data = temp_file['od_data']
        temp_tpc_totals = temp_file['tpc_totals']
        temp_od_totals = temp_file['od_totals']
        temp_meta = temp_file['meta']


        for key in tpc_data.keys():
            for i in range(len(temp_tpc_data['neutron_energy'])):
                tpc_data[key][current_tpc_data_size + i] = temp_tpc_data[key][i]

        for key in tpc_totals.keys():
            for i in range(len(temp_tpc_totals['neutron_energy'])):
                tpc_totals[key][current_tpc_tot_data_size + i] = temp_tpc_totals[key][i]
        
        for key in od_data.keys():
            for i in range(len(temp_od_data['neutron_energy'])):
                od_data[key][current_od_data_size + i] = temp_od_data[key][i]

        for key in od_totals.keys():
            for i in range(len(temp_od_totals['neutron_energy'])):
                od_totals[key][current_od_tot_data_size + i] = temp_od_totals[key][i]

        for key in meta.keys():
            meta[key][current_meta_size + meta_index] = temp_meta[key][0]

        meta_index += 1 # Goes up by one for every file
        current_od_data_size += len(temp_od_data['neutron_energy'])
        current_tpc_data_size += len(temp_tpc_data['neutron_energy'])
        current_od_tot_data_size += len(temp_od_totals['neutrons_counted'])
        current_tpc_tot_data_size += len(temp_tpc_totals['neutrons_counted'])
    
    for temp_file in file_list:
        temp_file.close()

    file.close()

    # Return the string name of the file
    return h5_output