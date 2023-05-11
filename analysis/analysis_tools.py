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

hdf5_structure = {
    ###---->  Meta data about the respective simulation
    'meta' : {   
        'seed':     {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'year':     {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'month':    {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'day':      {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'hour':     {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'minute':   {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
        'second':   {'shape' : (0,), 'dtype' : int, 'maxshape': (None,)},
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

#################################################
#             Analyses of H5 Files              }
#                                               }
#################################################


def get_hdf5_files(path, with_path = False)-> list:
    ''' Given a particular path, this returns a list of hdf5 files found at that location
    if with_path, the path to the file will be prepended '''

    from os import listdir
    from os.path import isfile, join

    h5_files = [f for f in listdir(path) if (isfile(join(path, f)) and f[-5:] == '.hdf5')]

    if with_path:
        h5_files = [path + file for file in h5_files]

        return h5_files
    else:
        return h5_files

def plot_e_vs_e(h5_file, tpc = False):
    ''' A simple scatter plot of muon vs neutron energies'''

    h5_file = h5.File(h5_file)

    if tpc:
        data = h5_file['tpc_data']
    else:
        data = h5.file['od_data']

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

def plot_neutron_energy_histogram(h5_file, bins = 100, tpc = False):
    ''' A simple histogram of the neutron energies from the file'''

    if tpc:
        data = h5_file['tpc_data']
    else:
        data = h5_file['od_data']

    neutron_energies = data['neutron_energy']
    logbins = np.logspace(-1, 1, bins)
    plt.hist(neutron_energies, bins = logbins)
    plt.yscale('log'); plt.xscale('log')
    plt.title('Neutron Energy Spectrum for N = ' + str(len(neutron_energies)))
    plt.xlabel('Energy [GeV]'); plt.ylabel('Count')
    plt.show()

def get_total_muons(h5_file):

    summary =  {'muons_simulated'           :       np.sum(h5_file['od_totals']['muons_simulated']),
                'neutrons_counted_od'        :       np.sum(h5_file['od_totals']['neutrons_counted']),
                'neutrons_counted_tpc'       :       np.sum(h5_file['tpc_totals']['neutrons_counted']),
                'muon_parents_od'           :       np.sum(h5_file['od_totals']['muon_parents']),
                'muon_parents_tpc'          :       np.sum(h5_file['tpc_totals']['muon_parents'])                          
                }
    
    return summary

def plot_impact_hist(h5_file_path, bins = 20, label = '', tpc = False):

    h5_file = h5.File(h5_file_path)

    if tpc:
        data = h5_file['tpc_data']
    else:
        data = h5_file['od_data']

    plt.hist(np.unique(data['muon_impact']), bins=bins, histtype='step', label = label)
    plt.title('Muon Impact Parameters')
    plt.xlabel('Impact Parameter [cm]'); plt.ylabel('Count')

def plot_both_impact_hist(h5_file, bins = 20, label = ''):


    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10,5))

    tpc_data = h5_file['tpc_data']

    od_data = h5_file['od_data']

    totals = get_total_muons(h5_file)

    ax1.hist(np.unique(od_data['muon_impact']), bins=bins, histtype='step', label = 'OD')
    ax1.set_title('OD Muon Impact Parameters N = ' + str(totals['muon_parents_od']) + ' muons \n making ' + str(totals['neutrons_counted_od']) + ' neutrons')
    ax1.set_xlabel('Impact Parameter [m]')

    ax2.hist(np.unique(tpc_data['muon_impact']), bins=bins, histtype='step', label = 'TPC')
    ax2.set_title('TPC Muon Impact Parameters N = ' + str(totals['muon_parents_tpc'])+ ' muons \n making ' + str(totals['neutrons_counted_tpc']) + ' neutrons')
    ax2.set_xlabel('Impact Parameter [m]')
    ax1.set_ylabel('Count'); ax2.set_ylabel('Count')

def plot_both_energy_hist(h5_file, bins = 20, label = ''):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10,5))

    tpc_data = np.array(h5_file['tpc_data']['neutron_energy']) 
    tpc_data = tpc_data - np.ones(len(tpc_data))*jtrack_rest_energies[8]/1000

    od_data = np.array(h5_file['od_data']['neutron_energy'])
    od_data = od_data - np.ones(len(od_data))*jtrack_rest_energies[8]/1000

    totals = get_total_muons(h5_file)

    logbins = np.logspace(-9, 1, bins)

    ax1.hist(np.unique(od_data), bins=logbins, histtype='step', label = 'OD')
    ax1.set_yscale('log'); ax1.set_xscale('log')
    ax1.set_title('OD Neutron Energies N = ' + str(totals['neutrons_counted_od']))
    ax1.set_xlabel('Neutron Energy [GeV]')

    ax2.hist(np.unique(tpc_data), bins=logbins, histtype='step', label = 'TPC')
    ax2.set_title('TPC Muon Impact Parameters N = ' + str(totals['neutrons_counted_tpc']))
    ax2.set_xlabel('Neutron Energy [GeV]')
    ax2.set_yscale('log'); ax2.set_xscale('log')
    ax1.set_ylabel('Count'); ax2.set_ylabel('Count')

def plot_coz_neutrons(h5_file, bins=50, tpc = False):
    h5_file = h5.File(h5_file)
    if tpc:
        data = h5_file['tpc_data']
    else:
        data = h5.file['od_data']
    plt.hist(data['neutron_direction'][:,2], bins=bins)
    plt.title('Neutron Zenith Angles')
    plt.xlabel(r'$\cos \theta$'); plt.ylabel('Count')

def initialize_h5_file(h5_filename) -> str:
    '''Creates an empty hdf5 file with the structure equivalent to the above dictionary. Renames the file if it already exists.'''

    # Some bare bones error correction
    if os.path.isfile(h5_filename):
        h5_filename = os.path.splitext(h5_filename)[0] + '_new.hdf5'
        
        # If this doesn't work...
        if os.path.isfile(h5_filename):
            print('CANNOT INITIALIZE HDF5 FILE: DUPLICATES!!!')
            return None
    
    file = h5.File(h5_filename, 'a')

    for group_name in hdf5_structure:
        group = file.create_group(group_name)
        group_dict = hdf5_structure[group_name]
        for dset in group_dict:
            data_set = group_dict[dset]
            group.create_dataset(dset, shape = data_set['shape'], dtype = data_set['dtype'], maxshape = data_set['maxshape'])

    file.close()

    return h5_filename

def different_seeds(file_paths) -> bool:
    ### PROBLEM
    seeds = []
    for path in file_paths:
        with h5.File(path, 'r') as file:
            for seed in file['meta']['seed']:
                seeds.append(seed)

    if len(np.unique(seeds)) == len(seeds):
        return True
    # The seeds are different

    else:
        return False

def merge_hdf5_files(file_paths, output_path):

    import h5py as h5

    # A quick safety check to make sure the seeds for each simulation are different
    # if not different_seeds(file_paths):
    #     return

    output_path = initialize_h5_file(output_path)

    # Open the output file in write mode
    with h5.File(output_path, 'a') as output:

        # Cycle through each of the given files
        for path in file_paths:
            file = h5.File(path, 'r')

            # Cycle through the groups of the file
            for groupname, group in file.items():
                
                output_group = output[groupname]

                # Cycle through the datasets in the file
                for dsetname, dset in group.items():
                    output_dset = output_group[dsetname]

                    current_size = output_dset.shape[0]

                    new_data = np.array(dset)

                    output_dset.resize(size = current_size + dset.shape[0], axis = 0)

                    output_dset[current_size:] = new_data

###
 # Not necessarily functional for this module yet _ salvaged from another notebook
 # Must be converted from tsv neutron file (fort.70) to pulling from hdf5
###

def event_pi_plot(neutron_list, title = '', creation = False):
    if creation: index = 11
    else: index = 0
    icode_list = [int(neutron[index]) for neutron in neutron_list]
    unique, counts = np.unique(icode_list, return_counts=True)
    labels = [icode_dictionary[icode] for icode in unique]

    fig, ax = plt.subplots()
    plt.title(title + ' N = ' + str(len(icode_list)))
    ax.pie(counts, labels=labels)

def plot_energy_spectrum(neutrons, title = '', bins = 0, logscale=True):
    if bins == 0:
        nbins = int(len(neutrons)/20)
    else:
        nbins = bins

    energies = [entry[4] - jtrack_rest_energies[8]/1000 for entry in neutrons]
    bins = np.linspace(np.min(energies),np.max(energies), nbins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),nbins)
    plt.hist(energies, log=logscale, histtype='step', bins = logbins)
    plt.title(title)
    plt.ylabel('Count'); plt.xlabel('Energy [GeV]')
    plt.xscale('log')
    plt.show()