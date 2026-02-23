#!/usr/bin/python

__version__ = 0.0
__author__ = 'Regan Ross'
## Last Edited June 7, 2024

'''
analysis.py

This part of the module is for performing analysis on data from the FLUKA simulations.
The data files have been created and organized, and are now ready to look at. This 
part of the module is certainly going to look the scrappiest when all is complete,
as there will be many files for manipulating the data into appropriate formats for 
plotting, and then also plotting functions.

'''
################################################################################
#                                                                              #
#                                   IMPORTS                                    #
#                                                                              #
################################################################################

# Import all the constants and dictionaries from the __init__ file
from . import *
import h5py as h5
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter

################################################################################
#                                                                              #
#                              DATA MANIPULATION                               #
#                                                                              #
################################################################################

def clean_resnuclei_data(resnuc_dict)-> None:
    '''Given a residual nuclei dictionary, this function sorts the data into a 3-column format: Z, A, and counts. It creates a new dictionary entry 'table' which will hold the newly sorted array. This can be used on residual nuclei dictionaries that have already been added together, OR individual ones.'''

    raw = resnuc_dict['raw']
    max_z = resnuc_dict['max_z']
    max_a = resnuc_dict['max_a']

    total = np.count_nonzero(raw)
    data = np.zeros((total, 3), dtype=float)

    new_raw = np.reshape(raw, (max_a, max_z))


    # Yes, the array parsing is FUNKY. Don't worry about it. This ought to be right
    # See this outdated (but still useful) man page: http://www.fluka.org/fluka.php?id=man_onl&sub=67
    count = 0
    k = -5
    for i in range(max_a):
        for j in range(max_z):
            Z = j + 1
            A = (i+1)+k+2*(j+1)
            if new_raw[i,j] > 0.0:
                data[count] = np.array([Z, A, new_raw[i,j]])
                count += 1

    resnuc_dict['table'] = data


def get_bootstrap_resnuc_data()-> np.ndarray:
    '''This function will be for producing an array of activation rates per year by sampling randomly
    and repeatedly from within a larger number of dictionaries.'''
    pass


################################################################################
#                                                                              #
#                              PLOTTING FUNCTIONS                              #
#                                                                              #
################################################################################

def plot_resnuclei_table_of_nuclides(nuclide_table : np.ndarray, title = 'Residual Nuclei', time = '') -> None:
    '''Anticipates an array whose rows are [Z, A, #]. Produces a plot of the table of nuclides
    as a sort of 2D histogram. One can hover the cursor over different boxes on the plot for 
    the particular information (Z, A, #)'''

    ## Turn the table into a 2D array 
    table = nuclide_table
    max_a = int(np.max(table[:,1]))
    max_z = int(np.max(table[:,0]))
    nuclides = np.zeros((max_z + 1, max_a + 1))

    for entry in table:
        x = int(entry[0])
        y = int(entry[1])
        val = entry[2]
        nuclides[x,y] = val

    nuclides[np.where(nuclides == 0)] = np.nan

    fig = px.imshow(nuclides, origin = 'lower')

    fig.update_layout(title = title, xaxis_title = 'Atomic Mass A', \
                    yaxis_title = 'Atomic Number Z', width = 1100, height = 600)

    fig.update_traces(hovertemplate='<b>Atomic Mass A:</b> %{x}<br>' + \
                    '<b>Atomic Number Z:</b> %{y}<br>' +
                    '<b>Count:</b> %{z}<br>' +
                    '<extra></extra>')
    
    fig.update_layout(
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True))
    fig.update_layout(coloraxis_colorbar = dict(title = 'Total Count ' + time, titleside = 'right', len = 0.6))
    fig.show()


def plot_neutron_parents_pie(neutrons : np.ndarray, title = 'Parents of Neutrons')-> None:
    '''Produces a plot of the immediate parent of neutrons as they are scored in the TPC.
    This simply tallies the parent JTRACK variables from the muon array.'''

    # Numbers representing parent particle species
    particle_species = neutrons[:,14]
    # Count occurrences of each particle species
    species_counts = Counter(particle_species)

    # Prepare labels and values for the pie chart
    labels = [JTRACK_LABELS[num] for num in species_counts.keys()]
    values = list(species_counts.values())

    # Create the pie chart
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, textinfo='none')])

    # Customize the layout
    fig.update_layout(
        title=title,
        title_x=0.5
    )

    # Show the plot
    fig.show()

def plot_neutron_icodes_pie(neutrons : np.ndarray, title = 'Neutron Birth Icodes')-> None:
    '''Produces a pie chart plot of the ICode values of the neutrons.
    In FLUKA lingo, this is basically their creation events.'''
    
    # Sample data: list of numbers representing particle species
    icodes = neutrons[:,0]
    # Count occurrences of each particle species
    icodes_counts = Counter(icodes)

    # Prepare labels and values for the pie chart
    labels = [ICODE_DICTIONARY[num] for num in icodes_counts.keys()]
    values = list(icodes_counts.values())

    # Create the pie chart
    fig = go.Figure(data=[go.Pie(labels=labels, values=values)])

    # Customize the layout
    fig.update_layout(
        title=title,
        title_x=0.5
    )
    # Show the plot
    fig.show()


################################################################################
#                                                                              #
#                         HDF5 FILE ANALYSIS FUNCTIONS                         #
#                                                                              #
################################################################################

def grab_file(file):
    '''If the file is only a string, this will open the file and return the handle'''
    if type(file) is str:
        return h5.File(file, 'r')
    else:
        return file

def plot_e_vs_e(h5_file, tpc = False):
    ''' A simple scatter plot of muon vs neutron energies'''

    h5_file = grab_file(h5_file)

    if tpc:
        data = h5_file['tpc_data']
    else:
        data = h5_file['od_data']

    meta = h5_file['meta']

    muon_energies = data['muon_energy']
    neutron_energies = data['neutron_energy']
    parents = meta['muon_parents'][0]

    plt.scatter(muon_energies, neutron_energies)
    plt.xscale('log'); plt.yscale('log')
    plt.title('Neutron vs. Muon Energies')
    plt.xlabel('Parent Muon Energy [GeV]')
    plt.ylabel('Neutron Energy [GeV]')

def plot_neutron_energy_histogram(h5_file, bins = 100, tpc = False):
    ''' A simple histogram of the neutron energies from the file'''
    h5_file = grab_file(h5_file)

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

def plot_impact_hist(h5_file, bins = 20, label = '', tpc = False):

    h5_file = grab_file(h5_file)

    if tpc:
        data = h5_file['tpc_data']
    else:
        data = h5_file['od_data']

    plt.hist(np.unique(data['muon_impact']), bins=bins, histtype='step', label = label)
    plt.title('Muon Impact Parameters')
    plt.xlabel('Impact Parameter [cm]'); plt.ylabel('Count')

def get_total_muons(h5_file):

    summary =  {'muons_simulated'           :       np.sum(h5_file['od_totals']['muons_simulated']),
                'neutrons_counted_od'        :       np.sum(h5_file['od_totals']['neutrons_counted']),
                'neutrons_counted_tpc'       :       np.sum(h5_file['tpc_totals']['neutrons_counted']),
                'muon_parents_od'           :       np.sum(h5_file['od_totals']['muon_parents']),
                'muon_parents_tpc'          :       np.sum(h5_file['tpc_totals']['muon_parents'])                          
                }
    
    return summary

def plot_both_impact_hist(h5_file, bins = 20, label = ''):

    h5_file = grab_file(h5_file)

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

    h5_file = grab_file(h5_file)

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
    ax2.set_title('TPC Neutron Energies N = ' + str(totals['neutrons_counted_tpc']))
    ax2.set_xlabel('Neutron Energy [GeV]')
    ax2.set_yscale('log'); ax2.set_xscale('log')
    ax1.set_ylabel('Count'); ax2.set_ylabel('Count')

def plot_coz_neutrons(h5_file, bins=50, tpc = False):

    h5_file = grab_file(h5_file)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10,5))

    tpc_data = h5_file['tpc_data']

    od_data = h5_file['od_data']

    totals = get_total_muons(h5_file)

    ax1.hist(np.unique(od_data['neutron_direction'][:,2]), bins=bins, histtype='step', label = 'OD')
    ax1.set_title('OD Neutron Zenith Angles N = ' + str(totals['neutrons_counted_od']) + ' neutrons')
    ax1.set_xlabel(r'$\cos \theta$'); ax1.set_ylabel('Count')

    ax2.hist(np.unique(tpc_data['neutron_direction'][:,2]), bins=bins, histtype='step', label = 'OD')
    ax2.set_title('TPC Neutron Zenith Angles N = ' + str(totals['neutrons_counted_tpc']) + ' neutrons')
    ax2.set_xlabel(r'$\cos \theta$'); ax2.set_ylabel('Count')
    plt.show()

def get_activation(h5_file, z, a, tpc = True):
    ''' Returns the absolute "count" of a residual nuclide given an hdf5 file'''

    h5_file = grab_file(h5_file)
    muons_simulated = h5_file['meta']['muons_simulated'][0]

    if tpc:
        resnuc_data = h5_file['resnuclei']['resnuclei']
    else:
        resnuc_data = h5_file['resnuclei']['resnuclei_cu']

    entries = []
    for entry in resnuc_data:
        if entry[1] == a and entry[0] == z:
            entries.append(entry)

    totals = []

    if len(entries) > 0:
        for entry in entries:
            totals.append(entry[2]*muons_simulated)

        return np.sum(totals)
    else:
        return 0

def get_xe137_activation(h5_file):
    ''' Returns the absolute "count" of Xe-137 tabled by the resnuclei FLUKA function given an hdf5 file'''

    h5_file = grab_file(h5_file)

    muons_simulated = h5_file['meta']['muons_simulated'][0]

    xe137_entries = []
    for entry in h5_file['resnuclei']['resnuclei']:
        if entry[1] == 137 and entry[0] == 54:
            xe137_entries.append(entry)

    totals = []

    if len(xe137_entries) > 0:
        for entry in xe137_entries:
            totals.append(entry[2]*muons_simulated)

        return np.sum(totals)
    else:
        return 0

def get_total_years(h5_file) -> float:
    h5_file = grab_file(h5_file)
    years = np.sum(h5_file['meta']['hours_simulated'])/8760
    return years 

def print_summary(h5_file) -> dict:
    '''Creates a summary dictionary of the file, prints a summary, and returns the dictionary'''
    h5_file = grab_file(h5_file)

    summary =  {'muons_simulated'           :       np.sum(h5_file['od_totals']['muons_simulated']),
                'neutrons_counted_od'       :       np.sum(h5_file['od_totals']['neutrons_counted']),
                'neutrons_counted_tpc'      :       np.sum(h5_file['tpc_totals']['neutrons_counted']),
                'muon_parents_od'           :       np.sum(h5_file['od_totals']['muon_parents']),
                'muon_parents_tpc'          :       np.sum(h5_file['tpc_totals']['muon_parents']),
                'roi_radius'                :       h5_file['meta']['roi_radius'][0],
                'roi_height'                :       h5_file['meta']['roi_height'][0],
                'xe137_activation'          :       get_xe137_activation(h5_file),
                'cu64_activation'           :       get_activation(h5_file, 29, 64, False),
                'cu66_activation'           :       get_activation(h5_file, 29, 66, False),
                'hours_simulated'           :       np.sum(h5_file['meta']['hours_simulated']),
                'years_simulated'           :       np.sum(h5_file['meta']['hours_simulated'])/8760,
                }
    
    print('Summary of Fluka File:\n')
    print('Muons Simulated: ' + str(summary['muons_simulated']))
    print('Intersecting ROI Radius: ' + str(summary['roi_radius']) + '[m] Height: ' + str(summary['roi_height']) + ' [m]')
    print(str(summary['muon_parents_od']) + ' muons creating ' + str(summary['neutrons_counted_od']) + ' neutrons in the OD')
    print(str(summary['muon_parents_tpc']) + ' muons creating ' + str(summary['neutrons_counted_tpc']) + ' neutrons in the TPC')
    print('Xenon-137 atoms counted: ' + str(summary['xe137_activation']))
    print('Copper-64 atoms counted in TPC shell: ' + str(summary['cu64_activation']))
    print('Copper-66 atoms counted in TPC shell: ' + str(summary['cu66_activation']))
    print('Time Simulated: ' + str(summary['hours_simulated']) + ' hours or ' + str(summary['years_simulated']) + ' years')
    
    return summary

def tabulate_resnuclei_data(resnuclei_dataset) -> np.ndarray:

    max_z = 0
    max_a = 0

    for entry in resnuclei_dataset:
        if entry[0] > max_z:
            max_z = int(entry[0])

        if entry[1] > max_a:
            max_a = int(entry[1])

    resnuclei_array = np.zeros((max_z + 1, max_a + 1))

    for entry in resnuclei_dataset:
        x , y, val = int(entry[0]), int(entry[1]), entry[2]

        resnuclei_array[x][y] += val

    return resnuclei_array

def activation_plot(h5_file):
    '''Plots a histogram of the table of nuclides counted in the detector'''
    h5_file = grab_file(h5_file)

    muons_per_run = h5_file['meta']['muons_simulated'][0]
    resnuclei = tabulate_resnuclei_data(h5_file['resnuclei']['resnuclei'])*muons_per_run
    max_z = np.shape(resnuclei)[0]; max_a = np.shape(resnuclei)[1]
    bins = np.arange(0.1, (max_z + 1.1), 1), np.arange(0.1, (max_a + 1.1), 1)

    resnuclei_cu = tabulate_resnuclei_data(h5_file['resnuclei']['resnuclei_cu'])*muons_per_run
    max_z_cu = np.shape(resnuclei_cu)[0]; max_a_cu = np.shape(resnuclei_cu)[1]
    bins_cu = np.arange(0.1, (max_z_cu + 1.1), 1), np.arange(0.1, (max_a_cu + 1.1), 1)

    resnuc_z, resnuc_a, resnuc_weights = [],[],[]
    for z in range(np.shape(resnuclei)[0]):
        for a in range(np.shape(resnuclei)[1]):
            resnuc_z.append(z); resnuc_a.append(a)
            resnuc_weights.append(resnuclei[z][a])
    
    resnuc_cu_z, resnuc_cu_a, resnuc_cu_weights = [],[],[]
    for z in range(np.shape(resnuclei_cu)[0]):
        for a in range(np.shape(resnuclei_cu)[1]):
            resnuc_cu_z.append(z); resnuc_cu_a.append(a)
            resnuc_cu_weights.append(resnuclei_cu[z][a])

    fig, ax = plt.subplots(1,2, figsize = (8, 8))
    hist1 = ax[0].hist2d(resnuc_z, resnuc_a, bins = bins, weights = resnuc_weights,  norm=colors.LogNorm(), cmap=mpl.colormaps['winter'])
    ax[0].grid(which = 'both')
    ax[0].set_title('Activation in TPC')
    ax[0].set_xticks(np.arange(0, np.shape(resnuclei)[0] + 1, 5))
    ax[0].set_yticks(np.arange(0, np.shape(resnuclei)[1] + 1, 5))

    hist2 = ax[1].hist2d(resnuc_cu_z, resnuc_cu_a, bins = bins_cu, weights = resnuc_cu_weights,  norm=colors.LogNorm(), cmap=mpl.colormaps['copper'])
    ax[1].grid(which = 'both')
    ax[1].set_title('Activation in TPC Copper')
    ax[1].set_xticks(np.arange(0, np.shape(resnuclei_cu)[0] + 1, 5))
    ax[1].set_yticks(np.arange(0, np.shape(resnuclei_cu)[1] + 1, 5))

    fig.colorbar(hist1[3], ax=ax[0])
    fig.colorbar(hist2[3], ax=ax[1])

    plt.show()

def activation_plot_per_year(h5_file):
    '''Plots a histogram of the table of nuclides counted in the detector, normalized per year'''
    h5_file = grab_file(h5_file)

    muons_per_run = h5_file['meta']['muons_simulated'][0]
    resnuclei = tabulate_resnuclei_data(h5_file['resnuclei']['resnuclei'])*muons_per_run/get_total_years(h5_file)
    max_z = np.shape(resnuclei)[0]; max_a = np.shape(resnuclei)[1]
    bins = np.arange(0.1, (max_z + 1.1), 1), np.arange(0.1, (max_a + 1.1), 1)

    resnuclei_cu = tabulate_resnuclei_data(h5_file['resnuclei']['resnuclei_cu'])*muons_per_run/get_total_years(h5_file)
    max_z_cu = np.shape(resnuclei_cu)[0]; max_a_cu = np.shape(resnuclei_cu)[1]
    bins_cu = np.arange(0.1, (max_z_cu + 1.1), 1), np.arange(0.1, (max_a_cu + 1.1), 1)

    resnuc_z, resnuc_a, resnuc_weights = [],[],[]
    for z in range(np.shape(resnuclei)[0]):
        for a in range(np.shape(resnuclei)[1]):
            resnuc_z.append(z); resnuc_a.append(a)
            resnuc_weights.append(resnuclei[z][a])
    
    resnuc_cu_z, resnuc_cu_a, resnuc_cu_weights = [],[],[]
    for z in range(np.shape(resnuclei_cu)[0]):
        for a in range(np.shape(resnuclei_cu)[1]):
            resnuc_cu_z.append(z); resnuc_cu_a.append(a)
            resnuc_cu_weights.append(resnuclei_cu[z][a])

    fig, ax = plt.subplots(1,2, figsize = (8, 8))
    hist1 = ax[0].hist2d(resnuc_z, resnuc_a, bins = bins, weights = resnuc_weights,  norm=colors.LogNorm(), cmap=mpl.colormaps['winter'])
    ax[0].grid(which = 'both')
    ax[0].set_title('Activation in TPC')
    ax[0].set_xticks(np.arange(0, np.shape(resnuclei)[0] + 1, 5))
    ax[0].set_xlabel('Atomic Number Z'); ax[0].set_ylabel('Atomic Mass A')
    ax[0].set_yticks(np.arange(0, np.shape(resnuclei)[1] + 1, 5))

    hist2 = ax[1].hist2d(resnuc_cu_z, resnuc_cu_a, bins = bins_cu, weights = resnuc_cu_weights,  norm=colors.LogNorm(), cmap=mpl.colormaps['copper'])
    ax[1].grid(which = 'both')
    ax[1].set_title('Activation in TPC Copper')
    ax[1].set_xticks(np.arange(0, np.shape(resnuclei_cu)[0] + 1, 5))
    ax[1].set_xlabel('Atomic Number Z')
    ax[1].set_yticks(np.arange(0, np.shape(resnuclei_cu)[1] + 1, 5))

    fig.colorbar(hist1[3], ax=ax[0])
    fig.colorbar(hist2[3], ax=ax[1])

    plt.show()

def get_hdf5_files(path, with_path = False)-> list:
    ''' Given a particular path, this returns a list of hdf5 files found at that location
    if with_path, the path to the file will be prepended '''

    from os import listdir
    from os.path import isfile, join

    h5_files = [f for f in listdir(path) if (isfile(join(path, f)) and f[-5:] == '.hdf5')]

    if with_path:
        h5_files = [path + '/' + file for file in h5_files]

        return h5_files
    else:
        return h5_files

def initialize_h5_file(h5_filename) -> str:
    '''Creates an empty hdf5 file with the structure equivalent to HDF5_STRUCTURE. Renames the file if it already exists.'''

    if os.path.isfile(h5_filename):
        h5_filename = os.path.splitext(h5_filename)[0] + '_new.hdf5'
        
        if os.path.isfile(h5_filename):
            print('CANNOT INITIALIZE HDF5 FILE: DUPLICATES!!!')
            return None
    
    file = h5.File(h5_filename, 'a')

    for group_name in HDF5_STRUCTURE:
        group = file.create_group(group_name)
        group_dict = HDF5_STRUCTURE[group_name]
        for dset in group_dict:
            data_set = group_dict[dset]
            group.create_dataset(dset, shape = data_set['shape'], dtype = data_set['dtype'], maxshape = data_set['maxshape'])

    file.close()

    return h5_filename

def different_seeds(file_paths) -> bool:
    seeds = []
    for path in file_paths:
        with h5.File(path, 'r') as file:
            for seed in file['meta']['seed']:
                seeds.append(seed)

    if len(np.unique(seeds)) == len(seeds):
        return True
    else:
        return False
    
def same_roi(file_paths) -> bool:

    roi_radii = []
    roi_heights = []
    for path in file_paths:
        with h5.File(path, 'r') as file:
            for rad in file['meta']['roi_radius']:
                roi_radii.append(rad)

            for height in file['meta']['roi_height']:
                roi_heights.append(height)

    if len(np.unique(roi_radii)) > 1:
        return False
    elif len(np.unique(roi_heights)) > 1:
        return False
    else:
        return True
    
def get_filenames(path):
    from os import listdir
    from os.path import isfile, join

    cwd = path
    onlyfiles = [os.path.join(cwd, f) for f in os.listdir(cwd) if 
    os.path.isfile(os.path.join(cwd, f))]

    return onlyfiles

def same_number_muons(file_paths)-> bool:

    numbers = []
    for path in file_paths:
        with h5.File(path, 'r') as file:
            for num in file['meta']['muons_simulated']:
                numbers.append(num)

    if len(np.unique(numbers)) == 1:
        return True
    else:
        return False

def check_all_runs_equal(file_paths) -> bool:
    '''determines whether or not all files have a unique seed, the same region of interest, and the same number of simulated muons'''
    if not different_seeds(file_paths):
        return False
    elif not same_roi(file_paths):
        return False
    elif not same_number_muons(file_paths):
        return False
    else:
        return True
    
def merge_hdf5_files(file_paths, output_path):

    if not check_all_runs_equal(file_paths):
        return

    output_path = initialize_h5_file(output_path)

    with h5.File(output_path, 'a') as output:

        for path in file_paths:
            file = h5.File(path, 'r')

            for groupname, group in file.items():
                
                output_group = output[groupname]

                for dsetname, dset in group.items():
                    output_dset = output_group[dsetname]

                    current_size = output_dset.shape[0]

                    new_data = np.array(dset)

                    output_dset.resize(size = current_size + dset.shape[0], axis = 0)

                    output_dset[current_size:] = new_data

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