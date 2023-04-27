#!/usr/bin/python

__version__ = 1.0
__author__ = 'Regan Ross'
## Last Edited April 26, 2023

'''
run_module.py


Contact:
Regan Ross
rross@laurentian.ca

A module for coordinating the FLUKA Simulation process on SDF.

'''
#################################################
#                    IMPORTS                    }
#                                               }
#################################################
import os
import re
import numpy as np
import h5py as h5
from datetime import datetime
from src import muon_functions as mf
import yaml

#################################################
#           CONSTANTS AND DICTIONARIES          }
#                                               }
#################################################


fluka_output_files = {
    'tpc_neutron_file' :    'nEXO_OD001_fort.72',
    'od_neutron_file' :     'nEXO_OD001_fort.70',
    'res_nuclei_file' :     'nEXO_OD001_fort.97',
    'res_nuclei_cu_file' :  'nEXO_OD001_fort.94',
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
        'parent':               {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'birth_icode':          {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
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
        'parent':               {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
        'birth_icode':          {'shape' :   (0,),   'dtype' : int,      'maxshape': (None,)},
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

def read_resnuclei_file(resnuclei_file):
    '''Scans the given ASCII resnuclei output file and returns an np array table 
        with the columns Z, A, and the number of stopping nuclei per primary'''
    
    with open(resnuclei_file) as file:
        header = file.read(500)

    # Create a temporary file to hold ONLY the array data
    os.system('sed -e \'/^[[:space:]]*[0-9].*[0-9]$/!d\' ' + resnuclei_file + ' > temp.tsv')
    # Read in the temporary data file
    data = np.loadtxt('temp.tsv')
    # Remove the temporary array-only file
    os.system('rm temp.tsv')

    # Collect the data that determine the iteration pattern.
    max_z = int(re.search('(?<=Z:).*(?=,)', header).group(0).strip())
    max_nz = int(re.search('(?<=-Z:).*(?=M)', header).group(0).strip())
    min_nz = int(re.search('(?<=Min. N-Z:).*(?=\n)', header).group(0).strip())
    k = -5
    max_a = max_nz - k

    # print_line = 'Max Z = {}, Max N-Z = {}, Min N-Z = {}'.format(max_z, max_nz, min_nz)
    # print(print_line)

    data = np.reshape(data, (max_a, max_z))

    arr_list = []

    for i in range(max_a):
        for j in range(max_z):
            Z = j + 1
            A = (i+1)+k+2*(j+1)
            if data[i,j] > 0.0:

                arr_list.append(np.array([Z, A, data[i,j]]))

    return np.array(arr_list)

def make_phase_space_file(num_muons, roi_radius = 0, roi_height = 0, filename = 'src/muon_file.txt'):
    '''A function to make a phase space file for a number of muons. 
        This file is used by the larger simulation as the source for each particle.
        NOTE: if the filename is changed, it must also be changed in the .inp file for the sim.'''

    from random import random

    if roi_radius == 0:
        roi_radius = mf.OD_RADIUS + 2
        roi_height = mf.OD_HEIGHT + 2

    roi = mf.OuterDetector(roi_radius, roi_height)

    file_stream = open(filename, 'w')

    muarray = mf.intersecting_muons(num_muons, roi)

    for muon in muarray:
        particle_code = 10
    
        if random() > 0.5:
            particle_code = 11
            
        file_stream.write(str(particle_code) + ' ' + str(muon) + '\n')

    file_stream.close()

def initialize_h5_file(h5_filename):
    '''Creates an empty hdf5 file with the structure equivalent to the above dictionary. Renames the file if it already exists.'''

    # Some bare bones error correction
    if os.path.isfile(h5_filename):
        h5_filename = os.path.splitext(h5_filename)[0] + '_new.hdf5'
        
        # If this doesn't work...
        if os.path.isfile(h5_filename):
            print('CANNOT INITIALIZE HDF5 FILE: DUPLICATES!!!')
            return None
    
    file = h5.File(h5_filename, 'a')

    for group in hdf5_structure:
        file.create_group(group)
        group = hdf5_structure[group]
        for dset in group:
            data_set = group[dset]
            file.create_dataset(dset, shape = data_set['shape'], dtype = data_set['dtype'], maxshape = data_set['maxshape'])

    file.close()

    return h5_filename
    
def remove_junk_files():
    junk_list = ['.o', '.exe', '.mod', '~']

def move_fluka_files(path):
    file_list = ['*fort*', '*lis*', '*tab*', '*.h5', '*.hdf5', '*fort*', '*.log', '*.inp', '*.err', '*.out', '*dump']

    for ext in file_list:
        os.system('mv ' + ext + path)

def read_neutron_file(neutron_filename) -> dict:
    '''Reads in a file of neutron entries with the anticipated format following:
    ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK, PICODE, PJTRACK

    Returns a dictionary of lists (one per type of entry)
    '''
    # Neutron attribute lists
    icode, ncase, jtrack, mreg, ltrack, etrack, xsco, ysco, zsco, cxtrck, cytrck, cztrck = [],[],[],[],[],[],[],[],[],[],[],[]
    # parent attribute lists
    picode, pjtrack, = [],[]

    
    with open(neutron_filename, 'r') as neutrons:
        for line in neutrons:
            temp = line.split()

            # Load the neutron lists
            icode.append(int(temp[0])); ncase.append(int(temp[1])); jtrack.append(int(temp[2]))
            mreg.append(int(temp[3])); ltrack.append(int(temp[4])); etrack.append(float(temp[5]))
            xsco.append(float(temp[6])); ysco.append(float(temp[7])); zsco.append(float(temp[8]))
            cxtrck.append(float(temp[9])); cytrck.append(float(temp[10]));  cztrck.append(float(temp[11]))

            # Load the parent lists
            picode.append(int(temp[12])); pjtrack.append(int(temp[14]))

    neutron_data = {'icode'     :   icode,
                    'ncase'     :   ncase,
                    'jtrack'    :   jtrack,
                    'mreg'      :   mreg,
                    'ltrack'    :   ltrack,
                    'etrack'    :   etrack,
                    'xsco'      :   xsco,
                    'ysco'      :   ysco,
                    'zsco'      :   zsco,
                    'cxtrck'    :   cxtrck,
                    'cytrck'    :   cytrck,
                    'cztrck'    :   cztrck,
                    'picode'    :   picode,
                    'pjtrack'   :   pjtrack
                    }
    
    return neutron_data

def resize_output_file(h5_filename, tpc_data_length, od_data_length, resnuc, resnuc_cu):

    ### Open the file
    file = h5.File(h5_filename, 'a')

    ### Retrieve current lengths
    tpc_data_current = int(len(file['tpc_data']['neutron_energy']))
    od_data_current = int(len(file['od_data']['neutron_energy']))
    meta_current = int(len(file['meta']['year']))
    tpc_tot_current = len(file['tpc_totals']['neutrons_counted'])
    od_tot_current = len(file['od_totals']['neutrons_counted'])

    ### Resize tpc and od data datasets by variable length
    for dset in file['tpc_data']:
        file['tpc_data'][dset].resize(tpc_data_current + tpc_data_length, axis = 0)
    for dset in file['od_data']:
        file['od_data'][dset].resize(od_data_current + od_data_length, axis = 0)
    
    ### Resize meta and totals datasets by one entry each

    for dset in file['meta']:
        file['meta'][dset].resize(meta_current + 1, axis = 0)
    for dset in file['od_totals']:
        file['od_totals'][dset].resize(od_tot_current + 1, axis = 0)
    for dset in file['tpc_totals']:
        file['tpc_totals'][dset].resize(tpc_tot_current + 1, axis = 0)

    ### For the residual nuclei in particular
    
    current_resnuc_len = len(file['resnuclei']['resnuclei'])
    current_resnuc_cu_len = len(file['resnuclei']['resnuclei_cu'])

    file['resnuclei']['resnuclei'].resize(current_resnuc_len + resnuc, axis = 0)
    file['resnuclei']['resnuclei_cu'].resize(current_resnuc_cu_len + resnuc_cu, axis = 0)


def store_data_in_h5(output_filename):
    '''Takes the data from the FLUKA output files and builds it into an hdf5 file.'''

    ### First we gather the data ###

    tpc_data = read_neutron_file(fluka_output_files['tpc_neutron_file'])
    tpc_length = len(tpc_data['etrack'])
    od_data = read_neutron_file(fluka_output_files['od_neutron_file'])
    od_length = len(od_data['etrack'])
    resnuclei_data = read_resnuclei_file(fluka_output_files['res_nuclei_file'])
    resnuc_length = len(resnuclei_data)
    resnuclei_cu_data = read_resnuclei_file(fluka_output_files['res_nuclei_cu_file'])
    resnuc_cu_length = len(resnuclei_cu_data)

    output_filename = initialize_h5_file(output_filename)
    resize_output_file(output_filename, tpc_length, od_length, resnuc_length, resnuc_cu_length)