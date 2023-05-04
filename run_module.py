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



## Loading in the parameters from the simconfig.yaml file

with open('simconfig.yaml') as yaml_file:

    input_yaml = yaml.safe_load(yaml_file)

    Simulation = input_yaml.get('Simulation')
    Input = input_yaml.get('Input')
    Output = input_yaml.get('Output')
    input_path = Input.get('InputPath')

    yaml_card = {
        # Simulation Parameters
        'num_muons' : Simulation.get('Muons'),
        'intersecting' : Simulation.get('Intersecting'),
        'make_new' : Simulation.get('MakeNewFile'),
        'reps' : Simulation.get('Repititions'),

        # Input Parameters
        'input_file' : input_path + Input.get('InputFile'),
        'source_routine' : input_path + Input.get('SourceFile'),
        'mgdraw_file' : input_path + Input.get('MGDrawFile'),

        # Output Parameters
        'output_dir' : Output.get('OutputDir'),
        'neutron_file' : Output.get('NeutronFile'),
        'progress_out' : Output.get('ProgressOut'),

        # Source Parameters
        'source_path' : input_yaml.get('Source').get('FlukaPath'),
    }

fluka_files = {
    'tpc_neutron_file'      :   yaml_card['input_file'][:-4] + '001_fort.72',
    'od_neutron_file'       :   yaml_card['input_file'][:-4] + '001_fort.70',
    'res_nuclei_file'       :   yaml_card['input_file'][:-4] + '001_fort.97',
    'res_nuclei_cu_file'    :   yaml_card['input_file'][:-4] + '001_fort.94',
    'muon_file'             :   'src/muon_file.txt',
    'mgdraw_file'           :   'mgdraw_neutron_count.f',
    'source_file'           :   'muon_from_file.f',  
    'input_file'            :   yaml_card['input_file'],
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

    for group_name in hdf5_structure:
        group = file.create_group(group_name)
        group_dict = hdf5_structure[group_name]
        for dset in group_dict:
            data_set = group_dict[dset]
            group.create_dataset(dset, shape = data_set['shape'], dtype = data_set['dtype'], maxshape = data_set['maxshape'])

    file.close()

    return h5_filename
    
def move_output_files(path):
    '''Moves simulation output files to a specified path'''

    try:
        os.system('mkdir ' + path)
    except: pass

    time_stamp = str(datetime.now())[10:16].replace(' ','').replace('-', '').replace(':', '')
    sub_dir = yaml_card['neutron_file'] + time_stamp
    last_dir = path + sub_dir + '/'

    try:
        os.system('mkdir ' + last_dir)
        os.system('mv *fort* ' + last_dir)
        os.system('mv *lis* *tab* ' + last_dir)
        os.system('mv *.hdf5 ' + path)
        os.system('mv *fort.97 ' + path)
        os.system('mv *.log *.err *.out *ran* *dump *fort* *.txt ' + last_dir)
    except: pass

    os.system('cp ' + fluka_files['input_file'] + ' ' + last_dir)

    ## STEP TEN: REMOVE COMPILED AND UNNECESSARY FILES
    os.system('rm *.o *.exe *.mod')

def move_fluka_files(path):
    file_list = ['*fort*', '*lis*', '*tab*', '*.h5', '*.hdf5', '*fort*', '*.log', '*.inp', '*.err', '*.out', '*dump']

    for ext in file_list:
        os.system('mv ' + ext + ' ' + path)

def change_seed(input_file = fluka_files['input_file']):
    '''Changes the seed to the simulation in a given input file'''

    seed = str(np.random.randint(0,999999))
    os.system('echo RR: Changing seed in input file to: ' + seed)
    space_string = '                      '
    num_spaces = len(space_string) - len(str(seed))
    num_string = 'RANDOMIZ' + space_string[:num_spaces] + seed
    os.system('sed -i \'s/^RANDOMIZ.*/' + num_string + '/\' ' + input_file)

    return seed

def link_and_compile(path_to_fluka, mgdraw_file = fluka_files['mgdraw_file'], 
                     source_routine = fluka_files['source_file'], 
                     progress_out = 'compile_summary.txt',
                     executable = 'nEXOsim.exe'):
    '''Links and compiles the fluka routines for the fluka executable'''

    if not path_to_fluka[-1] == '/':
        path_to_fluka = path_to_fluka + '/'
    
    compile_string = path_to_fluka + 'fff'

    os.system(compile_string + ' ' + mgdraw_file + ' >> '+ progress_out)
    os.system(compile_string + ' ' + source_routine + ' >> ' + progress_out)

    link_string = path_to_fluka + 'ldpmqmd -m fluka -o ' + executable + ' '
    mgd_compd = mgdraw_file[:-1] + 'o'
    source_compd = source_routine[:-1] + 'o'
    os.system(link_string + mgd_compd + ' ' + source_compd + '>> ' + progress_out)

def change_number_of_muons(num_muons, input_filename = fluka_files['input_file']):
    '''Changes the number of muons in the input file for a given simulation'''

    space_string = '               '
    num_spaces = len(space_string) - len(str(num_muons))
    num_string = 'START' + space_string[:num_spaces] + str(num_muons)
    os.system('sed -i \'s/^START.*/' + num_string + '/\' ' + input_filename)

def read_neutron_file(neutron_filename) -> dict:
    '''Reads in a file of neutron entries with the anticipated format following:
    ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK, PICODE, PJTRACK

    Returns a dictionary of lists (one per type of entry)
    '''

    # Make sure file exists!!! If it doesn't print file not found

    if not os.path.exists(neutron_filename):
        print('Neutron file not found.')
        return None

    # Neutron attribute lists
    icode, ncase, jtrack, mreg, ltrack, etrack, xsco, ysco, zsco, cxtrck, cytrck, cztrck = [],[],[],[],[],[],[],[],[],[],[],[]
    # parent attribute lists
    picode, pjtrack, = [],[]
    
    with open(neutron_filename, 'r') as neutrons:
        count = 0
        for line in neutrons:
            temp = line.split()
            count += 1

            # Load the neutron lists
            # ICODE, NCASE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK
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
                    'pjtrack'   :   pjtrack,
                    'total'     :   count,
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

    indices = {
        'tpc_data'      :   tpc_data_current,
        'od_data'       :   od_data_current,
        'od_totals'     :   od_tot_current,
        'tpc_totals'    :   tpc_tot_current,
        'meta'          :   meta_current,
        'resnuclei'     :   current_resnuc_len,
        'resnuclei_cu'  :   current_resnuc_cu_len
    }

    return indices

def retrieve_muons(muon_filename, ncase_list):
    '''Retrieves the list of muons that correspond to the list of ncase values. Don't forget, the ncase variable doesn't directly correspond to the muon number in the file. They are different by 1 as FORTRAN indices begin at 1, python indices begin at 0.'''

    with open(muon_filename) as muons: 
        muenergy, impact,  initx, inity, initz, mucosx, mucosy, mucosz, pos_neg = [],[],[],[],[],[],[],[],[]
        lines = muons.readlines()
        total = len(lines)
        for mu in ncase_list:
            # the number in the list "muon" is 1+ the index in the file.
            muon_array = lines[mu-1].split()

            pos_neg.append(int(muon_array[0])) # Whether the muon is positive or negative
            muenergy.append(float(muon_array[1])) # Muon Energy
            initx.append(float(muon_array[2]))
            inity.append(float(muon_array[3]))
            initz.append(float(muon_array[4]))
            mucosx.append(float(muon_array[5]))
            mucosy.append(float(muon_array[6]))
            mucosz.append(float(muon_array[7]))

            cos_z = float(muon_array[7])
            cos_x = float(muon_array[5])
            zenith = np.arccos(cos_z)
            azimuth = np.arccos(cos_x/(np.sin(zenith)))
            temp_muon = mf.Muon(zenith, azimuth, initial=(float(muon_array[2]), float(muon_array[3]), float(muon_array[4])))
            impact.append(temp_muon.impact_param)

    muon_dict = {
        'muon_energy'       :   muenergy,
        'muon_impact'       :   impact,
        'muon_initial'      :   np.array([initx, inity, initz]).transpose(),
        'muon_direction'    :   np.array([mucosx, mucosy, mucosz]).transpose(),
        'muon_pn'           :   pos_neg,
        'total'             :   total
    }

    return muon_dict

def store_data_in_h5(output_filename, seed) -> bool:
    '''Takes the data from the FLUKA output files and builds it into an hdf5 file.'''

    ### First we gather the data ###

    tpc_neutrons = read_neutron_file(fluka_files['tpc_neutron_file'])
    tpc_length = 0

    od_check = False
    tpc_check = False

    if tpc_neutrons is not None:
        tpc_length = len(tpc_neutrons['etrack'])
        tpc_check = True
        # Collect the muon data
        tpc_muons = retrieve_muons(fluka_files['muon_file'], tpc_neutrons['ncase'])
    
    od_neutrons = read_neutron_file(fluka_files['od_neutron_file'])
    od_length = 0

    if od_neutrons is not None:
        od_length = len(od_neutrons['etrack'])
        od_check = True
        od_muons = retrieve_muons(fluka_files['muon_file'], od_neutrons['ncase'])

    resnuclei_data = read_resnuclei_file(fluka_files['res_nuclei_file'])
    resnuc_length = len(resnuclei_data)
    resnuclei_cu_data = read_resnuclei_file(fluka_files['res_nuclei_cu_file'])
    resnuc_cu_length = len(resnuclei_cu_data)

    # Make sure we have a file to put the data into
    output_filename = initialize_h5_file(output_filename)

    # Resize the datasets of the newly initialized file
    indices = resize_output_file(output_filename, tpc_length, od_length, resnuc_length, resnuc_cu_length)

    if tpc_check or od_check:

        with h5.File(output_filename, 'a') as output_file:
            
            #   TPC DATASETS
            if tpc_check:

                output_file['tpc_data']['neutron_energy'][indices['tpc_data']:]       = tpc_neutrons['etrack']
                output_file['tpc_data']['neutron_generation'][indices['tpc_data']:]   = tpc_neutrons['ltrack']
                output_file['tpc_data']['neutron_icode'][indices['tpc_data']:]        = tpc_neutrons['icode']
                output_file['tpc_data']['neutron_region'][indices['tpc_data']:]       = tpc_neutrons['mreg']
                output_file['tpc_data']['neutron_xyz'][indices['tpc_data']:]          = np.array([tpc_neutrons['xsco'], tpc_neutrons['ysco'], tpc_neutrons['zsco']]).transpose()
                output_file['tpc_data']['neutron_direction'][indices['tpc_data']:]    = np.array([tpc_neutrons['cxtrck'], tpc_neutrons['cytrck'], tpc_neutrons['cztrck']]).transpose()
                output_file['tpc_data']['neutron_parent'][indices['tpc_data']:]       = tpc_neutrons['pjtrack']
                output_file['tpc_data']['neutron_birth_icode'][indices['tpc_data']:]  = tpc_neutrons['picode']

                output_file['tpc_data']['muon_pn'][indices['tpc_data']:]              = tpc_muons['muon_pn']
                output_file['tpc_data']['muon_impact'][indices['tpc_data']:]          = tpc_muons['muon_impact']
                output_file['tpc_data']['muon_energy'][indices['tpc_data']:]          = tpc_muons['muon_energy']
                output_file['tpc_data']['muon_initial'][indices['tpc_data']:]         = tpc_muons['muon_initial']
                output_file['tpc_data']['muon_direction'][indices['tpc_data']:]       = tpc_muons['muon_direction']

                output_file['tpc_totals']['neutrons_counted'][indices['tpc_totals']]   = tpc_neutrons['total']
                output_file['tpc_totals']['muons_simulated'][indices['tpc_totals']]    = tpc_muons['total']
                output_file['tpc_totals']['muon_parents'][indices['tpc_totals']]       = len(np.unique(tpc_muons['muon_energy']))


            #   OD DATASETS

            if od_check:

                output_file['od_data']['neutron_energy'][indices['od_data']:]       = od_neutrons['etrack']
                output_file['od_data']['neutron_generation'][indices['od_data']:]   = od_neutrons['ltrack']
                output_file['od_data']['neutron_icode'][indices['od_data']:]        = od_neutrons['icode']
                output_file['od_data']['neutron_region'][indices['od_data']:]       = od_neutrons['mreg']
                output_file['od_data']['neutron_xyz'][indices['od_data']:]          = np.array([od_neutrons['xsco'], od_neutrons['ysco'], od_neutrons['zsco']]).transpose()
                output_file['od_data']['neutron_direction'][indices['od_data']:]    = np.array([od_neutrons['cxtrck'], od_neutrons['cytrck'], od_neutrons['cztrck']]).transpose()
                output_file['od_data']['neutron_parent'][indices['od_data']:]       = od_neutrons['pjtrack']
                output_file['od_data']['neutron_birth_icode'][indices['od_data']:]  = od_neutrons['picode']

                output_file['od_data']['muon_pn'][indices['od_data']:]              = od_muons['muon_pn']
                output_file['od_data']['muon_impact'][indices['od_data']:]          = od_muons['muon_impact']
                output_file['od_data']['muon_energy'][indices['od_data']:]          = od_muons['muon_energy']
                output_file['od_data']['muon_initial'][indices['od_data']:]         = od_muons['muon_initial']
                output_file['od_data']['muon_direction'][indices['od_data']:]       = od_muons['muon_direction']

                output_file['od_totals']['neutrons_counted'][indices['od_totals']]   = od_neutrons['total']
                output_file['od_totals']['muons_simulated'][indices['od_totals']]    = od_muons['total']
                output_file['od_totals']['muon_parents'][indices['od_totals']]       = len(np.unique(od_muons['muon_energy']))


            #   RESNUCLEI DATASET

            output_file['resnuclei']['resnuclei'][indices['resnuclei']:]         = resnuclei_data
            output_file['resnuclei']['resnuclei_cu'][indices['resnuclei_cu']:]   = resnuclei_cu_data


            #   META DATASET

            now = datetime.now()

            output_file['meta']['hour'][indices['meta']] = int(now.strftime('%H'))
            output_file['meta']['minute'][indices['meta']] = int(now.strftime('%M'))
            output_file['meta']['second'][indices['meta']] = int(now.strftime('%S'))
            output_file['meta']['year'][indices['meta']] = int(now.strftime('%Y'))
            output_file['meta']['month'][indices['meta']] = int(now.strftime('%m'))
            output_file['meta']['day'][indices['meta']] = int(now.strftime('%d'))

            output_file['meta']['seed'][indices['meta']] = seed

        return True

    else:
        os.system('rm ' + output_filename)
        return False

def run_fluka():
    ''' Executes the command to run the simulation given everything else has been done'''

    run_string = yaml_card['source_path'] + 'rfluka -M 1 -e ./nEXOsim.exe ' + fluka_files['input_file'] 
    os.system(run_string)
    
def runsim():
    ''' The function for running the simulation from beginning to end'''
    
    ###     Step one: Make changes to the input file— Number of Muons

    change_number_of_muons(yaml_card['num_muons'], fluka_files['input_file'])

    ###     Step two: Compile and link the fluka files:

    link_and_compile(yaml_card['source_path'], fluka_files['mgdraw_file'], fluka_files['source_file'])

    ###     Step three: Looping from here on for however many repititions are demanded

    for i in range(yaml_card['reps']):

        ###     Make the phase space file

        make_phase_space_file(yaml_card['num_muons'])

        ###     Change the simulation seed

        seed = change_seed(fluka_files['input_file'])

        ###     Run the simulation

        run_fluka()

        ###     Deal with the output
        
        time_stamp = str(datetime.now())[10:16].replace(' ','').replace('-', '').replace(':', '')
        h5_filename = yaml_card['neutron_file'] + time_stamp + '.hdf5'

        if store_data_in_h5(h5_filename, seed):
            print('A neutron file was created')
        else:
            print('No neutron file was created')

        move_fluka_files(yaml_card['output_dir'])



runsim()








