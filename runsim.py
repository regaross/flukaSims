import os
import re
import numpy as np
import h5py as h5
from datetime import datetime
from src import muon_functions as mf
import yaml


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
        tpc_totals.create_dataset("resnuclei",(0,3), dtype=float, maxshape=(None,))

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

def store_events(tpc_filename, od_filename, resnuclei_filename, muon_filename, output_filename, meta_dict):
    ''' A function to collect the neutrons logged by the MGDRAW.f routine.
        First, the output file is parsed to collect the neutron data. These data are related to their source muons
        The data is appended to an hdf5 file. '''
    
    initialize_h5_file(output_filename)

    # Collect the time and date
    now = datetime.now()

    #                                           #
    #       TPC DATA COLLECTION AND FILING      #
    #                                           #
    
    # Order of output data in the simulation files:
    # ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK

    # Neutron attribute lists
    icode, ncase, jtrack, mreg, ltrack, etrack, xsco, ysco, zsco, cxtrck, cytrck, cztrck = [],[],[],[],[],[],[],[],[],[],[],[]
    # parent attribute lists
    picode, pjtrack, = [],[]

    with open(tpc_filename) as neutron_file:
        for line in neutron_file:
            # ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK
            temp = line.split()

            # Load the neutron lists
            icode.append(int(temp[0])); ncase.append(int(temp[1])); jtrack.append(int(temp[2]))
            mreg.append(int(temp[3])); ltrack.append(int(temp[4])); etrack.append(float(temp[5]))
            xsco.append(float(temp[6])); ysco.append(float(temp[7])); zsco.append(float(temp[8]))
            cxtrck.append(float(temp[9])); cytrck.append(float(temp[10]));  cztrck.append(float(temp[11]))

            # Load the parent lists
            picode.append(int(temp[12])); pjtrack.append(int(temp[14]))

    muenergy, impact,  initx, inity, initz, mucosx, mucosy, mucosz, pos_neg = [],[],[],[],[],[],[],[],[]
    with open(muon_filename) as muon_file:
        lines = muon_file.readlines()
        for mu in ncase:
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

    
    # Grab the file by the groups:
    file = h5.File(output_filename,'a')
    tpc_data = file['tpc_data']
    tpc_totals = file['tpc_totals']
    meta = file['meta']

    # Current length of the data in the file
    current_size = np.size(tpc_data['neutron_energy'])
    num_neutrons = len(ncase) # Number of elements to append
    current_meta_size = np.size(meta['year'])


    # resize the datasets

    for dset in tpc_data:
        tpc_data[dset].resize(current_size + num_neutrons, axis = 0)
    for dset in tpc_totals:
        if dset != 'resnuclei':
            tpc_totals[dset].resize(current_meta_size + 1, axis = 0)

    ## Residual nuclei data

    res_nuc_data = read_resnuclei_file(resnuclei_filename)
    current_res_nuc_size = np.size(tpc_totals['resnuclei'])
    res_nuc_entries = len(res_nuc_data)
    tpc_totals['resnuclei'].resize(current_res_nuc_size + res_nuc_entries, axis = 0)

    for i in range(res_nuc_entries):
        tpc_totals['resnuclei'][current_res_nuc_size + i] = res_nuc_data[i]


    meta['hour'][current_meta_size] = int(now.strftime('%H'))
    meta['minute'][current_meta_size] = int(now.strftime('%M'))
    meta['second'][current_meta_size] = int(now.strftime('%S'))
    meta['year'][current_meta_size] = int(now.strftime('%Y'))
    meta['month'][current_meta_size] = int(now.strftime('%m'))
    meta['day'][current_meta_size] = int(now.strftime('%d'))
    meta['seed'][current_meta_size] = meta_dict['seed']
    tpc_totals['neutrons_counted'][current_meta_size] = num_neutrons
    tpc_totals['muons_simulated'][current_meta_size] = meta_dict['number_of_muons']
    tpc_totals['muon_parents'][current_meta_size] = len(np.unique(ncase))

    # Now that the datasets have been resized, we must append the most recent sim data to the datasets.

    for i in range(num_neutrons):
        list_index = i
        dset_index = i + current_size

        tpc_data['muon_energy'][dset_index] = muenergy[list_index]
        tpc_data['muon_impact'][dset_index] = impact[list_index]
        tpc_data['muon_initial'][dset_index] = np.array([initx[list_index], inity[list_index], initz[list_index]])
        tpc_data['muon_direction'][dset_index] = np.array([mucosx[list_index], mucosy[list_index], mucosz[list_index]])
        tpc_data['muon_pn'][dset_index] = pos_neg[list_index]
        tpc_data['neutron_energy'][dset_index] = etrack[list_index]
        tpc_data['neutron_generation'][dset_index] = ltrack[list_index]
        tpc_data['neutron_xyz'][dset_index] = np.array([xsco[list_index], ysco[list_index], zsco[list_index]])
        tpc_data['neutron_direction'][dset_index] = np.array([cxtrck[list_index], cytrck[list_index], cztrck[list_index]])


    file.close()

    #                                           #
    #       OD DATA COLLECTION AND FILING       #
    #                                           #

    # Neutron attribute lists (empty from the TPC file)
    icode, ncase, jtrack, mreg, ltrack, etrack, xsco, ysco, zsco, cxtrck, cytrck, cztrck = [],[],[],[],[],[],[],[],[],[],[],[]
    # parent attribute lists
    picode, pjtrack, = [],[]

    
    with open(od_filename) as od_neutrons:
        for line in od_neutrons:
            temp = line.split()

            # Load the neutron lists
            icode.append(int(temp[0])); ncase.append(int(temp[1])); jtrack.append(int(temp[2]))
            mreg.append(int(temp[3])); ltrack.append(int(temp[4])); etrack.append(float(temp[5]))
            xsco.append(float(temp[6])); ysco.append(float(temp[7])); zsco.append(float(temp[8]))
            cxtrck.append(float(temp[9])); cytrck.append(float(temp[10]));  cztrck.append(float(temp[11]))

            # Load the parent lists
            picode.append(int(temp[12])); pjtrack.append(int(temp[14]))


    muenergy, impact,  initx, inity, initz, mucosx, mucosy, mucosz, pos_neg = [],[],[],[],[],[],[],[],[]
    with open(muon_filename) as muon_file:
        lines = muon_file.readlines()
        for mu in ncase:
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

    ## WE NOW HAVE THE DATA
    ## Must put it in an hdf5 file.
    

    # Open the h5 file in "append" mode so that it won't erase current data
    file = h5.File(output_filename,'a')

    # grabbing the file by the data groups
    od_data = file['od_data']
    od_totals = file['od_totals']

    # Current length of the data in the file
    current_size = np.size(od_data['neutron_energy'])
    num_neutrons = len(ncase) # Number of elements to append
    current_totals_size = np.size(od_totals['neutrons_counted'])

    for dset in tpc_data:
        # resize the datasets
        od_data[dset].resize(current_size + num_neutrons, axis = 0)
    for dset in tpc_totals:
        od_totals[dset].resize(current_totals_size + 1, axis = 0)

    od_totals['neutrons_counted'][current_totals_size] = num_neutrons
    od_totals['muons_simulated'][current_totals_size] = meta_dict['number_of_muons']
    od_totals['muon_parents'][current_totals_size] = len(np.unique(ncase))

    # Now that the datasets have been resized, we must append the most recent sim data to the datasets

    for i in range(num_neutrons):
        list_index = i
        dset_index = i + current_size

        od_data['muon_energy'][dset_index] = muenergy[list_index]
        od_data['muon_impact'][dset_index] = impact[list_index]
        od_data['muon_initial'][dset_index] = np.array([initx[list_index], inity[list_index], initz[list_index]])
        od_data['muon_direction'][dset_index] = np.array([mucosx[list_index], mucosy[list_index], mucosz[list_index]])
        od_data['muon_pn'][dset_index] = pos_neg[list_index]
        od_data['neutron_energy'][dset_index] = etrack[list_index]
        od_data['neutron_generation'][dset_index] = ltrack[list_index]
        od_data['neutron_xyz'][dset_index] = np.array([xsco[list_index], ysco[list_index], zsco[list_index]])
        od_data['neutron_direction'][dset_index] = np.array([cxtrck[list_index], cytrck[list_index], cztrck[list_index]])

    ## MUST CLOSE THE FILE
    file.close()

def main():
    yaml_file = open('simconfig.yaml')
    input_yaml = yaml.safe_load(yaml_file)

    os.system('echo RR: Loading in YAML Params')

    # Simulation Parameters
    Simulation = input_yaml.get('Simulation')
    num_muons = Simulation.get('Muons')
    make_new = Simulation.get('MakeNewFile')
    reps = Simulation.get('Repititions')
    muon_file = 'src/muon_file.txt'

    # Input Parameters
    Input = input_yaml.get('Input')
    input_path = Input.get('InputPath')
    input_file = input_path + Input.get('InputFile')
    source_routine = input_path + Input.get('SourceFile')
    mgdraw_file = input_path + Input.get('MGDrawFile')

    # Output Parameters
    Output = input_yaml.get('Output')
    output_dir = Output.get('OutputDir')
    neutron_file = Output.get('NeutronFile')
    progress_out = Output.get('ProgressOut')

    # Source Parameters
    source_path = input_yaml.get('Source').get('FlukaPath')

    yaml_file.close()

    ## MAKING A DICT TO CONTAIN META_DATA
    meta = dict()



    ## STEP ONE: CHANGE .inp FILE TO NUMBER OF MUONS
    os.system('echo RR: Changing number of muons in input file to ' + str(num_muons))
    space_string = '               '
    num_spaces = len(space_string) - len(str(num_muons))
    num_string = 'START' + space_string[:num_spaces] + str(num_muons)
    os.system('sed -i \'s/^START.*/' + num_string + '/\' ' + input_file)
    meta['number_of_muons'] = num_muons
    

    ## STEP THREE: COMPILE THE FILES
    compile_string = source_path + 'fff'
    os.system('echo RR: COMPILING FILES >> ' + progress_out)
    os.system(compile_string + ' ' + mgdraw_file + ' >> '+ progress_out)
    os.system(compile_string + ' ' + source_routine + ' >> ' + progress_out)

    ## STEP FOUR: LINK THE COMPILED FILES
    os.system('echo RR: Compiling the user routines')
    link_string = source_path + 'ldpmqmd -m fluka -o nEXOsim.exe '
    mgd_compd = mgdraw_file[:-1] + 'o'
    source_compd = source_routine[:-1] + 'o'
    os.system('echo RR: Linking the user routines')
    os.system(link_string + mgd_compd + ' ' + source_compd + '>> ' + progress_out)

    ## STEP FIVE (optional): MAKE THE PHASE SPACE FILE
    for i in range(reps):
        if make_new or not os.path.exists( muon_file):
            os.system('echo RR: Making the phase space file')
            make_phase_space_file(num_muons = num_muons, filename = muon_file)

        ## STEP SIX: CHANGE SEED IN .inp FILE
        seed = str(np.random.randint(0,999999))
        os.system('echo RR: Changing seed in input file to: ' + seed)
        space_string = '                      '
        num_spaces = len(space_string) - len(str(seed))
        num_string = 'RANDOMIZ' + space_string[:num_spaces] + seed
        os.system('sed -i \'s/^RANDOMIZ.*/' + num_string + '/\' ' + input_file)
        meta['seed'] = int(seed)
        
        ## STEP SEVEN: RUN THE SIM
        os.system('echo RR: Running the simulation')
        run_string = source_path + 'rfluka -M 1 -e ./nEXOsim.exe ' + input_file 
        os.system(run_string)

        ## STEP EIGHT: MOVE DATA TO H5 FILE
        os.system('echo RR: Copying neutron data to h5 file')
        tpc_filename = input_file[:-4] + "001_fort.72"
        od_filename = input_file[:-4] + "001_fort.70"
        res_nuc_filename = input_file[:-4] + "001_fort.97"
        

        time_stamp = str(datetime.now())[5:16].replace(' ','').replace('-', '').replace(':', '')

        neutron_filename = neutron_file + time_stamp + '.hdf5'
        try:
            store_events(tpc_filename, od_filename, muon_file, neutron_filename, meta)
        except:
            os.system('echo RR: The neutron hdf5 file WAS NOT CREATED\; probably because no neutrons were generated')
        

        ## STEP NINE: MOVING OUTPUT FILES TO SPECIFIED DIRECTORY
        os.system('echo RR: Moving output files to output directory: ' + output_dir)
        try:
            os.system('mkdir ' + output_dir)
        except: pass

        sub_dir = neutron_file + time_stamp
        last_dir = output_dir + sub_dir + '/'
        try:
            os.system('mkdir ' + last_dir)
            os.system('mv *fort* ' + last_dir)
            os.system('mv *lis* *tab* ' + last_dir)
            os.system('mv *.hdf5 ' + output_dir)
            os.system('mv *.log *.err *.out *ran* *dump *fort* *.txt ' + last_dir)
        except: pass
        os.system('echo RR: Copying input file to output dir for future reference')
        os.system('cp ' + input_file + ' ' + last_dir)

    ## STEP TEN: REMOVE COMPILED AND UNNECESSARY FILES
    os.system('echo RR: Removing unnecessary files')
    os.system('rm *.o *.exe *.mod')
    os.system('echo RR: DONE!')


main()