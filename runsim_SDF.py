import os
import numpy as np
import h5py as h5
from datetime import datetime
from src import muon_functions as mf
import yaml


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


def store_events(neutron_filename, muon_filename, output_filename):
    ''' A function to collect the neutrons logged by the MGDRAW.f routine.
        First, the output file is parsed to collect the neutron data. These data are related to their source muons
        The data is appended to an hdf5 file. '''

    now = datetime.now()

    # Neutron attribute lists
    muon, generation, energy, xsco, ysco, zsco, cosx, cosy, cosz = [],[],[],[],[],[],[],[],[]

    with open(neutron_filename) as neutron_file:
        for line in neutron_file:
            temp = line.split()
            muon.append(int(temp[0]))
            energy.append(float(temp[1]))
            generation.append(int(temp[2]))
            xsco.append(float(temp[4]))
            ysco.append(float(temp[5]))
            zsco.append(float(temp[6]))
            cosx.append(float(temp[7]))
            cosy.append(float(temp[8]))
            cosz.append(float(temp[9]))


    muenergy, impact,  initx, inity, initz, mucosx, mucosy, mucosz, pos_neg = [],[],[],[],[],[],[],[],[]
    with open(muon_filename) as muon_file:
        lines = muon_file.readlines()
        for mu in muon:
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
    # First, we check whether or not the file exists. If it does, we resize the datasets to accomodate the new data.

    # If the file does not exist
    if not os.path.isfile(output_filename): # The file must be created.
        file = h5.File(output_filename,'a')
        # Create a group for the data
        data = file.create_group('data')
        # Create a group for the meta data (to dilineate different simulation data sets)
        meta = file.create_group('meta')

        # So too must the datasets be created and instantiated with zero size.
        data.create_dataset("muon_energy", (0,), dtype=float, maxshape=(None,))
        data.create_dataset("muon_impact", (0,), dtype=float, maxshape=(None,))
        data.create_dataset("muon_initial", (0,3), dtype=float, maxshape=(None,3))
        data.create_dataset("muon_direction", (0,3), dtype=float, maxshape=(None,3))
        data.create_dataset("muon_pn", (0,), dtype=int, maxshape=(None,))

        data.create_dataset("neutron_energy", (0,), dtype=float, maxshape=(None,))
        data.create_dataset("neutron_generation", (0,), dtype=int, maxshape=(None,))
        data.create_dataset("neutron_xyz", (0,3), dtype=float, maxshape=(None,3))
        data.create_dataset("neutron_direction", (0,3), dtype=float, maxshape=(None,3))

        meta.create_dataset("year",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("month",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("day",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("hour",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("minute",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("second",(0,), dtype=int, maxshape=(None,))
        meta.create_dataset("how_many",(0,), dtype=int, maxshape=(None,))


        file.close()

    # Open the h5 file in "append" mode so that it won't erase current data
    file = h5.File(output_filename,'a')
    # grabbing the file by the data group
    data = file['data']
    meta = file['meta']

    # Current length of the data in the file
    current_size = np.size(data['neutron_energy'])
    n_points = len(muon) # Number of elements to append
    current_meta_size = np.size(meta['year'])

    for dset in data:
        # resize the datasets
        data[dset].resize(current_size + n_points, axis = 0)
    for dset in meta:
        meta[dset].resize(current_meta_size + 1, axis = 0)

    meta['hour'][current_meta_size] = int(now.strftime('%H'))
    meta['minute'][current_meta_size] = int(now.strftime('%M'))
    meta['second'][current_meta_size] = int(now.strftime('%S'))
    meta['year'][current_meta_size] = int(now.strftime('%Y'))
    meta['month'][current_meta_size] = int(now.strftime('%m'))
    meta['day'][current_meta_size] = int(now.strftime('%d'))
    meta['how_many'][current_meta_size] = n_points

    # Now that the datasets have been resized, we must append the most recent sim data to the datasets.


    for i in range(n_points):
        list_index = i
        dset_index = i + current_size

        data['muon_energy'][dset_index] = muenergy[list_index]
        data['muon_impact'][dset_index] = impact[list_index]
        data['muon_initial'][dset_index][:] = [initx[list_index], inity[list_index], initz[list_index]]
        data['muon_direction'][dset_index] = [mucosx[list_index], mucosy[list_index], mucosz[list_index]]
        data['muon_pn'][dset_index] = pos_neg[list_index]
        data['neutron_energy'][dset_index] = energy[list_index]
        data['neutron_generation'][dset_index] = generation[list_index]
        data['neutron_xyz'][dset_index][:] = [xsco[list_index], ysco[list_index], zsco[list_index]]
        data['neutron_direction'][dset_index][:] = [cosx[list_index], cosy[list_index], cosz[list_index]]

    ## MUST CLOSE THE FILE
    file.close()


def main():
    yaml_file = open('simconfig_SDF.yaml')
    input_yaml = yaml.safe_load(yaml_file)

    os.system('echo RR: Loading in YAML Params')

    # Simulation Parameters
    Simulation = input_yaml.get('Simulation')
    num_muons = Simulation.get('Muons')
    make_new = Simulation.get('MakeNewFile')
    reps = Simulation.get('Repititions')
    score_od = Simulation.get('ScoreOD')
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

    ## STEP ONE: CHANGE .inp FILE TO NUMBER OF MUONS
    os.system('echo RR: Changing number of muons in input file')
    space_string = '               '
    num_spaces = len(space_string) - len(str(num_muons))
    num_string = 'START' + space_string[:num_spaces] + str(num_muons)
    os.system('sed -i \'s/^START.*/' + num_string + '/\' ' + input_file)

    ## STEP TWO: MAKE SURE GEOMETRY SCORING GEOMETRY IS RIGHT FOR NEUTRON COUNT
    os.system('echo RR: Making sure scoring region is right')
    if score_od: 
        geo_num = 3  # Region number for OD water tank
        os.system('echo RR: Counting neutrons in OD Water Tank')

    else: 
        geo_num = 9         # Region number for TPC inside
        os.system('echo RR: Counting neutrons in TPC')

    mg_string = 'IF (MREG .EQ. ' + str(geo_num)
    os.system('sed -i \'s/IF (MREG .EQ. [0-9]/'+ mg_string + '/g\' ' + mgdraw_file)

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
        
        ## STEP SIX: RUN THE SIM
        os.system('echo RR: Running the simulation')
        run_string = source_path + 'rfluka -M 1 -e ./nEXOsim.exe ' + input_file 
        os.system(run_string)

        ## STEP SEVEN: MOVE DATA TO H5 FILE
        os.system('echo RR: Copying neutron data to h5 file')
        try:
            store_events('nEXO_OD001_fort.99', muon_file, neutron_file)
        except:
            os.system('echo RR: The neutron file WAS NOT CREATED\; probably because no neutrons were generated')
        

        ## STEP SIX: MOVING OUTPUT FILES TO SPECIFIED DIRECTORY
        os.system('echo RR: Moving output files to output directory: ' + output_dir)
        try:
            os.system('mkdir ' + output_dir)
        except: pass

        sub_dir = str(datetime.now())[:16].replace(' ', '_')
        last_dir = output_dir + sub_dir + '/'
        try:
            os.system('mkdir ' + last_dir)
            os.system('mv *fort* ' + last_dir)
            os.system('mv *.log *.err *.out *ran* *dump *fort* *.txt ' + last_dir)
        except: pass
        os.system('echo RR: Copying input file to output dir for future reference')
        os.system('cp ' + input_file + ' ' + last_dir)

    ## STEP SEVEN: REMOVE COMPILED AND UNNECESSARY FILES
    os.system('echo RR: Removing unnecessary files')
    os.system('rm *.o *.exe *.mod')
    os.system('echo RR: DONE!')


main()