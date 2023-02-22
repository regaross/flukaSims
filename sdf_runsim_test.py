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

def change_inp_muons(num_muons, input_file):
    ''' A small method to change the integer number of source muons in the input file '''

    os.system('echo RR: Changing number of muons in input file to ' + str(num_muons))
    space_string = '               '
    num_spaces = len(space_string) - len(str(num_muons))
    num_string = 'START' + space_string[:num_spaces] + str(num_muons)
    os.system('sed -i \'s/^START.*/' + num_string + '/\' ' + input_file)

def change_scoring_region(region_number, mgdraw_file):
    ''' Adjusts the scoring region integer in the input file'''

    os.system('echo RR: Making sure scoring region is right: 3 for OD, 9 for TPC')
    if region_number == 3 or region_number == 9:
        mg_string = 'IF (MREG .EQ. ' + str(region_number)
        os.system('sed -i \'s/IF (MREG .EQ. [0-9]/'+ mg_string + '/g\' ' + mgdraw_file)
    else:
        os.system('echo RR: WRONG REGION NUMBER!')

def compile_fluka_files(source_path, mgdraw_file, source_routine, output = 'compile_output.txt'):
    ''' Compiles the fluka routines '''

    compile_string = source_path + 'fff'
    os.system('echo RR: COMPILING FILES >> ' + output)
    os.system(compile_string + ' ' + mgdraw_file + ' >> '+ output)
    os.system(compile_string + ' ' + source_routine + ' >> ' + output)

def link_fluka_files(source_path, mgdraw_file, source_routine, output = 'compile_output.txt'):
    ''' Links fluka files for running the simulation'''
    
    os.system('echo RR: Compiling the user routines')
    link_string = source_path + 'ldpmqmd -m fluka -o nEXOsim.exe '
    mgd_compd = mgdraw_file[:-1] + 'o'
    source_compd = source_routine[:-1] + 'o'
    os.system('echo RR: Linking the user routines')
    os.system(link_string + mgd_compd + ' ' + source_compd + '>> ' + output)

def change_sim_seed(input_file):
    '''Changes the seed in the input file'''

    seed = str(np.random.randint(0,999999))
    os.system('echo RR: Changing seed in input file to: ' + seed)
    space_string = '                      '
    num_spaces = len(space_string) - len(str(seed))
    num_string = 'RANDOMIZ' + space_string[:num_spaces] + seed
    os.system('sed -i \'s/^RANDOMIZ.*/' + num_string + '/\' ' + input_file)

    return int(seed)

def run_sim(source_path, input_file):
    ''' Just runs the fluka simulation '''

    os.system('echo RR: Running the simulation')
    run_string = source_path + 'rfluka -M 1 -e ./nEXOsim.exe ' + input_file 
    os.system(run_string) 

def store_data(input_filename, muon_filename, output_filename, meta_dict):
    now = datetime.now()

    #  instantiate neutron attribute lists
    muon_numbers, generation, energy, xsco, ysco, zsco, cosx, cosy, cosz = [],[],[],[],[],[],[],[],[]
    fort_99_filename = input_filename[:-4] + "001_fort.99"
    with open(fort_99_filename) as neutron_file:
        for line in neutron_file:
            temp = line.split()
            muon_numbers.append(int(temp[0]))
            energy.append(float(temp[1]))
            generation.append(int(temp[2]))
            xsco.append(float(temp[4]))
            ysco.append(float(temp[5]))
            zsco.append(float(temp[6]))
            cosx.append(float(temp[7]))
            cosy.append(float(temp[8]))
            cosz.append(float(temp[9]))

            # THE DATA ARE HERE; These lists are full.

    # Instantiate muon attribute lists
    muenergy, impact,  initx, inity, initz, mucosx, mucosy, mucosz, pos_neg = [],[],[],[],[],[],[],[],[]
    with open(muon_filename) as muon_file:
        lines = muon_file.readlines()
        for mu in muon_numbers:
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
    ## Must put it in a file. 
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
    # number of neutrons created in file
    meta.create_dataset("neutrons_counted",(0,), dtype=int, maxshape=(None,))
    # number of muons simulated
    meta.create_dataset("muons_simulated",(0,), dtype=int, maxshape=(None,))
    # number of muons responsible for creating neutrons
    meta.create_dataset("muon_parents",(0,), dtype=int, maxshape=(None,))
    # the integer seed used in the simulation
    meta.create_dataset("seed",(0,), dtype=int, maxshape=(None,))
    # the region scored in the simulation
    meta.create_dataset("region",(0,), dtype=int, maxshape=(None,))

    # Current length of the data in the file
    current_size = len(data['neutron_energy'])
    num_neutrons = len(muon_numbers) # Number of elements to append
    current_meta_size = len(meta['year'])

    for dset in data:
        # resize the datasets
        data[dset].resize(current_size + num_neutrons, axis = 0)
    for dset in meta:
        meta[dset].resize(current_meta_size + 1, axis = 0)

    meta['hour'][current_meta_size] = int(now.strftime('%H'))
    meta['minute'][current_meta_size] = int(now.strftime('%M'))
    meta['second'][current_meta_size] = int(now.strftime('%S'))
    meta['year'][current_meta_size] = int(now.strftime('%Y'))
    meta['month'][current_meta_size] = int(now.strftime('%m'))
    meta['day'][current_meta_size] = int(now.strftime('%d'))
    meta['neutrons_counted'][current_meta_size] = num_neutrons
    meta['muons_simulated'] = int(meta_dict['number_of_muons'])
    meta['muon_parents'][current_meta_size] = len(np.unique(muon_numbers))

    meta['seed'] = meta_dict['seed']
    meta['region'] = meta_dict['scoring']


    for i in range(num_neutrons):
        list_index = i
        dset_index = i + current_size

        data['muon_energy'][dset_index] = muenergy[list_index]
        data['muon_impact'][dset_index] = impact[list_index]
        data['muon_initial'][dset_index] = [initx[list_index], inity[list_index], initz[list_index]]
        data['muon_direction'][dset_index] = [mucosx[list_index], mucosy[list_index], mucosz[list_index]]
        data['muon_pn'][dset_index] = pos_neg[list_index]
        data['neutron_energy'][dset_index] = float(energy[list_index])
        data['neutron_generation'][dset_index] = generation[list_index]
        data['neutron_xyz'][dset_index] = [xsco[list_index], ysco[list_index], zsco[list_index]]
        data['neutron_direction'][dset_index] = [float(cosx[list_index]), float(cosy[list_index]), float(cosz[list_index])]

    print(data['muon_direction'][0])

    file.close()


def main():

    # Collect the parameters from the yaml file

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

    # Define some of the meta parameters

    meta = dict.fromkeys(['number_of_muons', 'scoring', 'seed'])

    meta['number_of_muons'] = num_muons

    if score_od:
        meta['scoring'] = 3
    else: 
        meta['scoring'] = 9

    # Change the number of source muons in the input file
    change_inp_muons(num_muons, input_file)
    
    # Change the scoring region
    change_scoring_region(meta['scoring'], mgdraw_file)

    # Compile the files
    compile_fluka_files(source_path, mgdraw_file, source_routine)

    # Link the files
    link_fluka_files(source_path, mgdraw_file, source_routine)

    for i in range(reps):
        if make_new or not os.path.exists(muon_file):
            os.system('echo RR: Making the phase space file')
            make_phase_space_file(num_muons = num_muons, filename = muon_file)

        # Change seed in inp file
        meta['seed'] = change_sim_seed(input_file)

        # Run the simulation
        run_sim(source_path, input_file)

        # Move the data to the output file
        store_data(input_file, muon_file, 'output.h5', meta)

main()



