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

def store_data(neutron_fort_filename, muon_filename, output_filename, meta_dict):
    now = datetime.now()

    #  instantiate neutron attribute lists
    muon_numbers, generation, energy, xsco, ysco, zsco, cosx, cosy, cosz = [],[],[],[],[],[],[],[],[]

    with open(neutron_fort_filename) as neutron_file:
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

    print(muon_numbers)
    print(cosz)


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
        store_data(neutron_file, muon_file, 'output', meta)

main()



