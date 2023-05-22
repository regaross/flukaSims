import simtools
from simtools import run_tools as rt
from simtools import muon_functions as mf
import argparse
import numpy as np
import yaml

## Loading in the parameters from the simconfig.yaml file


parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seed', type=int, dest='first_seed')
args = parser.parse_args()

### Seed for random number generation from random number on system
STAMP = int(args.first_seed)
np.random.seed(STAMP)
stamp = STAMP

mf.set_seed(STAMP)
rt.set_seed(STAMP)

with open('simconfig.yaml') as yaml_file:

    input_yaml = yaml.safe_load(yaml_file)

    Simulation = input_yaml.get('Simulation')
    Input = input_yaml.get('Input')
    Output = input_yaml.get('Output')
    input_path = Input.get('InputPath')

    yaml_card = {
        # Simulation Parameters
        'num_muons'     : Simulation.get('Muons'),
        'intersecting'  : Simulation.get('Intersecting'),
        'make_new'      : Simulation.get('MakeNewFile'),
        'roi_radius'    : Simulation.get('ROI_Radius'),
        'roi_height'    : Simulation.get('ROI_Height'),


        # Input Parameters
        'input_file'        : input_path + Input.get('InputFile'),
        'source_routine'    : input_path + Input.get('SourceFile'),
        'mgdraw_file'       : input_path + Input.get('MGDrawFile'),

        # Output Parameters
        'output_dir'    : Output.get('OutputDir'),
        'neutron_file'  : Output.get('NeutronFile'),
        'progress_out'  : Output.get('ProgressOut'),

        # Source Parameters
        'source_path' : input_yaml.get('Source').get('FlukaPath'),
    }

fluka_files = {
    'tpc_neutron_file'      :   'input' + str(stamp) +  '001_fort.72',
    'od_neutron_file'       :   'input' + str(stamp) +  '001_fort.70',
    'res_nuclei_file'       :   'input' + str(stamp) +  '001_fort.97',
    'res_nuclei_cu_file'    :   'input' + str(stamp) +  '001_fort.94',
    'mgdraw_file'           :   'mgdrw' + str(stamp) + '.f',
    'source_file'           :   'musource' + str(stamp) + '.f',  
    'input_file'            :   'input' + str(stamp) + '.inp',
    'muon_file'             :   'muons' + str(stamp) + '.txt',
    'executable'            :   'exe' + str(stamp) + '.exe'
}


rt.runsim(STAMP, yaml_card, fluka_files)