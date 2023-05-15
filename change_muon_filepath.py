import os
import argparse
import time

### From SLURM environment variables
slurm_job_id = int(os.environ["SLURM_ARRAY_JOB_ID"])
slurm_task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
slurm_prefix = 'simrun-' + str(os.environ["SLURM_ARRAY_JOB_ID"]) + '-' + str(os.environ["SLURM_ARRAY_TASK_ID"])


parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seed', type=int, dest='first_seed')
args = parser.parse_args()

stamp = int(args.first_seed)

def change_muon_filepath(stamp):
    '''Changes the path to the muon_file in the provided fluka source file'''

    source_name = 'musource' + str(stamp) + '.f'

    if os.path.isfile(source_name):
        return 
    else:

        with open('muon_from_file.f', 'r') as source:
            lines = source.readlines()
        
        replace_string = '      call read_phase_space_file(\"'+ source_name + '\", \'GeV\', \'m\', phase_space_entry, .true. , nomore )'
        lines[527] = replace_string

        with open(source_name, 'w') as source:
            source.writelines(lines)

change_muon_filepath(stamp)