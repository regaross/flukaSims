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


def copy_input_files(stamp):
    stamp = str(stamp)
    if os.path.isfile('input' + str(stamp) + '.inp'):
        return 
    else:
        os.system('cp nEXO_OD.inp input' + stamp + '.inp')
        os.system('cp mgdraw_neutron_count.f mgdrw' + stamp + '.f')

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

def link_and_compile(path_to_fluka, stamp, path_to_sim = '/gpfs/slac/staas/fs1/g/exo/exo_data8/exo_data/users/rross/flukaSims/'):
    '''Links and compiles the fluka routines for the fluka executable'''

    if not path_to_fluka[-1] == '/':
        path_to_fluka = path_to_fluka + '/'
    
    compile_string = path_to_fluka + 'fff'

    os.system(compile_string + ' ' + path_to_sim + 'mgdrw' + stamp + '.f' )
    os.system(compile_string + ' ' + path_to_sim + 'musource' + stamp + '.f'  )

    link_string = path_to_fluka + 'ldpmqmd -m fluka -o exe'  + stamp + '.exe '
    mgd_compd = 'mgdrw' + stamp + '.o'
    source_compd = 'musource' + stamp + '.o'
    os.system(link_string + path_to_sim + mgd_compd + ' ' + path_to_sim + source_compd )

change_muon_filepath(stamp)

copy_input_files(stamp)

try:
    link_and_compile('/usr/local/fluka/bin/', stamp)
    time.sleep(2)
except:
    print('\n\nFailed at compiling')
