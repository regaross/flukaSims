# Import all the constants and dictionaries from the __init__ file
from . import *
import subprocess


def link_and_compile():
    '''This function links and compiles the fortran (.f) files into an executable file that the simulation will run'''

    path_to_fluka = YAML_PARAMS['source_path']

    # Just make sure the path is well defined
    if not path_to_fluka[-1] == '/':
        path_to_fluka = path_to_fluka + '/'
    
    compile_string = path_to_fluka + 'fff'

    subprocess.run(compile_string + ' ' + FLUKA_JOB_FILES['mgdraw'], check=True, shell=True)
    subprocess.run(compile_string + ' ' + FLUKA_JOB_FILES['source_routine'], check=True, shell=True)

    # If we have a residual nuclei file, compile it too
    resnuc_compd = ''
    if FLUKA_JOB_FILES['resnuc'] != '':
        subprocess.run(compile_string + ' ' + FLUKA_JOB_FILES['resnuc'], check=True, shell=True)
        resnuc_compd = FLUKA_JOB_FILES['resnuc'][:-2] + '.o'

    FLUKA_JOB_FILES['executable'] = PATHS['workdir'] +'execute' + str(SEED) + '.exe'


    link_string = path_to_fluka + 'ldpmqmd -m fluka -o ' + FLUKA_JOB_FILES['executable'] + ' '
    mgd_compd = FLUKA_JOB_FILES['mgdraw'][:-2] + '.o'
    source_compd = FLUKA_JOB_FILES['source_routine'][:-2] + '.o'
    subprocess.run(link_string + mgd_compd + ' ' + resnuc_compd + ' ' + source_compd, check=True, shell=True)

    # Move the .o and .mod files...
    subprocess.run('mv ' + PATHS['SIF'] + '*.o' + ' ' + PATHS['SIF'] + PATHS['workdir'], check=True, shell=True)
    subprocess.run('mv ' + PATHS['SIF'] + '*.mod' + ' ' + PATHS['SIF'] + PATHS['workdir'], check=True, shell=True)
    subprocess.run('chmod +x ' + PATHS['SIF'] + '*.exe', check=True, shell=True)
    subprocess.run('mv ' + PATHS['SIF'] + '*.exe' + ' ' + PATHS['SIF'] + PATHS['workdir'], check=True, shell=True)


def run_fluka():
    ''' Executes the command to run the simulation given everything else has been done'''
    source_path = YAML_PARAMS['source_path']

    run_string = source_path + 'rfluka -M 1 -e ./' + FLUKA_JOB_FILES['executable'] + ' ' + FLUKA_JOB_FILES['input']
    subprocess.run(run_string, check=True, shell=True)

