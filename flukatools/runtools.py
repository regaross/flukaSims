# Import all the constants and dictionaries from the __init__ file
from . import *
from os import system


def link_and_compile():
    '''This function links and compiles the fortran (.f) files into an executable file that the simulation will run'''

    path_to_fluka = YAML_PARAMS['source_path']

    # Just make sure the path is well defined
    if not path_to_fluka[-1] == '/':
        path_to_fluka = path_to_fluka + '/'
    
    compile_string = path_to_fluka + 'fff'

    system(compile_string + ' ' + FLUKA_FILES['mgdraw'] )
    system(compile_string + ' ' + FLUKA_FILES['source_routine']  )

    link_string = path_to_fluka + 'ldpmqmd -m fluka -o ' + FLUKA_FILES['executable'] + ' '
    mgd_compd = FLUKA_FILES['mgdraw'][:-2] + '.o'
    source_compd = FLUKA_FILES['source_routine'][:-2] + '.o'
    system(link_string + mgd_compd + ' ' + source_compd )

def run_fluka():
    ''' Executes the command to run the simulation given everything else has been done'''

    source_path = YAML_PARAMS['source_path']

    run_string = source_path + 'rfluka -M 1 -e ./' + FLUKA_FILES['executable'] + ' ' + FLUKA_FILES['input_file']
    system(run_string)

