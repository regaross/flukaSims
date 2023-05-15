import os
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seed', type=int, dest='first_seed')
args = parser.parse_args()

stamp = int(args.first_seed)

def change_muon_filepath(stamp):
    '''Changes the path to the muon_file in the provided fluka source file'''

    source_name = 'musource' + str(stamp) + '.f'


    with open('muon_from_file.f', 'r') as source:
        lines = source.readlines()
    
    replace_string = '      call read_phase_space_file(\"'+ source_name + '\", \'GeV\', \'m\', phase_space_entry, .true. , nomore )'
    lines[527] = replace_string

    with open(source_name, 'w') as source:
            source.writelines(lines)

change_muon_filepath(stamp)