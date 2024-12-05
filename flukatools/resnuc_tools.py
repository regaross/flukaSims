#!/usr/bin/python

__version__ = 0
__author__ = 'Regan Ross'
## Last Edited Sept 18, 2024

'''
resnuc_tools.py

Contact:
Regan Ross
regan.ross@mail.mcgill.ca

Tools for reading in and correlating events between the primary muon and the residual nuclei and fragments created

'''

# Need to import muons so that we can use the muon class
from .muons import Muon
import re
import numpy as np
from os import listdir


# class RNevent:
#     '''Residual Nucleus event : A class for storing the event information output by the mgdraw.f and usrrnc.f files and relating 
#     them to the primary muon'''

#     def __init__(self, seed : int, icode : int, region : int, jtrack : int, ltrack : int, secondaries : tuple, frags : tuple, mgrnc : tuple, usrrnc : tuple, muon : Muon):
#         '''Defines an instance of a Residual Nucleus event with all the attributes that are saved to the output files from mgdraw.f and usrrnc.f '''

#         self.seed           = seed
#         self.icode          = icode
#         self.region         = region
#         self.jtrack         = jtrack
#         self.ltrack         = ltrack
#         self.secondaries    = secondaries
#         self.frags          = frags
#         self.mgrnc          = mgrnc
#         self.usrrnc         = usrrnc
#         self.muon           = muon

#     def __init__(self, event_dict : dict, muon : 'Muon'):
#         '''Defines an instance of a Residual Nucleus event with all the attributes that are saved to the output files from mgdraw.f and usrrnc.f '''

#         self.seed           = event_dict['seed']
#         self.icode          = event_dict['icode']
#         self.region         = event_dict['region']
#         self.jtrack         = event_dict['jtrack']
#         self.ltrack         = event_dict['ltrack']
#         self.secondaries    = event_dict['secondaries']
#         self.frags          = event_dict['frags']
#         self.mgrnc          = event_dict['mgrnc']
#         self.usrrnc         = event_dict['usrrnc']
#         self.muon           = muon


def read_rn_event_sequence(sequence : str) -> dict:
    '''Reads a sequence from an mgdraw.f or usrrnc.f output file and returns a dictionary with the appropriate attributes for an RNevent.
    This function is only for parsing the event sequence string'''

    icode = int(re.search(r'ICODE\s*(\d{3})', sequence).group(1))
    primary = int(re.search(r'history:\s*(\d+)', sequence).group(1))
    region = int(re.search(r'region\s*(\d+)', sequence).group(1))
    jtrack = int(re.search(r'LTRACK:\s+(\d\d?)\s+(\d\d?)', sequence).group(1))
    ltrack = int(re.search(r'LTRACK:\s+(\d\d?)\s+(\d\d?)', sequence).group(2))
    secondaries = np.array(re.search(r'(?<=ID:)\s*(\d\d?\s*)*', sequence).group(0).split())
    n_fragments = int(re.search(r'\d(?=\sfragments)', sequence).group(0))

    masses = re.findall(r'(?<=A = )\s*\d*', sequence)
    atonums = re.findall(r'(?<=Z = )\s*\d*', sequence)

    if n_fragments > 0:
        frag_A = int(re.findall(r'(?<=A = )\s*\d*', sequence)[0].strip())
        frag_Z = int(re.findall(r'(?<=Z = )\s*\d*', sequence)[0].strip())
    else:
        frag_A, frag_Z = 0, 0

    frags = (frag_Z, frag_A)

    mgres_A = int(masses[1].strip())
    mgres_Z = int(atonums[1].strip())

    if type(re.search(r'USRRNC', sequence)) is not type(None):
        rncres_A = int(masses[2].strip())
        rncres_Z = int(atonums[2].strip())
    else:
        rncres_A, rncres_Z = 0, 0

    mgrnc = (mgres_Z, mgres_A)
    usrrnc = (rncres_Z, rncres_A)


    return {'icode'         : icode,
            'primary'       : primary,
            'region'        : region, 
            'jtrack'        : jtrack,
            'ltrack'        : ltrack,
            'secondaries'   : secondaries,
            'n_fragments'   : n_fragments, 
            'frags'         : frags,
            'mgrnc'         : mgrnc,
            'usrrnc'        : usrrnc}


def parse_rnc_file(filepath : str, only_residuals = True) -> list:
    '''Parses the given residual nuclei file and makes RNevents out of each entry
    Returns a list of these events pointing to their respective muon
    if only_residuals, then only the events that produce residual nuclei will be counted'''

    events = []
    seed = int(re.search(r'(\d*)(?=\.)', filepath).group(0))

    with open(filepath, 'r') as file:
        content = file.read()
            
        # Split the content by double newlines (indicating empty lines)
        sequences = content.strip().split('\n\n')
        for sequence in sequences:
            
            # Read in the sequence and retrieve the information for the event
            event_dict = read_rn_event_sequence(sequence)
            event_dict['seed'] = seed
            event_dict['muon_id'] = int(100000 * seed) + event_dict['primary']

            if only_residuals and event_dict['usrrnc'][1] == 0 and event_dict['mgrnc'][1] == 0 and event_dict['n_fragments'] == 0:
                continue

            events.append(event_dict)

    return events


def parse_rnc_dir(path_to_rnc_dir, only_residuals = True) -> tuple:


    events = []
    error_seeds = []

    event_files = listdir(path_to_rnc_dir)
    seeds = [int(re.search(r'(\d*)(?=\.)', event_file).group(0)) for event_file in event_files]

    for seed in seeds:
        resnuc_filename = path_to_rnc_dir + 'new_resnuclei_tpc' + str(seed) + '.asc'

        try:
            temp_events = parse_rnc_file(resnuc_filename, only_residuals)
            events.extend(temp_events)
        except:
            print('Error. Could not parse the following:')
            print(resnuc_filename)
            error_seeds.append(seed)
            continue

    return events, error_seeds


def get_muons_from_file(muon_file) -> list:
    '''Returns the list of muon objects from a particular muon file'''

    muarray = np.loadtxt(muon_file)
    muons = {}
    seed = int(re.search(r'(\d*)(?=\.)', muon_file).group(0))

    for i in range(len(muarray)):
        id = int(100000 * seed) + i + 1
        this_muon = {'pos_neg'      :   int(muarray[i][0]),
                     'energy'       :   muarray[i][1],
                     'init_pos'     :   muarray[i][2:5],
                     'cosines'      :   muarray[i][5:8],
                     'xe137'         :   0,
                     'cu62'         :   0,
                     'cu64'         :   0,
                     'cu66'         :   0,
                     'i130'         :   0,
                     'i132'         :   0,
                     'i134'         :   0,
                     'i135'         :   0,
                     'i136'         :   0,
                     'i137'         :   0,

                    }
        muons.update({id : this_muon})

    return muons

def get_muons_from_dir(path_to_muon_dir):

    muons = {}

    muon_files = listdir(path_to_muon_dir)

    for file in muon_files:
        muons.update(get_muons_from_file(path_to_muon_dir + file))

    return muons



    




            