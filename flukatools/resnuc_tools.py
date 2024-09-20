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


class RNevent:
    '''Residual Nucleus event : A class for storing the event information output by the mgdraw.f and usrrnc.f files and relating 
    them to the primary muon'''

    def __init__(self, seed : int, icode : int, region : int, jtrack : int, ltrack : int, secondaries : tuple, frags : tuple, mgrnc : tuple, usrrnc : tuple, muon : Muon):
        '''Defines an instance of a Residual Nucleus event with all the attributes that are saved to the output files from mgdraw.f and usrrnc.f '''

        self.seed           = seed
        self.icode          = icode
        self.region         = region
        self.jtrack         = jtrack
        self.ltrack         = ltrack
        self.secondaries    = secondaries
        self.frags          = frags
        self.mgrnc          = mgrnc
        self.usrrnc         = usrrnc
        self.muon           = muon

    def __init__(self, event_dict : dict, muon : 'Muon'):
        '''Defines an instance of a Residual Nucleus event with all the attributes that are saved to the output files from mgdraw.f and usrrnc.f '''

        self.seed           = event_dict['seed']
        self.icode          = event_dict['icode']
        self.region         = event_dict['region']
        self.jtrack         = event_dict['jtrack']
        self.ltrack         = event_dict['ltrack']
        self.secondaries    = event_dict['secondaries']
        self.frags          = event_dict['frags']
        self.mgrnc          = event_dict['mgrnc']
        self.rncresidual    = event_dict['rncresidual']
        self.muon           = muon


def read_rn_event_sequence(sequence : str) -> dict:
    '''Reads a sequence from an mgdraw.f or usrrnc.f output file and returns a dictionary with the appropriate attributes for an RNevent.
    This function is only for parsing the event sequence string'''

    icode = int(re.search(r'ICODE\s*(\d{3})', sequence).group(1))
    primary = int(re.search(r'history:\s*(\d+)', sequence).group(1))
    region = int(re.search(r'region\s*(\d+)', sequence).group(1))
    jtrack = int(re.search(r'LTRACK:\s+(\d\d?)\s+(\d\d?)', sequence).group(1))
    ltrack = int(re.search(r'LTRACK:\s+(\d\d?)\s+(\d\d?)', sequence).group(2))
    secondaries = tuple(re.search(r'(?<=ID:)\s*(\d\d?\s*)*', sequence).group(0).split())
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


def parse_rnc_file(filepath : str, muon_filepath :str, only_residuals = True) -> list:
    '''Parses the given residual nuclei file and makes RNevents out of each entry
    Returns a list of these events pointing to their respective muon
    if only_residuals, then only the events that produce residual nuclei will be counted'''

    events = []
    muons = []
    seed = int(re.search(r'(\d*)(?=\.)', filepath).group(0))
    prim = 0
    this_muon = Muon(seed, prim)

    with open(filepath, 'r') as file:
        content = file.read()
            
        # Split the content by double newlines (indicating empty lines)
        sequences = content.strip().split('\n\n')
        for sequence in sequences:
            
            # Read in the sequence and retrieve the information for the event
            event_dict = read_rn_event_sequence(sequence)

            if only_residuals and event_dict['usrrnc'][1] == 0 and event_dict['mgrnc'][1] == 0:
                continue

            # determine whether or not we have to instantiate a new muon
            # if we do, instantiate one and add it to the list.
            if event_dict['primary'] != prim:
                prim = event_dict['primary']

                this_muon = Muon(seed, prim)
                
                muons.append(this_muon)

            this_event = RNevent(event_dict, this_muon)
            this_muon.add_event(this_event)
            events.append(this_event)


    return (muons, events)



            