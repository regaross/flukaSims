#!/usr/bin/python

__version__ = 2.2
__author__ = 'Regan Ross'
## Last Edited March 8, 2023

'''
Contact:
Regan Ross
rross@laurentian.ca
'''
#################################################
#                    IMPORTS                    }
#                                               }
#################################################

import numpy as np
from itertools import islice
from scipy import constants as sc

#################################################
#                   CONSTANTS                    }
#               (And Dictionaries)               }
#################################################

icode_dictionary = {
    None: None,
    -1: 'event not completed',
    0: 'normal event termination',
    4: 'stack overflow', 

    # Icode = 1x: call from Kaskad
    10: 'elastic interaction recoil', 
    11: 'inelastic interaction recoil',
    12: 'stopping particle',
    13: 'pseudo-neutron deposition',
    14: 'escape', 
    15: 'time kill',

    # Icode = 2x: call from Emfsco  
    20: 'local energy deposition (i.e. photoelectric)',
    21: 'below threshold, iarg=1',
    22: 'below threshold, iarg=2',
    23: 'escape',
    24: 'time kill',

    # Icode = 3x: call from Kasneu
    30: 'target recoil',
    31: 'below threshold',
    32: 'escape',
    33: 'time kill',

    # Icode = 4x: call from Kashea
    40: 'escape',
    41: 'time kill',
    42: 'delta ray stack overflow',

    # Icode = 5x: call from Kasoph
    50: 'optical photon absorption',
    51: 'escape',
    52: 'time kill',

    # Icode = 10x: call from Kaskad 
    100: 'elastic interaction secondaries',
    101: 'inelastic interaction secondaries',
    102: 'particle decay  secondaries',
    103: 'delta ray  generation secondaries',
    104: 'pair production secondaries',
    105: 'bremsstrahlung  secondaries',
    110: 'decay products',

    # Icode = 20x: call from Emfsco
    208: 'bremsstrahlung secondaries',
    210: 'Moller secondaries',
    212: 'Bhabha secondaries',
    214: 'in-flight annihilation secondaries',
    215: 'annihilation at rest secondaries',
    217: 'pair production secondaries',
    219: 'Compton scattering secondaries',
    221: 'photoelectric secondaries',
    225: 'Rayleigh scattering secondaries',

    # Icode = 30x: call from Kasneu
    300: 'interaction secondaries',

    # Icode = 40x: call from Kashea
    400: 'delta ray generation secondaries',
}

jtrack_dictionary = {
    None:   None,
    -6:     'alpha, 4-HELIUM',
    -5:     'helium-3, 3-HELIUM',
    -4:     'triton, TRITON',
    -3:     'deuteron, DEUTERON',
    -2:     'generic heavy ion with Z > 2, HEAVYION',
    -1:     'optical photon, OPTIPHOT',
    0:      'ray: pseudo straight-line particle, RAY',
    1:      'proton, PROTON',
    2:      'antiproton, APROTON',
    3:      'electron, ELECTRON',
    4:      'positron, POSITRON',
    5:      'electron neutrino, NEUTRIE',
    6:      'electron antineutrino, ANEUTRIE',
    7:      'photon, PHOTON',
    8:      'neutron, NEUTRON',
    9:      'antineutron, ANEUTRON',
    10:     'positive muon, MUON+',
    11:     'negative muon, MUON-',
    12:     'Kaon-zero long, KAONLONG',
    13:     'Positive Pion, PION+',
    14:     'Negative Pion, PION-',
    15:     'Positive Kaon, KAON+',
    16:     'Negative Kaon, KAON-',

    23:     'Pion-zero, PIZERO',

    208:    'Energy Deposition',

    211:    'EM Energy Deposition',

    308:    'Low Energy Neutron Kerma'

}

jtrack_rest_energies = { # Expressed in MeV if known, 0 otherwise (for heavy ions for instance)
    None:   None,
    -6:     sc.physical_constants['alpha particle mass energy equivalent in MeV'][0],
    -5:     sc.physical_constants['helion mass energy equivalent in MeV'][0],
    -4:     sc.physical_constants['triton mass energy equivalent in MeV'][0],
    -3:     sc.physical_constants['deuteron mass energy equivalent in MeV'][0],
    -2:     0,
    -1:     0,
    0:      0,
    1:      sc.physical_constants['proton mass energy equivalent in MeV'][0],
    2:      sc.physical_constants['proton mass energy equivalent in MeV'][0],
    3:      sc.physical_constants['electron mass energy equivalent in MeV'][0],
    4:      sc.physical_constants['electron mass energy equivalent in MeV'][0],
    5:      0,
    6:      0,
    7:      0,
    8:      sc.physical_constants['neutron mass energy equivalent in MeV'][0],
    9:      sc.physical_constants['neutron mass energy equivalent in MeV'][0],
    10:     sc.physical_constants['muon mass energy equivalent in MeV'][0],
    11:     sc.physical_constants['muon mass energy equivalent in MeV'][0],
    12:     497.611, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    13:     139.57039, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    14:     139.57039, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    15:     493.677, # https://pdg.lbl.gov/2022/listings/contents_listings.html
    16:     493.677, # https://pdg.lbl.gov/2022/listings/contents_listings.html

    23:     134.9768, # https://pdg.lbl.gov/2022/listings/contents_listings.html

    208:    0,

    211:    0,

    308:    0,

}

jtrack_labels = {
    None:   None,
    -6:     'Alpha',
    -5:     'Helium-3',
    -4:     'Triton',
    -3:     'Deuteron',
    -2:     'Generic heavy ion with Z > 2',
    -1:     'Optical Photon',
    0:      'Ray: pseudo straight-line particle',
    1:      'Proton',
    2:      'Antiproton',
    3:      'Electron',
    4:      'Positron',
    5:      'Electron Neutrino',
    6:      'Electron Antineutrino',
    7:      'Photon',
    8:      'Neutron',
    9:      'Antineutron',
    10:     'Positive Muon',
    11:     'Negative Muon',
    12:     'Kaon-zero Long',
    13:     'Positive Pion',
    14:     'Negative Pion',
    15:     'Positive Kaon',
    16:     'Negative Kaon',

    23:     'Pion-zero',

    208:    'Energy Deposition',

    211:    'EM Energy Deposition',

    308:    'Low Energy Neutron Kerma'

}

particle_colour_dictionary = {
    None:   None,
    -6:     'darkorange',   #Alpha
    -5:     'deeppink',     #Helium-3
    -4:     'maroon',       #Triton
    -3:     'darkviolet',   #Deuteron
    -2:     'blue',         #Heavy Ion
    -1:     'ivory',        #Optical P
    0:      'white',        #Ray pseudo straight-line particle
    1:      'red',          #Proton
    2:      'white',        #Antiproton
    3:      'darkcyan',     #Electron
    4:      'magenta',      #Positron
    5:      'lightgrey',    #Electron neutrino
    6:      'slategrey',    #Electron antineutrino
    7:      'gold',         #Photon
    8:      'green',        #Neutron
    9:      'white',        #Antineutron
    10:     'turquoise',    #Positive muon
    11:     'turquoise',    #Negative muon
    12:     'mistyrose',    #Kaon-zero long
    13:     'sienna',       #Positive Pion
    14:     'sienna',       #Negative Pion
    15:     'olive',        #Positive Kaon
    16:     'olive',        #Negative Kaon
    23:     'white',        #Pion Zero
    208:    'aliceblue',    #Heavy recoil
    211:    'linen',        #EM Energy deposition
    308:    'honeydew',     #Low energy neutron kerma
}

particle_colors_capitalized = {
    -6:     'DarkOrange',   #Alpha
    -5:     'DeepPink',     #Helium-3
    -4:     'Maroon',       #Triton
    -3:     'DarkViolet',   #Deuteron
    -2:     'Blue',         #Heavy Ion
    -1:     'Ivory',        #Optical P
    0:      'White',        #Ray pseudo straight-line particle
    1:      'Red',          #Proton
    2:      'White',        #Antiproton
    3:      'Cyan',         #Electron
    4:      'Magenta',      #Positron
    5:      'LightGrey',    #Electron neutrino
    6:      'SlateGrey',    #Electron antineutrino
    7:      'Yellow',       #Photon
    8:      'Green',        #Neutron
    9:      'White',        #Antineutron
    10:     'Turquoise',    #Positive muon
    11:     'Turquoise',    #Negative muon
    12:     'MistyRose',    #Kaon-zero long
    13:     'Sienna',       #Positive Pion
    14:     'Sienna',       #Negative Pion
    15:     'Olive',        #Positive Kaon
    16:     'Olive',        #Negative Kaon
    23:     'White',        #Pion Zero
    208:    'AliceBlue',    #Heavy recoil
    211:    'Linen',        #EM Energy deposition
    308:    'Honeydew',     #Low energy neutron kerma
}

#################################################
#                    CLASSES                     }
#                                                }
#################################################
       
class Event:
    ''' A general class representing particle events in FLUKA from mgdraw output. 
    An event has input and output parts— particles, fragments, and residuals. An event also has a location which implies a region.
    Events also have types— what describes the event that took place? This is given by an integer from the FLUKA documentation.
    Obviously regions have materials; therefore there is an event substrate too. '''

    # ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK
    attribute_list = ['icode', 'jtrack', 'ltrack', 'mreg', 'etrack', 'loc', 'dir','n_seconds', 'seconds', 'frags', 'residual']
    parent_attribute_list = ['p_icode', 'p_jtrack', 'p_ltrack', 'p_mreg', 'p_etrack', 'p_loc', 'p_dir']

    def __init__(self, icode = -100, jtrack = -100, ltrack = -100, mreg = -100, etrack = -100.0,\
                  loc = None, dir = None, seconds = None, frags = None, residual = None) -> None:

        self.icode = icode
        self.jtrack = jtrack
        self.ltrack = ltrack
        self.mreg = mreg
        self.etrack = etrack
        self.loc = loc
        self.dir = dir
        self.n_seconds = 0
        self.seconds = seconds          # A list containing the particle types that are created in the event in the sequence N, P1, P2, ..., PN
        self.frags = frags      # A list of fragments in the sequence N, A_1, Z_1, ... A_N, Z_N
        self.residual = residual        # A tuple of the form (A, Z) to represent the residual nucleus created in the event (if any)

        # SECONDARY LIST to point to the secondary "events"
        self.secondaries = []

        # A pointer to be set to the parent
        self.parent = None

        # Parent event Attributes
        self.p_icode = -100
        self.p_jtrack = -100
        self.p_ltrack = -100
        self.p_mreg = -100
        self.p_etrack = -100.0
        self.p_loc = None
        self.p_dir = None
            
    def add_secondary(self, new_secondary) -> bool:

        if self.check_secondary(new_secondary):
            jtrack = new_secondary.jtrack
            # Add the secondary to the list
            self.secondaries.append(new_secondary)
            # remove one of the elements from the list
            self.seconds.remove(jtrack)
            return True
        else:
            return False
        
    def matches_parent_criteria(self, possible_parent) -> bool:

        for p_att in self.parent_attribute_list:
            reg_att = p_att[2:]
            
            if not (getattr(self,p_att) == getattr(possible_parent, reg_att)):
                return False
        return True
    

#################################################
#                   FUNCTIONS                    }
#                                                }
#################################################

def event_from_dict(event_dict):
    ''' Instantiates an event from a dictionary'''
    new_event = Event()

    for prop in new_event.attribute_list:
        setattr(new_event, prop, event_dict[prop])

    for prop in new_event.parent_attribute_list:
        setattr(new_event, prop, event_dict[prop])

    return new_event

def read_event_file(event_filename):
    '''Reads in the events from a particular event file and returns a list of proto-event string lists'''
    event_list = []
    with open(event_filename) as ef:
        # There is a new event every 6 lines
        while True:
            event_entry = list(islice(ef, 5))
            if not event_entry:
                break
            else:
                primary = event_entry[0].split()
                secondaries = event_entry[1].split()
                n_seconds = secondaries[0]
                seconds = secondaries[1:]
                fragments = event_entry[2].split()
                residual = event_entry[3].split()
                parent = event_entry[4].split()

                icode , jtrack, mreg, ltrack = int(primary[0]), int(primary[1]), int(primary[2]), int(primary[3])
                etrack, loc = float(primary[4]), (float(primary[5]), float(primary[6]), float(primary[7]))
                primary_direction = (float(primary[8]), float(primary[9]), float(primary[10]))

                # Make a dictionary with the events and append it
                # ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK
                event_dict = {'icode':       icode,
                             'jtrack':      jtrack,
                             'mreg':        mreg,
                             'ltrack':      ltrack,
                             'etrack':      etrack,
                             'loc':         loc,
                             'dir':         primary_direction,
                             'n_seconds':   int(n_seconds),
                             'seconds':     [int(sec) for sec in seconds],
                             'secondaries': [],
                             'parent':      None,
                             'frags':       [int(frag) for frag in fragments],
                             'residual':    residual,
                             'p_icode':     int(parent[0]),
                             'p_jtrack':    int(parent[1]),
                             'p_mreg':      int(parent[2]),
                             'p_ltrack':    int(parent[3]),
                             'p_etrack':    float(parent[4]),
                             'p_loc':       (float(parent[5]), float(parent[6]), float(parent[7])),
                             'p_dir':       (float(parent[8]), float(parent[9]), float(parent[10]))
                             }
                event_list.append(event_dict)

    return event_list

def make_tree(event_dict_list):
    events = []

    for event_dict in event_dict_list:

        if event_dict['p_ltrack'] == 0:
            # We are at the FIRST event
            first = event_from_dict(event_dict)
            # Set its parent to None
            first.parent = None
            # append the event to the list of events
            events.append(first)

        else: # it is not the first event and therefore has a parent.
            # We need to find the parent of this event and set the pointers correctly
            next = event_from_dict(event_dict)
            
            # Find the parent event of 'next'
            for event in events:
                if next.matches_parent_criteria(event):
                    event.secondaries.append(next)
                    next.parent = event

            events.append(next)

    return events

def plot_tree_mpl(event_tree_list):
    '''A function to plot the event tree elements in 3D to show the trajectories
    of all particles'''
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as axes3d
    
    fig = plt.figure(figsize = (10,8))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    
    xs, ys, zs = [],[],[]
    
    for event in event_tree_list:
        px, py, pz = event.p_loc
        x, y, z = event.loc
        xs.append(x); ys.append(y); zs.append(z)
        
        if event.p_loc != (0,0,0):
            color = particle_colour_dictionary[event.jtrack]
            ax.plot([x,px], [y, py], [z, pz], color = color)
        
    
    #ax.scatter(xs, ys, zs, color = 'blue', s = 1)
    ax.view_init(elev=0, azim=0)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    plt.show()
    
def get_dataframe(events):
    '''Creates a dataframe with the events and parent attributes'''
    import pandas as pd
    
    # The data columns to be stored in the dataframe
    
    data_dict = {'icode':       [],
                 'jtrack':      [],
                 'mreg':        [],
                 'ltrack':      [],
                 'etrack':      [],
                 'loc':         [],
                 'dir':         [],
                 'n_seconds':   [],
                 'p_icode':     [],
                 'p_jtrack':    [],
                 'p_mreg':      [],
                 'p_ltrack':    [],
                 'p_etrack':    [],
                 'p_loc':       [],
                 'p_dir':       []
                 }
    
    xs, ys, zs, cxs, cys, czs = [],[],[],[],[],[]
    pxs, pys, pzs, pcxs, pcys, pczs = [],[],[],[],[],[]
    
    
    if type(events[0]) is Event:
        
        
        for event in events:
            x, y, z = event.loc
            cx, cy, cz = event.dir
            
            px, py, pz = event.p_loc
            pcx, pcy, pcz = event.p_dir
            
            xs.append(x); ys.append(y); zs.append(z)
            cxs.append(cx); cys.append(cy); czs.append(cz)
    
            pxs.append(px); pys.append(py); pzs.append(pz)
            pcxs.append(pcx); pcys.append(pcy); pczs.append(pcz)
            
            for key in data_dict.keys():
                data_dict[key].append(getattr(event, key))
             
                
    elif type(events[0]) is dict :
        
        for event in events:
            
            x, y, z = event['loc']
            cx, cy, cz = event['dir']
            
            px, py, pz = event['p_loc']
            pcx, pcy, pcz = event['p_dir']
            
            xs.append(x); ys.append(y); zs.append(z)
            cxs.append(cx); cys.append(cy); czs.append(cz)
    
            pxs.append(px); pys.append(py); pzs.append(pz)
            pcxs.append(pcx); pcys.append(pcy); pczs.append(pcz)
    
            for key in data_dict.keys():
                data_dict[key].append(event[key])
                
    else:
        return None
    
    data_dict['x'] = xs; data_dict['y'] = ys; data_dict['z'] = zs
    data_dict['px'] = pxs; data_dict['py'] = pys; data_dict['pz'] = pzs
            
    dframe = pd.DataFrame(data_dict)
    
        
    return dframe

def cylinder_surface_z(radius, height, loc = (0,0,0)):

    u = np.linspace(0, np.pi*2, 60)
    t = np.linspace(0, radius, 60)
    z = np.linspace(loc[2], loc[2] + height, 2)

    #Defining the parametric space for making the top and bottom circles
    U, T = np.meshgrid(u, t)

    #For top and bottom
    X1 = T*np.cos(U)
    Y1 = T*np.sin(U)
    ZTop = height*np.ones(X1.shape)
    ZBottom = loc[2]*np.ones(X1.shape)

    #Defining the parametric space for making the cylindrical column
    u, z = np.meshgrid(u, z)

    x = radius * np.cos(u) + loc[0]
    y = radius * np.sin(u) + loc[1]

    return x, y, z

def circular_surface_z(radius, loc = (0,0,0)):
    u = np.linspace(0, np.pi*2, 60)
    t = np.linspace(0, radius, 60)
    U, T = np.meshgrid(u, t)
    
    #For top and bottom
    x = T*np.cos(U)
    y = T*np.sin(U)
    z = loc[2]*np.ones(x.shape)

    return x, y, z

def spherical_surface(radius, loc = (0,0,0)):
    u = np.linspace(0, np.pi*2, 60)
    v = np.linspace(0, np.pi, 60)
    u, v = np.meshgrid(u, v)

    x = radius*np.ones(v.shape)*np.sin(v)*np.cos(u) + loc[0]
    y = radius*np.ones(v.shape)*np.sin(v)*np.sin(u) + loc[1]
    z = radius*np.ones(v.shape)*np.cos(v) + loc[2]

    return x, y, z

def change_energy_unit(energy_gev) -> tuple:
    '''Changes energy unit to a more reasonable energy and returns the energy with the unit.'''
    # If the energy is too large for GeV:
    if energy_gev > 1000:
        return (energy_gev/1e3, 'TeV')
    elif energy_gev < 1:
        # MeV range
        if energy_gev < 1/1000:
            # keV range
            if energy_gev < 1/1e6:
                # eV range
                return (energy_gev*1e9, 'eV')
            return (energy_gev*1e6, 'keV')
        return (energy_gev*1e3, 'MeV')
    else:
        return (energy_gev, 'GeV')

def get_kinetic_energy(jtrack, etrack)-> float:
    '''Returns the kinetic energy of the particle in question as opposed to Etrack which is the rest + kinetic'''
    rest_energy_MeV = jtrack_rest_energies[jtrack]
    rest_energy_GeV = rest_energy_MeV/1000
    kinetic_energy_GeV = etrack - rest_energy_GeV
    return kinetic_energy_GeV

def plot_tree_plotly(event_dataframe, show_geometry = True, show_EM = False):
    '''A function to plot the 3D event tree from the simulation'''
    # import plotly.express as px
    import plotly.graph_objects as go
    
    xs, ys, zs, colors = [],[],[],[]
    tags = []
    
    for index, row in event_dataframe.iterrows():
        # Load the lists to create the scatter plot
        is_edep = row['jtrack'] == 211 or  row['jtrack'] == 308
        is_electron = row['jtrack'] == 3 or row['jtrack'] == 4
        has_electron_parent = row['p_jtrack'] == 3 or row['p_jtrack'] == 4
        is_photon = row['jtrack'] == 7


        if not show_EM and ( is_electron or has_electron_parent or is_photon):
            pass

        elif is_edep: 
            # Don't show energy deposition events 
            # These are for scoring dose and are identical to EM tracks above threshold
            pass

        else:

            # Must have entry, entry, None so the particle line segments are separated
            xs.append(row['x']); xs.append(row['px']);xs.append(None)
            ys.append(row['y']); ys.append(row['py']);ys.append(None)
            zs.append(row['z']); zs.append(row['pz']);zs.append(None)
            color = particle_colour_dictionary[row['jtrack']]
            colors.append(color);colors.append(color);colors.append('white')
            energy = get_kinetic_energy(row['jtrack'], row['etrack']) #format(row['etrack'], '.2E')
            energy_tuple = change_energy_unit(energy)
            tags.append((jtrack_labels[row['jtrack']], icode_dictionary[row['icode']], energy_tuple[0], energy_tuple[1]))
            tags.append((jtrack_labels[row['jtrack']], icode_dictionary[row['icode']], energy_tuple[0], energy_tuple[1]))
            tags.append((None, None, None, None))


        # Maybe add region numbers for clarity after superimposing the geometry
        
    # Renaming so that the axis labels make more sense
    event_dataframe.rename(columns={'x': 'X [cm]', 'y': 'Y [cm]', 'z': 'Z [cm]'}, inplace = True)
    
    tree = go.Scatter3d(x=xs[3:], y=ys[3:], z=zs[3:],
            marker=dict(size=2, opacity = 0.5),
            line=dict(color=colors[3:], width=2), hovertemplate = '%{text}',
            #text = ['{}<br>{} GeV</br>{}'.format(tag[0], tag[2], tag[1]) for tag in tags[3:]], name='Event Tree')
            text = ['{}<br>{} {}</br>{}'.format(tag[0], tag[2], tag[3], tag[1]) for tag in tags[3:]], name='Event Tree')
    

    if show_geometry:

        #OD Cylinder section
        odx, ody, odz = cylinder_surface_z(617.22, 1280, (0,0,-640))
        OD = go.Surface(x = odx, y = ody, z = odz, opacity = 0.3, showscale = False, hoverinfo='skip', name='OD Cylinder', colorscale='blues')
        # OD Top and Bottom
        odx, ody, odz = circular_surface_z(617.22, (0,0,-640))
        OD_Bottom = go.Surface(x = odx, y = ody, z = odz, opacity = 0.3, showscale = False, hoverinfo='skip', name='OD Cylinder', colorscale='blues')
        odx, ody, odz = circular_surface_z(617.22,(0,0,640))
        OD_Top = go.Surface(x = odx, y = ody, z = odz, opacity = 0.3, showscale = False, hoverinfo='skip', name='OD Cylinder', colorscale='blues')

        # TPC Cylinder
        tpc_x, tpc_y, tpc_z = cylinder_surface_z(64, 172, (0, 0, -86))
        TPC = go.Surface(x = tpc_x, y = tpc_y, z = tpc_z, opacity = 0.3, showscale = False, hoverinfo='skip', name='TPC Cylinder', colorscale='oranges')

        # Outer Cryostat
        cry_x, cry_y, cry_z = spherical_surface(227, (0,0,40))
        outer_cryostat = go.Surface(x = cry_x, y = cry_y, z = cry_z, opacity = 0.3, showscale = False, hoverinfo='skip', name='Outer Cryostat', colorscale='greens')

        fig = go.Figure(data=[tree, OD, OD_Top, OD_Bottom, TPC, outer_cryostat])
        fig.update_layout(scene = dict(xaxis = dict(range=[-1000,1000],),
                     yaxis = dict(range=[-1000,1000],),
                     zaxis = dict(range=[-1000,1000])),
                     scene_aspectmode='cube')
    else:
        fig = go.Figure(data=[tree])

    fig.update_layout(paper_bgcolor="slategrey", template = 'plotly_dark', showlegend = True)

    fig.show()