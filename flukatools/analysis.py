#!/usr/bin/python

__version__ = 0.0
__author__ = 'Regan Ross'
## Last Edited June 7, 2024

'''
analysis.py

This part of the module is for performing analysis on data from the FLUKA simulations.
The data files have been created and organized, and are now ready to look at. This 
part of the module is certainly going to look the scrappiest when all is complete,
as there will be many files for manipulating the data into appropriate formats for 
plotting, and then also plotting functions.

'''
################################################################################
#                                                                              #
#                                   IMPORTS                                    #
#                                                                              #
################################################################################

# Import all the constants and dictionaries from the __init__ file
from . import *
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter

################################################################################
#                                                                              #
#                              DATA MANIPULATION                               #
#                                                                              #
################################################################################

def clean_resnuclei_data(resnuc_dict)-> None:
    '''Given a residual nuclei dictionary, this function sorts the data into a 3-column format: Z, A, and counts. It creates a new dictionary entry 'data' which will hold the newly sorted array. This can be used on residual nuclei dictionaries that have already been added together, OR individual ones.'''

    raw = resnuc_dict['raw']
    max_z = resnuc_dict['max_z']
    max_a = resnuc_dict['max_a']

    total = np.count_nonzero(raw)
    data = np.zeros((total, 3), dtype=float)

    new_raw = np.reshape(raw, (max_a, max_z))


    # Yes, the array parsing is FUNKY. Don't worry about it. This ought to be right
    # See this outdated (but still useful) man page: http://www.fluka.org/fluka.php?id=man_onl&sub=67
    count = 0
    k = -5
    for i in range(max_a):
        for j in range(max_z):
            Z = j + 1
            A = (i+1)+k+2*(j+1)
            if new_raw[i,j] > 0.0:
                data[count] = np.array([Z, A, new_raw[i,j]])
                count += 1

    resnuc_dict['table'] = data


def get_bootstrap_resnuc_data()-> np.ndarray:
    '''This function will be for producing an array of activation rates per year by sampling randomly
    and repeatedly from within a larger number of dictionaries.'''
    pass


################################################################################
#                                                                              #
#                              PLOTTING FUNCTIONS                              #
#                                                                              #
################################################################################

def plot_resnuclei_table_of_nuclides(nuclide_table : np.ndarray, title = 'Residual Nuclei', time = '') -> None:
    '''Anticipates an array whose rows are [Z, A, #]. Produces a plot of the table of nuclides
    as a sort of 2D histogram. One can hover the cursor over different boxes on the plot for 
    the particular information (Z, A, #)'''

    ## Turn the table into a 2D array 
    table = nuclide_table
    max_a = int(np.max(table[:,1]))
    max_z = int(np.max(table[:,0]))
    nuclides = np.zeros((max_z + 1, max_a + 1))

    for entry in table:
        x = int(entry[0])
        y = int(entry[1])
        val = entry[2]
        nuclides[x,y] = val

    nuclides[np.where(nuclides == 0)] = np.nan

    fig = px.imshow(nuclides, origin = 'lower')

    fig.update_layout(title = title, xaxis_title = 'Atomic Mass A', \
                    yaxis_title = 'Atomic Number Z', width = 1100, height = 600)

    fig.update_traces(hovertemplate='<b>Atomic Mass A:</b> %{x}<br>' + \
                    '<b>Atomic Number Z:</b> %{y}<br>' +
                    '<b>Count:</b> %{z}<br>' +
                    '<extra></extra>')
    
    fig.update_layout(coloraxis_colorbar = dict(title = 'Total Count ' + time, titleside = 'right', len = 0.6))
    fig.show()


def plot_neutron_parents_pie(neutrons : np.ndarray, title = 'Parents of Neutrons')-> None:
    '''Produces a plot of the immediate parent of neutrons as they are scored in the TPC.
    This simply tallies the parent JTRACK variables from the muon array.'''

    # Numbers representing parent particle species
    particle_species = neutrons[:,14]
    # Count occurrences of each particle species
    species_counts = Counter(particle_species)

    # Prepare labels and values for the pie chart
    labels = [JTRACK_LABELS[num] for num in species_counts.keys()]
    values = list(species_counts.values())

    # Create the pie chart
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, textinfo='none')])

    # Customize the layout
    fig.update_layout(
        title=title,
        title_x=0.5
    )

    # Show the plot
    fig.show()

def plot_neutron_icodes_pie(neutrons : np.ndarray, title = 'Neutron Birth Icodes')-> None:
    '''Produces a pie chart plot of the ICode values of the neutrons.
    In FLUKA lingo, this is basically their creation events.'''
    
    # Sample data: list of numbers representing particle species
    icodes = neutrons[:,0]
    # Count occurrences of each particle species
    icodes_counts = Counter(icodes)

    # Prepare labels and values for the pie chart
    labels = [ICODE_DICTIONARY[num] for num in icodes_counts.keys()]
    values = list(icodes_counts.values())

    # Create the pie chart
    fig = go.Figure(data=[go.Pie(labels=labels, values=values)])

    # Customize the layout
    fig.update_layout(
        title=title,
        title_x=0.5
    )
    # Show the plot
    fig.show()