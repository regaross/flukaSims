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