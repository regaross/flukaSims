# Making a master file of FLUKA neutron events

# The goal of this step in the project is to write a script that takes simulation output and moves it to a larger file that will grow with each successive simulation.
# Ideally, this file will contain the following data:

# Neutron information:
# - Energy
# - Direction cosines
# - A point on its trajectory (instantiation point)
# - Generation number
# - Muon Ancestor (Index in Muon Dataset)
#     - Energy
#     - Direction cosines
#     - Starting point
#     - Impact Parameter

# WRITE(99,*) NCASE, ETRACK, LTRACK, WTRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK

import numpy as np
import h5py as h5
import muon_functions as mf
import os
from datetime import datetime
now = datetime.now()

neutron_filename = "NEWINPUT001_fort.99"
muon_filename = "muon_file.txt"

# Neutron attribute lists
muon, generation, energy, xsco, ysco, zsco, cosx, cosy, cosz = [],[],[],[],[],[],[],[],[]

with open(neutron_filename) as neutron_file:
    for line in neutron_file:
        temp = line.split()
        muon.append(int(temp[0]))
        energy.append(float(temp[1]))
        generation.append(int(temp[2]))
        xsco.append(float(temp[4]))
        ysco.append(float(temp[5]))
        zsco.append(float(temp[6]))
        cosx.append(float(temp[7]))
        cosy.append(float(temp[8]))
        cosz.append(float(temp[9]))


muenergy, impact,  initx, inity, initz, mucosx, mucosy, mucosz, pos_neg = [],[],[],[],[],[],[],[],[]
with open(muon_filename) as muon_file:
    lines = muon_file.readlines()
    for mu in muon:
        # the number in the list "muon" is 1+ the index in the file.
        muon_array = lines[mu-1].split()

        pos_neg.append(int(muon_array[0])) # Whether the muon is positive or negative
        muenergy.append(float(muon_array[1])) # Muon Energy
        initx.append(float(muon_array[2]))
        inity.append(float(muon_array[3]))
        initz.append(float(muon_array[4]))
        mucosx.append(float(muon_array[5]))
        mucosy.append(float(muon_array[6]))
        mucosz.append(float(muon_array[7]))

        cos_z = float(muon_array[7])
        cos_x = float(muon_array[5])
        zenith = np.arccos(cos_z)
        azimuth = np.arccos(cos_x/(np.sin(zenith)))
        temp_muon = mf.Muon(zenith, azimuth, initial=(float(muon_array[2]), float(muon_array[3]), float(muon_array[4])))
        impact.append(temp_muon.impact_param)

## WE NOW HAVE THE DATA
## Must put it in an hdf5 file.
# First, we check whether or not the file exists. If it does, we resize the datasets to accomodate the new data.

# If the file does not exist
if not os.path.isfile('./fluka_data.hdf5'): # The file must be created.
    file = h5.File('fluka_data.hdf5','a')
    # Create a group for the data
    data = file.create_group('data')
    # Create a group for the meta data (to dilineate different simulation data sets)
    meta = file.create_group('meta')

    # So too must the datasets be created and instantiated with zero size.
    data.create_dataset("muon_energy", (0,), dtype=float, maxshape=(None,))
    data.create_dataset("muon_impact", (0,), dtype=float, maxshape=(None,))
    data.create_dataset("muon_initial", (0,3), dtype=float, maxshape=(None,3))
    data.create_dataset("muon_direction", (0,3), dtype=float, maxshape=(None,3))
    data.create_dataset("muon_pn", (0,), dtype=int, maxshape=(None,))

    data.create_dataset("neutron_energy", (0,), dtype=float, maxshape=(None,))
    data.create_dataset("neutron_generation", (0,), dtype=int, maxshape=(None,))
    data.create_dataset("neutron_xyz", (0,3), dtype=float, maxshape=(None,3))
    data.create_dataset("neutron_direction", (0,3), dtype=float, maxshape=(None,3))

    meta.create_dataset("year",(0,), dtype=int, maxshape=(None,))
    meta.create_dataset("month",(0,), dtype=int, maxshape=(None,))
    meta.create_dataset("day",(0,), dtype=int, maxshape=(None,))
    meta.create_dataset("hour",(0,), dtype=int, maxshape=(None,))
    meta.create_dataset("minute",(0,), dtype=int, maxshape=(None,))
    meta.create_dataset("second",(0,), dtype=int, maxshape=(None,))
    meta.create_dataset("how_many",(0,), dtype=int, maxshape=(None,))


    file.close()

# Open the h5 file in "append" mode so that it won't erase current data
file = h5.File('fluka_data.hdf5','a')
# grabbing the file by the data group
data = file['data']
meta = file['meta']

# Current length of the data in the file
current_size = np.size(data['neutron_energy'])
n_points = len(muon) # Number of elements to append
current_meta_size = np.size(meta['year'])

for dset in data:
    # resize the datasets
    data[dset].resize(current_size + n_points, axis = 0)
for dset in meta:
    meta[dset].resize(current_meta_size + 1, axis = 0)

meta['hour'][current_meta_size] = int(now.strftime('%H'))
meta['minute'][current_meta_size] = int(now.strftime('%M'))
meta['second'][current_meta_size] = int(now.strftime('%S'))
meta['year'][current_meta_size] = int(now.strftime('%Y'))
meta['month'][current_meta_size] = int(now.strftime('%m'))
meta['day'][current_meta_size] = int(now.strftime('%d'))
meta['how_many'][current_meta_size] = n_points

# Now that the datasets have been resized, we must append the most recent sim data to the datasets.


for i in range(n_points):
    list_index = i
    dset_index = i + current_size

    data['muon_energy'][dset_index] = muenergy[list_index]
    data['muon_impact'][dset_index] = impact[list_index]
    data['muon_initial'][dset_index][:] = [initx[list_index], inity[list_index], initz[list_index]]
    data['muon_direction'][dset_index] = [mucosx[list_index], mucosy[list_index], mucosz[list_index]]
    data['muon_pn'][dset_index] = pos_neg[list_index]
    data['neutron_energy'][dset_index] = energy[list_index]
    data['neutron_generation'][dset_index] = generation[list_index]
    data['neutron_xyz'][dset_index][:] = [xsco[list_index], ysco[list_index], zsco[list_index]]
    data['neutron_direction'][dset_index][:] = [cosx[list_index], cosy[list_index], cosz[list_index]]

## MUST CLOSE THE FILE
file.close()