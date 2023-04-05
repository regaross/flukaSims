### Confirmed working on Apr 4, 2023


import argparse
import os
import re
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument(dest='resnuclei_file', type=str)
args = parser.parse_args()

resnuclei = args.resnuclei_file

with open(resnuclei) as file:
    header = file.read(500)

# Create a temporary file to hold ONLY the array data
os.system('sed -e \'/^[[:space:]]*[0-9].*[0-9]$/!d\' ' + resnuclei + ' > temp.tsv')
# Read in the temporary data file
data = np.loadtxt('temp.tsv')
# Remove the temporary array-only file
os.system('rm temp.tsv')

# Collect the data that determine the iteration pattern.
max_z = int(re.search('(?<=Z:).*(?=,)', header).group(0).strip())
max_nz = int(re.search('(?<=-Z:).*(?=M)', header).group(0).strip())
min_nz = int(re.search('(?<=Min. N-Z:).*(?=\n)', header).group(0).strip())
k = -5
max_a = max_nz - k

# print_line = 'Max Z = {}, Max N-Z = {}, Min N-Z = {}'.format(max_z, max_nz, min_nz)
# print(print_line)

data = np.reshape(data, (max_a, max_z))

print('Z\tA\tResult')

for i in range(max_a):
    for j in range(max_z):
        Z = j + 1
        A = (i+1)+k+2*(j+1)
        if data[i,j] > 0.0:
            print(str(Z) + '\t' + str(A) + '\t' + str(data[i,j]))