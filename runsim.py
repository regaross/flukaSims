from flukatools import filemanip, muons, runtools

# Make sure necessary files are present

# Copy files to working directory
filemanip.copy_to_workdir()

# Make the muon phase space file
muon_list, hours = muons.make_phase_space_file()

# Modify the source file