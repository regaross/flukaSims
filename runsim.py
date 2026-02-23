from flukatools import filemanip, muons, runtools

# 0. Make sure necessary files are present (future feature)

# 1. Import the parameters from the YAML file.
filemanip.read_in_config_yaml('simconfig.yaml')

# 2. Copy files to working directory
filemanip.copy_input_to_workdir()

# Make the output directory and copy the input to that directory
filemanip.make_output_directory()

# 3. Make the muon phase space file.
muons.make_phase_space_file()

# Copy the muon file to the output directory
filemanip.copy_muon_file()

# 4. Adjust the path to the phase space file in the source_routine
filemanip.change_muon_filepath()

# 5. Change the number of muons in the input file to the number specified in the YAML configuration file
filemanip.change_number_of_muons()

# 6. Link and Compile the .f files 
runtools.link_and_compile()

# 7. Run the Simulation
runtools.run_fluka()

# # 8. Deal with Output...
filemanip.manage_output_files()
