from flukatools import filemanip, muons, runtools

# 0. Make sure necessary files are present (future feature)

# 1. Import the parameters from the YAML file.
filemanip.read_in_config_yaml('simconfig.yaml')

# 2. Copy files to working directory
filemanip.copy_input_to_workdir()

# 3. Make the muon phase space file. Here we retain the time passed and the complete list of the muons for further analysis later on.
muon_list, hours = muons.make_phase_space_file()

# 4. Adjust the path to the phase space file in the source_routine
filemanip.change_muon_filepath()

# 5. Change the number of muons in the input file to the number specified in the YAML configuration file
filemanip.change_number_of_muons()

# 6. Link and Compile the .f files 
runtools.link_and_compile()

# 7. Run the Simulation
runtools.run_fluka()

# 8. Deal with Output...
