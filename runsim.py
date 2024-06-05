from flukatools import filemanip, muons, runtools
from time import sleep

# 0. Make sure necessary files are present (future feature)

# 1. Import the parameters from the YAML file.
filemanip.read_in_config_yaml('simconfig.yaml')

sleep(1)

# 2. Copy files to working directory
filemanip.copy_input_to_workdir()
# Make the output directory and copy the input to that directory
filemanip.make_output_directory()

sleep(1)

# 3. Make the muon phase space file. Here we retain the time passed and the complete list of the muons for further analysis later on.
muon_list, hours = muons.make_phase_space_file()

sleep(1)

# 4. Adjust the path to the phase space file in the source_routine
filemanip.change_muon_filepath()

sleep(1)

# 5. Change the number of muons in the input file to the number specified in the YAML configuration file
filemanip.change_number_of_muons()

sleep(1)

# 6. Link and Compile the .f files 
runtools.link_and_compile()

sleep(1)

# 7. Run the Simulation
runtools.run_fluka()

sleep(1)

# 8. Deal with Output...
filemanip.manage_output_files()
