#!/bin/bash
#
# SCRIPT: job_script.slurm
# AUTHOR: Regan_Ross
# DATE: May 7, 2024
#
# PURPOSE: For requesting resources on SDF and running simultaneous simulations
#

#################################################
#                   SLURM
#################################################

#SBATCH --partition=milano
#SBATCH --job-name=fluka_sims                   # a short name for your job
#SBATCH --output=simrun-%A-%a.out               # stdout file
#SBATCH --error=simrun-%A-%a.err                # stderr file
#SBATCH --nodes=1                               # node count
#SBATCH --ntasks=4                              # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=10G                       # memory per cpu-core (4G is default)
#SBATCH --time=03:30:00                         # total run time limit (HH:MM:SS)
#SBATCH --array=0-4                             # job array with index values 0, 1, 2, 3, 4
#SBATCH --mail-type=all                         # send email on job start, end and fault
#SBATCH --mail-user=regan.ross@mail.mcgill.ca   # This feature does work!


#################################################
#                  Variables
#################################################

output_dir="run"$SLURM_ARRAY_JOB_ID

rand=$RANDOM
path_to_sim='/sdf/home/r/rross/nEXO/flukaSims/'
image=$path_to_sim"FLUKA.sif"
job_num=$SLURM_ARRAY_JOB_ID         # An integer number representing this array job submission
task_num=$SLURM_ARRAY_TASK_ID       # An integer within the array range above representing which part of the parallelized run
fluka_path='/usr/local/fluka/bin/'

muon_source=$path_to_sim'muon_from_file.f'
mgdraw=$path_to_sim'mgdraw_neutron_count.f'
input=$path_to_sim'nEXO_2024.inp'

new_muon_source=$path_to_sim'musource'$rand'.f'
new_muon_source_o=$path_to_sim'musource'$rand'.o'
new_mgdraw=$path_to_sim'mgdrw'$rand'.f'
new_mgdraw_o=$path_to_sim'mgdrw'$rand'.o'
new_input=$path_to_sim'input'$rand'.inp'
exe_name=path_to_sim'exe'$rand'.exe'


#################################################
#    REMEMBER: These do not run sequentially.
################################################


mkdir $output_dir
mkdir $output_dir/subrun$SLURM_ARRAY_TASK_ID

#################################################
#         Printing meta data to stdout
################################################

cd $path_to_sim
date
echo "Executing on the machine:" $(hostname) 
echo "System random number used to seed the sim:" $rand
echo "Running Simulation!"
echo ; echo ; echo ; echo "COMPILING JARGON:"
echo




singularity exec -B /sdf $image python $path_to_sim/run_simulation.py -s $rand


mv *$rand* $output_dir/subrun$SLURM_ARRAY_TASK_ID
mv *.hdf5 $output_dir/