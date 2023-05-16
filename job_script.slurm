#!/bin/bash
#
# SCRIPT: run_fluka.sh
# AUTHOR: Regan_Ross
# DATE: May 12, 2023
#
# PURPOSE: For requesting resources on SDF and running simultaneous simulations
#

#################################################
#                   SLURM
#################################################

#SBATCH --partition=shared
#SBATCH --job-name=fluka_sims               # a short name for your job
#SBATCH --output=simrun-%A-%a.out            # stdout file
#SBATCH --error=simrun-%A-%a.err             # stderr file
#SBATCH --nodes=1                           # node count
#SBATCH --ntasks=4                          # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                   # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=20G                   # memory per cpu-core (4G is default)
#SBATCH --time=04:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --array=0-4                         # job array with index values 0, 1, 2, 3, 4
##SBATCH --mail-type=all                    # send email on job start, end and fault
##SBATCH --mail-user=rross@slac.stanford.edu


#################################################
#                  Variables
#################################################

rand=$RANDOM
path_to_sim='/gpfs/slac/staas/fs1/g/exo/exo_data8/exo_data/users/rross/flukaSims/'
image=$path_to_sim"fluka_nEXO.sif"
job_num=$SLURM_ARRAY_JOB_ID         # An integer number representing this array job submission
task_num=$SLURM_ARRAY_TASK_ID       # An integer within the array range above representing which part of the parallelized run
fluka_path='/usr/local/fluka/bin/'

muon_source=$path_to_sim'muon_from_file.f'
mgdraw=$path_to_sim'mgdraw_neutron_count.f'
input=$path_to_sim'nEXO_OD.inp'

new_muon_source=$path_to_sim'musource'$rand'.f'
new_muon_source_o=$path_to_sim'musource'$rand'.o'
new_mgdraw=$path_to_sim'mgdrw'$rand'.f'
new_mgdraw_o=$path_to_sim'mgdrw'$rand'.o'
new_input=$path_to_sim'input'$rand'.inp'
exe_name=path_to_sim'exe'$rand'.exe'



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

#################################################
#          Make Copies of the Files
#################################################

# echo "Making Copies of the Files"

# cp $muon_source $new_muon_source
# cp $mgdraw $new_mgdraw
# cp $input $new_input

#################################################
#   Changing the Muon filename in the source file
#################################################

# singularity shell -B /gpfs $image python $path_to_sim"change_muon_filepath.py" -s $rand

#################################################
#   Compile and link the simulation files
#################################################

# $image_exec $fluka_path"fff" $new_mgdraw
# $image_exec $fluka_path"fff" $new_muon_source
# $image_exec $fluka_path"ldpmqmd -m fluka -o " $exe_name $new_mgdraw_o $new_muon_source_o


#Only create the output paths if it's the first run— every try would be redundant

if (( "$SLURM_ARRAY_TASK_ID" < 1 )); then

    mkdir run$SLURM_ARRAY_JOB_ID
    mkdir "run$SLURM_ARRAY_JOB_ID"/subrun$SLURM_ARRAY_TASK_ID
fi

# executable=./exe$stamp".exe"

# while [ ! -f $executable ]
# do
#     echo "Trying to compile"
#     singularity exec -B /gpfs $image python $path_to_sim/make_executable.py -s $rand
# done

singularity exec -B /gpfs $image python $path_to_sim/run_module.py -s $rand
