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
image="$path_to_sim"fluka_nEXO.sif
job_num=$SLURM_ARRAY_JOB_ID         # An integer number representing this array job submission
task_num=$SLURM_ARRAY_TASK_ID       # An integer within the array range above representing which part of the parallelized run

#################################################
#         Printing meta data to stdout
################################################

date

echo "Executing on the machine:" $(hostname) 

echo "System random number used to seed the sim:" $rand

echo "Running Simulation!"

echo ; echo ; echo ; echo "COMPILING JARGON:"
echo

# Only create the output paths if it's the first run— every try would be redundant

if (( $SLURM_ARRAY_TASK_ID -lt 1 )); then

    mkdir run$SLURM_ARRAY_JOB_ID
    mkdir "run$SLURM_ARRAY_JOB_ID"/subrun$SLURM_ARRAY_TASK_ID
fi

while [ ! -f ./exe$stamp".exe" ]
do
    singularity exec -B /gpfs $image python $path_to_sim/make_executable.py -s $rand
done

singularity exec -B /gpfs $image python $path_to_sim/run_module.py -s $rand
