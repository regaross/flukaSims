#!/bin/bash
#
# SCRIPT: run_fluka.sh
# AUTHOR: Regan_Ross
# DATE: May 12, 2023
#
# PURPOSE: For requesting resources on SDF and running simultaneous simulations
#
#
#
#################################################
#                   SLURM
#################################################

#SBATCH --partition=shared
#SBATCH --job-name=fluka_sims               # a short name for your job
#SBATCH --output=slurm-%A.%a.out            # stdout file
#SBATCH --error=slurm-%A.%a.err             # stderr file
#SBATCH --nodes=1                           # node count
#SBATCH --ntasks=1                          # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                   # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=20G                    # memory per cpu-core (4G is default)
#SBATCH --time=04:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --array=0-4                         # job array with index values 0, 1, 2, 3, 4
##SBATCH --mail-type=all                     # send email on job start, end and fault
##SBATCH --mail-user=rross@slac.stanford.edu


#################################################
#                  Variables
#################################################

rand=$RANDOM
fluka_bin='/usr/local/fluka/bin/'
fff="$fluka_bin"fff
link="$fluka_bin""ldpmqmd -m fluka -o "
path_to_sim='/gpfs/slac/staas/fs1/g/exo/exo_data8/exo_data/users/rross/flukaSims/'
image="$path_to_sim"fluka_nEXO.sif
mgdraw='mgdraw_neutron_count'
source='muon_from_file'
exe="nEXO_OD"$rand".exe"

#################################################
#         Compile & Link FLUKA Routines
#################################################

singularity exec -B /gpfs $image $fff "$mgdraw"".f"
singularity exec -B /gpfs $image $fff "$source"".f"
singularity exec -B /gpfs $image $link $exe "$mgdraw"".o" "$source"".o"

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

singularity exec -B /gpfs $path_to_sim/fluka_nEXO.sif python $path_to_sim/run_module.py 