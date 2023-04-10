#!/bin/bash

    #SBATCH --partition=shared
    #
    #SBATCH --job-name=fluka_simulation
    #SBATCH --output=output-%j.txt
    #SBATCH --error=output-%j.txt
    #
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=10g
    #
    #SBATCH --time=0-10:00:00
    #
    #SBATCH --cpus 1

    cd /gpfs/slac/staas/fs1/g/exo/exo_data8/exo_data/users/rross/flukaSims

    singularity exec -B /gpfs fluka_nEXO.sif python runsim.py