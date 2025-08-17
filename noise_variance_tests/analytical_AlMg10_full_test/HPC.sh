#!/bin/bash

#SBATCH --job-name full500
#SBATCH --nodes 1
#SBATCH --cpus-per-task 48
#SBATCH --mem 30gb
#SBATCH --time 72:00:00

# variables
SANDBOXDIR="/home/${USER}/apptainer_2/MoDELib2_latest_arch.sandbox"

# change the working directory to the directory where the job is submitted
cd $SLURM_SUBMIT_DIR
# run DDD through a sandbox container
apptainer exec $SANDBOXDIR python main.py
