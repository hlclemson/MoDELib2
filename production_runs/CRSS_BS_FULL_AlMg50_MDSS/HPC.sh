#!/bin/bash

#SBATCH --job-name p615b_almg523_t-2-11
#SBATCH -C chip_manufacturer_intel,cpu_gen_haswell
#SBATCH --nodes 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 20gb
#SBATCH --time 72:00:00

# variables
SANDBOXDIR="/home/${USER}/apptainer_2/MoDELib2_latest_arch.sandbox"

# change the working directory to the directory where the job is submitted
cd $SLURM_SUBMIT_DIR
# run DDD through a sandbox container
apptainer exec $SANDBOXDIR python main.py initial_config.json
