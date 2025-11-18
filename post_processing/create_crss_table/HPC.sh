#!/bin/bash

#SBATCH --job-name crss_table
#SBATCH -C chip_manufacturer_intel,cpu_gen_haswell
#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 20gb
#SBATCH --time 22:00:00

# variables
SANDBOXDIR="/home/${USER}/apptainer_2/MoDELib2_latest_arch.sandbox"
VENV_PATH="/home/hyunsol/Environment/sandbox_env"

#source /home/hyunsol/Environment/sandbox_env/bin/activate

# change the working directory to the directory where the job is submitted
cd $SLURM_SUBMIT_DIR

# run DDD through a sandbox container
#apptainer exec $SANDBOXDIR python gen_crss_table.py

# run Python script inside container with bound virtual environment
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python gen_crss_table.py"
