#!/bin/bash

#SBATCH --job-name separation
#SBATCH -C chip_manufacturer_intel,cpu_gen_haswell
#SBATCH --nodes 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 4gb
#SBATCH --time 05:00:00

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
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e1_s01_mg5.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e1_s01_mg10.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e1_s01_mg15.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e1_s23_mg5.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e1_s23_mg10.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e1_s23_mg15.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e2_s01_mg5.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e2_s01_mg10.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e2_s01_mg15.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e2_s23_mg5.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e2_s23_mg10.json"
apptainer exec --bind $VENV_PATH --pwd $SLURM_SUBMIT_DIR $SANDBOXDIR \
    /bin/bash -c "source $VENV_PATH/bin/activate && python generate_separation_data.py config_e2_s23_mg15.json"
