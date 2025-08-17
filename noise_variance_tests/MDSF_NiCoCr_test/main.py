import os
import sys
import string
import subprocess
import numpy as np

# confg.py file
import config

# import matplotlib.pyplot as plt
sys.path.append("../../build/tools/pyMoDELib/")
import pyMoDELib

sys.path.append("../../python")
from pathlib import Path
from modlibUtils import *


def modify_parameter(target_config_fname: str, param: str, content) -> None:
    content_array = np.asarray(content)
    squeezed = np.squeeze(content_array)
    config_path = Path("inputFiles") / target_config_fname
    if squeezed.ndim > 1:  # matrix case
        setInputMatrix(config_path, param, squeezed)
    elif squeezed.ndim == 1:  # vector case
        setInputVector(config_path, param, squeezed, newCom="")
    else:  # scalar case
        setInputVariable(config_path, param, str(squeezed.item()))


def modify_dict_parameters(parameter_dict: dict, target_config_fname: str) -> None:
    """Batch process dictionary parameters"""
    for param, content in parameter_dict.items():
        try:
            modify_parameter(target_config_fname, param, content)
        except Exception as e:
            print(f"Error processing parameter '{param}': {str(e)}")
            continue


def set_initial_configuration() -> None:
    # Copy all necessary files from Library to inputFiles
    for key, src_path in config.files_to_copy_from_lib.items():
        dest = f"inputFiles/{src_path.name}"
        shutil.copy2(src_path.resolve(), dest)

    # set initial simulation configuration using the copied files
    ## DD.txt setup
    dd_fname = str(config.files_to_copy_from_lib["dd_file"]).split("/")[-1]
    modify_dict_parameters(config.DD_parameters, dd_fname)

    if config.enable_noise:
        ## noise_file.txt setup
        sf_noise_fname = str(config.files_to_copy_from_lib["noise_file_md_sf"]).split(
            "/"
        )[-1]
        modify_dict_parameters(config.MDStackingFault_parameters, sf_noise_fname)

    ## material file setup
    mat_fname = str(config.files_to_copy_from_lib["material_file"]).split("/")[-1]
    modify_dict_parameters(config.Material_parameters, mat_fname)

    ## microstructure file setup
    micro_fname = str(config.files_to_copy_from_lib["microstructure"]).split("/")[-1]
    modify_dict_parameters(config.Microstructure_parameters, micro_fname)
    with open("inputFiles/initialMicrostructure.txt", "w") as f:
        f.write(f"microstructureFile={micro_fname};\n")

    ## polycrystal file setup
    pf = PolyCrystalFile(mat_fname)
    mesh_fname = str(config.files_to_copy_from_lib["mesh"]).split("/")[-1]
    pf.meshFile = mesh_fname
    for param, value in config.Polycrystal_parameters.items():
        setattr(pf, param, value)
    pf.write("inputFiles")


def main() -> int:
    ################## Test Parameters #######################
    noise_seed_to_test = [1, 2, 3, 4, 5]
    glide_steps_to_test = np.array([
        [400, 405],
        [370, 405],
        [340, 405],
        [310, 405],
        [280, 405],
        [250, 405],
        [220, 405],
        [190, 405],
        [160, 405],
        [130, 405],
        [100, 405],
        [70, 405],
        [40, 405],
        [10, 405],
    ])

    # directory to store all the generated data
    dataStorageDir = Path("generatedData")
    os.makedirs(dataStorageDir, exist_ok=True)

    # remove old data
    # if os.path.exists(dataStorageDir):
    #    shutil.rmtree(dataStorageDir)

    for seed in noise_seed_to_test:
        for glide_step in glide_steps_to_test:
            # Preparing input files
            folders = ["evl", "F", "inputFiles"]
            for x in folders:
                # remove existing data
                if os.path.exists(x):
                    shutil.rmtree(x)
                # create necessary folder structure for the simulation
                os.makedirs(x)

            # set simulation parameters in inputFiles
            set_initial_configuration()

            # change glide step
            modify_parameter("periodicDipoleIndividual.txt", "glideSteps", np.array(glide_step))

            if config.enable_noise:
                # change the noise seed number
                modify_parameter("MDStackingFault.txt", "seed", seed)

            # boundary condition control file
            elasticDeformationFile = "ElasticDeformation.txt"

            # Voigt order is 11,22,33,12,23,13
            stressTensor = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            modify_parameter("ElasticDeformation.txt", "ExternalStress0", stressTensor)
            # generate microstructure
            subprocess.run(
                ["../../build/tools/MicrostructureGenerator/microstructureGenerator", "./"],
                capture_output=False,
                text=True,
            )
            exit()
            # run DDD
            subprocess.run(
                ["../../build/tools/DDomp/DDomp", "./"], capture_output=False, text=True
            )

            # store generated data
            seperation_dist = int(glide_step[-1]-glide_step[0])
            targetDir = Path(dataStorageDir) / f"seed{str(seed).strip()}" / f"sep{seperation_dist}b"
            os.makedirs(targetDir, exist_ok=True)
            for f in folders:
                shutil.move(f, targetDir)

    return 0


if __name__ == "__main__":
    main()
