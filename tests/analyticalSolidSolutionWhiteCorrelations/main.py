import os
import sys
import string
import subprocess
import numpy as np
import matplotlib.pyplot as plt

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

    ## noise_file.txt setup
    an_noise_fname = str(config.files_to_copy_from_lib["noise_file_an_ssw"]).split("/")[
        -1
    ]
    modify_dict_parameters(
        config.AnalyticalSolidSolutionWhite_parameters, an_noise_fname
    )

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


def circularShift(corrArray: np.ndarray, gridSize: np.ndarray) -> np.ndarray:
    NX, NY = np.squeeze(gridSize[0]), np.squeeze(gridSize[1])
    shiftY = NY / 2
    shiftX = NX / 2
    tempArr = np.zeros(NX * NY)
    for y in np.arange(NY):
        newY = int((y + shiftY) % NY)
        for x in np.arange(NX):
            newX = int((x + shiftX) % NX)
            tempArr[newY * NX + newX] += corrArray[y * NX + x]
    return tempArr


def main() -> int:
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

    simulationDir = os.path.abspath(".")
    ddBase = pyMoDELib.DislocationDynamicsBase(simulationDir)
    matFile = f"inputFiles/{str(config.files_to_copy_from_lib['material_file']).split('/')[-1]}"

    absoluteTemp = 1
    mat = pyMoDELib.PolycrystallineMaterialBase(matFile, absoluteTemp)
    b_SI = getValueInFile(f"{matFile}", "b_SI")
    mu0_SI = getValueInFile(f"{matFile}", "mu0_SI")
    rho_SI = getValueInFile(f"{matFile}", "rho_SI")

    tag = "1"
    seed = 0
    gridSize = np.array([256, 256, 64])
    gridSpacing = np.array([1e-10, 1e-10, 1e-10]) / b_SI
    latticeBasis = np.array([[1, 0], [0, 1]]).reshape(2, 2)

    msss_SI = 0.9e18 / mu0_SI**2
    # [Pa^2] Mean Square Shear Stress
    a = 1e-10 / b_SI
    a_cai_SI = 1e-09 / b_SI  # in meter
    anssNoise = pyMoDELib.AnalyticalSolidSolutionWhiteNoise(
        tag,
        seed,
        gridSize.reshape(3, 1),
        gridSpacing.reshape(3, 1),
        latticeBasis,
        a,
        a_cai_SI,
        msss_SI,
    )

    realizationNum = 1
    corr_xz, corr_yz = (
        np.array(anssNoise.averageNoiseCorrelation(realizationNum)) * mu0_SI**2
    )

    corr_xz = circularShift(corr_xz, gridSize)
    corr_yz = circularShift(corr_yz, gridSize)
    corr_xz = corr_xz.reshape(gridSize[0], gridSize[1])
    corr_yz = corr_yz.reshape(gridSize[0], gridSize[1])

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), dpi=200)
    # First plot with colorbar
    axes[0].set_title("corr_xz")
    im0 = axes[0].imshow(corr_xz, cmap="coolwarm", origin="lower")
    fig.colorbar(im0, ax=axes[0], label="Noise Amplitude")
    axes[0].set_xlabel("X (grid units)")
    axes[0].set_ylabel("Y (grid units)")

    # Second plot with colorbar
    axes[1].set_title("corr_yz")
    im1 = axes[1].imshow(corr_yz, cmap="coolwarm", origin="lower")
    fig.colorbar(im1, ax=axes[1], label="Noise Amplitude")
    axes[1].set_xlabel("X (grid units)")
    axes[1].set_ylabel("Y (grid units)")
    plt.tight_layout()
    fig.savefig(f"ensemCorr_R{realizationNum}.png", transparent=True)

    return 0


if __name__ == "__main__":
    main()
