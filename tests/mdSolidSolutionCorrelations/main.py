import os
import sys
import string
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# confg.py file
import config

#import matplotlib.pyplot as plt
sys.path.append("../../build/tools/pyMoDELib/")
import pyMoDELib

sys.path.append("../../python")
from pathlib import Path
from modlibUtils import *


def modify_parameter(target_config_fname: str, param: str, content) -> None:
    content_array = np.asarray(content)
    squeezed = np.squeeze(content_array)
    config_path = Path('inputFiles')/target_config_fname
    if squeezed.ndim > 1: # matrix case
        setInputMatrix(config_path, param, squeezed)
    elif squeezed.ndim == 1: # vector case
        setInputVector(config_path, param, squeezed, newCom='')
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
    dd_fname = str(config.files_to_copy_from_lib['dd_file']).split('/')[-1]
    modify_dict_parameters(config.DD_parameters, dd_fname)

    ## noise_file.txt setup
    ss_noise_fname = str(config.files_to_copy_from_lib['noise_file_md_ss']).split('/')[-1]
    modify_dict_parameters(config.MDSolidSolution_parameters, ss_noise_fname)

    ## material file setup
    mat_fname = str(config.files_to_copy_from_lib['material_file']).split('/')[-1]
    modify_dict_parameters(config.Material_parameters, mat_fname)

    ## microstructure file setup
    micro_fname = str(config.files_to_copy_from_lib['microstructure']).split('/')[-1]
    modify_dict_parameters(config.Microstructure_parameters, micro_fname)
    with open('inputFiles/initialMicrostructure.txt', 'w') as f:
        f.write(f"microstructureFile={micro_fname};\n")

    ## polycrystal file setup
    pf = PolyCrystalFile(mat_fname)
    mesh_fname = str(config.files_to_copy_from_lib['mesh']).split('/')[-1]
    pf.meshFile=mesh_fname
    for param, value in config.Polycrystal_parameters.items():
        setattr(pf, param, value)
    pf.write('inputFiles')


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
    folders=['evl','F', 'inputFiles']
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

    matFile = "inputFiles/FeCrAl_Fe.txt"
    absoluteTemp = 1
    mat = pyMoDELib.PolycrystallineMaterialBase(matFile, absoluteTemp)
    b_SI = getValueInFile(f"{matFile}", "b_SI")
    mu0_SI = getValueInFile(f"{matFile}", "mu0_SI")
    rho_SI = getValueInFile(f"{matFile}", "rho_SI")

    tag = "1"
    originalCorrFile_xz = "./inputFiles/MoDELCompatible_FeCr8_xz.vtk"
    originalCorrFile_yz = "./inputFiles/MoDELCompatible_FeCr8_yz.vtk"
    seed = 1
    gridSize = np.array([100, 100, 1])
    gridSpacing = np.array([1e-10, 1e-10, 1e-10]) / b_SI
    latticeBasis = np.array([[1, 0], [0, 1]]).reshape(2, 2)
    a_cai_SI = 10e-10 / b_SI
    #a_cai_SI = 0

    ssNoise = pyMoDELib.MDSolidSolutionNoise(
        mat,
        tag,
        originalCorrFile_xz,
        originalCorrFile_yz,
        seed,
        gridSize.reshape(3, 1),
        gridSpacing.reshape(3, 1),
        latticeBasis,
        a_cai_SI,
    )
    realizationNum = 100
    corr_xz, corr_yz = np.array(ssNoise.averageNoiseCorrelation(realizationNum)) * mu0_SI**2

    corr_xz = circularShift(corr_xz, gridSize)
    corr_yz = circularShift(corr_yz, gridSize)
    corr_xz = corr_xz.reshape(100, 100)
    corr_yz = corr_yz.reshape(100, 100)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), dpi=200)
    # First plot with colorbar
    axes[0].set_title("corr_xz")
    im0 = axes[0].imshow(corr_xz, cmap="coolwarm", origin="lower")
    fig.colorbar(im0, ax=axes[0], label="Pa2")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")

    # Second plot with colorbar
    axes[1].set_title("corr_yz")
    im1 = axes[1].imshow(corr_yz, cmap="coolwarm", origin="lower")
    fig.colorbar(im1, ax=axes[1], label="Pa2")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    plt.tight_layout()
    fig.savefig(f"ensemCorrSS_R{realizationNum}.png", transparent=True)

    return 0

if __name__ == "__main__":
    main()
