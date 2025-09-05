import os
import sys
import vtk
import string
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy as scp

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
    an_noise_fname = str(config.files_to_copy_from_lib['noise_file_an_ss']).split('/')[-1]
    modify_dict_parameters(config.AnalyticalSolidSolution_parameters, an_noise_fname)

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


def readVTKnoise(fname: str):
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(fname)  # declare the vtk filename
    reader.ReadAllVectorsOn()  # read all vector data
    reader.ReadAllScalarsOn()  # read all scalar data
    reader.Update()  # update to new file
    # read data dimensions
    dims = [0, 0, 0]
    reader.GetOutput().GetDimensions(dims)
    NX, NY, NZ = dims
    # NX, NY, NZ = reader.GetOutput().GetDimensions()
    # creates dynamic point data object
    pointData = reader.GetOutput().GetPointData()

    # Initialize a list to store scalar arrays
    scalarArrays = []
    # Iterate over the number of arrays in pointData
    for i in range(pointData.GetNumberOfArrays()):
        array = pointData.GetArray(i)
        # Check if it's a scalar array
        if array.GetNumberOfComponents() == 1:
            scalarArrays.append(np.array(array))
    noiseScalars = scalarArrays[0]
    noiseScalars = np.reshape(noiseScalars, shape=(NX, NY))
    return noiseScalars, dims


def main() -> int:
    original_noise_file_xz = "noise_patch_0.vtk"
    original_noise_file_yz = "noise_patch_1.vtk"
    input_noise_field_xz, i_dim1 = readVTKnoise(original_noise_file_xz)
    input_noise_field_yz, i_dim2 = readVTKnoise(original_noise_file_yz)

    matFile = "inputFiles/AlMg50.txt"
    #absoluteTemp = 1
    #mat = pyMoDELib.PolycrystallineMaterialBase(matFile, absoluteTemp)
    b_SI = getValueInFile(f"{matFile}", "b_SI")
    mu0_SI = getValueInFile(f"{matFile}", "mu0_SI")
    rho_SI = getValueInFile(f"{matFile}", "rho_SI")
    gridSize = 256*256

    input_noise_field_xz *= mu0_SI
    input_noise_field_yz *= mu0_SI

    input_noise_field_xz = scp.signal.correlate(input_noise_field_xz, input_noise_field_xz, mode='same')/gridSize
    input_noise_field_yz = scp.signal.correlate(input_noise_field_yz, input_noise_field_yz, mode='same')/gridSize
    print(np.max(input_noise_field_xz))
    print(np.max(input_noise_field_yz))

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), dpi=200)
    # First plot with colorbar
    axes[0].set_title("")
    im0 = axes[0].imshow(input_noise_field_xz, cmap="coolwarm", origin="lower")
    fig.colorbar(im0, ax=axes[0], label="Noise Amplitude")
    axes[0].set_xlabel("X (grid units)")
    axes[0].set_ylabel("Y (grid units)")

    # Second plot with colorbar
    axes[1].set_title("")
    im1 = axes[1].imshow(input_noise_field_yz, cmap="coolwarm", origin="lower")
    fig.colorbar(im1, ax=axes[1], label="Noise Amplitude")
    axes[1].set_xlabel("X (grid units)")
    axes[1].set_ylabel("Y (grid units)")
    plt.tight_layout()
    plt.show()
    fig.savefig(f"corr.png", transparent=True)

    return 0

if __name__ == "__main__":
    main()
