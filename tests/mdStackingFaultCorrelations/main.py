import os
import sys
import vtk
import string
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import norm
from matplotlib import cm
from matplotlib.colors import Normalize

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
    #ss_noise_fname = str(config.files_to_copy_from_lib['noise_file_md_ss']).split('/')[-1]
    sf_noise_fname = str(config.files_to_copy_from_lib['noise_file_md_sf']).split('/')[-1]
    #modify_dict_parameters(config.MDSolidSolution_parameters, ss_noise_fname)
    modify_dict_parameters(config.MDStackingFault_parameters, sf_noise_fname)

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


def initUnitVectors(originalCorrVTK: str):
    latVector1 = None
    latVector2 = None
    copyFlag1 = False
    copyFlag2 = False
    with open(originalCorrVTK) as o:
        for line in o:
            if copyFlag1:
                latVector1 = np.array(line.strip().split(), dtype=float)
                copyFlag1 = False
            elif copyFlag2:
                latVector2 = np.array(line.strip().split(), dtype=float)
                copyFlag2 = False
            elif "VECTORS lattice_basis1 double" in line:
                copyFlag1 = True
            elif "VECTORS lattice_basis2 double" in line:
                copyFlag2 = True

    unitVec1 = latVector1 / np.linalg.norm(latVector1)
    unitVec2 = latVector2 / np.linalg.norm(latVector2)
    unitVec1 = unitVec1[:2]  # drop the z axis value
    unitVec2 = unitVec2[:2]
    return unitVec1, unitVec2


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


def constructDotGrid(NX, NY, unitVec1, unitVec2):
    # Construct non-orthogonal dot grid based on the non-minimized atom structure
    firstNearNeighborDist = 2.8
    dotGrid = np.zeros((NX * NY, 2))  # each row contains x and y position of the dot
    vec1 = unitVec1 * firstNearNeighborDist
    vec2 = unitVec2 * firstNearNeighborDist
    for yIdx in range(NY):
        for xIdx in range(NX):
            dotGrid[yIdx * NX + xIdx] = vec1 * xIdx + vec2 * yIdx
    return dotGrid


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

    # setInputVariable(MDStackingFault_parameters['copy_to'], 'testNoiseSampling', str(1))
    simulationDir = os.path.abspath(".")
    ddBase = pyMoDELib.DislocationDynamicsBase(simulationDir)

    matFile = "inputFiles/AlMg5.txt"
    absoluteTemp = 1
    mat = pyMoDELib.PolycrystallineMaterialBase(matFile, absoluteTemp)
    b_SI = getValueInFile(f"{matFile}", "b_SI")
    mu0_SI = getValueInFile(f"{matFile}", "mu0_SI")
    rho_SI = getValueInFile(f"{matFile}", "rho_SI")

    tag = "0"
    originalCorrFile = "inputFiles/AlMg5_Cx_R100_ISF.vtk"
    seed = 1
    gridSize = np.array([100, 100, 1])
    gridSpacing = np.array([1e-10, 1e-10, 1e-10]) / b_SI
    latticeBasis = np.array([[1, 1], [1, 1]]).reshape(2, 2)
    sfNoise = pyMoDELib.MDStackingFaultNoise(
        mat,
        tag,
        originalCorrFile,
        seed,
        gridSize.reshape(3, 1),
        gridSpacing.reshape(3, 1),
        latticeBasis,
    )

    realizationNum = 100
    correlation = sfNoise.averageNoiseCorrelation(realizationNum)
    correlation = np.squeeze(np.array(correlation))
    # circular shift
    correlation = circularShift(correlation, gridSize)
    correlation = np.array(correlation).reshape(gridSize[0], gridSize[1])
    # convert unit from DDD to SI
    correlation *= (mu0_SI * b_SI) ** 2

    unitVec1, unitVec2 = initUnitVectors(originalCorrFile)

    inputPath = "./inputFiles/"
    # reshape the correlation data to 1D array
    NX_s, NY_s, NZ_s = gridSize

    # create callable color map
    color_map = mpl.colormaps["coolwarm"]

    # normalize the correlation coefficient array for the color assignment
    dotGrid_s = constructDotGrid(NX_s, NY_s, unitVec1, unitVec2)
    correlation = correlation.reshape(NY_s * NX_s, 1)
    c_normalized = (correlation - np.min(correlation.flatten())) / (
        np.max(correlation.flatten()) - np.min(correlation.flatten())
    )
    color = color_map(c_normalized)

    # load input correlation values from the correlation vtk files
    inputCorrelationField, i_dim = readVTKnoise(originalCorrFile)
    NX_i, NY_i, NZ_i = i_dim
    dotGrid_i = constructDotGrid(NX_i, NY_i, unitVec1, unitVec2)
    inputCorrelationField = inputCorrelationField.reshape(NY_i * NX_i, 1)
    c_normalized_i = (
        inputCorrelationField - np.min(inputCorrelationField.flatten())
    ) / (
        np.max(inputCorrelationField.flatten())
        - np.min(inputCorrelationField.flatten())
    )
    color_i = color_map(c_normalized_i)

    # plot correlation data on the dot grid
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 8.6), dpi=500)

    axes[0].set_title("Input Correlation")
    axes[0].scatter(dotGrid_i[:, 0], dotGrid_i[:, 1], c=color_i, s=30)
    axes[1].set_title(f"Sampled R")
    axes[1].scatter(dotGrid_s[:, 0], dotGrid_s[:, 1], c=color, s=10)
    fig.colorbar(
        plt.cm.ScalarMappable(
            norm=Normalize(
                np.min(inputCorrelationField), np.max(inputCorrelationField)
            ),
            cmap=color_map,
        ),
        ax=axes[0], label="J/m2"
    )
    fig.colorbar(
        plt.cm.ScalarMappable(
            norm=Normalize(np.min(correlation), np.max(correlation)), cmap=color_map
        ),
        ax=axes[1], label="J/m2"
    )
    plt.tight_layout()
    fig.savefig(f"ensemCorr_R{realizationNum}.png", transparent=True)
    plt.close()

    return 0

if __name__ == "__main__":
    main()
