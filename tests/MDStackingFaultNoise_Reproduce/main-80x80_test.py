import os
import sys
import vtk
import glob
import string
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib as mpl
import matplotlib.transforms as mtransforms
sys.path.append("../../build/tools/pyMoDELib/")
import pyMoDELib
from scipy.signal import correlate2d
from scipy.stats import norm
from matplotlib import cm
from matplotlib.colors import Normalize
sys.path.append("../../python")
from modlibUtils import *
from pathlib import Path

# Define all file templates (using pathlib for cross-platform paths)
file_templates = {
    'dd_file': Path('../../Library/DislocationDynamics/DD.txt'),
    'noise_file': Path('../../Library/GlidePlaneNoise/MDStackingFault.txt'),
    'material_file': Path('../../Library/Materials/AlMg5.txt'),
    'elastic_deformation_file': Path('../../Library/ElasticDeformation/ElasticDeformation.txt'),
    'mesh': Path('../../Library/Meshes/unitCube24.msh'),
    'microstructure': Path('../../Library/Microstructures/periodicDipoleIndividual.txt')
}

MDStackingFault_parameters = {
    'variables': {
        # Scalar parameters (use setInputVariable)
        'type': 'MDStackingFaultNoise',
        'tag': 1,
        'seed': 2,
        'outputNoise': 1,
        'noiseFile': 'AlMg5_ISF.vtk',
        # File paths (convert to absolute paths)
        'correlationFile': Path('../../Library/GlidePlaneNoise/AlMg5_Cx_R100_ISF.vtk').resolve(),
    },
    'vectors': {
        # Vector parameters with comments (use setInputVector)
        'gridSize': {
            'value': np.array([80, 80, 1]),
            'comment': 'number of grid points in each direction'
        },
        'gridSpacing_SI': {
            'value': np.array([1.0e-10, 1.0e-10, 1e-10]),
            'comment': 'grid spacing in each direction'
        }
    },
    'copy_to': 'inputFiles/MDStackingFault.txt'
}

Material_parameters = {
    'variables': {
        'enabledSlipSystems': 'Shockley',
        #'glidePlaneNoise': ['MDSolidSolution.txt', 'MDStackingFault.txt'],
        'glidePlaneNoise': 'MDStackingFault.txt',
        'atomsPerUnitCell': '1',
        'dislocationMobilityType': 'default'
    },
    'copy_to': 'inputFiles/AlMg5.txt'
}

Polycrystal_parameters= {
    'parameters': {
        'absoluteTemperature': 1,
        'grain1globalX1': np.array([0, 1, 1]),
        'grain1globalX3': np.array([-1, 1, -1]),
        'boxEdges': np.array([[0,1,1], [-2,-1,1], [-1,1,-1]]),
        'boxScaling': np.array([200, 200, 200]),
        #'X0': np.array([0.5, 0.5, 0.5]),
        'X0': np.array([0., 0., 0.]),
        'periodicFaceIDs': np.array([0,1,2,3,4,5]),
        'gridSpacing_SI': np.array([1.0e-10, 1.0e-10])
    }
}

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
    #NX, NY, NZ = reader.GetOutput().GetDimensions()
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

def processInitialConfigurations() -> None:
    # Copy all template files
    print("\033[1;32mCopying template files...\033[0m")
    for key, src_path in file_templates.items():
        dest = f"inputFiles/{src_path.name}"
        shutil.copy2(src_path.resolve(), dest)
        #shutil.copy2(src_path, dest)
        print(f"Created {dest}")

    #print("\033[1;32mCreating ddFile\033[0m")
    ## Apply all DD.txt parameters from the dictionary
    #for param, value in DD_parameters['variables'].items():
    #    setInputVariable(DD_parameters['copy_to'], param, str(value))

    print("\033[1;32mCreating  noiseFile\033[0m")
    # copy noise file into inputFiles
    shutil.copy2(MDStackingFault_parameters['variables']['correlationFile'], 'inputFiles/')
    # copy declared variables to the configuration txt file
    for param, value in MDStackingFault_parameters['variables'].items():
        if 'correlationFile' in param:
            relativePath = f"{str(value).split('/')[-1]}"
            setInputVariable(MDStackingFault_parameters['copy_to'], param, relativePath)
        else:
            setInputVariable(MDStackingFault_parameters['copy_to'], param, str(value))
    for param, data in MDStackingFault_parameters['vectors'].items():
        setInputVector(
            MDStackingFault_parameters['copy_to'],
            param,
            data['value'],
            data['comment']
        )

    # Process material file
    print("\033[1;32mCreating materialFile...\033[0m")
    #print(Material_parameters['variables'].items())
    matParams = Material_parameters['variables']
    setInputVariable(Material_parameters['copy_to'], 'enabledSlipSystems', matParams['enabledSlipSystems'])

    #for gPlaneNoise in matParams['glidePlaneNoise']:
    #    setInputVariable(Material_parameters['copy_to'],'glidePlaneNoise',gPlaneNoise)
    setInputVariable(Material_parameters['copy_to'],'glidePlaneNoise',matParams['glidePlaneNoise'])

    setInputVariable(Material_parameters['copy_to'],'atomsPerUnitCell',matParams['atomsPerUnitCell'])
    setInputVariable(Material_parameters['copy_to'],'dislocationMobilityType',matParams['dislocationMobilityType'])

    # Process polycrystal
    #pf = PolyCrystalFile('FeCrAl_Fe.txt')
    pf = PolyCrystalFile(Material_parameters['copy_to'].split('/')[-1])
    #pf.meshFile='unitCube24.msh'
    pf.meshFile=f'{str(file_templates["mesh"]).split("/")[-1]}'
    for param, value in Polycrystal_parameters['parameters'].items():
        setattr(pf, param, value)
    pf.write('inputFiles')

def constructDotGrid(NX, NY, unitVec1, unitVec2):
    # Construct non-orthogonal dot grid based on the non-minimized atom structure
    firstNearNeighborDist = 2.8
    dotGrid = np.zeros((NX*NY, 2)) # each row contains x and y position of the dot
    vec1 = unitVec1*firstNearNeighborDist;
    vec2 = unitVec2*firstNearNeighborDist;
    for yIdx in range(NY):
        for xIdx in range(NX):
            dotGrid[yIdx*NX + xIdx] = vec1*xIdx + vec2*yIdx
    return dotGrid

def plotComparisonCorrelation(realizationNum: int):
    matFile = str(file_templates['material_file']).split('/')[-1]
    b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
    mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')
    rho_SI = getValueInFile(f'inputFiles/{matFile}', 'rho_SI')

    latVector1 = None
    latVector2 = None
    copyFlag1 = False
    copyFlag2 = False
    with open('inputFiles/AlMg5_Cx_R100_ISF.vtk') as o:
        for line in o:
            if copyFlag1:
                latVector1 = np.array(line.strip().split(), dtype=float)
                copyFlag1 = False
            elif copyFlag2:
                latVector2 = np.array(line.strip().split(), dtype=float)
                copyFlag2 = False
            elif 'VECTORS lattice_basis1 double' in line:
                copyFlag1 = True
            elif 'VECTORS lattice_basis2 double' in line:
                copyFlag2 = True

    unitVec1 = latVector1/np.linalg.norm(latVector1)
    unitVec2 = latVector2/np.linalg.norm(latVector2)
    unitVec1 = unitVec1[:2] # drop the z axis value
    unitVec2 = unitVec2[:2]

    inputCorrelationField, i_dim = readVTKnoise('inputFiles/AlMg5_Cx_R100_ISF.vtk')
    NX_i, NY_i, NZ_i = i_dim

    # Construct non-orthogonal dot grid based on the non-minimized atom structure
    dotGrid_i = constructDotGrid(NX_i, NY_i, unitVec1, unitVec2)

    inputPath = "./inputFiles/"
    # reshape the correlation data to 1D array
    NX_s, NY_s, NZ_s = MDStackingFault_parameters['vectors']['gridSize']['value']
    correlationData = np.loadtxt(Path(f'{inputPath}/ensembledCorrelationR{realizationNum}.txt')).reshape(NX_s, NY_s)
    # convert unit from DDD to SI
    correlationData *= (mu0_SI*b_SI)**2

    # create callable color map
    color_map = mpl.colormaps["coolwarm"]

    # normalize the correlation coefficient array for the color assignment
    dotGrid_s = constructDotGrid(NX_s, NY_s, unitVec1, unitVec2)
    correlation = correlationData.reshape(NY_s*NX_s,1)
    c_normalized = (correlation - np.min(correlation.flatten())) / (np.max(correlation.flatten()) - np.min(correlation.flatten()))
    color = color_map(c_normalized)

    inputCorrelationField= inputCorrelationField.reshape(NY_i*NX_i,1)
    c_normalized_i = (inputCorrelationField - np.min(inputCorrelationField.flatten())) / (np.max(inputCorrelationField.flatten()) - np.min(inputCorrelationField.flatten()))
    color_i = color_map(c_normalized_i)

    # plot correlation data on the dot grid
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 8.6), dpi=500)

    axes[0].set_title("Input Correlation")
    axes[0].scatter(dotGrid_i[:,0], dotGrid_i[:,1], c=color_i, s=30)
    axes[1].set_title(f"Sampled R{realizationNum}")
    axes[1].scatter(dotGrid_s[:,0], dotGrid_s[:,1], c=color, s=10)
    fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(np.min(inputCorrelationField), np.max(inputCorrelationField)), cmap=color_map), ax=axes[0])
    fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(np.min(correlation), np.max(correlation)), cmap=color_map), ax=axes[1])
    plt.tight_layout() 
    fig.savefig(f'ensembleCorrelationR{realizationNum}.png', transparent=True)
    plt.close()

def plotComparisonHistogram(realizationNumbers: list):
    # Plotting on the first subplot
    fig, axs = plt.subplots(1, 1, figsize=(8,6), dpi=200)  # Adjusted for two subplots

    matFile = str(file_templates['material_file']).split('/')[-1]
    b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
    mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')

    inputPath = "./inputFiles/"
    for i, realizationNum in enumerate(realizationNumbers):
        data = np.loadtxt(f'{inputPath}/noiseDistributionR{realizationNum}.txt')
        data *= mu0_SI*b_SI
        axs.hist(data, density=True, alpha=0.3, bins='auto', label=f"R{realizationNum}")

        std = np.std(data)
        # Annotate in the top right, offsetting each annotation
        axs.text(
            0.98, 0.95 - 0.05*i,
            f'R{realizationNum} = STD {std:.3g} J/m2',
            transform=axs.transAxes,
            ha='right',
            va='top'
        )

    axs.grid(True)
    axs.set_title(f"")
    axs.set_xlabel(f"")
    axs.set_ylabel(f"")
    axs.legend()
    fig.savefig(f'sampledNoiseDistribution.png')
    plt.close()

def main() -> int:
    # Preparing input files
    folders=['inputFiles']
    for x in folders:
        # remove existing data
        if os.path.exists(x):
            shutil.rmtree(x)
        # create necessary folder structure for the simulation
        os.makedirs(x)

    # set simulation parameters in inputFiles
    processInitialConfigurations()

    setInputVariable(MDStackingFault_parameters['copy_to'], 'testNoiseSampling', str(1))
    simulationDir=os.path.abspath(".")
    ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)

    realizationNumbers = [10, 100, 1000]
    for realizationNum in realizationNumbers:
        plotComparisonCorrelation(realizationNum)

    # draw noise histogram
    plotComparisonHistogram(realizationNumbers)

    return 0

if __name__ == "__main__":
    main()
