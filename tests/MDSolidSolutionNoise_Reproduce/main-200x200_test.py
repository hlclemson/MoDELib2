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
    'noise_file': Path('../../Library/GlidePlaneNoise/MDSolidSolution.txt'),
    'material_file': Path('../../Library/Materials/FeCrAl_Fe.txt'),
    'elastic_deformation_file': Path('../../Library/ElasticDeformation/ElasticDeformation.txt'),
    'mesh': Path('../../Library/Meshes/unitCube24.msh'),
    'microstructure': Path('../../Library/Microstructures/periodicDipoleIndividual.txt')
}

MDSolidSolution_parameters = {
    'variables': {
        # Scalar parameters (use setInputVariable)
        'type': 'MDSolidSolutionNoise',
        'tag': 1,
        'seed': 2,
        'outputNoise': 1,
        'noiseFile_L': 'FeCr8_L_noise.vtk',
        'noiseFile_T': 'FeCr8_T_noise.vtk',
        # File paths (convert to absolute paths)
        'correlationFile_L': Path('../../Library/GlidePlaneNoise/MoDELCompatible_FeCr8_xz.vtk').resolve(),
        'correlationFile_T': Path('../../Library/GlidePlaneNoise/MoDELCompatible_FeCr8_yz.vtk').resolve(),
    },
    'vectors': {
        # Vector parameters with comments (use setInputVector)
        'gridSize': {
            'value': np.array([200, 200, 1]),
            'comment': 'number of grid points in each direction'
        },
        'gridSpacing_SI': {
            'value': np.array([1.0e-10, 1.0e-10, 1e-10]),
            'comment': 'grid spacing in each direction'
        }
    },
    'copy_to': 'inputFiles/MDSolidSolution.txt'
}

Material_parameters = {
    'variables': {
        'enabledSlipSystems': 'Shockley',
        #'glidePlaneNoise': ['MDSolidSolution.txt', 'MDStackingFault.txt'],
        'glidePlaneNoise': 'MDSolidSolution.txt',
        'atomsPerUnitCell': '1',
        'dislocationMobilityType': 'default'
    },
    'copy_to': 'inputFiles/FeCrAl_Fe.txt'
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
    shutil.copy2(MDSolidSolution_parameters['variables']['correlationFile_L'], 'inputFiles/')
    shutil.copy2(MDSolidSolution_parameters['variables']['correlationFile_T'], 'inputFiles/')
    # copy declared variables to the configuration txt file
    for param, value in MDSolidSolution_parameters['variables'].items():
        if 'correlationFile' in param:
            relativePath = f"{str(value).split('/')[-1]}"
            setInputVariable(MDSolidSolution_parameters['copy_to'], param, relativePath)
        else:
            setInputVariable(MDSolidSolution_parameters['copy_to'], param, str(value))
    for param, data in MDSolidSolution_parameters['vectors'].items():
        setInputVector(
            MDSolidSolution_parameters['copy_to'],
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

    inputCorrelationField_xz, i_dim_xz = readVTKnoise('./inputFiles/MoDELCompatible_FeCr8_xz.vtk')
    NX_i_xz, NY_i_xz, NZ_i_xz = i_dim_xz
    inputCorrelationField_yz, i_dim_yz = readVTKnoise('./inputFiles/MoDELCompatible_FeCr8_yz.vtk')
    NX_i_yz, NY_i_yz, NZ_i_yz = i_dim_yz

    inputPath = "./inputFiles/"
    # reshape the correlation data to 1D array
    NX_s, NY_s, NZ_s = MDSolidSolution_parameters['vectors']['gridSize']['value']
    correlationData_xz = np.loadtxt(Path(f'{inputPath}/ensembledCorrelation_xz_R{realizationNum}.txt')).reshape(NX_s, NY_s)
    correlationData_yz = np.loadtxt(Path(f'{inputPath}/ensembledCorrelation_yz_R{realizationNum}.txt')).reshape(NX_s, NY_s)
    # convert unit from DDD to SI
    correlationData_xz *= mu0_SI**2
    correlationData_yz *= mu0_SI**2

    # create callable color map
    #color_map = mpl.colormaps["coolwarm"]
    # set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 5}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)

    fig, axs = plt.subplots(1, 2, figsize=(16,6), dpi=200)  # Adjusted for two subplots
    x = np.arange(NY_i_xz)
    y = np.arange(NX_i_xz)
    xTickPos = np.arange(len(x))
    yTickPos = np.arange(len(y))
    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
    axs[0].set_xticks(xTickPos)
    axs[0].set_xticklabels(labels=xTicks, rotation=45, ha="right")
    axs[0].set_yticks(yTickPos)
    axs[0].set_yticklabels(labels=yTicks)
    axs[0].xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[0].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[0].set_title('original')  # Set the title for each axes
    X = x
    Y = y
    X, Y = np.meshgrid(X, Y)
    Z = inputCorrelationField_xz
    c = axs[0].pcolormesh(X, Y, Z, cmap=cm.coolwarm)
    fig.colorbar(c, ax=axs[0])

    x = np.arange(NY_s)
    y = np.arange(NX_s)
    xTickPos = np.arange(len(x))
    yTickPos = np.arange(len(y))
    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
    axs[1].set_xticks(xTickPos)
    axs[1].set_xticklabels(labels=xTicks, rotation=45, ha="right")
    axs[1].set_yticks(yTickPos)
    axs[1].set_yticklabels(labels=yTicks)
    axs[1].xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[1].set_title('sampled')  # Set the title for each axes
    X = x
    Y = y
    X, Y = np.meshgrid(X, Y)
    Z = correlationData_xz
    c = axs[1].pcolormesh(X, Y, Z, cmap=cm.coolwarm)
    fig.colorbar(c, ax=axs[1])

    # Save figure 
    plt.savefig(f'ensembleCorrelation_xz_R{realizationNum}.png', transparent=True)
    plt.close()

    fig, axs = plt.subplots(1, 2, figsize=(16,6), dpi=200)  # Adjusted for two subplots
    x = np.arange(NY_i_yz)
    y = np.arange(NX_i_yz)
    xTickPos = np.arange(len(x))
    yTickPos = np.arange(len(y))
    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
    axs[0].set_xticks(xTickPos)
    axs[0].set_xticklabels(labels=xTicks, rotation=45, ha="right")
    axs[0].set_yticks(yTickPos)
    axs[0].set_yticklabels(labels=yTicks)
    axs[0].xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[0].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[0].set_title('original')  # Set the title for each axes
    X = x
    Y = y
    X, Y = np.meshgrid(X, Y)
    Z = inputCorrelationField_yz
    c = axs[0].pcolormesh(X, Y, Z, cmap=cm.coolwarm)
    fig.colorbar(c, ax=axs[0])

    x = np.arange(NY_s)
    y = np.arange(NX_s)
    xTickPos = np.arange(len(x))
    yTickPos = np.arange(len(y))
    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
    axs[1].set_xticks(xTickPos)
    axs[1].set_xticklabels(labels=xTicks, rotation=45, ha="right")
    axs[1].set_yticks(yTickPos)
    axs[1].set_yticklabels(labels=yTicks)
    axs[1].xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[1].set_title('sampled')  # Set the title for each axes
    X = x
    Y = y
    X, Y = np.meshgrid(X, Y)
    Z = correlationData_yz
    c = axs[1].pcolormesh(X, Y, Z, cmap=cm.coolwarm)
    fig.colorbar(c, ax=axs[1])

    # Save figure 
    plt.savefig(f'ensembleCorrelation_yz_R{realizationNum}.png', transparent=True)
    plt.close()

def plotComparisonHistogram(realizationNumbers: list):
    # Plotting on the first subplot
    fig, axs = plt.subplots(1, 2, figsize=(16,6), dpi=200)  # Adjusted for two subplots

    matFile = str(file_templates['material_file']).split('/')[-1]
    b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
    mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')

    inputPath = "./inputFiles/"
    for i, realizationNum in enumerate(realizationNumbers):
        data_xz = np.loadtxt(f'{inputPath}/noiseDistribution_xz_R{realizationNum}.txt')
        data_yz = np.loadtxt(f'{inputPath}/noiseDistribution_yz_R{realizationNum}.txt')
        data_xz *= mu0_SI*b_SI
        data_yz *= mu0_SI*b_SI
        axs[0].hist(data_xz, density=True, alpha=0.3, bins='auto', label=f"R{realizationNum}")
        axs[1].hist(data_yz, density=True, alpha=0.3, bins='auto', label=f"R{realizationNum}")

        std_xz = np.std(data_xz)
        std_yz = np.std(data_yz)
        # Annotate in the top right, offsetting each annotation
        axs[0].text(
            0.98, 0.95 - 0.05*i,
            f'R{realizationNum} = STD {std_xz:.3g} Pa2',
            transform=axs[0].transAxes,
            ha='right',
            va='top'
        )
        axs[1].text(
            0.98, 0.95 - 0.05*i,
            f'R{realizationNum} = STD {std_yz:.3g} Pa2',
            transform=axs[1].transAxes,
            ha='right',
            va='top'
        )

    axs[0].grid(True)
    axs[0].set_title(f"noise xz")
    axs[0].set_xlabel(f"")
    axs[0].set_ylabel(f"")
    axs[0].legend()
    axs[1].grid(True)
    axs[1].set_title(f"noise yz")
    axs[1].set_xlabel(f"")
    axs[1].set_ylabel(f"")
    axs[1].legend()
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

    setInputVariable(MDSolidSolution_parameters['copy_to'], 'testNoiseSampling', str(1))
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
