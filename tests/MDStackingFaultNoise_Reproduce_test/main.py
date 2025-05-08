import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
from pathlib import Path
import matplotlib as mpl
sys.path.append("../../build/tools/pyMoDELib/")
import pyMoDELib
from scipy.signal import correlate2d
from scipy.stats import norm
import vtk
import matplotlib.transforms as mtransforms
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import Axes, HostAxes, angle_helper
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axisartist.grid_finder import DictFormatter

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
            'value': np.array([100, 100, 1]),
            #'value': np.array([60, 60, 1]),
            #'value': np.array([30, 27, 1]),
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

def process_configurations() -> None:
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

def drawDotCorrelationPlot(dotGrid: np.ndarray, basis1, basis2, correlationData, labels, fileName, subDataDir):
    # create plot object with preconfigured specification
    fig, ax = createPlotObjectsNonOrtho()
    title = f"{labels['title']}"
    xlabel = f"{labels['xlabel']}"
    ylabel = f"{labels['ylabel']}"
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # reshape the correlation data to 1D array
    xGridSize = correlationData.shape[1]
    yGridSize = correlationData.shape[0]
    correlationData = correlationData.reshape(xGridSize*yGridSize,1)

    # create callable color map
    color_map = mpl.colormaps["viridis"]

    # normalize the correlation coefficient array for the color assignment
    c_normalized = (correlationData - np.min(correlationData.flatten())) / (np.max(correlationData.flatten()) - np.min(correlationData.flatten()))
    color = color_map(c_normalized)

    # plot correlation data on the dot grid
    ax.scatter(dotGrid[:,0], dotGrid[:,1], c=color)
    fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(np.min(correlationData), np.max(correlationData)), cmap=color_map), ax=ax)

    #ax.set_aspect(1)
    #ax.set_xlim(0, 10)
    #ax.set_ylim(0, 10)

    #ax.axis["t"] = ax.new_floating_axis(0, 0, axis_direction='bottom')
    ##ax.axis["t"].set_axisline_style("->", size=1.0)
    ##ax.axis["t"].set_ticklabel_direction('+')
    #ax.axis["t"].invert_ticklabel_direction()
    #ax.axis["t"].set_axislabel_direction('-')
    #ax.axis["t"].set_axis_direction('left')
    #ax.axis["t"].set(transform=ax.transAxes)

    #ax.axis["t2"] = ax.new_floating_axis(0, 20)
    ax.grid(True, zorder=0)

    # save figure 
    figFolder = './figures'
    if not os.path.exists(f'{figFolder}/{subDataDir}'):
        os.system(f'mkdir {figFolder}/{subDataDir}')
    plt.savefig(f'{figFolder}/{subDataDir}/{fileName}')
    # free up memory
    plt.close()

def mapPlot_AA(x, y, x_AA, y_AA, dat, labels, fileName):
    # set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 10}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)
    fig, ax = plt.figure(figsize=(8,6), dpi=200), plt.subplot(1,1,1)

    title = f"{labels['title']}"
    xlabel = f"{labels['xlabel']}"
    ylabel = f"{labels['ylabel']}"
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    xTickPos = x_AA
    yTickPos = y_AA
    xTicks = x_AA
    yTicks = y_AA
    ax.set_xticks(xTickPos, labels=x_AA)
    ax.set_yticks(yTickPos, labels=y_AA)

    ax.xaxis.set_major_formatter('{x:.2f}')
    ax.yaxis.set_major_formatter('{x:.2f}')
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

    X = x_AA
    Y = y_AA
    X, Y = np.meshgrid(X, Y)
    Z = dat
    c = ax.pcolormesh(X, Y, Z, cmap=cm.viridis)
    fig.colorbar(c, ax=ax)

    # save figure 
    figFolder = './figures'
    if not os.path.exists(figFolder):
        os.system(f'mkdir {figFolder}')
    plt.savefig(f'{figFolder}/{fileName}')
    # free up memory
    plt.close()

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

def main() -> int:
    # boundary condition control file
    #elasticDeformationFile=str(file_templates['elastic_deformation_file']).split('/')[-1]

    # extract material info
    #matFile = str(file_templates['material_file']).split('/')[-1]
    #b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
    #mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')
    #rho_SI = getValueInFile(f'inputFiles/{matFile}', 'rho_SI')
    #cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
    #convertTimeUnit = b_SI / cs  # [sec]
    #convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))

    # Preparing input files
    folders=['sampledNoises']
    for x in folders:
        # remove existing data
        if os.path.exists(x):
            shutil.rmtree(x)
        # create necessary folder structure for the simulation
        os.makedirs(x)


    realizations = 10
    sampledNoises = []
    for iterNum in np.arange(1,realizations+1,1):
        # Preparing input files
        folders=['inputFiles']
        for x in folders:
            # remove existing data
            if os.path.exists(x):
                shutil.rmtree(x)
            # create necessary folder structure for the simulation
            os.makedirs(x)

        # set simulation parameters in inputFiles
        process_configurations()

        #VECTORS lattice_basis1 double
        #89.853 0 0
        #VECTORS lattice_basis2 double
        #40.63565057065618 70.38301138699174 0
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

        unitVec1_i = latVector1/np.linalg.norm(latVector1)
        unitVec2_i = latVector2/np.linalg.norm(latVector2)
        unitVec1_i = unitVec1_i[:2] # drop the z axis value
        unitVec2_i = unitVec2_i[:2] 

        inputCorrelationField, i_dim = readVTKnoise('inputFiles/AlMg5_Cx_R100_ISF.vtk')
        NX_i, NY_i, NZ_i = i_dim
        ## Iterate over configurations and run one step
        simulationDir=os.path.abspath(".")
        #
        setInputVariable(MDStackingFault_parameters['copy_to'], 'seed', str(iterNum))
        ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)
        #print(dir(ddBase))
        shutil.copy2(f'inputFiles/{MDStackingFault_parameters['variables']['noiseFile']}', f'sampledNoises/sampled_noise_{iterNum}.vtk')
        if iterNum==1:
            sampledNoiseField, s_dim = readVTKnoise(f'sampledNoises/sampled_noise_{iterNum}.vtk')
        else:
            sampledNoiseField, _ = readVTKnoise(f'sampledNoises/sampled_noise_{iterNum}.vtk')
        sampledNoises.append(sampledNoiseField)


    sampledNoises = np.array(sampledNoises)
    sampledNoises = np.mean(sampledNoises, axis=0)
    print(sampledNoises.shape)
    print(f"mean = {np.mean(sampledNoises)}")
    print(f"std = {np.std(sampledNoises)}")

    #print(f"mu0_SI = {mu0_SI}")
    #print(f"b_SI = {b_SI}")
    #print(f"mu0_SI*b_SI = {mu0_SI*b_SI}")
    NX_s, NY_s, NZ_s = s_dim


    #correlation = correlate2d(sampledNoiseField, sampledNoiseField, mode='same')
    correlation = correlate2d(sampledNoises, sampledNoises, mode='same')

    # Construct non-orthogonal dot grid based on the non-minimized atom structure
    firstNearNeighborDist = 2.8
    dotGrid = np.zeros((NX_s*NY_s, 2)) # each row contains x and y position of the dot
    vec1 = unitVec1*firstNearNeighborDist;
    vec2 = unitVec2*firstNearNeighborDist;
    for yIdx in range(NY_s):
        for xIdx in range(NX_s):
            dotGrid[yIdx*NX_s + xIdx] = vec1*xIdx + vec2*yIdx
    dotGrid_i = constructDotGrid(NX_i, NY_i, unitVec1_i, unitVec2_i)

    dotGrid_s = constructDotGrid(NX_s, NY_s, unitVec1, unitVec2)

    # reshape the correlation data to 1D array
    #xGridSize = correlationData.shape[1]
    #yGridSize = correlationData.shape[0]
    #correlationData = correlationData.reshape(xGridSize*yGridSize,1)

    # create callable color map
    color_map = mpl.colormaps["coolwarm"]

    # normalize the correlation coefficient array for the color assignment
    correlation = correlation.reshape(NY_s*NX_s,1)
    print(correlation.shape)
    c_normalized = (correlation - np.min(correlation.flatten())) / (np.max(correlation.flatten()) - np.min(correlation.flatten()))
    color = color_map(c_normalized)
    #color = color_map(correlation)
    #color = color_map(correlation)

    inputCorrelationField= inputCorrelationField.reshape(NY_i*NX_i,1)
    c_normalized_i = (inputCorrelationField - np.min(inputCorrelationField.flatten())) / (np.max(inputCorrelationField.flatten()) - np.min(inputCorrelationField.flatten()))
    color_i = color_map(c_normalized_i)

    # plot correlation data on the dot grid
    factorFor60degree = 0.577350269189626
    def tr(x, y): return x+y*factorFor60degree , y
    def inv_tr(x, y): return x-y*factorFor60degree , y
    # remove x and y ticks
    xticks, yticks = [], []
    tick_formatter1 = DictFormatter(dict(xticks))
    tick_formatter2 = DictFormatter(dict(yticks))
    grid_helper = GridHelperCurveLinear((tr, inv_tr), tick_formatter1=tick_formatter1, tick_formatter2=tick_formatter2)
    #fig, ax = plt.figure(figsize=(10, 4.3), dpi=500), plt.subplot(1, 1, 1, axes_class=Axes, grid_helper=grid_helper)
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 8.6), dpi=500)

    axes[0].scatter(dotGrid_i[:,0], dotGrid_i[:,1], c=color_i, s=10)
    axes[1].scatter(dotGrid_s[:,0], dotGrid_s[:,1], c=color, s=10)
    fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(np.min(inputCorrelationField), np.max(inputCorrelationField)), cmap=color_map), ax=axes[0])
    fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(np.min(correlation), np.max(correlation)), cmap=color_map), ax=axes[1])
    plt.tight_layout() 
    fig.savefig('test.png', transparent=True)
    plt.close()

    # Plot the tiled heatmap with coolwarm colormap and colorbar
    #plt.figure(figsize=(10, 10), dpi=200)
    #im = plt.imshow(sampledNoiseField, cmap='coolwarm')
    #plt.colorbar(im, label='Noise Intensity')  # Add colorbar
    #plt.xlabel('X-axis')
    #plt.ylabel('Z-axis')
    #plt.savefig("noise.png")

    ## Plot the tiled heatmap with coolwarm colormap and colorbar
    #plt.figure(figsize=(10, 10), dpi=200)
    #correlation = correlate2d(sampledNoiseField, sampledNoiseField, mode='same')
    #im = plt.imshow(correlation, cmap='coolwarm')
    #plt.colorbar(im, label='Correlation Intensity')  # Add colorbar
    #plt.xlabel('X-axis')
    #plt.ylabel('Z-axis')
    #plt.savefig("correlationFromNoise.png")
    #plt.close()

    return 0

if __name__ == "__main__":
    main()
