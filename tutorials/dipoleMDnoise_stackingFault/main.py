import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
from pathlib import Path
sys.path.append("../../build/tools/pyMoDELib/")
import pyMoDELib

# Define all file templates, these files get copied to the inputFiles directory
file_templates = {
    'dd_file': Path('../../Library/DislocationDynamics/DD.txt'),
    'noise_file': Path('../../Library/GlidePlaneNoise/MDStackingFault.txt'),
    'correlationFile': Path(f'../../Library/GlidePlaneNoise/AlMg5_Cx_R100_ISF.vtk'),
    'material_file': Path('../../Library/Materials/AlMg5.txt'),
    'elastic_deformation_file': Path('../../Library/ElasticDeformation/ElasticDeformation.txt'),
    'mesh': Path('../../Library/Meshes/unitCube24.msh'),
    'microstructure': Path('../../Library/Microstructures/periodicDipoleIndividual.txt')
}

# initial config
DD_parameters = {
    'variables': {
        'setNodalVelocityToZeroInX': '0',
        'setNodalVelocityToZeroInY': '1',
        'setNodalVelocityToZeroInZ': '0',
        'useFEM': '0',
        'useDislocations': '1',
        'useInclusions': '0',
        'useElasticDeformation': '1',
        'useClusterDynamics': '0',
        'maxJunctionIterations': '1',
        'EwaldLengthFactor': '1',
        'timeSteppingMethod': 'fixed',  # 'adaptive' or 'fixed'
        'dtMax': '0.5',
        'dxMax': '1',  # max nodal displacement for adaptive time stepping
        'use_velocityFilter': '0',  # disable if noise is enabled
        'use_stochasticForce': '0',  # Langevin thermal noise
        'alphaLineTension': '1',  # dimensionless line tension scale factor
        'remeshFrequency': '0',  # node redistribution
        'quadPerLength': '3',
        'coreSize': '0.5',
        'Lmin': '5',  # min segment length (Burgers vector units)
        'Lmax': '20',  # max segment length (Burgers vector units)
        'outputFrequency': '10',  # output frequency
        'outputQuadraturePoints': '1',  # output quadrature data
        'computeElasticEnergyPerLength': '1',  # output energy data
        'glideSolverType': 'Galerkin',  # solver type ('Galerkin' or 'none')
        'climbSolverType': 'none', # type of clim solver, or none 
        'Nsteps': '100' # number of simulation steps 
    },
    'copy_to': 'inputFiles/DD.txt'
}

# initial config
MDStackingFault_parameters = {
    'variables': {
        # Scalar parameters (use setInputVariable)
        'type': 'MDStackingFaultNoise',
        'tag': 1,
        'seed': 1,
        # correlation file path in MDStackingFault.txt
        'correlationFile': Path('./AlMg5_Cx_R100_ISF.vtk'),
    },
    'vectors': {
        # Vector parameters with comments (use setInputVector)
        'gridSize': {
            'value': np.array([100, 100, 1]),
            'comment': 'number of grid points in each direction'
        },
        'gridSpacing_SI': {
            'value': np.array([1.0e-10, 1.0e-10, 1e-10]),
            'comment': 'grid spacing in each direction'
        }
    },
    'copy_to': 'inputFiles/MDStackingFault.txt'
}

# initial config
Material_parameters = {
    'variables': {
        'enabledSlipSystems': 'Shockley',
        'glidePlaneNoise': 'MDStackingFault.txt',
        'atomsPerUnitCell': '1',
        'dislocationMobilityType': 'default'
    },
    'copy_to': 'inputFiles/AlMg5.txt'
}

# initial config
Microstructure_parameters = {
    'vectors': {
        'slipSystemIDs': {
            'value': np.array([0,1]),
            'comment': 'slip system IDs for each dipole'
        },
        'exitFaceIDs': {
            'value': np.array([1,1]),
            'comment': '1 is for edge, 0 for screw'
        },
        'nodesPerLine': {
            'value': np.array([10,10]),
            'comment': 'number of extra nodes on each dipole'
        },
        'dipoleHeights': {
            'value': np.array([25,25]),
            'comment': 'height of each dipole, in number of planes'
        },
        'glideSteps': {
            'value': np.array([400,405]),
            'comment': 'step of each dipole in the glide plane'
        }
    },
    'matrix': {
        'dipoleCenters': {
            'value': np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]]),
        }
    },
    'initial_microstructure_file': 'inputFiles/initialMicrostructure.txt',
    'copy_to': 'inputFiles/periodicDipoleIndividual.txt'
}


# initial config
Polycrystal_parameters= {
    'parameters': {
        'absoluteTemperature': 1,
        'grain1globalX1': np.array([0, 1, 1]),
        'grain1globalX3': np.array([-1, 1, -1]),
        'boxEdges': np.array([[0,1,1], [-2,-1,1], [-1,1,-1]]),
        'boxScaling': np.array([1020, 36, 5800]),
        #'X0': np.array([0.5, 0.5, 0.5]),
        'X0': np.array([0., 0., 0.]),
        'periodicFaceIDs': np.array([0,1,2,3,4,5]),
        'gridSpacing_SI': np.array([2.86e-10, 2.86e-10])
    }
}

def setInitialConfiguration() -> None:

    # Copy all template files
    print("\033[1;32mCopying template files...\033[0m")
    for key, src_path in file_templates.items():
        dest = f"inputFiles/{src_path.name}"
        shutil.copy2(src_path.resolve(), dest)
        #shutil.copy2(src_path, dest)
        print(f"Created {dest}")

    print("\033[1;32mCreating ddFile\033[0m")
    # Apply all DD.txt parameters from the dictionary
    for param, value in DD_parameters['variables'].items():
        setInputVariable(DD_parameters['copy_to'], param, str(value))

    print("\033[1;32mCreating  noiseFile\033[0m")
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

    #for param, value in Material_parameters['variables'].items():
    #    setInputVariable(Material_parameters['copy_to'], param, str(value))

    # Process elastic deformation
    #print("\033[1;32mCreating elastic deformation file...\033[0m")
    #for param, data in Elastic_deformation_parameters['vectors'].items():
    #    setInputVector(
    #        Elastic_deformation_parameters['copy_to'],
    #        param,
    #        data['value'],
    #        data['comment']
    #    )

    # Process microstructure
    print("\033[1;32mCreating microstructureFile...\033[0m")
    for param, data in Microstructure_parameters['vectors'].items():
        setInputVector(
            Microstructure_parameters['copy_to'],
            param,
            data['value'],
            data['comment']
        )
    for param, data in Microstructure_parameters['matrix'].items():
        setInputMatrix(
            Microstructure_parameters['copy_to'],
            param,
            data['value']
        )

    # Write initial microstructure file
    with open(Microstructure_parameters['initial_microstructure_file'], 'w') as f:
        f.write(f"microstructureFile={Microstructure_parameters['copy_to'].split('/')[-1]};\n")

    # Process polycrystal
    #pf = PolyCrystalFile('FeCrAl_Fe.txt')
    pf = PolyCrystalFile(Material_parameters['copy_to'].split('/')[-1])
    #pf.meshFile='unitCube24.msh'
    pf.meshFile=f'{str(file_templates['mesh']).split('/')[-1]}'
    for param, value in Polycrystal_parameters['parameters'].items():
        setattr(pf, param, value)
    pf.write('inputFiles')

def main() -> int:
    # directory to store all the generated data
    dataStorageDir = Path("generatedData")

    # remove old data
    if os.path.exists(dataStorageDir):
        shutil.rmtree(dataStorageDir)

    # boundary condition control file
    elasticDeformationFile=str(file_templates['elastic_deformation_file']).split('/')[-1]

    # initial stresses to test in MPa
    stressToTests = [3000]
    strSign = -1
    for stress in stressToTests:
        # change sign
        stress *= strSign
        # Preparing input files
        folders=['evl','F', 'inputFiles']
        for x in folders:
            # remove existing data
            if os.path.exists(x):
                shutil.rmtree(x)
            # create necessary folder structure for the simulation
            os.makedirs(x)

        # set simulation parameters in inputFiles
        setInitialConfiguration()

        # extract material info
        matFile = str(file_templates['material_file']).split('/')[-1]
        b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
        mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')
        rho_SI = getValueInFile(f'inputFiles/{matFile}', 'rho_SI')
        cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
        convertTimeUnit = b_SI / cs  # [sec]
        convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))

        s13 = stress*convertMPaToMu
        #s13 = stress*convertMPaToMu
        # Voigt order is 11,22,33,12,23,13
        stressTensor = np.array([0.0, 0.0, 0.0, 0.0, 0.0, s13])
        print("\033[1;32mSetting elasticDeformatinoFile\033[0m")
        print(f"{str(Path("inputFiles") / elasticDeformationFile)}")
        setInputVector(Path("inputFiles") / elasticDeformationFile,
                    'ExternalStress0', stressTensor, 'applied stress')

        ## Iterate over configurations and run one step
        simulationDir=os.path.abspath(".")
        fFile=file_path = "F/F_0.txt"
        if os.path.exists(fFile):
            os.remove(fFile)
        fLabes=file_path = "F/F_labels.txt"
        if os.path.exists(fLabes):
            os.remove(fLabes)
        #periodicImageCases=np.array([[4,4,4,2]])
        periodicImageSize=np.array([1,1,1,1]).astype(int)
        #EwaldLengthFactor=[1]
        #setInputVariable('inputFiles/'+DDfile,'EwaldLengthFactor',str(EwaldLengthFactor))
        #setInputVector('inputFiles/'+DDfile,'periodicImageSize',periodicImageSize,'number of periodic images along each period vector')
        ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)
        mesh=ddBase.mesh
        xMax=mesh.xMax()
        xMin=mesh.xMin()
        xCenter=mesh.xCenter()

        # set simulation microstructure parameters
        spec=pyMoDELib.PeriodicDipoleIndividualSpecification()
        spec.slipSystemIDs=[0,1]
        spec.exitFaceIDs=[1,1]
        spec.dipoleCenters=[[0,0,0],[0,0,0]]
        spec.dipoleHeights=[2500, 2500]
        spec.nodesPerLine=[10,10]
        spec.glideSteps=[400, 405]

        microstructureGenerator=pyMoDELib.MicrostructureGenerator(ddBase)
        # write evl_0.txt (optional)
        microstructureGenerator.writeConfigFiles(0)

        microstructureGenerator.addPeriodicDipoleIndividual(spec)
        defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)
        defectiveCrystal.initializeConfiguration(microstructureGenerator.configIO)
        defectiveCrystal.runSteps()

        f_File = readFfile(f'{simulationDir}/F')
        print(f_File)

        #if not os.path.exists(dataStorageDir):
        #    os.makedirs(dataStorageDir)

        #targetDir = Path(dataStorageDir)/str(stress)
        #os.makedirs(targetDir, exist_ok=True)
        #for f in folders:
        #    shutil.move(f, targetDir)

    return 0

if __name__ == "__main__":
    main()
