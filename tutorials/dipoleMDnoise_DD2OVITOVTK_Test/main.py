# /opt/local/bin/python3.12 test.py
import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
from pathlib import Path
sys.path.append("../../build/tools/pyMoDELib/")
import pyMoDELib


# Define all file templates (using pathlib for cross-platform paths)
file_templates = {
    'dd_file': Path('../../Library/DislocationDynamics/DD.txt'),
    'noise_file': Path('../../Library/GlidePlaneNoise/MDSolidSolution.txt'),
    'material_file': Path('../../Library/Materials/FeCrAl_Fe.txt'),
    'elastic_deformation_file': Path('../../Library/ElasticDeformation/ElasticDeformation.txt'),
    'mesh': Path('../../Library/Meshes/unitCube24.msh'),
    'microstructure': Path('../../Library/Microstructures/periodicDipoleIndividual.txt')
}


DD_parameters = {
    'variables': {
        'setNodalVelocityBaseX': '0',
        'setNodalVelocityBaseY': '0',
        'setNodalVelocityBaseZ': '1',
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


MDSolidSolution_parameters = {
    'variables': {
        # Scalar parameters (use setInputVariable)
        'type': 'MDSolidSolutionNoise',
        'tag': 1,
        'seed': 1,
        'a_cai_SI': 0.0,
        # File paths (convert to absolute paths)
        'correlationFile_L': Path('../../Library/GlidePlaneNoise/MoDELCompatible_FeCr8_xz.vtk').resolve(),
        'correlationFile_T': Path('../../Library/GlidePlaneNoise/MoDELCompatible_FeCr8_yz.vtk').resolve(),
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
    'copy_to': 'inputFiles/MDSolidSolution.txt'
}


Material_parameters = {
    'variables': {
        'enabledSlipSystems': 'full<111>{110}',
        'glidePlaneNoise': 'MDSolidSolution.txt',
        'atomsPerUnitCell': '1',
        'dislocationMobilityType': 'default'
    },
    'copy_to': 'inputFiles/FeCrAl_Fe.txt'
}


#Elastic_deformation_parameters= {
#    'vectors': {
#        'ExternalStress0': {
#            'value': np.array([0.0, 0.0, 0.0, 0.001, 0.0, 0.0]),
#            'comment': 'applied stress'
#        }
#    },
#    'copy_to': 'inputFiles/ElasticDeformation.txt'
#}


Microstructure_parameters = {
    'vectors': {
        'slipSystemIDs': {
            'value': np.array([0]),
            'comment': 'slip system IDs for each dipole'
        },
        'exitFaceIDs': {
            'value': np.array([0]),
            'comment': '4 is for edge, 2 for screw'
        },
        'nodesPerLine': {
            'value': np.array([10]),
            'comment': 'number of extra nodes on each dipole'
        },
        'dipoleHeights': {
            'value': np.array([800]),
            'comment': 'height of each dipole, in number of planes'
        },
        'glideSteps': {
            'value': np.array([10.0]),
            'comment': 'step of each dipole in the glide plane'
        }
    },
    'initial_microstructure_file': 'inputFiles/initialMicrostructure.txt',
    'copy_to': 'inputFiles/periodicDipoleIndividual.txt'
}


Polycrystal_parameters= {
    'parameters': {
        'absoluteTemperature': 1,
        'grain1globalX1': np.array([1, 1, -1]),
        'grain1globalX3': np.array([1, -2, -1]),
        'boxEdges': np.array([[1, 1, -1], [1, 0, 1], [1, -2, -1]]),
        'boxScaling': np.array([1200, 1000, 110]),
        'X0': np.array([0, 0, 0]),
        'periodicFaceIDs': np.array([-1]),
        'gridSpacing_SI': np.array([1.0e-10, 1.0e-10])
    }
}


def process_configurations() -> None:
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
    for param, value in Material_parameters['variables'].items():
        setInputVariable(Material_parameters['copy_to'], param, str(value))

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
    #stressToTests = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900]
    stressToTests = [50]
    for stress in stressToTests:
        # Preparing input files
        folders=['evl','F', 'inputFiles']
        for x in folders:
            # remove existing data
            if os.path.exists(x):
                shutil.rmtree(x)
            # create necessary folder structure for the simulation
            os.makedirs(x)

        # set simulation parameters in inputFiles
        process_configurations()

        # extract material info
        matFile = str(file_templates['material_file']).split('/')[-1]
        b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
        mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')
        rho_SI = getValueInFile(f'inputFiles/{matFile}', 'rho_SI')
        cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
        convertTimeUnit = b_SI / cs  # [sec]
        convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))

        s12 = stress*convertMPaToMu
        # Voigt order is 11,22,33,12,23,13
        stressTensor = np.array([0.0, 0.0, 0.0, s12, 0.0, 0.0])
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
        spec=pyMoDELib.PeriodicDipoleIndividualSpecification()
        #spec=pyMoDELib.PeriodicDipoleIndividualSpecification()
        spec.slipSystemIDs=[0]
        #spec.slipSystemIDs=[0]
        spec.exitFaceIDs=[0]
        #spec.dipoleCenters=np.array(xCenter)
        spec.dipoleCenters=np.array(xCenter)
        spec.dipoleHeights=[600]
        spec.nodesPerLine=[80]
        spec.glideSteps=[10]
        microstructureGenerator=pyMoDELib.MicrostructureGenerator(ddBase)
        microstructureGenerator.writeConfigFiles(0) # write evl_0.txt (optional)

        microstructureGenerator.addPeriodicDipoleIndividual(spec)
        defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)
        defectiveCrystal.initializeConfiguration(microstructureGenerator.configIO)
        defectiveCrystal.runSteps()

        pyMoDELib.DD2OvitoVtk(simulationDir)

        #if not os.path.exists(dataStorageDir):
        #    os.makedirs(dataStorageDir)

        #targetDir = Path(dataStorageDir)/str(stress)
        #os.makedirs(targetDir, exist_ok=True)
        #for f in folders:
        #    shutil.move(f, targetDir)

    #X=np.linspace(0, 101, num=101)
    #print(f"x.shape = {X.shape}")
    #F,Flabels=readFfile('./F')
    #E=getFarray(F,Flabels,'dislocation elastic energy [mu b^3]')
    #data=np.empty([np.size(X),8])
    #data[:,0]=X
    #data[:,1]=E
    #data[:,2]=S11
    #data[:,3]=S22
    #data[:,4]=S33
    #data[:,5]=S23
    #data[:,6]=S13
    #data[:,7]=S12

    return 0

if __name__ == "__main__":
    main()
