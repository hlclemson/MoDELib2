import numpy as np
from pathlib import Path

# Define all file templates, these files get copied to the inputFiles directory
files_to_copy_from_lib = {
    'dd_file': Path('../../Library/DislocationDynamics/DD.txt'),
    #'noise_file_md_sf': Path('../../Library/GlidePlaneNoise/MDStackingFault.txt'),
    'noise_file_md_ss': Path('../../Library/GlidePlaneNoise/MDSolidSolution.txt'),
    #"noise_file_an_ss": Path("../../Library/GlidePlaneNoise/AnalyticalSolidSolutionNoise.txt"),
    #'correlation_file_md_sf': Path(f'../../Library/GlidePlaneNoise/AlMg5_Cx_R100_ISF.vtk'),
    'correlation_file_md_ss_xz': Path(f'../../Library/GlidePlaneNoise/MoDELCompatible_FeCr8_xz.vtk'),
    'correlation_file_md_ss_yz': Path(f'../../Library/GlidePlaneNoise/MoDELCompatible_FeCr8_yz.vtk'),
    'material_file': Path('../../Library/Materials/FeCrAl_Fe.txt'),
    'elastic_deformation_file': Path('../../Library/ElasticDeformation/ElasticDeformation.txt'),
    'mesh': Path('../../Library/Meshes/unitCube24.msh'),
    'microstructure': Path('../../Library/Microstructures/periodicDipoleIndividual.txt')
}

# initial config
DD_parameters = {
    'useFEM': '0',
    'useDislocations': '1',
    'useInclusions': '0',
    'useElasticDeformation': '1',
    'useClusterDynamics': '0',
    'maxJunctionIterations': '0',
    'EwaldLengthFactor': '1',
    'timeSteppingMethod': 'fixed',  # 'adaptive' or 'fixed'
    'dtMax': '0.06',
    'dxMax': '1',  # max nodal displacement for adaptive time stepping
    'use_velocityFilter': '0',  # disable if noise is enabled
    'use_stochasticForce': '0',  # Langevin thermal noise
    'alphaLineTension': '1',  # dimensionless line tension scale factor
    'remeshFrequency': '0',  # node redistribution
    'quadPerLength': '6',
    'coreSize': '0',
    'Lmin': '1',  # min segment length (Burgers vector units)
    'Lmax': '5',  # max segment length (Burgers vector units)
    'outputFrequency': '500',  # output frequency
    #'outputFrequency': '10',  # output frequency
    'outputQuadraturePoints': '1',  # output quadrature data
    'computeElasticEnergyPerLength': '1',  # output energy data
    'glideSolverType': 'Galerkin',  # solver type ('Galerkin' or 'none')
    'climbSolverType': 'none', # type of clim solver, or none 
    'Nsteps': '10001', # number of simulation steps
    #'Nsteps': '101', # number of simulation steps
    'nodalVelocityConstraints': np.array([0,0,1])
}

# initial config
MDStackingFault_parameters = {
    # Scalar parameters (use setInputVariable)
    'type': 'MDStackingFaultNoise',
    'tag': 1,
    'seed': 1,
    # correlation file path in MDStackingFault.txt
    'correlationFile': Path('./AlMg5_Cx_R100_ISF.vtk'),
    # Vector parameters with comments (use setInputVector)
    'gridSize': np.array([2048, 2048, 1]),
    'gridSpacing_SI': np.array([2.86e-10, 2.86e-10, 1e-10]),
}

# initial config
MDSolidSolution_parameters = {
    # Scalar parameters (use setInputVariable)
    'type': 'MDSolidSolutionNoise',
    'tag': 1,
    'seed': 1,
    "a_cai_SI": 0,
    # correlation file path in MDSolidSolution.txt
    'correlationFile_xz': Path('./MoDELCompatible_FeCr8_xz.vtk'),
    'correlationFile_yz': Path('./MoDELCompatible_FeCr8_yz.vtk'),
    # Vector parameters with comments (use setInputVector)
    'gridSize': np.array([2000, 2000, 1]),
    'gridSpacing_SI': np.array([2.86e-10, 2.86e-10, 1e-10]),
}

AnalyticalSolidSolution_parameters = {
    # Scalar parameters (use setInputVariable)
    "type": "AnalyticalSolidSolutionNoise",
    "tag": 1,
    "seed": 1,
    "a_cai_SI": 0,
    "gridSize": np.array([256, 256, 64]),
    "gridSpacing_SI": np.array([1.0e-10, 1.0e-10, 1e-10]),
    "MSSS_SI": 0.09e18
}

# initial config
Material_parameters = {
    #'enabledSlipSystems': 'Shockley',
    'enabledSlipSystems': 'full<111>{110}',
    #'glidePlaneNoise': 'MDStackingFault.txt',
    'glidePlaneNoise': 'MDSolidSolution.txt',
    #'glidePlaneNoise': 'AnalyticalSolidSolutionNoise.txt',
    'atomsPerUnitCell': '1',
    'dislocationMobilityType': 'default'
}

# initial config
Microstructure_parameters = {
    'slipSystemIDs': np.array([0]),
    'exitFaceIDs': np.array([0]),
    #'nodesPerLine': np.array([80]),
    'nodesPerLine': np.array([120]),
    'dipoleHeights': np.array([1200]),
    'glideSteps': np.array([400]),
    #'dipoleHeights': np.array([10]),
    #'glideSteps': np.array([10]),
    'dipoleCenters': np.array([[0.0,0.0,0.0]]),
}

# initial config
Polycrystal_parameters= {
    'absoluteTemperature': 1,
    'grain1globalX1': np.array([1, 1, -1]),
    'grain1globalX3': np.array([1, -2, -1]),
    'boxEdges': np.array([[1, 1, -1], [1, 0, 1], [1, -2, -1]]),
    'boxScaling': np.array([1000, 2000, 109]),
    'X0': np.array([0, 0, 0]),
    'periodicFaceIDs': np.array([0,1,2,3,4,5]),
}

