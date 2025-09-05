import numpy as np
from pathlib import Path

# Define all file templates, these files get copied to the inputFiles directory
files_to_copy_from_lib = {
    'dd_file': Path('../../Library/DislocationDynamics/DD.txt'),
    "noise_file_an_ssw": Path("../../Library/GlidePlaneNoise/AnalyticalSolidSolutionWhiteNoise.txt"),
    'material_file': Path('../../Library/Materials/AlMg50.txt'),
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
    'alphaLineTension': '1.356',  # dimensionless line tension scale factor
    'remeshFrequency': '0',  # node redistribution
    'quadPerLength': '6',
    'coreSize': '1.1713',
    #'Lmin': '5',  # min segment length (Burgers vector units)
    #'Lmax': '20',  # max segment length (Burgers vector units)
    'outputFrequency': '100',  # output frequency
    'outputQuadraturePoints': '1',  # output quadrature data
    'computeElasticEnergyPerLength': '1',  # output energy data
    'glideSolverType': 'Galerkin',  # solver type ('Galerkin' or 'none')
    'climbSolverType': 'none', # type of clim solver, or none 
    'Nsteps': '12001', # number of simulation steps
    #'Nsteps': '1001', # number of simulation steps
    #'nodalVelocityConstraints': np.array([0,1,0])
    'nodalVelocityConstraints': np.array([1,0,0]),
    'computeDDinteractions': 0
}

AnalyticalSolidSolutionWhite_parameters = {
    # Scalar parameters (use setInputVariable)
    "type": "AnalyticalSolidSolutionWhiteNoise",
    "tag": 1,
    "seed": 1,
    #"spreadLstress_SI": 1.0355e-10, # eigenstrain spreading distance (a)
    #"a_cai_SI": 5.6e-10, # non-singular param (a_n), edge
    "a_cai_SI": 3.34e-10, # non-singular param (a_n), screw
    #"gridSize": np.array([2560, 2560, 64]),
    "gridSize": np.array([256, 256, 64]),
    "gridSpacing_SI": np.array([1.0e-10, 1.0e-10, 1e-10]),
    "dislocation_length_SI": 1430,
    "MSSS_SI": 0.6114e18
}

# initial config
Material_parameters = {
    'enabledSlipSystems': 'full',
    #'enabledSlipSystems': 'Shockley',
    #'glidePlaneNoise': 'MDStackingFault.txt',
    #'glidePlaneNoise': 'MDSolidSolution.txt',
    'glidePlaneNoise': 'AnalyticalSolidSolutionWhiteNoise.txt',
    'atomsPerUnitCell': '1',
    'dislocationMobilityType': 'default'
}

# initial config
Microstructure_parameters = {
    'slipSystemIDs': np.array([0]),
    'exitFaceIDs': np.array([2]),
    'nodesPerLine': np.array([10]),
    'dipoleHeights': np.array([2000]),
    'glideSteps': np.array([400]),
    'dipoleCenters': np.array([0.0,0.0,0.0]),
    #'slipSystemIDs' : np.array([0,1]),
    #'exitFaceIDs': np.array([1,1]), # 1 edge, 2 screw
    #'nodesPerLine': np.array([10,10]),
    #'dipoleHeights': np.array([2000,2000]),
    #'glideSteps': np.array([400,405]),
    #'dipoleCenters': np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]]),
}

# initial config
Polycrystal_parameters= {
    'absoluteTemperature': 1,
    'grain1globalX1': np.array([0, 1, 1]),
    'grain1globalX3': np.array([-1, 1, -1]),
    'boxEdges': np.array([[0,1,1], [-2,-1,1], [-1,1,-1]]),
    #'boxScaling': np.array([1020, 500, 3800]),
    'boxScaling': np.array([500, 1020, 3800]),
    'X0': np.array([0., 0., 0.]),
    'periodicFaceIDs': np.array([0,1,2,3,4,5]),
}
