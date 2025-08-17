import numpy as np
from pathlib import Path

# Define all file templates, these files get copied to the inputFiles directory
files_to_copy_from_lib = {
    'dd_file': Path('../../Library/DislocationDynamics/DD.txt'),
    'noise_file_md_sf': Path('../../Library/GlidePlaneNoise/MDStackingFault.txt'),
    "noise_file_an_ss": Path("../../Library/GlidePlaneNoise/AnalyticalSolidSolutionNoise.txt"),
    'material_file': Path('../../Library/Materials/AlMg5.txt'),
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
    'alphaLineTension': '0.422',  # dimensionless line tension scale factor
    'remeshFrequency': '0',  # node redistribution
    'quadPerLength': '6',
    'coreSize': '0',
    'Lmin': '5',  # min segment length (Burgers vector units)
    'Lmax': '20',  # max segment length (Burgers vector units)
    'outputFrequency': '500',  # output frequency
    'outputQuadraturePoints': '1',  # output quadrature data
    'computeElasticEnergyPerLength': '1',  # output energy data
    'glideSolverType': 'Galerkin',  # solver type ('Galerkin' or 'none')
    'climbSolverType': 'none', # type of clim solver, or none 
    'Nsteps': '1001', # number of simulation steps
    'nodalVelocityConstraints': np.array([0,1,0])
}

AnalyticalSolidSolution_parameters = {
    # Scalar parameters (use setInputVariable)
    "type": "AnalyticalSolidSolutionNoise",
    "tag": 1,
    "seed": 1,
    "spreadLstress_SI": 1e-10, # eigenstrain spreading distance (a)
    "a_cai_SI": 0, # non-singular param (a_n), edge
    "gridSize": np.array([2560, 2560, 64]),
    "gridSpacing_SI": np.array([1.0e-10, 1.0e-10, 1e-10]),
    "MSSS_SI": 0.16e18
}

# initial config
Material_parameters = {
    'enabledSlipSystems': 'full',
    'glidePlaneNoise': 'AnalyticalSolidSolutionNoise.txt',
    'atomsPerUnitCell': '1',
    'dislocationMobilityType': 'default'
}

# initial config
Microstructure_parameters = {
    'slipSystemIDs': np.array([0]),
    'exitFaceIDs': np.array([1]),
    'nodesPerLine': np.array([10]),
    'dipoleHeights': np.array([2000]),
    'glideSteps': np.array([400]),
    'dipoleCenters': np.array([0.0,0.0,0.0]),
}

# initial config
Polycrystal_parameters= {
    'absoluteTemperature': 1,
    'grain1globalX1': np.array([0, 1, 1]),
    'grain1globalX3': np.array([-1, 1, -1]),
    'boxEdges': np.array([[0,1,1], [-2,-1,1], [-1,1,-1]]),
    'boxScaling': np.array([1020, 500, 3800]),
    'X0': np.array([0., 0., 0.]),
    'periodicFaceIDs': np.array([0,1,2,3,4,5]),
}
