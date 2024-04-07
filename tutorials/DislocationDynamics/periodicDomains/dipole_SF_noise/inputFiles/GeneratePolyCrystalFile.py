import sys
sys.path.append("../../../../../python/")
from modlibUtils import *

# Make a local copy of simulation parameters file and modify that copy if necessary
DDfile='DD.txt'
# shutil.copy2('../../'+DDfile, '.') 
# setInputVariable(DDfi0le,'timeIntegrationMethod','0') # 0 for constant time integration
# setInputVariable(DDfile,'timeStep','20') # 0 for constant time integration
# setInputVariable(DDfile,'dxMax','1') # max time step if timeIntegrationMethod=1
# setInputVariable(DDfile,'use_velocityFilter','0') # don't filter velocity if noise is enabled
# setInputVariable(DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
# setInputVariable(DDfile,'alphaLineTension','0.1') # dimensionless scale factor in for line tension forces
# setInputVariable(DDfile,'Lmin','5')  # min segment length (in Burgers vector units)
# setInputVariable(DDfile,'Lmax','20')  # max segment length (in Burgers vector units)
# setInputVariable(DDfile,'outputFrequency','1')  # output frequency
# setInputVariable(DDfile,'outputQuadraturePoints','1')  # output quadrature data

# make a local copy of material file, and modify that copy if necessary
materialTemplate='Ni.txt';
shutil.copy2('../../../MaterialsLibrary/'+materialTemplate, '.') # target filename is /dst/dir/file.ext
setInputVariable(materialTemplate,'enabledSlipSystems','Shockley')
setInputVariable(materialTemplate,'ISF_SI','-27.0e-3')
setInputVariable(materialTemplate,'USF_SI','312e-3')
setInputVariable(materialTemplate,'MSF_SI','789e-3')


# Create polycrystal.txt using local material file
pf=PolyCrystalFile(materialTemplate);
pf.absoluteTemperature=300;
pf.meshFile='../../../MeshLibrary/unitCube.msh'
pf.grain1globalX1=np.array([0,1,1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
pf.grain1globalX3=np.array([-1,1,-1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.boxEdges=np.array([[0,1,1],[2,1,-1],[-1,1,-1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([200,200,200]) # must be a vector of integers
pf.X0=np.array([0.5,0.5,0.5]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([0,1,2,3,4,5])
pf.solidSolutionNoiseMode=0
pf.stackingFaultNoiseMode=1; # 0=no noise
pf.gridSize=np.array([27,30])
pf.gridSpacing_SI=np.array([2.515e-10,2.515e-10])

pf.write()


# Edit dipoles partials dislocations
microstructureTemplate='periodicDipole.txt';
shutil.copy2('../../../MicrostructureLibrary/'+microstructureTemplate, '.') # target filename is /dst/dir/file.ext
setInputVector(microstructureTemplate,'periodicDipoleSlipSystemIDs',np.array([0,1]),'slip system IDs for each dipole')
setInputVector(microstructureTemplate,'periodicDipoleExitFaceIDs',np.array([1,1]),'1 is for edge, 0 for screw')
setInputMatrix(microstructureTemplate,'periodicDipolePoints',np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]]))
setInputVector(microstructureTemplate,'periodicDipoleNodes',np.array([10,10]),'number of extra nodes on each dipole')
setInputVector(microstructureTemplate,'periodicDipoleHeights',np.array([50,50]),'height of each dipole, in number of planes')
setInputVector(microstructureTemplate,'periodicDipoleGlideSteps',np.array([10.0,30.0]),'step of each dipole in the glide plane')

# Edit stress
uniformExternalLoadControllerFile='uniformExternalLoadController.txt'
loadAxis=np.array([0,1,1])@pf.C2G.transpose()
loadAxis=loadAxis/np.linalg.norm(loadAxis)
stress=0.01
ExternalStress0=stress*np.outer(loadAxis,loadAxis)
setInputMatrix(uniformExternalLoadControllerFile,'ExternalStress0',ExternalStress0)

