import os
import sys
import string
import subprocess
import numpy as np

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
    #sf_noise_fname = str(config.files_to_copy_from_lib['noise_file_md_sf']).split('/')[-1]
    an_noise_fname = str(config.files_to_copy_from_lib['noise_file_an_ss']).split('/')[-1]
    #modify_dict_parameters(config.MDSolidSolution_parameters, ss_noise_fname)
    #modify_dict_parameters(config.MDStackingFault_parameters, sf_noise_fname)
    modify_dict_parameters(config.AnalyticalSolidSolution_parameters, an_noise_fname)

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


def main() -> int:
    ################## Test Parameters #######################
    glidePlasticStrain = 1e-10 # abstract value found through trial and error
    #noise_seed_to_test = [3, 4]
    #noise_seed_to_test = [7,8]
    noise_seed_to_test = [1,2]
    node_num_to_test = [ 300 ]

    search_start_stress = 1 # in MPa
    #search_start_stress = 500 # in MPa

    # directory to store all the generated data
    dataStorageDir = Path("generatedData")
    os.makedirs(dataStorageDir, exist_ok=True)

    # remove old data
    #if os.path.exists(dataStorageDir):
    #    shutil.rmtree(dataStorageDir)

    zero_str_tensor = np.array([0, 0, 0, 0, 0, 0])
    #min_steps = 201 # total minimization steps
    min_steps = int(config.DD_parameters['outputFrequency']) + 1
    strSign = -1 # +1 or -1
    for seed in noise_seed_to_test:
        for node in node_num_to_test:
            stress = search_start_stress
            searchStepSize = 100 # initial search step size
            searchResolution = 10
            dislocMovingSlowly = False
            dislocationMoving = False

            while not dislocMovingSlowly:
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

                # change the node number on dislocation line
                modify_parameter('periodicDipoleIndividual.txt', 'nodesPerLine', node)

                # change the noise seed number
                #modify_parameter('MDSolidSolution.txt', 'seed', seed)
                #modify_parameter('MDStackingFault.txt', 'seed', seed)
                modify_parameter('AnalyticalSolidSolutionNoise.txt', 'seed', seed)

                # boundary condition control file
                elasticDeformationFile='ElasticDeformation.txt'

                # minimize the system before the run
                modify_parameter('DD.txt', 'Nsteps', min_steps)
                modify_parameter(elasticDeformationFile, 'ExternalStress0', zero_str_tensor)

                subprocess.run(["../../build/tools/MicrostructureGenerator/microstructureGenerator", "./"], capture_output=False, text=True)
                exit()
                subprocess.run(["../../build/tools/DDomp/DDomp", "./"], capture_output=False, text=False)

                #except subprocess.CalledProcessError as e:
                #    print(f"DDomp failed: {e}")

                # change the target run back to original
                original_target_run = config.DD_parameters['Nsteps']
                modify_parameter('DD.txt', 'Nsteps', original_target_run)

                # extract material info
                matFile = str(config.files_to_copy_from_lib['material_file']).split('/')[-1]
                b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
                mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')
                convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))

                # apply stress and run simulation
                # Voigt order is 11,22,33,12,23,13
                #s12 = stress*convertMPaToMu
                s13 = strSign*stress*convertMPaToMu
                #stressTensor = np.array([0.0, 0.0, 0.0, s12, 0, 0])
                stressTensor = np.array([0.0, 0.0, 0.0, 0, 0, s13])
                modify_parameter('ElasticDeformation.txt', 'ExternalStress0', stressTensor)
                subprocess.run(["../../build/tools/DDomp/DDomp", "./"], capture_output=False, text=True)

                # determine whether crss or not
                f,fLabels=readFfile('./F')
                runIDs=getFarray(f,fLabels,'runID').astype(int)
                #betaP_12=getFarray(f,fLabels,'betaP_12').astype(float)
                betaP_13=getFarray(f,fLabels,'betaP_13').astype(float)

                # store generated data
                #target_dir_temp = dataStorageDir/f"seed{str(seed).strip()}"/f"node{str(node[0])}"
                target_dir_temp = dataStorageDir/f"seed{str(seed).strip()}"/f"node{str(node)}"
                #targetDir = Path(dataStorageDir)/f"seed{str(seed).strip()}"/f"node{str(node[0])}"/f"{str(stress).strip('-')}MPa"
                targetDir = target_dir_temp/f"{str(stress).strip('-')}MPa"
                os.makedirs(targetDir, exist_ok=True)
                for f in folders:
                    shutil.move(f, targetDir)

                # filter the duplicate data (minimization creates duplicates on runID)
                # Find duplicates and their indices
                unique_values, counts = np.unique(runIDs, return_counts=True)
                duplicates = unique_values[counts > 1]
                # remove the first duplicate data entry
                dup_idx = np.where(runIDs==duplicates)[0][0]
                runIDs = np.delete(runIDs, dup_idx)
                #betaP_12 = np.delete(betaP_12, dup_idx)
                betaP_13 = np.delete(betaP_13, dup_idx)

                # calculate the absolute plastic distortion
                #abs_betaP = np.abs(betaP_12 - betaP_12[0])
                abs_betaP = np.abs(betaP_13 - betaP_13[0])
                #print("Duplicate values:", duplicates)

                dydx = np.gradient(abs_betaP, runIDs)
                dydx_mean = np.mean(dydx[-3:])

                #with open(Path(targetDir)/'search_log.txt', 'a') as log:
                #f = open(Path(dataStorageDir)/f"seed{str(seed).strip()}"/"search_log.txt", "a")
                f = open(target_dir_temp/"search_log.txt", "a")
                f.write(f"{stress}MPa, seed{seed}\n")

                if dydx_mean > glidePlasticStrain:
                    # while loop stop condition
                    #if searchStepSize < 10 and dislocationMoving:
                    if searchStepSize < searchResolution:
                        dislocMovingSlowly = True
                        f.write(f"found CRSS, {stress-searchStepSize}MPa, seed{seed}")
                    else:
                        dislocationMoving = True
                        searchStepSize -= searchStepSize//2
                        stress -= searchStepSize
                else:
                    if searchStepSize < searchResolution and dislocationMoving:
                        dislocMovingSlowly = True
                        f.write(f"found CRSS, {stress-searchStepSize}MPa, seed{seed}")
                    elif dislocationMoving:
                        searchStepSize -= searchStepSize//2
                        stress += searchStepSize
                    else:
                        stress += searchStepSize

                # stop the run if the stress becomes less than zero
                if stress < 0:
                    break

                f.close()

    return 0

if __name__ == "__main__":
    main()
