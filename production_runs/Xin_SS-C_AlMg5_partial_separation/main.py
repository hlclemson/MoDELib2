import re
import os
import sys
import json
import string
import subprocess
import numpy as np

# import matplotlib.pyplot as plt
sys.path.append("../../build/tools/pyMoDELib/")
import pyMoDELib

sys.path.append("../../python")
from pathlib import Path
from modlibUtils import *

# read initial configuration file
config_fname = "initial_config.json"
with open(config_fname, "r") as f:
    config = json.load(f)

files_to_copy_from_lib = {
    key: src_path for key, src_path in config["files_to_copy_from_lib"].items()
}


def _set_noise_types(fileName: Path, config_params: dict) -> None:
    """
    Compile the glide-plane noise list and assign it in the material file
    """
    # compile the glideplane noise list based on the config
    noise_list = []
    for k,v in config_params.items():
        match = re.fullmatch(r'glidePlaneNoise\d*', k)
        if match:
            # is it a list? else, make it as a list
            for entry in (v if isinstance(v, list) else [v]):
                tmp_string = f"glidePlaneNoise ={entry};"
                noise_list.append(tmp_string)

    fileName = Path("inputFiles") / fileName
    with open(fileName) as f:
        lines = f.readlines()
    pattern = re.compile(r'glidePlaneNoise')
    # ensure each list entry ends with \n
    noise_list = [s + '\n' for s in noise_list]
    for idx, line in enumerate(lines):
        if pattern.search(line):
            # replace the existing glidePlaneNoise string with the new one
            # swap 1 line → N lines
            lines[idx:idx+1] = noise_list
            break
    else: # no match at all, return
        return
    # replace the existing file
    with open(fileName, 'w') as f:
        f.writelines(lines)


def modify_parameter(target_config_fname: str, param: str, content) -> None:
    content_array = np.asarray(content)
    squeezed = np.squeeze(content_array)
    config_path = Path("inputFiles") / target_config_fname
    if squeezed.ndim > 1:  # matrix case
        setInputMatrix(config_path, param, squeezed)
    elif squeezed.ndim == 1:  # vector case
        setInputVector(config_path, param, squeezed, newCom="")
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


def copy_run_files() -> None:
    """Copy all necessary files from Library to inputFiles"""
    dest_dir = Path("inputFiles")
    dest_dir.mkdir(exist_ok=True)

    for key, src in files_to_copy_from_lib.items():
        # loop through all the lists if it's list object
        for path in (src if isinstance(src, list) else [src]):
            # if it's empty, continue
            if not path:
                continue
            path = Path(path)
            shutil.copy2(path.resolve(), dest_dir / path.name)


def create_polycrystal(mat_fname: str, **params) -> None:
    # Set defaults while allowing overrides
    settings = {
        "absoluteTemperature": 1,
        "grain1globalX1": [0, 1, 1],
        "grain1globalX3": [-1, 1, -1],
        "boxEdges": [[0, 1, 1], [-2, -1, 1], [-1, 1, -1]],
        "boxScaling": [100, 100, 100],
        "X0": [0, 0, 0],
        "periodicFaceIDs": [-1],
        **params,  # Overrides defaults with passed values
    }
    pf = PolyCrystalFile(mat_fname)
    mesh_fname = Path(files_to_copy_from_lib["mesh"]).name
    pf.meshFile = mesh_fname
    for param, value in settings.items():
        setattr(pf, param, np.array(value))
    pf.write("inputFiles")


def getVectorInFile(fname: str, param: str, dtype=int):
    return list(map(dtype, getStringInFile(fname, param).split()))


def _remove_dup_F() -> None:
    """
    remove duplicate entry introudced by minimization in F_0.txt file
    """
    # output sanitization, remove duplicate entry on F
    f, fLabels = readFfile("./F")
    runIDs = getFarray(f, fLabels, "runID").astype(np.int32)
    # unique_values, counts = np.unique(runIDs, return_counts=True)
    _, idx, counts = np.unique(runIDs, return_inverse=True, return_counts=True)
    # construct boolean array that is true for every duplicate row
    duplicate_mask = counts[idx] > 1
    # get idxes of unique values
    first_of_each = np.unique(runIDs, return_index=True)[1]
    # set the duplicate value to value
    duplicate_mask[first_of_each] = False
    # find duplicate entry row idx
    duplicate_row_idx = np.squeeze(np.where(duplicate_mask))

    f_file = Path("./F/F_0.txt")
    f_file_copy = Path("./F/F_0_copy.txt")
    with f_file.open() as fin, f_file_copy.open("w") as fout:
        for idx, line in enumerate(fin):
            if idx == duplicate_row_idx:
                continue
            fout.write(line)
    # replace the old one with new file
    f_file.unlink()
    f_file_copy.rename(f_file)


def main() -> int:
    ####### Binaries #########
    MICROSTUCTEXE = "../../build/tools/MicrostructureGenerator/microstructureGenerator"
    DDOMP = "../../build/tools/DDomp/DDomp"

    # directory to store all the generated data
    dataStorageDir = Path("generatedData")
    os.makedirs(dataStorageDir, exist_ok=True)

    ####### Test Parameters #######
    glidePlasticStrain = 1e-10  # abstract value found through trial and error
    #noise_seed_to_test = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # test 100 pair
    noise_seed_to_test = range(1,101)
    ref_glide_pos = 400
    glide_inc_step = 4
    glide_steps_to_test = [ [ref_glide_pos, ref_glide_pos+(glide_inc_step*i)] for i in range(1,101) ]

    stress = 0  # in MPa
    crss_search_resolution = 10
    zero_str_tensor = np.array([0, 0, 0, 0, 0, 0])
    min_steps = int(config["DD_parameters"]["outputFrequency"]) + 1
    strSign = -1  # +1 or -1
    for seed, glide_steps in zip(noise_seed_to_test, glide_steps_to_test):
        # Preparing input files
        folders = ["evl", "F", "inputFiles"]
        for x in folders:
            # remove existing data
            if os.path.exists(x):
                shutil.rmtree(x)
            # create necessary folder structure for the simulation
            os.makedirs(x)

        # copy simulation files to inputFiles
        copy_run_files()

        # DD.txt setup
        dd_fname = Path(files_to_copy_from_lib["dd_file"]).name
        modify_dict_parameters(config["DD_parameters"], dd_fname)

        # material file setup
        mat_fname = Path(files_to_copy_from_lib["material_file"]).name
        modify_dict_parameters(config["Material_parameters"], mat_fname)

        # check if the configuration includes noise
        has_noise = any(re.fullmatch(r'glidePlaneNoise\d*', k) for k in config["Material_parameters"].keys())
        # compile the glideplane noise list and
        # assign it in the material file
        if has_noise:
            _set_noise_types(mat_fname, config["Material_parameters"])

        # microstructure file setup
        micro_fname = Path(files_to_copy_from_lib["microstructure"]).name
        modify_dict_parameters(config["Microstructure_parameters"], micro_fname)
        with open("inputFiles/initialMicrostructure.txt", "w") as f:
            f.write(f"microstructureFile={micro_fname};\n")

        # change glide step for partials
        modify_parameter(micro_fname, "glideSteps", glide_steps)

        # noise_file.txt setup
        # correletion_fname = Path( files_to_copy_from_lib["correlation_file_md_sf"]).name
        noise_path = files_to_copy_from_lib["noise_file"]
        for noise in (noise_path if isinstance(noise_path, list) else [noise_path]):
            noise_fname = Path(noise).name
            noise_setting = config["noise_settings"][noise_fname]
            modify_dict_parameters(noise_setting, noise_fname)
            if "MDStackingFault.txt" == noise_fname:
                correletion_fname = Path(
                    files_to_copy_from_lib["correlation_file_md_sf"]
                ).name
                modify_parameter(noise_fname, "correlationFile", correletion_fname)
            if "MDSolidSolution.txt" == noise_fname:
                correletion_fname_xz = Path(
                    files_to_copy_from_lib["correlation_file_md_ss_xz"]
                ).name
                correletion_fname_yz = Path(
                    files_to_copy_from_lib["correlation_file_md_ss_yz"]
                ).name
                modify_parameter(
                    noise_fname, "correlationFile_xz", correletion_fname_xz
                )
                modify_parameter(
                    noise_fname, "correlationFile_yz", correletion_fname_yz
                )
            # change the noise seed number
            modify_parameter(noise_fname, "seed", seed)

        # polycrystal setting
        poly_settings = config["polycrystal"]
        create_polycrystal(
            mat_fname,
            absoluteTemperature=poly_settings["absoluteTemperature"],
            grain1globalX1=poly_settings["grain1globalX1"],
            grain1globalX3=poly_settings["grain1globalX3"],
            boxEdges=poly_settings["boxEdges"],
            boxScaling=poly_settings["boxScaling"],
            X0=poly_settings["X0"],
            periodicFaceIDs=poly_settings["periodicFaceIDs"],
        )

        # boundary condition control file
        elasticDeformationFile = "ElasticDeformation.txt"

        # minimize the system before the run
        modify_parameter("DD.txt", "Nsteps", min_steps)
        modify_parameter(elasticDeformationFile, "ExternalStress0", zero_str_tensor)
        subprocess.run(
            [MICROSTUCTEXE, "./"],
            capture_output=False,
            text=True,
        )
        subprocess.run(
            [DDOMP, "./"],
            capture_output=False,
            text=False,
        )

        # change the target run back to original
        original_target_run = config["DD_parameters"]["Nsteps"]
        modify_parameter("DD.txt", "Nsteps", original_target_run)

        # extract material info
        # b_SI = getValueInFile(f"inputFiles/{matFile}", "b_SI")
        mu0_SI = getValueInFile(f"inputFiles/{mat_fname}", "mu0_SI")
        convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))

        # apply stress and run simulation
        # Voigt order is 11,22,33,12,23,13
        # s12 = strSign * stress * convertMPaToMu
        s13 = strSign * stress * convertMPaToMu
        stressTensor = np.array([0, 0, 0, 0, 0, s13])
        modify_parameter("ElasticDeformation.txt", "ExternalStress0", stressTensor)
        subprocess.run(
            [DDOMP, "./"],
            capture_output=False,
            text=False,
        )

        # output sanitization, remove duplicate entry on F
        _remove_dup_F()

        # build the storage directory path
        target_dir = (
            dataStorageDir / f"seed{seed}_steps{''.join(map(str, glide_steps))}"
        )
        # remove old data and create a new empty dir
        shutil.rmtree(target_dir, ignore_errors=True)
        target_dir.mkdir(parents=True)
        # move the generated folders in
        for f in folders:
            shutil.move(f, target_dir)

    return 0


if __name__ == "__main__":
    main()
