# standard lib
from os import path
import re
import sys
import json
import shutil
import tarfile
import tempfile

# 3rd party lib
import numpy as np
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt

# ----- MoDELib / utils paths -----
sys.path.append("../../python")
#sys.path.append("../../python/lib")
from visUtils import *
from modlibUtils import *

sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib

def readValFromMaterialFile(matDir: str, alloy: str, var: str) -> float:
    with open(matDir / alloy, "r") as mFile:
        for line in mFile:
            # strip down comments from the data
            line = line.split(";")[0]
            if line.startswith(f"{var}"):
                # Split the line by '=' and take the second part, then remove leading/trailing whitespace
                value = float(line.split("=")[1].strip())
                break
    return value


def main():
    try:
        config_fname = sys.argv[1]
    except:
        exit("usage) generate_separation_data.py config.json")

    box_axis_idx = {"x": 0, "y": 1, "z": 2}

    with open(config_fname, "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])
    glide_axis = box_axis_idx[config["glide_axis"]]
    line_tangent_axis = box_axis_idx[config["line_tangent_axis"]]

    #data_paths = [p for pat in ("*MPa.tar.gz", "MPa") for p in simulationDir.rglob(pat)]
    # sort by seed number
    #data_paths = sorted(
    #    data_paths, key=lambda p: int(re.search(r"seed(\d+)", str(p)).group(1))
    #)
    data_path = simulationDir

    # sort each paths based on seed
    #paths_per_seed = defaultdict(list)
    #for i in data_paths:
    #    seed_num = re.search(r"seed(\d+)", str(i.parent)).group(1)
    #    paths_per_seed[int(seed_num)].append(i)

    # sort stress value
    # extract stress value
    stress = int(re.search(r"(\d+)MPa", str(data_path.name)).group(1))
    # open data
    src = data_path 
    tmp = None  # will hold TemporaryDirectory handle
    work_dir = None  # will point to the folder we actually use
    if src.is_dir():  # in uncompressed format
        work_dir = src
    elif src.suffixes == [".tar", ".gz"]:  # if .tar.gz format
        tmp = tempfile.TemporaryDirectory()
        with tarfile.open(src, "r:gz") as tf:
            tf.extractall(tmp.name, filter="data")
        work_dir = Path(tmp.name) / src.name.removesuffix(".tar.gz")
    else:
        raise FileNotFoundError("neither directory nor .tar.gz found")

    ddBase = pyMoDELib.DislocationDynamicsBase(str(work_dir))
    ddio = pyMoDELib.DDconfigIO(str(work_dir / "evl"))
    # Box geometry
    xMin = np.array(ddBase.mesh.xMin(), dtype=float)
    xMax = np.array(ddBase.mesh.xMax(), dtype=float)
    xCenter = np.array(ddBase.mesh.xCenter(), dtype=float)

    # read labels
    with open(work_dir/"F/F_labels.txt", "r") as label:
        fLabels = label.read()
    # remove empty element and store it as a list
    fLabels = np.array([x for x in fLabels.split("\n") if x], dtype=str)
    # open F file
    fData = np.loadtxt(work_dir/"F/F_0.txt")
    # find the last step
    runIDs = getFarray(fData,fLabels,'runID').astype(int)
    #last_runID=getFarray(fData,fLabels,'runID').astype(int)[-1]
    stable_runID = runIDs[len(runIDs)//3]
    # read evl file
    ddio.readTxt(stable_runID)

    # First pass: find active planes over ALL EVLs
    active_plane_keys_global = set()
    defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
    defectiveCrystal.initializeConfiguration(ddio)
    DN = defectiveCrystal.dislocationNetwork()
    planes = ddBase.glidePlanes()

    # find the number of glide planes
    plane_keys = [stable_plane_key(gp) for gp in planes]  # Not used
    print(f" ... Total glide planes: {len(plane_keys)} ...")

    node_pos_per_loop = defaultdict(list)
    # collect network IDs for each loop number
    for key in DN.networkNodes().keys():
        node = DN.networkNodes().getRef(key)
        node_pos = np.array(node.position(), dtype=np.float64)
        node_pos_per_loop[tuple(node.loopIDs())].append(node_pos)
    # filter the loopIDs that has more than 1 element (connecting nodes)
    node_pos_per_loop = {k[0]: np.array(v) for k, v in node_pos_per_loop.items() if v and len(k)==1}
    print(node_pos_per_loop[5].shape)
    test_pos = np.asarray(node_pos_per_loop[5][0], dtype=np.float64).reshape(3,1)
    print(test_pos.shape)# Should be (3, 1)
    print(test_pos.dtype)# Should be float64
    print(test_pos.flags)# Check if C_CONTIGUOUS is True

    print(planes)
    gppy = pyMoDELib.GlidePlane
    test_plane = planes[0]
    test_loopid = 5
    localPoses = np.asarray([ gppy.localPosition(test_plane, x) for x in node_pos_per_loop[test_loopid] ])
    print(localPoses)
    plt.scatter(localPoses[:,0], localPoses[:,1])
    plt.show()
    #print(dir(gppy.meshSegments))
    #print(gppy.meshSegments(planes[0]))
    #print(gppy.key_tuple(planes[0]))
    #print(gppy.localPosition(planes[0], test_pos))
    exit()

    # find loop pair
    #glide_plane_loops = defaultdict(list)
    #for loop_id, v in node_pos_per_loop.items():
    #    unique_z_pos = np.unique(v[:,2]).item()
    #    glide_plane_loops[unique_z_pos].append(loop_id)
    #partial_pairs = list(glide_plane_loops.values())

    ## calc separation distance
    #separation_bt_nodes_all = []
    #for partial_pair in partial_pairs:
    #    glide_pos_per_id = []
    #    for partial_loop_id in partial_pair:
    #        node_pos = node_pos_per_loop[partial_loop_id]
    #        line_glide_pos = node_pos[:, glide_axis]
    #        glide_pos_per_id.append(line_glide_pos)
    #    # skip the partial pair if node nums are not the same
    #    if glide_pos_per_id[0].shape != glide_pos_per_id[1].shape:
    #        continue
    #    glide_pos_per_id = np.array(glide_pos_per_id)
    #    separation_bt_nodes_all.append(np.abs(glide_pos_per_id[1,:]-glide_pos_per_id[0,:]))
    ## reshape to 1D
    #separation_bt_nodes_all = [item for row in separation_bt_nodes_all for item in row]
    #separation_bt_nodes_all = np.array(separation_bt_nodes_all)

    ##if len(separation_bt_nodes_all)>1:
    ##    separation_bt_nodes_all = np.reshape(separation_bt_nodes_all, -1)
    #separation_out_fname = data_path.parent/f"{re.search(r"(\d+)MPa", str(data_path.name)).group(0)}.txt"
    #np.savetxt(separation_out_fname, separation_bt_nodes_all)

    ## automatically deletes the tmp tree if there is any
    #if tmp:
    #    tmp.cleanup()

if __name__ == "__main__":
    main()
