# standard lib
#from os import path, sep
import os
import re
import sys
import json
import shutil
import tarfile
import tempfile

# 3rd party lib
import numpy as np
#import pandas as pd
#from scipy import stats
from pathlib import Path
from collections import defaultdict
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

# ----- MoDELib / utils paths -----
sys.path.append("../../python")
#sys.path.append("../../python/lib")
from visUtils import *
from modlibUtils import *

sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib

#def readValFromMaterialFile(matDir: str, alloy: str, var: str) -> float:
#    with open(matDir / alloy, "r") as mFile:
#        for line in mFile:
#            # strip down comments from the data
#            line = line.split(";")[0]
#            if line.startswith(f"{var}"):
#                # Split the line by '=' and take the second part, then remove leading/trailing whitespace
#                value = float(line.split("=")[1].strip())
#                break
#    return value


def main():
    try:
        config_fname = sys.argv[1]
    except:
        exit("usage) generate_separation_data.py config.json")

    box_axis_idx = {"x": 0, "y": 1, "z": 2}
    with open(config_fname, "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])
    out_file =  Path(config["out_file"])
    glide_axis = box_axis_idx[config["glide_axis"]]
    line_tangent_axis = box_axis_idx[config["line_tangent_axis"]]
    os.makedirs(out_file.parent, exist_ok=True)

    data_paths = [p for pat in ("seed*_steps*.tar.gz", "MPa") for p in simulationDir.rglob(pat)]
    # sort by seed number
    data_paths = sorted(
        data_paths, key=lambda p: int(re.search(r"seed(\d+)", str(p)).group(1))
    )

    total_separation_dist_data = []
    for data_path in data_paths:
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
        ddgp = pyMoDELib.GlidePlane
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
        #stable_runID = runIDs[len(runIDs)//3]
        stable_runID = runIDs[-1]
        # read evl file
        ddio.readTxt(stable_runID)

        # create defective crystal C object
        defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
        defectiveCrystal.initializeConfiguration(ddio)
        DN = defectiveCrystal.dislocationNetwork()
        planes = ddBase.glidePlanes()

        node_pos = {}
        for key in DN.networkNodes().keys():
            node = DN.networkNodes().getRef(key)
            node_pos[node.networkID()] = np.asarray(node.position(), dtype=np.float64)

        plane_geometry = {ddgp.key_tuple(plane): (ddgp.unitNormal(plane), ddgp.planeOrigin(plane)) for plane in planes}
        plane_to_segments = defaultdict(list)
        tol = 1e-4
        for lkey in DN.networkLinks().keys():
            link = DN.networkLinks().getRef(lkey)
            # filter segment that has zero Burgers vector
            if hasattr(link, "hasZeroBurgers") and link.hasZeroBurgers():
                continue
            a = node_pos[link.source.networkID()]
            b = node_pos[link.sink.networkID()]
            # convert dict -> list, return last element
            loop_id = list(link.loopIDs()).pop()
            best_plane_key = None
            best_dist = float("inf")
            for pk, (n_hat, p0) in plane_geometry.items():
                d = segment_plane_distance(a, b, n_hat, p0)
                if d < best_dist:
                    best_dist = d
                    best_plane_key = pk
            if best_plane_key is not None and best_dist < tol:
                #plane_to_segments[best_plane_key].append((a, b))
                plane_to_segments[(best_plane_key, loop_id)].append((a, b))

        xyz_col_num = 3
        for k, node_pos in plane_to_segments.items():
            # combine source and sink pos dimension
            node_pos = np.reshape(node_pos, (-1, xyz_col_num))
            # remove duplicate node pos (combine overlapping source and sink pos)
            node_pos = np.unique(node_pos, axis=0)
            plane_to_segments[k] = node_pos

        # find partial loop pair (in case of partial dislocation)
        glide_plane_loops = defaultdict(list)
        for k, node_pos in plane_to_segments.items():
            plane_key, loop_id = k
            z_pos = node_pos[:,2]
            unique_z_pos = np.unique(z_pos).item()
            glide_plane_loops[unique_z_pos].append(loop_id)
        partial_pairs = list(glide_plane_loops.values())

        separation_dists = []
        for partial_pair in partial_pairs:
            #pair_separation_dist = 0
            pair_separation_dist = []
            for loop_id_pair in partial_pair:
                for k, node_pos in plane_to_segments.items():
                    plane_key, loop_id = k
                    if loop_id != loop_id_pair:
                        continue
                    for plane in planes:
                        tmp_plane_key = ddgp.key_tuple(plane)
                        if tmp_plane_key == plane_key:
                            # project global pos (3D) to glide plane pos (2D)
                            local_node_pos = np.asarray([ddgp.localPosition(plane, x) for x in node_pos])
                            avg_pos = np.mean(local_node_pos[:, line_tangent_axis])
                            pair_separation_dist.append(avg_pos)
            separation_dist = np.abs(np.diff(pair_separation_dist))
            #unreasonable_separation = 10
            #if separation_dist > unreasonable_separation:
            #    continue
            total_separation_dist_data.append(separation_dist)
        seed = int(re.search(r"seed(\d+)", str(data_path)).group(1))
        sep_match = re.search(r"steps(\d+)", str(data_path)).group(1)
        #split the string into two num
        midpoint = len(sep_match)//2
        g_step1 = int(sep_match[:midpoint])
        g_step2 = int(sep_match[midpoint:])
        init_sep = np.abs(np.diff([g_step1, g_step2]))
        #match = re.search(r"seed(\d+)", str(data_path))
        #init_sep = re.search("steps*", data_path)
        json_data = {"seed":seed , "init_sep_dist": init_sep.item(), "relaxed_sep_dist": separation_dist.item()}
        with open(out_file, "a") as json_f:
            json_f.write(f"{json.dumps(json_data)}\n")

    #plt.hist(total_separation_dist_data)
    #fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=200)
    #n, bins, patches = ax.hist(total_separation_dist_data, bins=30, density=True, alpha=0.7, color='skyblue', edgecolor='black', linewidth=0.5)
    ## Calculate Gaussian fit parameters
    #mu, sigma = np.mean(total_separation_dist_data), np.std(total_separation_dist_data)
    ## Generate x values for the Gaussian curve
    #x = np.linspace(total_separation_dist_data.min(), total_separation_dist_data.max(), 1000)
    ## Calculate Gaussian probability density function
    #gaussian_fit = stats.norm.pdf(x, mu, sigma)
    ## Plot the Gaussian fit
    ##ax.plot(x, gaussian_fit, 'r-', linewidth=2, label=f'Gaussian fit (μ={mu:.2f}, σ={sigma:.2f})')

    ## Add labels and legend
    #ax.set_xlabel('Value')
    #ax.set_ylabel('Density')
    #ax.set_title('Histogram with Gaussian Fit')
    #ax.legend()
    #plt.tight_layout()
    #os.makedirs(out_file.parent, exist_ok=True)
    #fig.savefig(out_file, transparent=True)
    ##plt.show()

if __name__ == "__main__":
    main()
