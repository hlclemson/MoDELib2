#import pudb; pu.db

# standard lib
import sys
import json
import math
import shutil
import tarfile
import pathlib
import tempfile

# 3rd party lib
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
from collections import defaultdict
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# Configure global plot settings (applies to all figures)
rcParams.update(
    {
        "figure.dpi": 200,
        #"figure.autolayout": True,  # Prevent label clipping
        #"axes.grid": True,
        #"grid.alpha": 0.6,
        "text.usetex": False,
        "font.size": 10,  # Default font size for text
        "mathtext.fontset": "stix",  # Use STIX font for math text
        "font.family": "serif",  # Use serif font (matches LaTeX default)
    }
)

# ----- MoDELib / utils paths -----
sys.path.append("../../python")
#sys.path.append("../../python/lib")
from visUtils import *
from modlibUtils import *

sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib


def main():
    box_axis_idx = {"x": 0, "y": 1, "z": 2}

    with open("config.json", "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])
    evl_step = config["evl_step"]
    glide_axis = box_axis_idx[config["glide_axis"]]
    line_tangent_axis = box_axis_idx[config["line_tangent_axis"]]
    num_samples = config["number_of_samples"]


    # test 100 pair
    noise_seed_to_test = range(1, num_samples+1)
    glide_inc_step = 4
    ref_glide_pos = 400
    glide_steps_to_test = [
        [ref_glide_pos, ref_glide_pos + (glide_inc_step * i)] for i in range(1, num_samples+1)
    ]


    separation_dists_total = defaultdict(list)
    for seed, glide_steps in zip(noise_seed_to_test, glide_steps_to_test):
        work_dir = simulationDir / f"seed{seed}_steps{''.join(map(str, glide_steps))}"
        ddBase = pyMoDELib.DislocationDynamicsBase(str(work_dir))
        ddio = pyMoDELib.DDconfigIO(str(work_dir / "evl"))

        # Box geometry
        xMin = np.array(ddBase.mesh.xMin(), dtype=float)
        xMax = np.array(ddBase.mesh.xMax(), dtype=float)
        xCenter = np.array(ddBase.mesh.xCenter(), dtype=float)

        # read evl file
        ddio.readTxt(evl_step)

        # ------------- First pass: find active planes over ALL EVLs -------------
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

        # find loop pair
        glide_plane_loops = defaultdict(list)
        for loop_id, v in node_pos_per_loop.items():
            unique_z_pos = np.unique(v[:,2]).item()
            glide_plane_loops[unique_z_pos].append(loop_id)
        partial_pairs = list(glide_plane_loops.values())
        print(partial_pairs)

        # calc separation distance
        #separation_dists = {}
        for partial_pair in partial_pairs:
            avg_line_pos = []
            for partial_loop_id in partial_pair:
                node_pos = node_pos_per_loop[partial_loop_id]
                mean_line_pos = np.mean(node_pos[:, glide_axis])
                avg_line_pos.append(mean_line_pos)
            #separation_dists[tuple(partial_pair)] = np.abs(np.diff(avg_line_pos)) 
            separation_dists_total[tuple(partial_pair)].append(np.abs(np.diff(avg_line_pos)))

    separation_dists_total = {k: np.array(v) for k, v in separation_dists_total.items()}

    print(separation_dists_total)
    exit()

    # automatically deletes the tmp tree if there is any
    if tmp:
        tmp.cleanup()

if __name__ == "__main__":
    main()
