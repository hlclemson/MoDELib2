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
from visUtils import *
from modlibUtils import *

sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib


def main():
    with open("config.json", "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])

    # ------------- Setup MoDELib objects -------------
    src = simulationDir
    tmp = None # will hold TemporaryDirectory handle
    work_dir = None # will point to the folder we actually use
    if src.is_dir(): # in uncompressed format
        work_dir = src
    elif src.suffixes == ['.tar', '.gz']: # if .tar.gz format
        tmp = tempfile.TemporaryDirectory()
        with tarfile.open(src, 'r:gz') as tf:
            tf.extractall(tmp.name)
        work_dir = pathlib.Path(tmp.name) / src.name.removesuffix('.tar.gz')
    else:
        raise FileNotFoundError('neither directory nor .tar.gz found')

    # ----- use work_dir here -----
    ddBase = pyMoDELib.DislocationDynamicsBase(str(work_dir))
    ddio = pyMoDELib.DDconfigIO(str(work_dir / "evl"))
    ddauxio = pyMoDELib.DDauxIO(str(work_dir / "evl"))

    # --- Box geometry ---
    xMin = np.array(ddBase.mesh.xMin(), dtype=float)
    xMax = np.array(ddBase.mesh.xMax(), dtype=float)
    xCenter = np.array(ddBase.mesh.xCenter(), dtype=float)

    # ------------- EVL range & discovery -------------
    start_evl = config["evl_start_step"]
    end_evl = config["evl_end_step"]
    evl_indices = []
    for i in range(start_evl, end_evl + 1):
        try:
            ddio.readTxt(i)
            evl_indices.append(i)
        except Exception:
            continue  # missing evl file, skip
    if not evl_indices:
        raise RuntimeError("No EVL files found in requested range.")

    # ------------- First pass: find active planes over ALL EVLs -------------
    active_plane_keys_global = set()
    for evl_idx in evl_indices:
        ddio.readTxt(evl_idx)
        ddauxio.readTxt(evl_idx)

        defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
        defectiveCrystal.initializeConfiguration(ddio)
        DN = defectiveCrystal.dislocationNetwork()
        #glidePlanes = ddBase.glidePlanes()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')


        qps = []
        for qp in ddauxio.quadraturePoints:
            pos = qp.r
            qps.append(pos)
            qp_velocity = qp.velocity
            ax.quiver(*pos, *qp_velocity, length=100, normalize=False, arrow_length_ratio=.2, color='r')

            #print(qp.stackingFaultForce)
            #print(qp.lineTensionForce)
            #print(qp.velocity)
            #print(qp.pkForce)
        qps = np.array(qps)
        ax.scatter(qps[:, 0], qps[:, 1], qps[:, 2])
        plt.show()

        exit()


if __name__ == "__main__":
    main()
