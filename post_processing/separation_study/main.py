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
    with open("config.json", "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])
    evl_step = config["evl_step"]

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

    # test 100 pair
    num_pairs = 100
    noise_seed_to_test = range(1, num_pairs+1)
    glide_inc_step = 4
    ref_glide_pos = 400
    glide_steps_to_test = [
        [ref_glide_pos, ref_glide_pos + (glide_inc_step * i)] for i in range(1, num_pairs+1)
    ]

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
    glidePlanes = ddBase.glidePlanes()

    planes = getattr(ddBase, "glidePlanes", lambda: [])()
    plane_list = list(planes) if planes else []
    plane_keys = [stable_plane_key(gp) for gp in plane_list]  # Not used
    print(f" ... Total glide planes: {len(plane_keys)} ...")

    node_pos = {}
    for key in DN.networkNodes().keys():
        node = DN.networkNodes().getRef(key)
        node_pos[node.networkID()] = np.array(node.position(), float)

    plane_segments_by_key = {}
    for gp in plane_list:
        key = stable_plane_key(gp)
        segs = gp.meshSegments() if hasattr(gp, "meshSegments") else []
        plane_segments_by_key[key] = [
            (np.array(P0, float), np.array(P1, float)) for (P0, P1) in segs
        ]

    plane_geometry = build_plane_geometry_from_segments(plane_segments_by_key)
    plane_to_segments, _ = build_plane_to_segments_geometric(
        DN, node_pos, plane_geometry, tol=1e-4
    )

    for pk, segs in plane_to_segments.items():
        if len(segs) > 0:
            active_plane_keys_global.add(pk)

    active_plane_keys_global = sorted(active_plane_keys_global)
    print("Planes that EVER have segments:", active_plane_keys_global)
    if not active_plane_keys_global:
        raise RuntimeError("No active planes with segments found over the EVL range.")

    # Build static 2D frames (basis + edges + global x/y)
    plane_frames = {}  # pk -> (origin, u, v, n)
    plane_proj_edges = {}  # pk -> list of ((x0,y0), (x1,y1))
    all_xy = []
    for pk in active_plane_keys_global:
        segs = plane_segments_by_key.get(pk, [])
        if len(segs) < 1:
            continue
        pts = []
        for P0, P1 in segs:
            pts += [P0, P1]
        if len(pts) < 3:
            continue
        origin, u, v, n = plane_basis_from_points(np.asarray(pts, float))
        plane_frames[pk] = (origin, u, v, n)

        edges2d = []
        for P0, P1 in segs:
            x0, y0 = project(P0, origin, u, v)
            x1, y1 = project(P1, origin, u, v)
            edges2d.append(((x0, y0), (x1, y1)))
            all_xy += [(x0, y0), (x1, y1)]
        plane_proj_edges[pk] = edges2d

    if not all_xy:
        raise RuntimeError("No projected edge points; cannot build 2D frames.")
    all_xy = np.asarray(all_xy, float)

    node_pos = {}
    for key in DN.networkNodes().keys():
        node = DN.networkNodes().getRef(key)
        node_pos[node.networkID()] = np.array(node.position(), float)

    plane_segments_by_key = {}
    for gp in plane_list:
        key = stable_plane_key(gp)
        segs = gp.meshSegments() if hasattr(gp, "meshSegments") else []
        plane_segments_by_key[key] = [
            (np.array(P0, float), np.array(P1, float)) for (P0, P1) in segs
        ]

    plane_geometry = build_plane_geometry_from_segments(plane_segments_by_key)
    plane_to_segments, _ = build_plane_to_segments_geometric(
        DN, node_pos, plane_geometry, tol=1e-4
    )

    # convert value to nparray and remove empty list
    plane_to_segments = {k: np.array(v) for k,v in plane_to_segments.items() if v}

    # convert source_pos and sink_pos pair to plain node pos info
    for k, v in plane_to_segments.items():
        # stack source and sink positions (N, 2, 3) -> (N*2, 3)
        tmp_val = v.reshape(-1, 3)
        # Remove exact duplicate rows (source sink overlap)
        tmp_val = np.unique(tmp_val, axis=0)
        # replace the dictionary entry
        plane_to_segments[k] = tmp_val

    # for debugging
    # visualizes the node positional information
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Scatter plot
    for k, v in plane_to_segments.items():
        ax.scatter(v[:,0], v[:,1], v[:,2])
    plt.show()

    # automatically deletes the tmp tree if there is any
    if tmp:
        tmp.cleanup()

if __name__ == "__main__":
    main()
