# standard lib
from math import e
from os import path
import re
import sys
import json
import shutil
import tarfile
import tempfile

# 3rd party lib
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
from matplotlib import rcParams
from collections import defaultdict

# ----- MoDELib / utils paths -----
sys.path.append("../../python")
from visUtils import *
from modlibUtils import *

sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib

# Configure global plot settings (applies to all figures)
rcParams.update(
    {
        "figure.dpi": 400,
        # "figure.autolayout": True,  # Prevent label clipping
        # "axes.grid": True,
        # "grid.alpha": 0.6,
        "text.usetex": False,
        "font.size": 10,  # Default font size for text
        "mathtext.fontset": "stix",  # Use STIX font for math text
        "font.family": "serif",  # Use serif font (matches LaTeX default)
    }
)


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
    data_path = Path(config["data_path"])
    glide_axis = box_axis_idx[config["glide_axis"]]
    line_tangent_axis = box_axis_idx[config["line_tangent_axis"]]
    runID = config["target_runID"]

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
            tf.extractall(tmp.name)
        work_dir = Path(tmp.name) / src.name.removesuffix(".tar.gz")
    else:
        raise FileNotFoundError("neither directory nor .tar.gz found")

    # read labels
    with open(work_dir / "F/F_labels.txt", "r") as label:
        fLabels = label.read()
    # remove empty element and store it as a list
    fLabels = np.array([x for x in fLabels.split("\n") if x], dtype=str)
    # open F file
    fData = np.loadtxt(work_dir / "F/F_0.txt")
    # find the last step
    runIDs = getFarray(fData, fLabels, "runID").astype(int)

    ddBase = pyMoDELib.DislocationDynamicsBase(str(work_dir))
    ddio = pyMoDELib.DDconfigIO(str(work_dir / "evl"))
    # ddseg = pyMoDELib.DislocationSegmentIO()
    # Box geometry
    xMin = np.array(ddBase.mesh.xMin(), dtype=float)
    xMax = np.array(ddBase.mesh.xMax(), dtype=float)
    xCenter = np.array(ddBase.mesh.xCenter(), dtype=float)
    # read evl file
    ddio.readTxt(runID)

    defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
    defectiveCrystal.initializeConfiguration(ddio)
    DN = defectiveCrystal.dislocationNetwork()
    planes = ddBase.glidePlanes()

    node_pos = {}
    node_velocity = {}
    node_pos_per_loop = defaultdict(list)
    for key in DN.networkNodes().keys():
        node = DN.networkNodes().getRef(key)
        node_pos[node.networkID()] = np.asarray(node.position(), dtype=np.float32)
        node_velocity[node.networkID()] = np.asarray(node.velocity(), dtype=np.float32)

    plane_segments_by_key = {}
    for gp in planes:
        key = stable_plane_key(gp)
        segs = gp.meshSegments() if hasattr(gp, "meshSegments") else []
        plane_segments_by_key[key] = [
            (np.array(P0, float), np.array(P1, float)) for (P0, P1) in segs
        ]

    plane_geometry = build_plane_geometry_from_segments(plane_segments_by_key)
    (
        plane_to_segments,
        plane_to_segments_velocity,
        plane_to_segments_b,
        plane_to_segments_b_tail,
    ) = build_plane_to_segments_geo_b(
        DN, node_pos, node_velocity, plane_geometry, tol=1e-4
    )

    # combine source and sink pos dimension
    xyz_colnum = 3
    for (k1, v1), (k2, v2) in zip(
        plane_to_segments.items(), plane_to_segments_velocity.items()
    ):
        # combine source and sink pos dimension
        v1 = np.reshape(v1, (-1, xyz_colnum))
        v2 = np.reshape(v2, (-1, xyz_colnum))
        # remove duplicate (source and sink)
        v1 = np.unique(v1, axis=0)
        v2 = np.unique(v2, axis=0)
        plane_to_segments[k1] = v1
        plane_to_segments_velocity[k2] = v2

    # find active planes over ALL EVLs
    active_plane_keys_global = set()
    for pk, segs in plane_to_segments.items():
        if len(segs) > 0:
            active_plane_keys_global.add(pk)

    active_plane_keys_global = sorted(active_plane_keys_global)
    print("Planes that EVER have segments:", active_plane_keys_global)
    if not active_plane_keys_global:
        raise RuntimeError("No active planes with segments found over the EVL range.")

    # cColormap for active planes
    cmap = plt.cm.get_cmap("tab20", len(active_plane_keys_global))
    key_to_color = {pk: cmap(i % 20) for i, pk in enumerate(active_plane_keys_global)}

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
        pts = np.asarray(pts, dtype=np.float32)
        origin, u, v, n = plane_basis_from_points(pts)
        # Return a consistently-oriented (u,v,n) frame
        #u, v, n = orient_basis(pts, origin, u, v, n)
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

    all_xy = np.asarray(all_xy, dtype=np.float32)
    xmin, ymin = all_xy.min(axis=0)
    xmax, ymax = all_xy.max(axis=0)
    dx, dy = (xmax - xmin), (ymax - ymin)
    pad = config["pad_glide_edge"] * (max(dx, dy) if max(dx, dy) > 0 else 1.0)
    xlim = (xmin - pad, xmax + pad)
    ylim = (ymin - pad, ymax + pad)

    # Build 2D subplot grid (one panel per active plane)
    occupied_pks = list(plane_frames.keys())
    # nplanes2d = len(occupied_pks)
    # cols = min(2, nplanes2d)
    # rows = int(math.ceil(nplanes2d / cols))
    cols = 2
    rows = (len(active_plane_keys_global) + cols - 1) // cols  # ceiling division
    fig, axes = plt.subplots(rows, cols, figsize=(8, 4 * rows), squeeze=False)

    for ax in axes.ravel():
        ax.axis("off")

    ax_for_pk = {}
    segment_artists = {pk: [] for pk in occupied_pks}
    for i, pk in enumerate(occupied_pks):
        r, c = divmod(i, cols)
        ax = axes[r][c]
        ax_for_pk[pk] = ax
        ax.axis("on")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.grid(True, ls=":", alpha=0.35)

        # static plane outline, glide plane edges
        for (x0, y0), (x1, y1) in plane_proj_edges.get(pk, []):
            ax.plot(
                [x0, x1],
                [y0, y1],
                lw=config["glide_plane_line_width"],
                alpha=0.6,
                color="k",
            )
        #ax.set_title(f"Plane key={pk}")

    # draw new segments
    for pk in occupied_pks:
        ax = ax_for_pk[pk]
        col = key_to_color.get(pk, (0.1, 0.1, 0.1, 1.0))
        if pk not in plane_frames:
            continue
        origin, u, v, n = plane_frames[pk]
        pos_v = [project(x, origin, u, v) for x in plane_to_segments.get(pk, [])]
        pos_v = np.asarray(pos_v)

        pos_v_x = pos_v[:, 0]
        pos_v_y = pos_v[:, 1]
        # draw dislocation line
        ax.scatter(pos_v_x, pos_v_y, s=10, color=col)

        # make the plot follow dislocation line
        if config["track_dislocation"]:
            # if config["glide_axis"]=="x":
            dis_buffer = config["dislocation_buffer"]
            avg_x = np.mean(pos_v_x)
            xlim = (avg_x - dis_buffer, avg_x + dis_buffer)
            ax.set_xlim(*xlim)
            # elif config["glide_axis"]=="y":
            #    avg_y = np.mean(pos_v_y)
            #    ylim = (avg_y - 10, avg_y + 10)
            #    ax.set_ylim(*ylim)

        # draw burgers vector
        burgers_v = [project(x, origin, u, v) for x in plane_to_segments_b.get(pk, [])]
        burgers_v_tail = [
            project(x, origin, u, v) for x in plane_to_segments_b_tail.get(pk, [])
        ]
        burgers_v = np.asarray(burgers_v)
        burgers_v_tail = np.asarray(burgers_v_tail)
        ax.quiver(
            burgers_v_tail[:, 0], burgers_v_tail[:, 1], burgers_v[:, 0], burgers_v[:, 1], label="Burgers"
        )
    ax.legend(fontsize=5)
    fig.savefig(config["figure_name"], transparent=True)

    # automatically deletes the tmp tree if there is any
    if tmp:
        tmp.cleanup()


if __name__ == "__main__":
    main()
