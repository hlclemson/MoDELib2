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
from collections import defaultdict

# ----- MoDELib / utils paths -----
sys.path.append("../../python")
#sys.path.append("../../python/lib")
#from python.visUtils import build_plane_to_segments_geometric_vel
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
    glide_axis = box_axis_idx[config["glide_axis"]]
    line_tangent_axis = box_axis_idx[config["line_tangent_axis"]]
    data_path = Path(config["data_path"])
    runID = 1000

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

    ddBase = pyMoDELib.DislocationDynamicsBase(str(work_dir))
    ddio = pyMoDELib.DDconfigIO(str(work_dir / "evl"))
    #ddseg = pyMoDELib.DislocationSegmentIO()
    # Box geometry
    xMin = np.array(ddBase.mesh.xMin(), dtype=float)
    xMax = np.array(ddBase.mesh.xMax(), dtype=float)
    xCenter = np.array(ddBase.mesh.xCenter(), dtype=float)
    # read evl file
    ddio.readTxt(runID)

    # First pass: find active planes over ALL EVLs
    active_plane_keys_global = set()
    defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
    defectiveCrystal.initializeConfiguration(ddio)
    DN = defectiveCrystal.dislocationNetwork()
    #planes = ddBase.glidePlanes()
    glide_planes = ddBase.glidePlanes()
    # ------------- First pass: find active planes over ALL EVLs -------------
    active_plane_keys_global = set()

    planes = getattr(ddBase, "glidePlanes", lambda: [])()
    plane_list = list(planes) if planes else []
    #plane_keys = [stable_plane_key(gp) for gp in plane_list]  # Not used
    #print(f" ... Total glide planes: {len(plane_keys)} ...")

    ### this is how you return burgers vector
    #for lkey in DN.networkLinks().keys():
    #    link = DN.networkLinks().getRef(lkey)
    #    print(link)
    #    print(dir(link))
    #    print(link.burgers())
    #    exit()

    node_pos = {}
    #node_velocity = {}
    for key in DN.networkNodes().keys():
        node = DN.networkNodes().getRef(key)
        node_pos[node.networkID()] = np.asarray(node.position(), dtype=np.float32)
        #node_velocity[node.networkID()] = np.asarray(node.velocity(), dtype=np.float32)

    plane_segments_by_key = {}
    for gp in plane_list:
        key = stable_plane_key(gp)
        segs = gp.meshSegments() if hasattr(gp, "meshSegments") else []
        plane_segments_by_key[key] = [
            (np.array(P0, float), np.array(P1, float)) for (P0, P1) in segs
        ]

    plane_geometry = build_plane_geometry_from_segments(plane_segments_by_key)
    #plane_to_segments, unassigned = build_plane_to_segments_geometric(
    plane_to_segments = build_plane_to_segments_geometric(
        DN, node_pos, plane_geometry, tol=1e-4
    )
    #plane_to_segments, plane_to_segments_velocity = build_plane_to_segments_geometric(
        #DN, node_pos, plane_geometry, tol=1e-4
        #DN, node_pos, node_velocity, plane_geometry, tol=1e-4
    #)

    for pk, segs in plane_to_segments.items():
        if len(segs) > 0:
            active_plane_keys_global.add(pk)

    active_plane_keys_global = sorted(active_plane_keys_global)
    print("Planes that EVER have segments:", active_plane_keys_global)
    if not active_plane_keys_global:
        raise RuntimeError("No active planes with segments found over the EVL range.")

    # ---cColormap for active planes ---
    cmap = plt.cm.get_cmap("tab20", len(active_plane_keys_global))
    key_to_color = {pk: cmap(i % 20) for i, pk in enumerate(active_plane_keys_global)}

    # ------------- Build static 2D frames (basis + edges + global x/y) -------------
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
        origin, u, v, n = plane_basis_from_points(np.asarray(pts, dtype=np.float32))
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
    #pad = 0.05 * (max(dx, dy) if max(dx, dy) > 0 else 1.0)
    pad = 0.00 * (max(dx, dy) if max(dx, dy) > 0 else 1.0)
    xlim = (xmin - pad, xmax + pad)
    ylim = (ymin - pad, ymax + pad)

    # ------------- Build 2D subplot grid (one panel per active plane) -------------
    occupied_pks = list(plane_frames.keys())
    nplanes2d = len(occupied_pks)
    #cols = min(6, nplanes2d)
    cols = min(2, nplanes2d)
    rows = int(math.ceil(nplanes2d / cols))

    #fig2, axes2 = plt.subplots(
    #    #rows, cols, figsize=(4.2 * cols, 4.2 * rows), squeeze=False
    #    rows, cols, figsize=(4.2 * cols, 8 * rows), squeeze=False
    #)
    cols = 2
    rows = (len(active_plane_keys_global) + cols - 1) // cols  # ceiling division
    fig2, axes2 = plt.subplots(rows, cols, figsize=(8, 4*rows), squeeze=False)

    for ax in axes2.ravel():
        ax.axis("off")

    ax_for_pk = {}
    segment_artists = {pk: [] for pk in occupied_pks}

    for i, pk in enumerate(occupied_pks):
        r, c = divmod(i, cols)
        ax = axes2[r][c]
        ax_for_pk[pk] = ax
        ax.axis("on")
        #ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.grid(True, ls=":", alpha=0.35)

        # static plane outline, glide plane edges
        for (x0, y0), (x1, y1) in plane_proj_edges.get(pk, []):
            #ax.plot([x0, x1], [y0, y1], lw=1.0, alpha=0.6, color="k")
            ax.plot([x0, x1], [y0, y1], lw=config["glide_plane_line_width"], alpha=0.6, color="k")
        ax.set_title(f"Plane key={pk}")

    # hide unused axes
    for j in range(i + 1, rows * cols):
        r, c = divmod(j, cols)
        axes2[r][c].set_visible(False)
    fig2.suptitle("Glide planes – evolution")
    plt.tight_layout()

    # ------------- Animation: update per EVL -------------
    def draw_frame_for_evl(evl_idx):
        # clear old segment lines
        for pk in occupied_pks:
            for ln in segment_artists[pk]:
                ln.remove()
            segment_artists[pk].clear()

        ddio.readTxt(evl_idx)
        defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
        defectiveCrystal.initializeConfiguration(ddio)
        DN = defectiveCrystal.dislocationNetwork()
        glidePlanes = ddBase.glidePlanes()

        planes = getattr(ddBase, "glidePlanes", lambda: [])()
        plane_list = list(planes) if planes else []
        plane_keys = [stable_plane_key(gp) for gp in plane_list]  # Not used
        print(f" ... Total glide planes: {len(plane_keys)} ...")

        node_pos = {}
        node_velocity = {}
        for key in DN.networkNodes().keys():
            node = DN.networkNodes().getRef(key)
            node_pos[node.networkID()] = np.asarray(node.position(), dtype=np.float32)
            node_velocity[node.networkID()] = np.asarray(node.velocity(), dtype=np.float32)

        plane_segments_by_key = {}
        for gp in plane_list:
            key = stable_plane_key(gp)
            segs = gp.meshSegments() if hasattr(gp, "meshSegments") else []
            plane_segments_by_key[key] = [
                (np.array(P0, float), np.array(P1, float)) for (P0, P1) in segs
            ]

        plane_geometry = build_plane_geometry_from_segments(plane_segments_by_key)
        plane_to_segments, plane_to_segments_velocity, plane_to_segments_burger = build_plane_to_segments_geometric_vel_burger(
            DN, node_pos, node_velocity, plane_geometry, tol=1e-4
        )

        # convert value to nparray and remove empty list
        plane_to_segments = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments.items() if v}
        plane_to_segments_velocity = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_velocity.items() if v}
        plane_to_segments_burger = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_burger.items() if v}
        # the problem is that source and sink shares the same burgers vector
        # how do I assign correct position for the burgers vector arrows?
        print(plane_to_segments_burger)
        exit()

        # draw new segments
        for pk in occupied_pks:
            ax = ax_for_pk[pk]
            col = key_to_color.get(pk, (0.1, 0.1, 0.1, 1.0))
            if pk not in plane_frames:
                continue
            origin, u, v, n = plane_frames[pk]
            avg_x = []
            avg_y = []
            for a, b in plane_to_segments.get(pk, []):
                x0, y0 = project(a, origin, u, v)
                x1, y1 = project(b, origin, u, v)
                avg_x.append(x0)
                avg_y.append(y0)
                (ln,) = ax.plot([x0, x1], [y0, y1], lw=config["dislocation_line_width"], alpha=1.0, color=col)
                #(ln,) = ax.plot([x0, x1], [y0, y1], lw=config["dislocation_line_width"], alpha=1.0, color="k")
                segment_artists[pk].append(ln)
            avg_x = np.mean(avg_x)
            avg_y = np.mean(avg_y)
            xlim = (avg_x - 20, avg_x + 20)
            #ylim = (avg_y - 10, avg_y + 10)
            ax.set_xlim(*xlim)
            #ax.set_ylim(*ylim)
        #exit()

        # Show EVL index in a super-title or text
        fig2.suptitle(f"Glide planes segment – (ID: {evl_idx})")

        return [ln for lines in segment_artists.values() for ln in lines]

    def animate(frame_idx):
        runID = runIDs[frame_idx]
        print(f"Rendering EVL {runID}")
        return draw_frame_for_evl(runID)

    ani = animation.FuncAnimation(
        fig2, animate, frames=len(runIDs), interval=200, blit=False
    )

    # --------- Save movie (requires ffmpeg) ---------
    writer = animation.FFMpegWriter(fps=config["FPS"], bitrate=1800)
    ani.save("glide_planes_evolution.mp4", writer=writer, dpi=200)
    #ani.save("glide_planes_evolution.mp4", writer=writer)
    plt.close(fig2)

    # automatically deletes the tmp tree if there is any
    if tmp:
        tmp.cleanup()

if __name__ == "__main__":
    main()
