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
sys.path.append("../../python/lib")
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

    #ddBase.simulationParameters.traitsIO
    #print(ddBase.simulationParameters.traitsIO)

    # ------------- First pass: find active planes over ALL EVLs -------------
    active_plane_keys_global = set()
    for evl_idx in evl_indices:
        ddio.readTxt(evl_idx)
        ddauxio.readTxt(evl_idx)

        print(pyMoDELib.DislocationNodeIO)
        print(dir(pyMoDELib.DislocationNodeIO))
        print(np.array(pyMoDELib.DislocationNodeIO.V))
        print(dir(pyMoDELib.DislocationNodeIO.V))
        print(pyMoDELib.DislocationNodeIO.P)

        defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
        defectiveCrystal.initializeConfiguration(ddio)
        DN = defectiveCrystal.dislocationNetwork()
        glidePlanes = ddBase.glidePlanes()

        for qp in ddauxio.quadraturePoints:
            print(qp.r)
            print(qp.stackingFaultForce)
            print(qp.lineTensionForce)
            print(qp.velocity)
            print(qp.pkForce)
            exit()

        exit()

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
        plane_to_segments, unassigned = build_plane_to_segments_geometric(
            DN, node_pos, plane_geometry, tol=1e-4
        )

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

    fig2, axes2 = plt.subplots(
        #rows, cols, figsize=(4.2 * cols, 4.2 * rows), squeeze=False
        rows, cols, figsize=(4.2 * cols, 8 * rows), squeeze=False
    )
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

        #ax.set_title(f"Plane pk={pk}", fontsize=9)
        ax.set_title(f"Plane pk={pk}")

    # hide unused axes
    for j in range(i + 1, rows * cols):
        r, c = divmod(j, cols)
        axes2[r][c].set_visible(False)

    #fig2.suptitle("Glide planes – evolution", fontsize=14)
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
        plane_to_segments, unassigned = build_plane_to_segments_geometric(
            DN, node_pos, plane_geometry, tol=1e-4
        )

        # convert value to nparray and remove empty list
        plane_to_segments = {k: np.array(v) for k,v in plane_to_segments.items() if v}

        # draw new segments
        for pk in occupied_pks:
            ax = ax_for_pk[pk]
            col = key_to_color.get(pk, (0.1, 0.1, 0.1, 1.0))
            if pk not in plane_frames:
                continue
            origin, u, v, n = plane_frames[pk]
            for a, b in plane_to_segments.get(pk, []):
                x0, y0 = project(a, origin, u, v)
                x1, y1 = project(b, origin, u, v)
                (ln,) = ax.plot([x0, x1], [y0, y1], lw=config["dislocation_line_width"], alpha=1.0, color=col)
                segment_artists[pk].append(ln)

        # Show EVL index in a super-title or text
        fig2.suptitle(f"Glide planes segment – evolution ({evl_idx})")

        return [ln for lines in segment_artists.values() for ln in lines]

    def animate(frame_idx):
        evl_idx = evl_indices[frame_idx]
        print(f"Rendering EVL {evl_idx}")
        return draw_frame_for_evl(evl_idx)

    ani = animation.FuncAnimation(
        fig2, animate, frames=len(evl_indices), interval=200, blit=False
    )

    # --------- Save movie (requires ffmpeg) ---------
    writer = animation.FFMpegWriter(fps=config["FPS"], bitrate=1800)
    #ani.save("glide_planes_evolution.mp4", writer=writer, dpi=200)
    ani.save("glide_planes_evolution.mp4", writer=writer)
    plt.close(fig2)


    if tmp:               # automatically deletes the tmp tree
        tmp.cleanup()

if __name__ == "__main__":
    main()
