#!/opt/local/bin/python3.12
import sys, os, math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import matplotlib.cm as cm
from matplotlib import animation
import json
from pathlib import Path

plt.rcParams["text.usetex"] = True
# ----- MoDELib / utils paths -----
sys.path.append("../../python")
sys.path.append("../../python/lib")
from visUtils import *
#from readFile import *
from modlibUtils import *

sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib


# ================== MAIN SCRIPT ==================
if __name__ == "__main__":
    # ------------- Setup MoDELib objects -------------
    #simulationDir = os.path.join("../..", "tutorials", "prismaticDensityFCC")

    with open("config.json", "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])
    #ddBase = pyMoDELib.DislocationDynamicsBase(simulationDir)
    ddBase = pyMoDELib.DislocationDynamicsBase(str(simulationDir))
    #ddio = pyMoDELib.DDconfigIO(simulationDir + "/evl")
    ddio = pyMoDELib.DDconfigIO(str(simulationDir / "evl"))
    evl_file = 500


    # --- Box geometry ---
    xMin = np.array(ddBase.mesh.xMin(), dtype=float)
    xMax = np.array(ddBase.mesh.xMax(), dtype=float)
    xCenter = np.array(ddBase.mesh.xCenter(), dtype=float)

    # ------------- EVL range & discovery -------------
    start_evl = 0
    end_evl = 10_000
    evl_indices = []
    for i in range(start_evl, end_evl + 1):
        try:
            ddio.readTxt(i)
            evl_indices.append(i)
        except Exception:
            continue  # missing evl file, skip
    print("Found EVL files:", evl_indices)
    if not evl_indices:
        raise RuntimeError("No EVL files found in requested range.")

    # ------------- First pass: find active planes over ALL EVLs -------------
    active_plane_keys_global = set()
    for evl_idx in evl_indices:

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
    pad = 0.05 * (max(dx, dy) if max(dx, dy) > 0 else 1.0)
    xlim = (xmin - pad, xmax + pad)
    ylim = (ymin - pad, ymax + pad)

    # ------------- Build 2D subplot grid (one panel per active plane) -------------
    occupied_pks = list(plane_frames.keys())
    nplanes2d = len(occupied_pks)
    cols = min(6, nplanes2d)
    rows = int(math.ceil(nplanes2d / cols))

    fig2, axes2 = plt.subplots(
        rows, cols, figsize=(4.2 * cols, 4.2 * rows), squeeze=False
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
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.grid(True, ls=":", alpha=0.35)

        # static plane outline
        for (x0, y0), (x1, y1) in plane_proj_edges.get(pk, []):
            ax.plot([x0, x1], [y0, y1], lw=1.2, alpha=0.6, color="k")

        ax.set_title(f"Plane pk={pk}", fontsize=9)

    # hide unused axes
    for j in range(i + 1, rows * cols):
        r, c = divmod(j, cols)
        axes2[r][c].set_visible(False)

    fig2.suptitle("Glide planes – evolution", fontsize=14)
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
                (ln,) = ax.plot([x0, x1], [y0, y1], lw=2.0, alpha=1.0, color=col)
                segment_artists[pk].append(ln)

        # Show EVL index in a super-title or text
        fig2.suptitle(f"Glide planes segment – evolution ({evl_idx})", fontsize=14)

        return [ln for lines in segment_artists.values() for ln in lines]

    def animate(frame_idx):
        evl_idx = evl_indices[frame_idx]
        print(f"Rendering EVL {evl_idx}")
        return draw_frame_for_evl(evl_idx)

    ani = animation.FuncAnimation(
        fig2, animate, frames=len(evl_indices), interval=200, blit=False
    )

    # --------- Save movie (requires ffmpeg) ---------
    writer = animation.FFMpegWriter(fps=5, bitrate=1800)
    ani.save("glide_planes_evolution.mp4", writer=writer, dpi=200)
    plt.close(fig2)
