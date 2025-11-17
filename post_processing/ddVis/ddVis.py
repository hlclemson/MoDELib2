# /opt/local/bin/python3.12 test.py
import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict
import matplotlib.cm as cm
import json

plt.rcParams["text.usetex"] = True
sys.path.append("../../python")
sys.path.append("../../python/lib")
# from readFile import *
from modlibUtils import *

sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib
import networkx as nx
import math
from pathlib import Path


def build_plane_geometry_from_segments(plane_segments_by_key, eps=1e-10):
    plane_geometry = {}
    for pk, segs in plane_segments_by_key.items():
        if len(segs) < 1:
            continue
        # first segment, point + direction
        P0a, P1a = segs[0]
        v1 = P1a - P0a
        p0 = 0.5 * (P0a + P1a)
        # find a second non-parallel direction
        v2 = None
        for Q0, Q1 in segs[1:]:
            w = Q1 - Q0
            cross = np.cross(v1, w)
            if np.linalg.norm(cross) > eps:
                v2 = w
                break
        if v2 is None:
            # degenerate box intersection: fabricate v2 in-plane
            axis = np.eye(3)[np.argmin(np.abs(v1))]
            v2 = axis - np.dot(axis, v1) * v1 / (np.dot(v1, v1) + eps)
        n = np.cross(v1, v2)
        n_norm = np.linalg.norm(n)
        if n_norm < eps:
            continue
        n_hat = n / n_norm
        plane_geometry[pk] = (n_hat, p0)
    return plane_geometry


def build_plane_to_segments_geometric(DN, node_pos, plane_geometry, tol=1e-4):
    plane_to_segments = {pk: [] for pk in plane_geometry.keys()}
    unassigned = []  # for debugging
    for lkey in DN.networkLinks().keys():
        link = DN.networkLinks().getRef(lkey)
        if hasattr(link, "hasZeroBurgers") and link.hasZeroBurgers():
            continue
        a = node_pos[link.source.networkID()]
        b = node_pos[link.sink.networkID()]
        best_pk = None
        best_dist = float("inf")
        for pk, (n_hat, p0) in plane_geometry.items():
            d = segment_plane_distance(a, b, n_hat, p0)
            if d < best_dist:
                best_dist = d
                best_pk = pk
        if best_pk is not None and best_dist < tol:
            plane_to_segments[best_pk].append((a, b))
        else:
            unassigned.append((a, b))
    return plane_to_segments, unassigned


def segment_plane_distance(a, b, n_hat, p0):
    da = np.dot(a - p0, n_hat)
    db = np.dot(b - p0, n_hat)
    return max(abs(da), abs(db))


def draw_box(ax, xmin, xmax, lw=1.0, alpha=0.25):
    X0, X1 = np.array(xmin), np.array(xmax)
    xs = [X0[0], X1[0]]
    ys = [X0[1], X1[1]]
    zs = [X0[2], X1[2]]
    edges = [  # 12 edges
        ([xs[0], xs[1]], [ys[0], ys[0]], [zs[0], zs[0]]),
        ([xs[0], xs[1]], [ys[1], ys[1]], [zs[0], zs[0]]),
        ([xs[0], xs[1]], [ys[0], ys[0]], [zs[1], zs[1]]),
        ([xs[0], xs[1]], [ys[1], ys[1]], [zs[1], zs[1]]),
        ([xs[0], xs[0]], [ys[0], ys[1]], [zs[0], zs[0]]),
        ([xs[1], xs[1]], [ys[0], ys[1]], [zs[0], zs[0]]),
        ([xs[0], xs[0]], [ys[0], ys[1]], [zs[1], zs[1]]),
        ([xs[1], xs[1]], [ys[0], ys[1]], [zs[1], zs[1]]),
        ([xs[0], xs[0]], [ys[0], ys[0]], [zs[0], zs[1]]),
        ([xs[1], xs[1]], [ys[0], ys[0]], [zs[0], zs[1]]),
        ([xs[0], xs[0]], [ys[1], ys[1]], [zs[0], zs[1]]),
        ([xs[1], xs[1]], [ys[1], ys[1]], [zs[0], zs[1]]),
    ]
    for ex, ey, ez in edges:
        ax.plot(ex, ey, ez, linewidth=lw, alpha=alpha, color="k")


############################################################################

# ------------- Create data objects ------------- #
#simulationDir = Path(
#    "../partial_orient_615AA_crss_2-3_almg5/generatedData/seed1/seed1/101MPa"
#)
# read initial configuration file
#config_fname = sys.argv[1]
with open("config.json", "r") as f:
    config = json.load(f)

simulationDir = Path(config["data_path"])
ddBase = pyMoDELib.DislocationDynamicsBase(str(simulationDir))
ddio = pyMoDELib.DDconfigIO(str(simulationDir / "evl"))
evl_file = 500
ddio.readTxt(evl_file)
defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
defectiveCrystal.initializeConfiguration(ddio)
DN = defectiveCrystal.dislocationNetwork()
glidePlanes = ddBase.glidePlanes()
# ----------------------------------------------------- #


# --- Box  ---
xMin = np.array(ddBase.mesh.xMin(), dtype=float)
xMax = np.array(ddBase.mesh.xMax(), dtype=float)
xCenter = np.array(ddBase.mesh.xCenter(), dtype=float)


# --- Function for Getting key to Glide Plane ---
def stable_plane_key(gp):
    h, k, l, layer, latID = gp.key_tuple()
    return (h, k, l, layer, latID)


# --- All planes & a color per plane ---
planes = getattr(ddBase, "glidePlanes", lambda: [])()
plane_list = list(planes) if planes else []
plane_keys = [stable_plane_key(gp) for gp in plane_list]  # Not used
print(f" ... Total glide planes: {len(plane_keys)} ...")

# --- Collect glide-plane box intersections ---
plane_segments_by_key = {}
for gp in plane_list:
    key = stable_plane_key(gp)
    segs = gp.meshSegments() if hasattr(gp, "meshSegments") else []
    plane_segments_by_key[key] = [
        (np.array(P0, float), np.array(P1, float)) for (P0, P1) in segs
    ]

# --- Build node lookup from DN (network nodes) ---
node_pos = {}
for key in DN.networkNodes().keys():
    node = DN.networkNodes().getRef(key)
    node_pos[node.networkID()] = np.array(node.position(), float)

# --- Assign segments to planes geometrically ---
plane_geometry = build_plane_geometry_from_segments(plane_segments_by_key)
plane_to_segments, unassigned = build_plane_to_segments_geometric(
    DN, node_pos, plane_geometry, tol=1e-4
)
print(
    f"Active planes (with segments, geometric): {sum(len(v)>0 for v in plane_to_segments.values())}"
)
print(f"Unassigned segments (didn't fit any plane within tol): {len(unassigned)}")
geom_keys = [
    pk for pk in plane_geometry.keys()
]  # build colors from the geometric plane keys:
active_plane_keys = [pk for pk, segs in plane_to_segments.items() if len(segs) > 0]
print(f"Active planes (with segments): {len(active_plane_keys)}")
cmap = plt.cm.get_cmap("tab20", len(active_plane_keys) if active_plane_keys else 1)
key_to_color = {pk: cmap(i % 20) for i, pk in enumerate(active_plane_keys)}


# --- Plotting ---
fig = plt.figure(figsize=(9, 8))
ax = fig.add_subplot(111, projection="3d")
draw_box(ax, xMin, xMax, lw=3, alpha=0.6)

# --- Plane outlines: only for planes with segments ---
for pk in active_plane_keys:
    col = key_to_color[pk]
    for P0, P1 in plane_segments_by_key.get(pk, []):
        ax.plot(
            [P0[0], P1[0]],
            [P0[1], P1[1]],
            [P0[2], P1[2]],
            linestyle="-",
            linewidth=1.5,
            alpha=0.35,
            color=col,
        )

# -- Plot segments, by geometric plane
for pk in active_plane_keys:
    col = key_to_color[pk]
    for a, b in plane_to_segments[pk]:
        ax.plot(
            [a[0], b[0]],
            [a[1], b[1]],
            [a[2], b[2]],
            lw=2.0,
            alpha=1.0,
            color=col,
        )

ax.set_axis_off()
ax.set_box_aspect((xMax - xMin))
# plt.show()

# ========= Per-plane 2D plots (all to the same scale) =========
def _plane_basis_from_points(P):
    """Orthonormal (u,v,n) via PCA on points P (Nx3)."""
    P = np.asarray(P, float)
    origin = P.mean(axis=0)
    Q = P - origin
    if Q.shape[0] < 3:
        # fabricate a tiny in-plane basis if degenerate
        u = np.array([1.0, 0.0, 0.0])
        v = np.array([0.0, 1.0, 0.0])
        n = np.cross(u, v)
        return origin, u, v, n / np.linalg.norm(n)
    _, _, Vt = np.linalg.svd(Q, full_matrices=False)
    u = Vt[0]
    v = Vt[1]
    u /= np.linalg.norm(u)
    v /= np.linalg.norm(v)
    n = np.cross(u, v)
    n /= np.linalg.norm(n)
    return origin, u, v, n


def _project(P, origin, u, v):
    d = np.asarray(P, float) - origin
    return np.dot(d, u), np.dot(d, v)


# 1) Build a plane frame from either box edges (preferred) or plane's segments
plane_frames = {}  # pk -> (origin,u,v,n)
plane_proj_segments = {}  # pk -> list of ((x0,y0),(x1,y1))
plane_proj_edges = {}  # pk -> list of ((x0,y0),(x1,y1))

for pk in plane_to_segments.keys():
    pts = []
    # Prefer box intersections (span the plane footprint)
    for P0, P1 in plane_segments_by_key.get(pk, []):
        pts += [P0, P1]
    # Fallback to dislocation segments if needed
    if len(pts) < 3:
        for a, b in plane_to_segments.get(pk, []):
            pts += [a, b]
    if len(pts) < 3:
        continue  # nothing we can do for this plane

    origin, u, v, n = _plane_basis_from_points(np.asarray(pts, float))
    plane_frames[pk] = (origin, u, v, n)

# 2) Project everything and collect global bounds
all_xy = []
for pk, (origin, u, v, n) in plane_frames.items():
    # Dislocation segments
    segs2d = []
    for a, b in plane_to_segments.get(pk, []):
        x0, y0 = _project(a, origin, u, v)
        x1, y1 = _project(b, origin, u, v)
        segs2d.append(((x0, y0), (x1, y1)))
        all_xy += [(x0, y0), (x1, y1)]
    plane_proj_segments[pk] = segs2d

    # Plane outline (box intersections)
    edges2d = []
    for P0, P1 in plane_segments_by_key.get(pk, []):
        x0, y0 = _project(P0, origin, u, v)
        x1, y1 = _project(P1, origin, u, v)
        edges2d.append(((x0, y0), (x1, y1)))
        all_xy += [(x0, y0), (x1, y1)]
    plane_proj_edges[pk] = edges2d

if not all_xy:
    print("No projected points; skipping per-plane figure.")
else:
    all_xy = np.asarray(all_xy, float)
    xmin, ymin = all_xy.min(axis=0)
    xmax, ymax = all_xy.max(axis=0)
    dx, dy = (xmax - xmin), (ymax - ymin)
    pad = 0.05 * (max(dx, dy) if max(dx, dy) > 0 else 1.0)
    xlim = (xmin - pad, xmax + pad)
    ylim = (ymin - pad, ymax + pad)

    # 3) Grid layout and plot
    occupied_pks = list(plane_frames.keys())
    nplanes2d = len(occupied_pks)
    cols = min(4, nplanes2d)  # up to 4 columns; adjust as you like
    rows = int(math.ceil(nplanes2d / cols))

    fig2, axes2 = plt.subplots(
        rows, cols, figsize=(4.2 * cols, 4.2 * rows), squeeze=False
    )
    for ax in axes2.ravel():
        ax.axis("off")

    for i, pk in enumerate(occupied_pks):
        r, c = divmod(i, cols)
        ax = axes2[r][c]
        ax.axis("on")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.grid(True, ls=":", alpha=0.35)

        # outline (faint)
        for (x0, y0), (x1, y1) in plane_proj_edges.get(pk, []):
            ax.plot([x0, x1], [y0, y1], lw=1.2, alpha=0.6, color="k")

        # dislocation segments (colored)
        col = key_to_color.get(pk, (0.1, 0.1, 0.1, 1.0))
        for (x0, y0), (x1, y1) in plane_proj_segments.get(pk, []):
            ax.plot([x0, x1], [y0, y1], lw=2.0, alpha=1.0, color=col)

        ax.set_title(
            f"Plane pk={pk} • segs: {len(plane_proj_segments[pk])}", fontsize=9
        )

    # hide unused
    for j in range(i + 1, rows * cols):
        r, c = divmod(j, cols)
        axes2[r][c].set_visible(False)

    fig2.suptitle("Glide planes (2D view, common scale)", fontsize=14)
    plt.tight_layout()

plt.show()
