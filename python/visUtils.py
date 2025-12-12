from collections import defaultdict
import sys, os, math
import numpy as np

# ================== Geometry helpers ==================
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
        for (Q0, Q1) in segs[1:]:
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

#def build_plane_to_segments_geometric(DN, node_pos, plane_geometry, tol=1e-4):
def build_plane_to_segments_geometric(DN, node_pos, plane_geometry, tol=1e-4):
    plane_to_segments = {pk: [] for pk in plane_geometry.keys()}
    unassigned = []   # for debugging
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

    #return plane_to_segments, unassigned
    return plane_to_segments


def build_plane_to_segments_geo_b(DN, node_pos, node_velocity, plane_geometry, tol=1e-4):
    plane_to_segments = defaultdict(list)
    plane_to_segments_velocity = defaultdict(list)
    plane_to_segments_b = defaultdict(list)
    plane_to_segments_b_tail = defaultdict(list)
    unassigned = []   # for debugging
    unassigned_velocity = []   # for debugging
    for lkey in DN.networkLinks().keys():
        link = DN.networkLinks().getRef(lkey)
        if hasattr(link, "hasZeroBurgers") and link.hasZeroBurgers():
            continue
        a = node_pos[link.source.networkID()]
        b = node_pos[link.sink.networkID()]
        # midpoint
        mean_pos = (a + b) / 2
        av = node_velocity[link.source.networkID()]
        bv = node_velocity[link.sink.networkID()]
        link_loopID = tuple(link.loopIDs())
        link_loopID = link_loopID[0]
        burgers = DN.loops().getRef(link_loopID).burgers()
        #print(burgers)
        #exit()
        #if link.source.networkID()>link.sink.networkID():
        #    link_burger = link.burgers()
        #else:
        #    link_burger = -link.burgers()
        #link_burger = link.burgers()
        best_pk = None
        best_dist = float("inf")
        for pk, (n_hat, p0) in plane_geometry.items():
            d = segment_plane_distance(a, b, n_hat, p0)
            if d < best_dist:
                best_dist = d
                best_pk = pk
        if best_pk is not None and best_dist < tol:
            plane_to_segments[best_pk].append((a, b))
            plane_to_segments_velocity[best_pk].append((av, bv))
            #plane_to_segments_b[best_pk].append(link_burger)
            plane_to_segments_b[best_pk].append(burgers)
            plane_to_segments_b_tail[best_pk].append(mean_pos)
        else:
            continue
    # convert value to nparray and remove empty list
    plane_to_segments = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments.items() if v}
    plane_to_segments_velocity = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_velocity.items() if v}
    plane_to_segments_b = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_b.items() if v}
    plane_to_segments_b_tail = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_b_tail.items() if v}

    return plane_to_segments, plane_to_segments_velocity, plane_to_segments_b, plane_to_segments_b_tail


def build_plane_to_segments_geo_b_testing(DN, node_pos, node_velocity, plane_geometry, tol=1e-4):
    #plane_to_segments = {pk: [] for pk in plane_geometry.keys()}
    #plane_to_segments_velocity = {pk: [] for pk in plane_geometry.keys()}
    #plane_to_segments_b = {pk: [] for pk in plane_geometry.keys()}
    #plane_to_segments_b_tail = {pk: [] for pk in plane_geometry.keys()}
    plane_to_segments = defaultdict(list)
    plane_to_segments_velocity = defaultdict(list)
    plane_to_segments_b = defaultdict(list)
    plane_to_segments_b_tail = defaultdict(list)
    unassigned = []   # for debugging
    unassigned_velocity = []   # for debugging
    for lkey in DN.networkLinks().keys():
        link = DN.networkLinks().getRef(lkey)
        if hasattr(link, "hasZeroBurgers") and link.hasZeroBurgers():
            continue
        # first entry is loopID
        loopID_sc = node_pos[link.source.networkID()][0]
        loopID_si = node_pos[link.sink.networkID()][0]
        a = node_pos[link.source.networkID()][1:]
        b = node_pos[link.sink.networkID()][1:]
        # midpoint
        mean_pos = (a + b) / 2
        # first entry is loopID
        av = node_velocity[link.source.networkID()][1:]
        bv = node_velocity[link.sink.networkID()][1:]
        link_burger = link.burgers()
        best_pk = None
        best_dist = float("inf")
        for pk, (n_hat, p0) in plane_geometry.items():
            d = segment_plane_distance(a, b, n_hat, p0)
            if d < best_dist:
                best_dist = d
                best_pk = pk
        if best_pk is not None and best_dist < tol:
            #plane_to_segments[best_pk].append((a, b))
            #plane_to_segments_velocity[best_pk].append((av, bv))
            #plane_to_segments_b[best_pk].append(link_burger)
            #plane_to_segments_b_tail[best_pk].append(mean_pos)
            entry = np.asarray(loopID_sc, loopID_si, a, b, av, bv)
            plane_to_segments[best_pk].append(entry)
            plane_to_segments_b[best_pk].append(link_burger)
            plane_to_segments_b_tail[best_pk].append(mean_pos)
        else:
            continue
            #unassigned.append((a, b))
            #unassigned_velocity.append((av, bv))

    print(plane_to_segments)
    exit()
    # convert value to nparray and remove empty list
    plane_to_segments = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments.items() if v}
    plane_to_segments_velocity = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_velocity.items() if v}
    plane_to_segments_b = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_b.items() if v}
    plane_to_segments_b_tail = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_b_tail.items() if v}

    #return plane_to_segments, unassigned
    return plane_to_segments, plane_to_segments_velocity, plane_to_segments_b, plane_to_segments_b_tail


def build_plane_to_burger(DN, node_pos, plane_geometry, tol=1e-4):
    plane_to_segments_b = {pk: [] for pk in plane_geometry.keys()}
    plane_to_segments_b_tail = {pk: [] for pk in plane_geometry.keys()}
    for lkey in DN.networkLinks().keys():
        link = DN.networkLinks().getRef(lkey)
        if hasattr(link, "hasZeroBurgers") and link.hasZeroBurgers():
            continue
        a = node_pos[link.source.networkID()]
        b = node_pos[link.sink.networkID()]
        # midpoint
        mean_pos = (a + b) / 2
        link_burger = link.burgers()
        best_pk = None
        best_dist = float("inf")
        for pk, (n_hat, p0) in plane_geometry.items():
            d = segment_plane_distance(a, b, n_hat, p0)
            if d < best_dist:
                best_dist = d
                best_pk = pk
        if best_pk is not None and best_dist < tol:
            plane_to_segments_b[best_pk].append(link_burger)
            plane_to_segments_b_tail[best_pk].append(mean_pos)
        else:
            continue

    # convert value to nparray and remove empty list
    plane_to_segments_b = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_b.items() if v}
    plane_to_segments_b_tail = {k: np.asarray(v, dtype=np.float32) for k,v in plane_to_segments_b_tail.items() if v}

    return plane_to_segments_b, plane_to_segments_b_tail


def segment_plane_distance(a, b, n_hat, p0):
    da = np.dot(a - p0, n_hat)
    db = np.dot(b - p0, n_hat)
    return max(abs(da), abs(db))

def stable_plane_key(gp):
    h, k, l, layer, latID = gp.key_tuple()
    return (h, k, l, layer, latID)

def draw_box(ax, xmin, xmax, lw=1.0, alpha=0.25):
    X0, X1 = np.array(xmin), np.array(xmax)
    xs = [X0[0], X1[0]]
    ys = [X0[1], X1[1]]
    zs = [X0[2], X1[2]]
    edges = [    # 12 edges
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
        ax.plot(ex, ey, ez, linewidth=lw, alpha=alpha, color='k')


# ========= 2D plane projection helpers =========
def plane_basis_from_points(P):
    """Orthonormal (u,v,n) via PCA on points P (Nx3)."""
    P = np.asarray(P, float)
    origin = P.mean(axis=0)
    Q = P - origin
    if Q.shape[0] < 3:
        u = np.array([1.0, 0.0, 0.0])
        v = np.array([0.0, 1.0, 0.0])
        n = np.cross(u, v)
        return origin, u, v, n / np.linalg.norm(n)
    _, _, Vt = np.linalg.svd(Q, full_matrices=False)
    u = Vt[0]; v = Vt[1]
    u /= np.linalg.norm(u); v /= np.linalg.norm(v)
    n = np.cross(u, v); n /= np.linalg.norm(n)
    return origin, u, v, n


#def orient_basis(P, origin, u, v, n):
#    """
#    Return a consistently-oriented (u,v,n) frame.
#    Sometimes the projected vector gets fliped 180 degrees because
#    the SVD gives the two dominant directions, but their sign is 
#    arbitrary—u or v can come out pointing the opposite way on different runs.
#    this function fixes the issue by pinning the orientation to something 
#    deterministic, 
#    e.g. always keep the first vertex on the “positive” side of each axis
#    """
#    # pick a reference point (first vertex after centroid)
#    ref = np.asarray(P[0], float) - origin
#    # flip u if ref is on the negative side
#    if np.dot(ref, u) < 0:
#        u = -u
#    # flip v if ref is on the negative side
#    if np.dot(ref, v) < 0:
#        v = -v
#    # rebuild right-handed normal
#    n = np.cross(u, v)
#    n /= np.linalg.norm(n)
#    return u, v, n


def project(P, origin, u, v):
    d = np.asarray(P, np.float32) - origin
    return np.dot(d, u), np.dot(d, v)


def project_b(v_burgers, origin, u, v):
    d_burger = np.asarray(v_burgers, np.float32) - origin
    return np.dot(d_burger, u)

