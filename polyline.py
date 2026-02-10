#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Route-B (Feng 2021) for closed STL contacts, with robust loop reconstruction.

Pipeline (per pair A,B):
1) Triangle-triangle intersections -> raw segment soup (p0,p1)
2) Optional coarse split into contact components (greedy segment distance BFS)
3) Within each component:
   a) Endpoint snap: merge endpoints within tol_snap
   b) T-junction split: if a snapped vertex lies on interior of an edge within tol_snap,
      split that edge at the projected point (iterative single pass)
   c) Loop tracing: build undirected graph; extract closed polylines by walking edges
      (prefer vertices with degree==2; still works if extra branches exist)
4) Compute Sn, Gn, nA/nB, xc(0) from each closed polyline loop
5) Output VTK:
   - A_transformed.vtk, B_transformed.vtk : transformed meshes
   - routeB_intersection.vtk : reconstructed loops as polyline LINES
     plus contact points (xc0) as VERTICES
     with CELL_DATA contact_id for all cells

Deps: numpy only (C-friendly).
"""

import os
import struct
import argparse
import numpy as np

EPS = 1e-12


# ----------------------------
# Quaternion / transform
# ----------------------------
def quat_normalize(q):
    q = np.asarray(q, dtype=float)
    n = np.linalg.norm(q)
    if n < 1e-30:
        return np.array([1.0, 0.0, 0.0, 0.0])
    return q / n

def quat_rotate(q, v):
    q = quat_normalize(q)
    w, x, y, z = q
    v = np.asarray(v, dtype=float)
    qv = np.array([x, y, z], dtype=float)
    t = 2.0 * np.cross(qv, v)
    return v + w * t + np.cross(qv, t)

def transform_points(pts, t, q):
    pts = np.asarray(pts, dtype=float)
    return quat_rotate(q, pts) + np.asarray(t, dtype=float)


# ----------------------------
# STL loading (binary/ascii)
# ----------------------------
def load_stl(path):
    with open(path, "rb") as f:
        _ = f.read(80)
        n_tri_bytes = f.read(4)
        if len(n_tri_bytes) < 4:
            raise ValueError("Invalid STL.")
        n_tri = struct.unpack("<I", n_tri_bytes)[0]
        file_size = os.path.getsize(path)
        expected = 84 + n_tri * 50
        if expected == file_size:
            f.seek(84)
            tris = np.empty((n_tri, 3, 3), dtype=float)
            for i in range(n_tri):
                rec = f.read(50)
                if len(rec) != 50:
                    raise ValueError("Truncated binary STL.")
                data = struct.unpack("<12fH", rec)
                tris[i, 0] = data[3:6]
                tris[i, 1] = data[6:9]
                tris[i, 2] = data[9:12]
            return tris

    # ASCII fallback
    tris = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        vs = []
        for line in f:
            s = line.strip().lower()
            if s.startswith("vertex"):
                parts = s.split()
                if len(parts) >= 4:
                    vs.append([float(parts[1]), float(parts[2]), float(parts[3])])
            if s.startswith("endloop"):
                if len(vs) >= 3:
                    tris.append(vs[-3:])
    if not tris:
        raise ValueError("Failed to parse STL (neither binary nor ascii).")
    return np.array(tris, dtype=float)


def tri_normals(tris):
    e1 = tris[:, 1] - tris[:, 0]
    e2 = tris[:, 2] - tris[:, 0]
    n = np.cross(e1, e2)
    ln = np.linalg.norm(n, axis=1)
    good = ln > 1e-30
    n[good] /= ln[good][:, None]
    return n


# ----------------------------
# Geometry helpers
# ----------------------------
def plane_from_tri(tri):
    n = np.cross(tri[1] - tri[0], tri[2] - tri[0])
    ln = np.linalg.norm(n)
    if ln < 1e-30:
        return None, None
    n = n / ln
    d = -np.dot(n, tri[0])
    return n, d

def point_in_tri(p, tri, n):
    a, b, c = tri
    ab = b - a; bc = c - b; ca = a - c
    ap = p - a; bp = p - b; cp = p - c
    c1 = np.dot(np.cross(ab, ap), n)
    c2 = np.dot(np.cross(bc, bp), n)
    c3 = np.dot(np.cross(ca, cp), n)
    return (c1 >= -1e-10) and (c2 >= -1e-10) and (c3 >= -1e-10)

def seg_plane_intersection(p0, p1, n, d):
    u = p1 - p0
    denom = np.dot(n, u)
    num = -(np.dot(n, p0) + d)
    if abs(denom) < 1e-30:
        return False, None
    t = num / denom
    if t < -1e-10 or t > 1.0 + 1e-10:
        return False, None
    t = min(1.0, max(0.0, t))
    return True, p0 + t * u

def unique_points(pts, tol=1e-9):
    out = []
    for p in pts:
        keep = True
        for q in out:
            if np.linalg.norm(p - q) <= tol:
                keep = False
                break
        if keep:
            out.append(p)
    return out

def tri_tri_intersection_segment(triA, nA, dA, triB, nB, dB):
    tau = np.cross(nA, nB)
    lt = np.linalg.norm(tau)
    if lt < 1e-20:
        return None
    tau = tau / lt

    candidates = []
    for (i, j) in [(0,1),(1,2),(2,0)]:
        hit, p = seg_plane_intersection(triA[i], triA[j], nB, dB)
        if hit and point_in_tri(p, triB, nB):
            candidates.append(p)

    for (i, j) in [(0,1),(1,2),(2,0)]:
        hit, p = seg_plane_intersection(triB[i], triB[j], nA, dA)
        if hit and point_in_tri(p, triA, nA):
            candidates.append(p)

    candidates = unique_points(candidates, tol=1e-9)
    if len(candidates) < 2:
        return None

    ts = [np.dot(p, tau) for p in candidates]
    i0 = int(np.argmin(ts))
    i1 = int(np.argmax(ts))
    if i0 == i1:
        return None
    p0 = candidates[i0]
    p1 = candidates[i1]
    if np.linalg.norm(p1 - p0) < 1e-18:
        return None
    if np.dot(p1 - p0, tau) < 0:
        p0, p1 = p1, p0
    return p0, p1


# ----------------------------
# Feng accumulators from closed polyline
# ----------------------------
def accumulate_Sn_Gn_from_polyline(loop_pts):
    """
    loop_pts: list/array of shape (M,3), closed (first!=last is ok; we'll wrap)
    Uses same discrete formulas on each edge (p_i, p_{i+1}).
    """
    P = np.asarray(loop_pts, float)
    m = len(P)
    Sn = np.zeros(3, float)
    Gn = np.zeros(3, float)
    for i in range(m):
        p0 = P[i]
        p1 = P[(i+1) % m]
        Sn += 0.5 * np.cross(p0, p1)
        dx = p1 - p0
        s = np.dot(p0, p1) + (1.0/3.0) * np.dot(dx, dx)
        Gn += -(1.0/3.0) * s * dx
    return Sn, Gn

def contact_point_xc0(Sn, Gn):
    a = float(np.linalg.norm(Sn))
    if a < 1e-14:
        return None, np.zeros(3), np.zeros(3), 0.0
    nA = -Sn / a
    nB = -nA
    xc0 = np.cross(nA, Gn) / a
    return xc0, nA, nB, a


# ----------------------------
# Coarse split into contact components (greedy BFS using segment-segment distance)
# ----------------------------
def _cell_index(p, h):
    return (int(np.floor(p[0]/h)), int(np.floor(p[1]/h)), int(np.floor(p[2]/h)))

def seg_seg_dist2(p1, q1, p2, q2):
    p1 = np.asarray(p1, float); q1 = np.asarray(q1, float)
    p2 = np.asarray(p2, float); q2 = np.asarray(q2, float)
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    a = float(np.dot(d1, d1))
    e = float(np.dot(d2, d2))
    f = float(np.dot(d2, r))

    if a <= EPS and e <= EPS:
        return float(np.dot(p1 - p2, p1 - p2))
    if a <= EPS:
        t = np.clip(f / e, 0.0, 1.0)
        c = p2 + t * d2
        return float(np.dot(p1 - c, p1 - c))
    if e <= EPS:
        s = np.clip(-np.dot(d1, r) / a, 0.0, 1.0)
        c = p1 + s * d1
        return float(np.dot(c - p2, c - p2))

    b = float(np.dot(d1, d2))
    c = float(np.dot(d1, r))
    denom = a * e - b * b
    if abs(denom) > 0.0:
        s = np.clip((b * f - c * e) / denom, 0.0, 1.0)
    else:
        s = 0.0
    t = (b * s + f) / e
    if t < 0.0:
        t = 0.0
        s = np.clip(-c / a, 0.0, 1.0)
    elif t > 1.0:
        t = 1.0
        s = np.clip((b - c) / a, 0.0, 1.0)
    cp1 = p1 + s * d1
    cp2 = p2 + t * d2
    diff = cp1 - cp2
    return float(np.dot(diff, diff))

def label_segments_by_contact_greedy(segments, tol=1e-6):
    h = float(tol)
    nseg = len(segments)
    if nseg == 0:
        return [], 0, []

    grid = {}
    seg_aabb = []
    rep_points = []

    for i, (p0, p1) in enumerate(segments):
        p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
        mid = 0.5 * (p0 + p1)
        rep_points.append((p0, p1, mid))
        mn = np.minimum(p0, p1) - h
        mx = np.maximum(p0, p1) + h
        seg_aabb.append((mn, mx))
        for rp in (p0, p1, mid):
            c = _cell_index(rp, h)
            grid.setdefault(c, []).append(i)

    def nearby_candidates(i):
        cand = set()
        for rp in rep_points[i]:
            c = _cell_index(rp, h)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        cc = (c[0]+dx, c[1]+dy, c[2]+dz)
                        for j in grid.get(cc, []):
                            if j != i:
                                cand.add(j)
        return cand

    def aabb_overlap(i, j):
        mn1, mx1 = seg_aabb[i]
        mn2, mx2 = seg_aabb[j]
        return not (mx1[0] < mn2[0] or mn1[0] > mx2[0] or
                    mx1[1] < mn2[1] or mn1[1] > mx2[1] or
                    mx1[2] < mn2[2] or mn1[2] > mx2[2])

    tol2 = h * h
    visited = [False] * nseg
    comp_ids = [-1] * nseg
    comps = []

    for seed in range(nseg):
        if visited[seed]:
            continue
        cid = len(comps)
        queue = [seed]
        visited[seed] = True
        comp_ids[seed] = cid
        comp = [seed]
        while queue:
            i = queue.pop()
            p1, q1 = segments[i]
            for j in nearby_candidates(i):
                if visited[j]:
                    continue
                if not aabb_overlap(i, j):
                    continue
                p2, q2 = segments[j]
                if seg_seg_dist2(p1, q1, p2, q2) <= tol2:
                    visited[j] = True
                    comp_ids[j] = cid
                    queue.append(j)
                    comp.append(j)
        comps.append(comp)

    order = sorted(range(len(comps)), key=lambda k: len(comps[k]), reverse=True)
    remap = {old: new for new, old in enumerate(order)}
    comps_sorted = [comps[old] for old in order]
    comp_ids = [remap[c] for c in comp_ids]
    return comp_ids, len(comps_sorted), comps_sorted


# ----------------------------
# Snap endpoints + T-junction split + loop tracing
# ----------------------------
class UnionFind:
    def __init__(self, n=0):
        self.parent = list(range(n))
        self.rank = [0]*n

    def find(self, x):
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a, b):
        ra = self.find(a); rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1

def snap_endpoints(segments, tol):
    """Return (V, edges) where V are snapped vertex positions, edges are pairs of vertex indices."""
    h = float(tol)
    nseg = len(segments)
    pts = np.empty((2*nseg, 3), float)
    for i, (p0, p1) in enumerate(segments):
        pts[2*i] = p0
        pts[2*i+1] = p1

    uf = UnionFind(2*nseg)
    grid = {}  # cell -> list[point_index]

    for i in range(2*nseg):
        p = pts[i]
        c = _cell_index(p, h)
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    cc = (c[0]+dx, c[1]+dy, c[2]+dz)
                    for j in grid.get(cc, []):
                        if np.linalg.norm(p - pts[j]) <= h:
                            uf.union(i, j)
        grid.setdefault(c, []).append(i)

    # map roots to vertex id, with averaged position
    root_to_vid = {}
    sums = {}
    counts = {}
    for i in range(2*nseg):
        r = uf.find(i)
        if r not in root_to_vid:
            root_to_vid[r] = len(root_to_vid)
            sums[r] = pts[i].copy()
            counts[r] = 1
        else:
            sums[r] += pts[i]
            counts[r] += 1

    V = np.zeros((len(root_to_vid), 3), float)
    for r, vid in root_to_vid.items():
        V[vid] = sums[r] / counts[r]

    edges = []
    for s in range(nseg):
        a = root_to_vid[uf.find(2*s)]
        b = root_to_vid[uf.find(2*s+1)]
        if a != b:
            edges.append((a, b))
    return V, edges

def point_segment_dist2_and_t(p, a, b):
    """Return (dist2, t) where proj = a + t*(b-a), t in [0,1]."""
    p = np.asarray(p, float); a = np.asarray(a, float); b = np.asarray(b, float)
    ab = b - a
    denom = float(np.dot(ab, ab))
    if denom <= EPS:
        return float(np.dot(p-a, p-a)), 0.0
    t = float(np.dot(p-a, ab) / denom)
    t_clamped = min(1.0, max(0.0, t))
    proj = a + t_clamped * ab
    d2 = float(np.dot(p-proj, p-proj))
    return d2, t_clamped

def split_t_junctions(V, edges, tol):
    """
    If any vertex lies on interior of an edge within tol, split the edge.
    Returns updated (V, edges).
    """
    h = float(tol)
    tol2 = h*h
    # spatial hash for vertices
    vgrid = {}
    for vid, p in enumerate(V):
        c = _cell_index(p, h)
        vgrid.setdefault(c, []).append(vid)

    new_edges = []
    for (a, b) in edges:
        pa = V[a]; pb = V[b]
        mn = np.minimum(pa, pb) - h
        mx = np.maximum(pa, pb) + h
        # collect candidate vertices in bbox via grid
        cmin = _cell_index(mn, h)
        cmax = _cell_index(mx, h)
        splits = []  # (t, vid)
        for ix in range(cmin[0], cmax[0]+1):
            for iy in range(cmin[1], cmax[1]+1):
                for iz in range(cmin[2], cmax[2]+1):
                    for vid in vgrid.get((ix,iy,iz), []):
                        if vid == a or vid == b:
                            continue
                        d2, t = point_segment_dist2_and_t(V[vid], pa, pb)
                        # interior only (avoid snapping to endpoints)
                        if d2 <= tol2 and t > 1e-6 and t < 1.0-1e-6:
                            splits.append((t, vid))
        if not splits:
            new_edges.append((a, b))
            continue
        # sort by t, unique by vid
        splits.sort(key=lambda x: x[0])
        uniq = []
        seen = set()
        for t, vid in splits:
            if vid in seen:
                continue
            seen.add(vid)
            uniq.append((t, vid))
        chain = [a] + [vid for _, vid in uniq] + [b]
        for i in range(len(chain)-1):
            u = chain[i]; v = chain[i+1]
            if u != v:
                new_edges.append((u, v))
    return V, new_edges

def build_adjacency(nv, edges):
    adj = [[] for _ in range(nv)]
    for ei, (a,b) in enumerate(edges):
        adj[a].append((b, ei))
        adj[b].append((a, ei))
    return adj

def extract_loops(V, edges):
    """
    Extract closed polylines (loops) from an undirected graph.
    We walk unused edges; if we return to the start vertex, we accept as loop.
    """
    nv = len(V)
    adj = build_adjacency(nv, edges)
    used = [False]*len(edges)
    loops = []

    # try starting from degree==2 vertices first
    start_vertices = list(range(nv))
    start_vertices.sort(key=lambda i: (len(adj[i]) != 2, -len(adj[i])))

    for sv in start_vertices:
        for nb, ei in adj[sv]:
            if used[ei]:
                continue
            # start walking from sv -> nb along ei
            path_vids = [sv, nb]
            used[ei] = True
            prev = sv
            cur = nb
            steps = 0
            while steps < 100000:
                steps += 1
                if cur == sv:
                    # closed; drop last repeated sv
                    loops.append(path_vids[:-1])
                    break
                # choose next unused edge; prefer not going back to prev
                cand = [(n2, e2) for (n2, e2) in adj[cur] if not used[e2]]
                if not cand:
                    # open chain; discard (or keep if you want)
                    break
                # heuristic: pick edge whose neighbor != prev if possible
                pick = None
                for n2, e2 in cand:
                    if n2 != prev:
                        pick = (n2, e2)
                        break
                if pick is None:
                    pick = cand[0]
                n2, e2 = pick
                used[e2] = True
                path_vids.append(n2)
                prev, cur = cur, n2
    return loops


# ----------------------------
# VTK writers
# ----------------------------
def write_vtk_polydata_tris(path, tris, comment="mesh"):
    pts = []
    faces = []

    def add_point(p):
        p = np.asarray(p, float)
        for idx, q in enumerate(pts):
            if np.linalg.norm(p - q) < 1e-9:
                return idx
        pts.append(p.copy())
        return len(pts) - 1

    for tri in tris:
        ids = [add_point(tri[0]), add_point(tri[1]), add_point(tri[2])]
        faces.append(ids)

    with open(path, "w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{comment}\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write(f"POINTS {len(pts)} float\n")
        for p in pts:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")
        f.write(f"POLYGONS {len(faces)} {len(faces)*4}\n")
        for a, b, c in faces:
            f.write(f"3 {a} {b} {c}\n")

def write_vtk_polydata_polylines(path, loops_pts, contact_points, loop_contact_ids):
    """
    loops_pts: list of loops, each loop is list of 3D points (M,3) without repeated last point
    contact_points: list of 3D points (xc0), length = n_loops (or subset)
    loop_contact_ids: list[int], length = n_loops

    VTK:
      - Each loop is one LINE cell with M+1 points (repeat first point at end for closure)
      - Each xc0 is a VERTEX cell
      - CELL_DATA contact_id is provided for all cells: n_lines + n_vertices
    """
    pts = []
    # map point to index (snap exact)
    def add_point(p):
        p = np.asarray(p, float)
        for idx, q in enumerate(pts):
            if np.linalg.norm(p - q) < 1e-9:
                return idx
        pts.append(p.copy())
        return len(pts)-1

    line_cells = []
    for loop in loops_pts:
        ids = [add_point(p) for p in loop]
        if len(ids) < 2:
            continue
        # close for display
        ids.append(ids[0])
        line_cells.append(ids)

    vtx_ids = []
    for cp in contact_points:
        if cp is None:
            continue
        vtx_ids.append(add_point(cp))

    n_lines = len(line_cells)
    n_verts = len(vtx_ids)
    total_cells = n_lines + n_verts

    # Build CELL_DATA contact_id for all cells
    # For vertices: use same id as their loop (assume 1:1, else -1)
    cell_ids = []
    for cid in loop_contact_ids[:n_lines]:
        cell_ids.append(int(cid))
    # vertex ids
    if len(contact_points) == n_lines:
        for cid in loop_contact_ids[:n_lines]:
            cell_ids.append(int(cid))
    else:
        for _ in range(n_verts):
            cell_ids.append(-1)

    with open(path, "w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("reconstructed_contact_loops\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write(f"POINTS {len(pts)} float\n")
        for p in pts:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")

        # LINES section
        total_ints = sum(1 + len(ids) for ids in line_cells)
        f.write(f"LINES {n_lines} {total_ints}\n")
        for ids in line_cells:
            f.write(str(len(ids)) + " " + " ".join(str(i) for i in ids) + "\n")

        if n_verts > 0:
            f.write(f"VERTICES {n_verts} {n_verts*2}\n")
            for vid in vtx_ids:
                f.write(f"1 {vid}\n")

        f.write(f"\nCELL_DATA {total_cells}\n")
        f.write("SCALARS contact_id int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for cid in cell_ids:
            f.write(f"{cid}\n")


# ----------------------------
# DEM mapping
# ----------------------------
def build_contacts_from_loops(idA, idB, loops_pts, contact_ids, k_stiff=1.0):
    """
    Build DEM-style ContactManifold from reconstructed closed loops.

    Paper-consistent *minimal energy* force model (normal only):
      - (Given by Feng 2021 framework)  f = - dPsi/dλ * n
      - (Choice, minimal)              Psi(λ) = (k/2) * Ω(λ)
      - (Geometry from Route-B)        dΩ/dλ = A_c  with  A_c = |S_n|
    => Normal force magnitude:
         Fn = (k/2) * |S_n|   (per contact loop)
    Direction:
         Force on particle A:  F_A = Fn * nA
         Force on particle B:  F_B = -F_A
    """
    contacts = []
    k_stiff = float(k_stiff)

    for k, loop in enumerate(loops_pts):
        Sn, Gn = accumulate_Sn_Gn_from_polyline(loop)
        xc0, nA, nB, area = contact_point_xc0(Sn, Gn)

        # Minimal energy model: Fn = 0.5 * k * area
        Fn = 0.5 * k_stiff * float(area)
        FA = (Fn * np.asarray(nA, float)).tolist()
        FB = (-np.asarray(FA, float)).tolist()

        contacts.append({
            "idA": idA, "idB": idB,
            "cid": int(contact_ids[k]),
            "nVert": int(len(loop)),

            # geometry (Route-B)
            "Sn": Sn.tolist(),
            "Gn": Gn.tolist(),
            "area": float(area),
            "nA": nA.tolist(),
            "nB": nB.tolist(),
            "xc0": None if xc0 is None else xc0.tolist(),

            # mechanics (minimal, paper-consistent)
            "k": k_stiff,
            "Fn": float(Fn),
            "F_A": FA,   # force on A
            "F_B": FB,   # force on B
        })

    return {"idA": idA, "idB": idB, "contacts": contacts}


# ----------------------------
# Main
# ----------------------------
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("stlA")
    ap.add_argument("stlB")
    ap.add_argument("--tA", nargs=3, type=float, default=[0.0, 0.0, 0.0])
    ap.add_argument("--qA", nargs=4, type=float, default=[1.0, 0.0, 0.0, 0.0])
    ap.add_argument("--tB", nargs=3, type=float, default=[0.0, 0.0, 0.0])
    ap.add_argument("--qB", nargs=4, type=float, default=[1.0, 0.0, 0.0, 0.0])

    ap.add_argument("--split-contacts", action="store_true",
                    help="coarsely split into components before loop reconstruction")
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="snap / adjacency tolerance (try 1e-5 if needed)")

    ap.add_argument("--k", type=float, default=1.0,
                    help="contact stiffness parameter k for minimal energy model (paper-consistent): Fn = 0.5*k*|Sn|")

    ap.add_argument("--vtk", default="routeB_intersection.vtk")
    ap.add_argument("--outA", default="A_transformed.vtk")
    ap.add_argument("--outB", default="B_transformed.vtk")
    return ap.parse_args()


def main():
    args = parse_args()

    trisA0 = load_stl(args.stlA)
    trisB0 = load_stl(args.stlB)

    trisA = transform_points(trisA0.reshape(-1, 3), args.tA, args.qA).reshape(-1, 3, 3)
    trisB = transform_points(trisB0.reshape(-1, 3), args.tB, args.qB).reshape(-1, 3, 3)

    write_vtk_polydata_tris(args.outA, trisA, comment="Mesh A (transformed)")
    write_vtk_polydata_tris(args.outB, trisB, comment="Mesh B (transformed)")
    print(f"Transformed meshes written: {args.outA}, {args.outB}")

    nA_all = tri_normals(trisA)
    nB_all = tri_normals(trisB)
    planesA = [plane_from_tri(trisA[i]) for i in range(len(trisA))]
    planesB = [plane_from_tri(trisB[i]) for i in range(len(trisB))]

    aabbA_min = trisA.min(axis=1); aabbA_max = trisA.max(axis=1)
    aabbB_min = trisB.min(axis=1); aabbB_max = trisB.max(axis=1)

    segments = []
    for i in range(len(trisA)):
        nAi, dAi = planesA[i]
        if nAi is None:
            continue
        minA = aabbA_min[i]; maxA = aabbA_max[i]
        for j in range(len(trisB)):
            if (maxA[0] < aabbB_min[j,0] or minA[0] > aabbB_max[j,0] or
                maxA[1] < aabbB_min[j,1] or minA[1] > aabbB_max[j,1] or
                maxA[2] < aabbB_min[j,2] or minA[2] > aabbB_max[j,2]):
                continue
            nBj, dBj = planesB[j]
            if nBj is None:
                continue
            seg = tri_tri_intersection_segment(trisA[i], nA_all[i], dAi, trisB[j], nB_all[j], dBj)
            if seg is not None:
                segments.append(seg)

    if not segments:
        print("No intersection segments found.")
        return

    print(f"Raw intersection segments: {len(segments)}")
    tol = float(args.tol)

    # coarse components
    if args.split_contacts:
        _, ncomp, comps = label_segments_by_contact_greedy(segments, tol=tol)
        print(f"Coarse components: {ncomp} (tol={tol})")
    else:
        comps = [list(range(len(segments)))]
        print("Coarse components: 1 (no split)")

    all_loops_pts = []
    all_loop_contact_ids = []
    all_xc0 = []

    next_cid = 0
    for comp_idx, seg_idx_list in enumerate(comps):
        segs = [segments[i] for i in seg_idx_list]
        # snap endpoints
        V, edges = snap_endpoints(segs, tol=tol)
        # split T-junctions
        V2, edges2 = split_t_junctions(V, edges, tol=tol)
        # extract loops
        loops_vids = extract_loops(V2, edges2)

        # convert to point loops, filter tiny ones
        loops_pts = []
        for lv in loops_vids:
            if len(lv) < 3:
                continue
            pts_loop = [V2[i] for i in lv]
            # reject degenerate tiny loops
            if np.linalg.norm(np.max(pts_loop, axis=0) - np.min(pts_loop, axis=0)) < tol*2:
                continue
            loops_pts.append(pts_loop)

        print(f"\n[Component {comp_idx}] rawSeg={len(segs)} snappedV={len(V)} edges={len(edges2)} loops={len(loops_pts)}")

        # each loop becomes a DEM "Contact" (multi-connected => multiple loops)
        for lp in loops_pts:
            Sn, Gn = accumulate_Sn_Gn_from_polyline(lp)
            xc0, nA, nB, area = contact_point_xc0(Sn, Gn)
            print(f"  [Contact {next_cid}] nVert={len(lp)} |Sn|={area}")
            print(f"    nA={nA.tolist()}")
            print(f"    xc(0)={None if xc0 is None else xc0.tolist()}")

            all_loops_pts.append(lp)
            all_loop_contact_ids.append(next_cid)
            all_xc0.append(xc0)
            next_cid += 1

    # output VTK polylines
    write_vtk_polydata_polylines(args.vtk, all_loops_pts, all_xc0, all_loop_contact_ids)
    print(f"\nReconstructed-loop VTK written: {args.vtk}")
    print('ParaView: color by CELL "contact_id". Each loop is one polyline cell.')

    manifold = build_contacts_from_loops(0, 1, all_loops_pts, all_loop_contact_ids, k_stiff=args.k)
    print("\n=== DEM ContactManifold ===")
    print(f"idA={manifold['idA']} idB={manifold['idB']} nContact={len(manifold['contacts'])}")
    for c in manifold["contacts"]:
        print(f"  [cid={c['cid']}] nVert={c['nVert']} area={c['area']} Fn={c['Fn']}")
        print(f"    nA={c['nA']} xc0={c['xc0']} F_A={c['F_A']}")


if __name__ == "__main__":
    main()
