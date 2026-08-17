# gpu_compare.py -- per-frame CPU/GPU VTK parity comparison (GPU port Task 4+)
#
# usage:
#   python gpu_compare.py <dir_cpu> <dir_gpu> <n_particles> \
#        [--pos-tol 1e-6] [--energy-tol 1e-9] [--g 9.81] [--mass 0.422] [--iavg 7e-4]
#
# For every frame present in BOTH dirs (matched by the step number in the
# file name, prefix_000010.vtk):
#   - parse POINTS / POLYGONS / SCALARS id / VECTORS velocity,omega
#   - per particle: COM from the mesh points of its cells, velocity and omega
#     from its first cell
#   - compare max |dCOM| over particles and energy difference
#       E = sum( 0.5*M*|v|^2 + 0.5*IAVG*|w|^2 + M*g*z_com )
#     (same constants as pile_analyze.py: M=0.422, IAVG=7e-4 per banana)
# Prints the per-frame table, the max deviations over all frames and
# PASS/FAIL vs the tolerances. Exit code 1 on FAIL or on any parse/mismatch
# error (missing frames, particle count mismatch, differing frame sets).

import argparse
import math
import os
import re
import sys

STEP_RE = re.compile(r'_(\d+)\.vtk$')


def frame_steps(d):
    out = {}
    if not os.path.isdir(d):
        sys.exit("ERROR: not a directory: %s" % d)
    for f in os.listdir(d):
        m = STEP_RE.search(f)
        if m and f.endswith('.vtk'):
            out[int(m.group(1))] = os.path.join(d, f)
    return out


class Frame:
    __slots__ = ('points', 'cells', 'cid', 'vel', 'omega')


def parse_vtk(path):
    fr = Frame()
    fr.points = []
    fr.cells = []
    fr.cid = []
    fr.vel = []
    fr.omega = []
    ncell = None
    with open(path) as f:
        lines = f.readlines()
    i = 0
    n = len(lines)
    while i < n:
        ln = lines[i].strip()
        if ln.startswith('POINTS '):
            cnt = int(ln.split()[1])
            i += 1
            for j in range(cnt):
                p = lines[i + j].split()
                fr.points.append((float(p[0]), float(p[1]), float(p[2])))
            i += cnt
            continue
        if ln.startswith('POLYGONS '):
            cnt = int(ln.split()[1])
            i += 1
            for j in range(cnt):
                c = lines[i + j].split()
                fr.cells.append((int(c[1]), int(c[2]), int(c[3])))
            i += cnt
            ncell = cnt
            continue
        if ln.startswith('SCALARS id'):
            i += 2  # skip LOOKUP_TABLE default
            for j in range(ncell):
                fr.cid.append(int(lines[i + j]))
            i += ncell
            continue
        if ln.startswith('VECTORS velocity'):
            i += 1
            for j in range(ncell):
                t = lines[i + j].split()
                fr.vel.append((float(t[0]), float(t[1]), float(t[2])))
            i += ncell
            continue
        if ln.startswith('VECTORS omega'):
            i += 1
            for j in range(ncell):
                t = lines[i + j].split()
                fr.omega.append((float(t[0]), float(t[1]), float(t[2])))
            i += ncell
            continue
        i += 1
    if not fr.points or not fr.cells or not fr.cid or not fr.vel or not fr.omega:
        sys.exit("ERROR: incomplete VTK (missing section): %s" % path)
    return fr


def particle_states(fr, npart):
    # Segment cells by the id SCALAR (robust to mixed meshes), gather each
    # particle's point set from its polygon indices.
    per_p = [[] for _ in range(npart)]
    for c, pid in zip(fr.cells, fr.cid):
        if not (0 <= pid < npart):
            sys.exit("ERROR: cell id %d out of range (n=%d)" % (pid, npart))
        per_p[pid].append(c)
    states = []
    for pid in range(npart):
        cells = per_p[pid]
        if not cells:
            sys.exit("ERROR: particle %d has no cells" % pid)
        idxs = set()
        for c in cells:
            idxs.update(c)
        pts = [fr.points[k] for k in sorted(idxs)]
        com = tuple(sum(p[q] for p in pts) / len(pts) for q in range(3))
        first = next(k for k, c_id in enumerate(fr.cid) if c_id == pid)
        states.append((com, fr.vel[first], fr.omega[first]))
    return states


def energy(states, mass, iavg, g):
    e = 0.0
    for com, v, w in states:
        e += 0.5 * mass * sum(c * c for c in v)
        e += 0.5 * iavg * sum(c * c for c in w)
        e += mass * g * com[2]
    return e


def dist(a, b):
    return math.sqrt(sum((a[q] - b[q]) ** 2 for q in range(3)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir_cpu')
    ap.add_argument('dir_gpu')
    ap.add_argument('n_particles', type=int)
    ap.add_argument('--pos-tol', type=float, default=1e-6)
    ap.add_argument('--energy-tol', type=float, default=1e-9)
    ap.add_argument('--g', type=float, default=9.81)
    ap.add_argument('--mass', type=float, default=0.422)
    ap.add_argument('--iavg', type=float, default=7.0e-4)
    args = ap.parse_args()

    cpu = frame_steps(args.dir_cpu)
    gpu = frame_steps(args.dir_gpu)
    if not cpu:
        sys.exit("ERROR: no VTK frames in %s" % args.dir_cpu)
    if not gpu:
        sys.exit("ERROR: no VTK frames in %s (GPU run missing/failed)" % args.dir_gpu)
    common = sorted(set(cpu) & set(gpu))
    if not common:
        sys.exit("ERROR: no common frame steps (cpu=%s gpu=%s)"
                 % (sorted(cpu), sorted(gpu)))
    only_cpu = sorted(set(cpu) - set(gpu))
    only_gpu = sorted(set(gpu) - set(cpu))
    if only_cpu or only_gpu:
        print("FAIL: frame set mismatch: only_cpu=%s only_gpu=%s"
              % (only_cpu, only_gpu))
        sys.exit(1)

    print("frames: cpu=%d gpu=%d common=%d%s%s"
          % (len(cpu), len(gpu), len(common),
             (" only_cpu=%s" % only_cpu) if only_cpu else "",
             (" only_gpu=%s" % only_gpu) if only_gpu else ""))
    print("tolerances: pos=%.3g energy=%.3g" % (args.pos_tol, args.energy_tol))
    print("%8s %12s %12s %12s %12s" % ("step", "max|dpos|", "dE", "E_cpu", "E_gpu"))

    max_dpos = 0.0
    max_dE = 0.0
    for s in common:
        sc = particle_states(parse_vtk(cpu[s]), args.n_particles)
        sg = particle_states(parse_vtk(gpu[s]), args.n_particles)
        dpos = max(dist(a[0], b[0]) for a, b in zip(sc, sg))
        Ec = energy(sc, args.mass, args.iavg, args.g)
        Eg = energy(sg, args.mass, args.iavg, args.g)
        dE = abs(Ec - Eg)
        max_dpos = max(max_dpos, dpos)
        max_dE = max(max_dE, dE)
        print("%8d %12.3e %12.3e %12.6f %12.6f" % (s, dpos, dE, Ec, Eg))

    print("\nmax over %d frames: |dpos|=%.3e  |dE|=%.3e" % (len(common), max_dpos, max_dE))
    ok = max_dpos <= args.pos_tol and max_dE <= args.energy_tol
    print("RESULT: %s (pos %.3e <= %.3g, energy %.3e <= %.3g)"
          % ("PASS" if ok else "FAIL", max_dpos, args.pos_tol,
             max_dE, args.energy_tol))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
