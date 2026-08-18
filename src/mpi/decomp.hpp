#pragma once

#include <mpi.h>

#include "core/vec3.hpp"
#include "host/config_io.hpp"
#include "host/sim_build.hpp"

// Static 3D brick decomposition of the global simulation box.
// All ranks call make_decomp() with identical inputs; every rank gets its own
// brick [sub_lo, sub_hi) and its cartesian coordinates in the process grid.
struct Decomp {
    double box_lo[3], box_hi[3];     // global box
    int dims[3];                     // px py pz
    int coords[3];                   // cartesian coords of this rank
    MPI_Comm cart = MPI_COMM_NULL;   // 3D topology (reorder=false, ranks match COMM_WORLD)
    int rank = 0, nprocs = 1;
    double sub_lo[3], sub_hi[3];     // this rank's brick
    double ghost_depth = 0.0;        // 2*max_radius
};

// cfg.mpi_box wins; otherwise the initial particle bounding box padded on every
// side by margin (default 5*max_radius, overridable via cfg.mpi_margin).
// dims by greedy splitting: start from {1,1,1}, repeatedly multiply the
// remaining process count onto the axis whose dims[a]*box_len[a] is largest so
// brick footprints stay balanced; then validate per-axis brick thickness >=
// ghost_depth -- a violating axis gets its dims halved and the split retries
// (still violating at dims==1 -> fatal: box too small / too many ranks).
Decomp make_decomp(const SimConfig& cfg, const SimBuild& sim);
int owner_rank(const Decomp& d, const Vec3& pos);   // pos -> global rank (floor into brick; out-of-box clamping is detected by the caller)
bool in_sub(const Decomp& d, const Vec3& pos);      // owned by this rank (left-closed right-open, same floor semantics as owner_rank); positions on/above the outer high boundary clamp into the last brick -- the caller-side out-of-box error catches true escapes
