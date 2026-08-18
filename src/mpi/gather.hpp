#pragma once

#include <mpi.h>

#include <vector>

#include "core/vec3.hpp"
#include "host/config_io.hpp"
#include "host/sim_build.hpp"
#include "host/vtk_io.hpp"
#include "mpi/decomp.hpp"

// One particle's output record, shipped with MPI_Gatherv as raw bytes
// (MPI_BYTE, homogeneous-cluster assumption). Particle is POD
// (Transform/Vec3/Mat3 only, no strings), so sizeof(FrameSnapshot) is the
// fixed byte stride of the packed array.
struct FrameSnapshot {
    Particle p;
    int gid;
    int cc;
    Vec3 force;
    Vec3 torque;
};

// Gather every rank's local particles on rank 0, sort by gid (== config
// order, so the VTK id field matches the CPU driver), and write one frame
// with the CPU file naming (vtk_prefix_%06d.vtk under cfg.output_dir).
// Non-root ranks only ship their data.
void gather_write_frame(const Decomp& d, const SimConfig& cfg, const SimBuild& sim,
                        const std::vector<Particle>& local, const std::vector<int>& gids,
                        const std::vector<int>& cc, const std::vector<Vec3>& forces,
                        const std::vector<Vec3>& torques, int step);

// MPI_Reduce SUM over d.cart. Rank 0 returns the total, other ranks 0.
int reduce_add(const Decomp& d, int local_val);

// Rank 0 prints "rank_n=a,b,c": the entries of every rank's
// local_sizes_hint collected with Allgatherv (consumed by Task 5+).
void log_rank_n(const Decomp& d, const std::vector<int>& local_sizes_hint);
