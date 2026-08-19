#pragma once

#include "core/vec3.hpp"
#include "geometry/mesh.hpp"
#include "host/vtk_io.hpp"

using namespace zdem;

// One pair: v8 containment detection + per-vertex penalty force. Returns the
// total contained-vertex count n_inc (0 = no contact).
// f_i/t_i/f_j/t_j are outputs (zeroed then accumulated): the owner=i half and
// the other=j half. Caller's choice: CPU adds both halves; MPI Newton-off
// takes only the own half.
// Note: the per-vertex order inside the blocks matches the original main.cpp
// interleaving (all incA first, then all incB, normal then tangential per
// vertex), but adding the blocks back to forces[] regroups the FP evaluation
// vs the original per-vertex direct accumulation -- gated by the
// scratch/mpi_v0_* and scratch/pp_v0_* bit-parity baselines.
int pp_contact_pair(const Particle& pa, const Particle& pb,
                    const Mesh& ma, const Mesh& mb,
                    Vec3& f_i, Vec3& t_i, Vec3& f_j, Vec3& t_j,
                    double tangential_damping);

// Cached-triangles overload: trisA/trisB are transform_tris(ma, pa.tf) /
// transform_tris(mb, pb.tf) computed ONCE per particle per step by the
// caller. transform_tris is a pure function of (mesh, tf) and tf is frozen
// within a step, so the cached values are bit-identical to the per-pair
// recomputation -- per-pair numerics and accumulation order unchanged.
// The Mesh args stay: the containment scan walks ma.vertices/mb.vertices
// in mesh order, which is NOT reproducible from the triangle list.
int pp_contact_pair(const Particle& pa, const Particle& pb,
                    const Mesh& ma, const Mesh& mb,
                    const std::vector<std::array<Vec3, 3>>& trisA,
                    const std::vector<std::array<Vec3, 3>>& trisB,
                    Vec3& f_i, Vec3& t_i, Vec3& f_j, Vec3& t_j,
                    double tangential_damping);
