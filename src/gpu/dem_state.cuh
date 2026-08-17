// src/gpu/dem_state.cuh -- device state layout + GpuSim (GPU port Task 3)
//
// Device buffers for the three static-after-upload registries:
//   DeviceMeshes    body-frame mesh registry, flattened (vertices + tris as
//                   GLOBAL vertex indices; the CPU Mesh stores triangles as
//                   positions, the conversion happens at upload)
//   DeviceWalls     static world-space coplanar wall groups (one-time
//                   grouping at upload; full footprints, no AABB filter --
//                   the kernel filters geometrically via point_in_tri)
//   DeviceParticles SoA particle state
// Host double -> real conversion happens exactly at the upload copy;
// downloads convert back to double.
#pragma once
#include <vector>
#include <cuda_runtime.h>

#include "gpu/real.hpp"
#include "host/config_io.hpp"
#include "host/vtk_io.hpp"  // Particle (+ Mesh via geometry/mesh.hpp)

struct DeviceMeshes {          // flattened body-frame mesh registry
    // vertices: 3*total_v floats, tris: 3*total_t ints (vertex indices, GLOBAL)
    real* d_verts = nullptr; int* d_tris = nullptr;
    int* d_voffset = nullptr;            // [n_mesh+1] vertex offset per mesh
    int* d_toffset = nullptr;            // [n_mesh+1] tri offset per mesh
    real* d_mean_edge = nullptr;         // [n_mesh]
    int n_mesh = 0;
};

struct DeviceWalls {           // static world-space coplanar groups
    real* d_gn = nullptr;      // [3*ng] group normal (push side = particle side)
    real* d_gd = nullptr;      // [ng] s = dot(n,x)+d; kernel must handle EITHER sign (push along n*sign(s))
    real* d_fp = nullptr;      // [9*n_fp] footprints (3 verts each), row-major
    int* d_fp_start = nullptr; // [ng+1] footprint offset per group
    real* d_mu = nullptr;      // [ng] friction (min particle/wall applied in kernel)
    int n_groups = 0;
};

struct DeviceParticles {       // SoA state
    real *pos = nullptr, *vel = nullptr, *omega = nullptr, *quat = nullptr, *L = nullptr;  // [3N],[3N],[3N],[4N],[3N]
    real *mass = nullptr, *inv_mass = nullptr, *radius = nullptr, *equiv_radius = nullptr;
    real *inertia_body_inv = nullptr;    // [9N] row-major body frame
    real *young = nullptr, *poisson = nullptr, *mu = nullptr, *restitution = nullptr;
    int *mesh_index = nullptr;
    real *force = nullptr, *torque = nullptr;  // [3N] zeroed each step
    int *contact_count = nullptr;              // [N]
    int *contacts = nullptr;                   // [1] per-step contact counter (wall groups now, pp pairs in Task 7)
    int n = 0;
};

class GpuSim {
public:
    DeviceParticles P; DeviceMeshes M; DeviceWalls W;
    real dt = 0; real gravity[3] = {0, 0, 0}; real tangential_damping = 0;
    real* d_gravity = nullptr;          // [3] device copy of gravity
    // Host registry copies kept from upload() for periodic VTK output: the
    // mesh registry (body frame; write_vtk_particles rotates to world with
    // the per-particle quaternion) and the initial Particle array (static
    // fields mass/radius/mesh_index; dynamic fields overwritten per output
    // from download_frame before writing).
    std::vector<Mesh> host_meshes;
    std::vector<Particle> host_particles;
    // Task-3 checksum echo: sums of the exact staged (real-typed) values that
    // were uploaded, computed on the host right before the device copy. The
    // --check-sums mode compares these against a device roundtrip; identical
    // summation order makes the pair bit-exact in both precision modes.
    double host_pos_sum = 0.0;
    double host_mass_sum = 0.0;
    void upload(const SimConfig& cfg);           // builds meshes/particles on host (same code path as CPU main), uploads
    void free_all();
    void download_soa(std::vector<double>& pos, std::vector<double>& vel,
                      std::vector<double>& omega) const;  // for output/verification
    // Full frame download for VTK output: pos/vel/omega/quat/force/torque
    // (real -> double) + contact_count.
    void download_frame(std::vector<double>& pos, std::vector<double>& vel,
                        std::vector<double>& omega, std::vector<double>& quat,
                        std::vector<double>& force, std::vector<double>& torque,
                        std::vector<int>& contact_count) const;
};
