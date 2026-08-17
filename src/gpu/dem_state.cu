// src/gpu/dem_state.cu -- host-side build + device upload (GPU port Task 3)
//
// upload() replicates the CPU driver (src/main.cpp) construction path:
//   load_stl -> scale -> build_mesh(center_mesh) (compute_mass_properties
//   inside) -> inertia_body / inertia_body_inv (mat3_inverse) ->
//   radius = mesh.radius, equiv_radius = cbrt(3V/4pi) ->
//   L = inertia_world(p) * omega.
// Mesh registry dedupes STLs by "path|scale", same key as the CPU driver.
//
// Walls: one-time coplanar grouping at upload, same quantized (n,d) merge as
// the CPU per-step grouping (pq = 1e6), but per WALL (groups never merge
// across walls, matching the CPU loop nesting) and WITHOUT the per-particle
// AABB pre-filter -- full footprints are uploaded and the kernel filters
// geometrically (point_in_tri). Group normals are oriented toward the FIRST
// particle's position side at upload time. NOTE for Task 5: the wall kernel
// must treat any-sign penetration (s = dot(n,wv)+d, push direction
// n*sign(s)) so a particle that crosses to the other side of the plane still
// gets pushed away from the plane -- equivalent to the CPU per-particle
// orientation for static walls.
//
// Host double -> real conversion happens only at the staging buffers;
// downloads convert back to double. All comments ASCII-only (nvcc + MSVC
// code page safety).

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

#include <cuda_runtime.h>

#include "core/mat3.hpp"
#include "core/quat.hpp"
#include "core/transform.hpp"
#include "core/vec3.hpp"
#include "geometry/mesh.hpp"
#include "geometry/mesh_build.hpp"
#include "geometry/stl_io.hpp"
#include "gpu/cuda_check.hpp"
#include "gpu/dem_state.cuh"
#include "host/vtk_io.hpp"  // Particle

namespace {

// Upload a host array to a fresh device buffer. Returns nullptr for count 0.
template <typename T>
T* upload_array(const T* host_data, std::size_t count) {
    if (count == 0) return nullptr;
    T* d = nullptr;
    CUDA_CHECK(cudaMalloc(&d, count * sizeof(T)));
    CUDA_CHECK(cudaMemcpy(d, host_data, count * sizeof(T), cudaMemcpyHostToDevice));
    return d;
}

// Allocate a zeroed device buffer (force/torque/contact_count accumulators).
template <typename T>
T* upload_zeroed(std::size_t count) {
    if (count == 0) return nullptr;
    T* d = nullptr;
    CUDA_CHECK(cudaMalloc(&d, count * sizeof(T)));
    CUDA_CHECK(cudaMemset(d, 0, count * sizeof(T)));
    return d;
}

// Map a triangle-corner position to its index in mesh.vertices using the
// same first-match tolerance as unique_vertices (dist2 <= 1e-10), so the
// position-based CPU triangles become global vertex indices.
int vertex_index_of(const Mesh& mesh, const Vec3& v) {
    const double tol2 = 1e-10;
    for (std::size_t k = 0; k < mesh.vertices.size(); ++k) {
        Vec3 d = v - mesh.vertices[k];
        if (dot(d, d) <= tol2) return static_cast<int>(k);
    }
    return -1;  // unreachable: every tri corner is a unique_vertices output
}

// World-frame inertia tensor R * I_body * R^T (same as CPU inertia_world).
Mat3 inertia_world(const Particle& p) {
    Mat3 R = quat_to_mat3(p.tf.rot);
    Mat3 Rt = mat3_transpose(R);
    return mat3_mul(mat3_mul(R, p.inertia_body), Rt);
}

}  // namespace

void GpuSim::upload(const SimConfig& cfg) {
    // ================= host build: same code path as CPU main =================
    std::unordered_map<std::string, std::vector<Triangle>> tri_cache;
    std::unordered_map<std::string, int> mesh_cache;
    std::vector<Mesh> meshes;

    auto load_mesh_for = [&](const std::string& path, double scale) -> int {
        std::string key = path + "|" + std::to_string(scale);
        auto it_cache = mesh_cache.find(key);
        if (it_cache != mesh_cache.end()) {
            return it_cache->second;
        }
        auto it = tri_cache.find(path);
        if (it == tri_cache.end()) {
            std::vector<Triangle> tris;
            if (!load_stl(path, tris)) {
                return -1;
            }
            it = tri_cache.emplace(path, std::move(tris)).first;
        }
        std::vector<Triangle> tris = it->second;
        for (auto& tri : tris) {
            tri[0] *= scale;
            tri[1] *= scale;
            tri[2] *= scale;
        }
        Mesh m = build_mesh(tris, cfg.center_mesh);
        meshes.push_back(m);
        int idx = static_cast<int>(meshes.size() - 1);
        mesh_cache[key] = idx;
        return idx;
    };

    std::vector<ParticleInit> inits;
    if (!cfg.particle_inits.empty()) {
        inits = cfg.particle_inits;
    } else {
        inits.resize(cfg.n);
        for (int i = 0; i < cfg.n; ++i) {
            inits[i].stl_path = cfg.stl_path;
            inits[i].pos = Vec3{cfg.spacing * i, 0.0, 0.0};
            inits[i].vel = (i == 0) ? cfg.v0 : (i == 1 ? cfg.v1 : Vec3{0.0, 0.0, 0.0});
            inits[i].rot = Quat{};
            inits[i].omega = Vec3{0.0, 0.0, 0.0};
            inits[i].scale = 1.0;
            inits[i].density = 1.0;
            inits[i].young = 1e7;
            inits[i].poisson = 0.25;
            inits[i].mu = 0.5;
            inits[i].restitution = 0.5;
        }
    }

    std::vector<Particle> particles(inits.size());
    for (std::size_t i = 0; i < inits.size(); ++i) {
        const auto& init = inits[i];
        std::string stl_path = init.stl_path.empty() ? cfg.stl_path : init.stl_path;
        if (stl_path.empty()) {
            std::fprintf(stderr, "Missing stl for particle %zu\n", i);
            std::exit(1);
        }
        double scale = init.scale > 0.0 ? init.scale : 1.0;
        int mesh_index = load_mesh_for(stl_path, scale);
        if (mesh_index < 0) {
            std::fprintf(stderr, "Failed to load STL: %s\n", stl_path.c_str());
            std::exit(1);
        }
        const Mesh& mesh = meshes[mesh_index];
        if (mesh.vertices.empty()) {
            std::fprintf(stderr, "Mesh has no vertices.\n");
            std::exit(1);
        }
        Particle p;
        p.tf.pos = init.pos;
        p.tf.rot = init.rot;
        p.vel = init.vel;
        p.omega = init.omega;
        if (mesh.volume > 0.0) {
            double density = init.density > 0.0 ? init.density : 1.0;
            p.mass = density * mesh.volume;
            p.inv_mass = 1.0 / p.mass;
            p.inertia_body = mat3_scale(mesh.inertia_unit, density);
            p.inertia_body_inv = mat3_inverse(p.inertia_body);
        } else {
            double density = init.density > 0.0 ? init.density : 1.0;
            p.mass = density;
            p.inv_mass = 1.0 / p.mass;
            double I = 0.4 * p.mass * mesh.radius * mesh.radius;
            Mat3 Ibody = mat3_identity();
            Ibody.m[0][0] = I;
            Ibody.m[1][1] = I;
            Ibody.m[2][2] = I;
            p.inertia_body = Ibody;
            p.inertia_body_inv = mat3_inverse(Ibody);
        }
        p.radius = mesh.radius;
        p.equiv_radius = mesh.volume > 0.0
            ? std::cbrt(3.0 * mesh.volume / (4.0 * 3.141592653589793))
            : mesh.radius;
        p.young = init.young;
        p.poisson = init.poisson;
        p.mu = init.mu;
        p.restitution = init.restitution;
        p.mesh_index = mesh_index;
        p.L = mat3_mul_vec3(inertia_world(p), p.omega);
        particles[i] = p;
    }

    // ================= walls: static world-space coplanar groups ==============
    // One-time grouping (CPU regroups per step and per particle; walls are
    // static so the geometry never changes). Grouping granularity = whole
    // wall's triangles, no per-particle AABB filter: every triangle of every
    // wall is either merged into a coplanar group or starts one.
    struct WallGroup {
        Vec3 n; double d; double mu;
        std::vector<Triangle> footprints;
    };
    std::vector<WallGroup> wall_groups;
    for (const auto& wi : cfg.wall_inits) {
        if (wi.stl_path.empty()) {
            std::fprintf(stderr, "Wall missing stl path.\n");
            continue;
        }
        std::vector<Triangle> wall_tris;
        if (!load_stl(wi.stl_path, wall_tris)) {
            std::fprintf(stderr, "Failed to load wall STL: %s\n", wi.stl_path.c_str());
            continue;
        }
        double scale = wi.scale > 0.0 ? wi.scale : 1.0;
        for (auto& tri : wall_tris) {
            tri[0] = tri[0] * scale;
            tri[1] = tri[1] * scale;
            tri[2] = tri[2] * scale;
        }
        Mesh wall_mesh = build_mesh(wall_tris, cfg.center_mesh);
        Transform tf;
        tf.pos = wi.pos;
        tf.rot = wi.rot;
        std::vector<Triangle> world_tris = transform_tris(wall_mesh, tf);
        std::printf("Loaded wall: %s (tris=%zu)\n", wi.stl_path.c_str(),
                    world_tris.size());

        // Coplanar groups local to this wall (never merged across walls,
        // same nesting as the CPU loop).
        std::vector<WallGroup> local;
        const double pq = 1e6;
        for (const auto& wtri : world_tris) {
            Vec3 pn = tri_normal(wtri);
            // Orient toward the FIRST particle's position side: push side =
            // particle side (walls static, pile sits on one side). Task 5's
            // kernel treats any-sign penetration so this is only the sign
            // convention of the stored group normal.
            if (!particles.empty() && dot(pn, particles[0].tf.pos - wtri[0]) < 0.0) {
                pn = pn * -1.0;
            }
            double pd = -dot(pn, wtri[0]);
            bool merged = false;
            for (auto& g : local) {
                if (std::llround(g.n.x * pq) == std::llround(pn.x * pq) &&
                    std::llround(g.n.y * pq) == std::llround(pn.y * pq) &&
                    std::llround(g.n.z * pq) == std::llround(pn.z * pq) &&
                    std::llround(g.d * pq) == std::llround(pd * pq)) {
                    g.footprints.push_back(wtri);
                    merged = true;
                    break;
                }
            }
            if (!merged) {
                WallGroup g;
                g.n = pn; g.d = pd; g.mu = wi.mu;
                g.footprints.push_back(wtri);
                local.push_back(g);
            }
        }
        for (auto& g : local) {
            wall_groups.push_back(std::move(g));
        }
    }

    // ================= device upload: meshes (flattened registry) ============
    M.n_mesh = static_cast<int>(meshes.size());
    std::vector<real> st_verts;
    std::vector<int> st_tris;          // GLOBAL vertex indices
    std::vector<int> st_voff(1, 0);
    std::vector<int> st_toff(1, 0);
    std::vector<real> st_mean_edge;
    for (const Mesh& m : meshes) {
        const int vbase = st_voff.back();
        for (const Vec3& v : m.vertices) {
            st_verts.push_back(static_cast<real>(v.x));
            st_verts.push_back(static_cast<real>(v.y));
            st_verts.push_back(static_cast<real>(v.z));
        }
        for (const Triangle& t : m.tris) {
            for (int k = 0; k < 3; ++k) {
                int vi = vertex_index_of(m, t[k]);
                if (vi < 0) {
                    std::fprintf(stderr, "Mesh triangle corner not in vertex list\n");
                    std::exit(1);
                }
                st_tris.push_back(vbase + vi);
            }
        }
        st_voff.push_back(st_voff.back() + static_cast<int>(m.vertices.size()));
        st_toff.push_back(st_toff.back() + static_cast<int>(m.tris.size()));
        st_mean_edge.push_back(static_cast<real>(m.mean_edge));
    }
    M.d_verts = upload_array(st_verts.data(), st_verts.size());
    M.d_tris = upload_array(st_tris.data(), st_tris.size());
    M.d_voffset = upload_array(st_voff.data(), st_voff.size());
    M.d_toffset = upload_array(st_toff.data(), st_toff.size());
    M.d_mean_edge = upload_array(st_mean_edge.data(), st_mean_edge.size());

    // ================= device upload: wall groups ============================
    W.n_groups = static_cast<int>(wall_groups.size());
    std::vector<real> st_gn;          // [3*ng]
    std::vector<real> st_gd;          // [ng]
    std::vector<real> st_fp;          // [9*n_fp] row-major
    std::vector<int> st_fp_start(1, 0);
    std::vector<real> st_gmu;         // [ng]
    for (const auto& g : wall_groups) {
        st_gn.push_back(static_cast<real>(g.n.x));
        st_gn.push_back(static_cast<real>(g.n.y));
        st_gn.push_back(static_cast<real>(g.n.z));
        st_gd.push_back(static_cast<real>(g.d));
        st_gmu.push_back(static_cast<real>(g.mu));
        for (const Triangle& fp : g.footprints) {
            for (int k = 0; k < 3; ++k) {
                st_fp.push_back(static_cast<real>(fp[k].x));
                st_fp.push_back(static_cast<real>(fp[k].y));
                st_fp.push_back(static_cast<real>(fp[k].z));
            }
        }
        st_fp_start.push_back(static_cast<int>(st_fp.size() / 9));
    }
    W.d_gn = upload_array(st_gn.data(), st_gn.size());
    W.d_gd = upload_array(st_gd.data(), st_gd.size());
    W.d_fp = upload_array(st_fp.data(), st_fp.size());
    W.d_fp_start = upload_array(st_fp_start.data(), st_fp_start.size());
    W.d_mu = upload_array(st_gmu.data(), st_gmu.size());

    // ================= device upload: particles (SoA) ========================
    const int n = static_cast<int>(particles.size());
    P.n = n;
    std::vector<real> s_pos(3 * n), s_vel(3 * n), s_omega(3 * n), s_quat(4 * n), s_L(3 * n);
    std::vector<real> s_mass(n), s_inv_mass(n), s_radius(n), s_equiv_radius(n);
    std::vector<real> s_ibinv(9 * n);
    std::vector<real> s_young(n), s_poisson(n), s_mu(n), s_restitution(n);
    std::vector<int> s_mesh_index(n);
    for (int i = 0; i < n; ++i) {
        const Particle& p = particles[i];
        s_pos[3 * i + 0] = static_cast<real>(p.tf.pos.x);
        s_pos[3 * i + 1] = static_cast<real>(p.tf.pos.y);
        s_pos[3 * i + 2] = static_cast<real>(p.tf.pos.z);
        s_vel[3 * i + 0] = static_cast<real>(p.vel.x);
        s_vel[3 * i + 1] = static_cast<real>(p.vel.y);
        s_vel[3 * i + 2] = static_cast<real>(p.vel.z);
        s_omega[3 * i + 0] = static_cast<real>(p.omega.x);
        s_omega[3 * i + 1] = static_cast<real>(p.omega.y);
        s_omega[3 * i + 2] = static_cast<real>(p.omega.z);
        s_quat[4 * i + 0] = static_cast<real>(p.tf.rot.w);
        s_quat[4 * i + 1] = static_cast<real>(p.tf.rot.x);
        s_quat[4 * i + 2] = static_cast<real>(p.tf.rot.y);
        s_quat[4 * i + 3] = static_cast<real>(p.tf.rot.z);
        s_L[3 * i + 0] = static_cast<real>(p.L.x);
        s_L[3 * i + 1] = static_cast<real>(p.L.y);
        s_L[3 * i + 2] = static_cast<real>(p.L.z);
        s_mass[i] = static_cast<real>(p.mass);
        s_inv_mass[i] = static_cast<real>(p.inv_mass);
        s_radius[i] = static_cast<real>(p.radius);
        s_equiv_radius[i] = static_cast<real>(p.equiv_radius);
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                s_ibinv[9 * i + 3 * r + c] = static_cast<real>(p.inertia_body_inv.m[r][c]);
            }
        }
        s_young[i] = static_cast<real>(p.young);
        s_poisson[i] = static_cast<real>(p.poisson);
        s_mu[i] = static_cast<real>(p.mu);
        s_restitution[i] = static_cast<real>(p.restitution);
        s_mesh_index[i] = p.mesh_index;
    }
    P.pos = upload_array(s_pos.data(), s_pos.size());
    P.vel = upload_array(s_vel.data(), s_vel.size());
    P.omega = upload_array(s_omega.data(), s_omega.size());
    P.quat = upload_array(s_quat.data(), s_quat.size());
    P.L = upload_array(s_L.data(), s_L.size());
    P.mass = upload_array(s_mass.data(), s_mass.size());
    P.inv_mass = upload_array(s_inv_mass.data(), s_inv_mass.size());
    P.radius = upload_array(s_radius.data(), s_radius.size());
    P.equiv_radius = upload_array(s_equiv_radius.data(), s_equiv_radius.size());
    P.inertia_body_inv = upload_array(s_ibinv.data(), s_ibinv.size());
    P.young = upload_array(s_young.data(), s_young.size());
    P.poisson = upload_array(s_poisson.data(), s_poisson.size());
    P.mu = upload_array(s_mu.data(), s_mu.size());
    P.restitution = upload_array(s_restitution.data(), s_restitution.size());
    P.mesh_index = upload_array(s_mesh_index.data(), s_mesh_index.size());
    P.force = upload_zeroed<real>(3 * static_cast<std::size_t>(n));
    P.torque = upload_zeroed<real>(3 * static_cast<std::size_t>(n));
    P.contact_count = upload_zeroed<int>(static_cast<std::size_t>(n));

    // Scalar sim parameters (converted to real like everything else).
    dt = static_cast<real>(cfg.dt);
    gravity[0] = static_cast<real>(cfg.gravity.x);
    gravity[1] = static_cast<real>(cfg.gravity.y);
    gravity[2] = static_cast<real>(cfg.gravity.z);
    tangential_damping = static_cast<real>(cfg.tangential_damping);

    // Checksum echo: sum the exact staged values (real -> double) in buffer
    // order, right before/after the device copy. Same summation order as the
    // --check-sums device roundtrip, so equality is bit-exact in both
    // precision modes and verifies the upload path itself.
    host_pos_sum = 0.0;
    for (real v : s_pos) host_pos_sum += static_cast<double>(v);
    host_mass_sum = 0.0;
    for (real v : s_mass) host_mass_sum += static_cast<double>(v);
}

void GpuSim::free_all() {
    auto free_buf = [](void* p) {
        if (p) CUDA_CHECK(cudaFree(p));
    };
    free_buf(M.d_verts); free_buf(M.d_tris);
    free_buf(M.d_voffset); free_buf(M.d_toffset); free_buf(M.d_mean_edge);
    M.n_mesh = 0;
    free_buf(W.d_gn); free_buf(W.d_gd); free_buf(W.d_fp);
    free_buf(W.d_fp_start); free_buf(W.d_mu);
    W.n_groups = 0;
    free_buf(P.pos); free_buf(P.vel); free_buf(P.omega);
    free_buf(P.quat); free_buf(P.L);
    free_buf(P.mass); free_buf(P.inv_mass);
    free_buf(P.radius); free_buf(P.equiv_radius);
    free_buf(P.inertia_body_inv);
    free_buf(P.young); free_buf(P.poisson); free_buf(P.mu); free_buf(P.restitution);
    free_buf(P.mesh_index);
    free_buf(P.force); free_buf(P.torque); free_buf(P.contact_count);
    P.n = 0;
}

void GpuSim::download_soa(std::vector<double>& pos, std::vector<double>& vel,
                          std::vector<double>& omega) const {
    const std::size_t n3 = 3 * static_cast<std::size_t>(P.n);
    pos.assign(n3, 0.0);
    vel.assign(n3, 0.0);
    omega.assign(n3, 0.0);
    if (P.n == 0) return;
    std::vector<real> tmp(n3);
    CUDA_CHECK(cudaMemcpy(tmp.data(), P.pos, n3 * sizeof(real), cudaMemcpyDeviceToHost));
    for (std::size_t k = 0; k < n3; ++k) pos[k] = static_cast<double>(tmp[k]);
    CUDA_CHECK(cudaMemcpy(tmp.data(), P.vel, n3 * sizeof(real), cudaMemcpyDeviceToHost));
    for (std::size_t k = 0; k < n3; ++k) vel[k] = static_cast<double>(tmp[k]);
    CUDA_CHECK(cudaMemcpy(tmp.data(), P.omega, n3 * sizeof(real), cudaMemcpyDeviceToHost));
    for (std::size_t k = 0; k < n3; ++k) omega[k] = static_cast<double>(tmp[k]);
}
