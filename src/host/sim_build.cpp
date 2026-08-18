#include "host/sim_build.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>

#include "core/mat3.hpp"
#include "geometry/mesh_build.hpp"
#include "geometry/stl_io.hpp"

using namespace zdem;

// Moved verbatim from src/main.cpp: used only by the particle build below.
static Mat3 inertia_world(const Particle& p) {
    Mat3 R = quat_to_mat3(p.tf.rot);
    Mat3 Rt = mat3_transpose(R);
    return mat3_mul(mat3_mul(R, p.inertia_body), Rt);
}

bool build_sim(const SimConfig& cfg, SimBuild& out, std::string& err) {
    std::unordered_map<std::string, std::vector<std::array<Vec3, 3>>> tri_cache;
    std::unordered_map<std::string, int> mesh_cache;
    std::vector<Mesh>& meshes = out.meshes;

    auto load_mesh_for = [&](const std::string& path, double scale) -> int {
        std::string key = path + "|" + std::to_string(scale);
        auto it_cache = mesh_cache.find(key);
        if (it_cache != mesh_cache.end()) {
            return it_cache->second;
        }

        auto it = tri_cache.find(path);
        if (it == tri_cache.end()) {
            std::vector<std::array<Vec3, 3>> tris;
            if (!load_stl(path, tris)) {
                return -1;
            }
            tri_cache[path] = tris;
            it = tri_cache.find(path);
        }
        std::vector<std::array<Vec3, 3>> tris = it->second;
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

    std::vector<ParticleInit>& inits = out.inits;
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

    std::vector<Particle>& particles = out.particles;
    particles.resize(inits.size());
    for (std::size_t i = 0; i < inits.size(); ++i) {
        const auto& init = inits[i];
        std::string stl_path = init.stl_path.empty() ? cfg.stl_path : init.stl_path;
        if (stl_path.empty()) {
            err = "Missing stl for particle " + std::to_string(i);
            return false;
        }
        double scale = init.scale > 0.0 ? init.scale : 1.0;
        int mesh_index = load_mesh_for(stl_path, scale);
        if (mesh_index < 0) {
            err = "Failed to load STL: " + stl_path;
            return false;
        }
        const Mesh& mesh = meshes[mesh_index];
        if (mesh.vertices.empty()) {
            err = "Mesh has no vertices.";
            return false;
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
        p.equiv_radius = mesh.volume > 0.0 ? std::cbrt(3.0 * mesh.volume / (4.0 * 3.141592653589793)) : mesh.radius;
        p.young = init.young;
        p.poisson = init.poisson;
        p.mu = init.mu;
        p.restitution = init.restitution;
        p.mesh_index = mesh_index;
        p.L = mat3_mul_vec3(inertia_world(p), p.omega);
        particles[i] = p;
    }

    // Initialize walls
    std::vector<Wall>& walls = out.walls;
    for (const auto& wi : cfg.wall_inits) {
        if (wi.stl_path.empty()) {
            std::cerr << "Wall missing stl path.\n";
            continue;
        }
        std::vector<std::array<Vec3, 3>> wall_tris;
        if (!load_stl(wi.stl_path, wall_tris)) {
            std::cerr << "Failed to load wall STL: " << wi.stl_path << "\n";
            continue;
        }
        double scale = wi.scale > 0.0 ? wi.scale : 1.0;
        for (auto& tri : wall_tris) {
            tri[0] = tri[0] * scale;
            tri[1] = tri[1] * scale;
            tri[2] = tri[2] * scale;
        }
        Wall wall;
        wall.mesh = build_mesh(wall_tris, cfg.center_mesh);
        wall.tf.pos = wi.pos;
        wall.tf.rot = wi.rot;
        wall.mu = wi.mu;
        wall.restitution = wi.restitution;

        // Precompute wall triangle normals in world space
        std::vector<std::array<Vec3, 3>> world_tris = transform_tris(wall.mesh, wall.tf);
        wall.tri_normals.resize(world_tris.size());
        for (std::size_t t = 0; t < world_tris.size(); ++t) {
            wall.tri_normals[t] = tri_normal(world_tris[t]);
        }
        walls.push_back(wall);
        std::cout << "Loaded wall: " << wi.stl_path << " (tris=" << wall.mesh.tris.size() << ")\n";
    }

    out.max_radius = 0.0;
    for (const auto& p : out.particles) {
        out.max_radius = std::max(out.max_radius, p.radius);
    }
    return true;
}
