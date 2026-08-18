#pragma once

#include <string>
#include <vector>

#include "core/transform.hpp"
#include "core/vec3.hpp"
#include "geometry/mesh.hpp"
#include "host/config_io.hpp"
#include "host/vtk_io.hpp"

using namespace zdem;

struct Wall {
    Mesh mesh;
    Transform tf;
    double mu = 0.5;
    double restitution = 0.5;
    std::vector<Vec3> tri_normals;
};

struct SimBuild {
    std::vector<Mesh> meshes;          // dedup mesh registry (by stl path + scale, same as the old main.cpp logic)
    std::vector<Particle> particles;   // config order; gid == index
    std::vector<ParticleInit> inits;   // particle init table (main's state-txt output / future MPI restart)
    std::vector<Wall> walls;
    double max_radius = 0.0;           // max(p.radius), for ghost depth / margin
};

// Build meshes/particles/walls from SimConfig (extracted verbatim from
// main.cpp; error paths now report through err). On false, err holds the
// message.
bool build_sim(const SimConfig& cfg, SimBuild& out, std::string& err);
