#pragma once

#include <string>
#include <vector>

#include "core/mat3.hpp"
#include "core/transform.hpp"
#include "core/vec3.hpp"
#include "geometry/mesh.hpp"

using namespace zdem;

struct Particle {
    Transform tf;
    Vec3 vel;
    Vec3 omega;
    Vec3 L;
    double mass = 1.0;
    double inv_mass = 1.0;
    Mat3 inertia_body;
    Mat3 inertia_body_inv;
    double radius = 0.0;
    double equiv_radius = 0.0;
    double young = 0.0;
    double poisson = 0.0;
    double mu = 0.0;
    double restitution = 0.0;
    int mesh_index = 0;
};

void write_vtk_particles(const std::string& path,
                         const std::vector<Mesh>& meshes,
                         const std::vector<Particle>& particles,
                         const std::vector<Vec3>& forces,
                         const std::vector<Vec3>& torques,
                         const std::vector<int>& contact_counts);
