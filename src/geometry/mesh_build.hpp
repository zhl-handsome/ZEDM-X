#pragma once

#include <vector>

#include "core/transform.hpp"
#include "geometry/mesh.hpp"

namespace zdem {

std::vector<Vec3> unique_vertices(const std::vector<Triangle>& tris);
void compute_mass_properties(const std::vector<Triangle>& tris,
                             double& out_volume,
                             Vec3& out_centroid,
                             Mat3& out_inertia);
Mesh build_mesh(const std::vector<Triangle>& tris, bool center_mesh);
std::vector<Triangle> transform_tris(const Mesh& mesh, const Transform& tf);
Vec3 tri_normal(const Triangle& tri);

}  // namespace zdem
