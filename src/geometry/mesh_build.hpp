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
// Buffer form of transform_tris: resizes out to mesh.tris.size() and fills it
// with the same per-element expression. resize never shrinks capacity, so a
// caller holding out across steps reuses the allocation (zero steady-state
// heap traffic); the by-value form above delegates here.
void transform_tris_into(const Mesh& mesh, const Transform& tf,
                         std::vector<Triangle>& out);
Vec3 tri_normal(const Triangle& tri);

// Point-in-closed-mesh test (ray parity along +x). Returns false for points
// outside or exactly on the surface; used with point_mesh_distance to detect
// a vertex that has slipped through a face of a concave mesh BEFORE the
// surfaces start crossing (the segment-intersection pipeline is blind then).
bool point_inside_mesh(const std::vector<Triangle>& tris, const Vec3& p);

// Unsigned distance from p to the closest point on the mesh surface, plus the
// closest point itself and the surface normal of the closest face (outward-
// facing as wound in the STL). The (point, closest-point) pair is the contact
// element for the per-vertex penalty force.
double point_mesh_distance(const std::vector<Triangle>& tris, const Vec3& p,
                           Vec3& out_closest_point, Vec3& out_closest_normal);

}  // namespace zdem
