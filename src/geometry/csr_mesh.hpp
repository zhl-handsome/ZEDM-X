#pragma once

#include <cstdint>
#include <vector>

#include "core/vec3.hpp"
#include "geometry/mesh.hpp"

namespace zdem {

struct CsrMesh {
    std::vector<Vec3> vertices;
    std::vector<std::uint32_t> face_offsets;
    std::vector<std::uint32_t> face_indices;
    std::vector<Vec3> face_normals;
    double radius = 0.0;
    double mean_edge = 0.0;
};

inline CsrMesh build_csr_mesh_from_triangles(const Mesh& mesh) {
    CsrMesh out;
    out.vertices = mesh.vertices;
    out.face_offsets.reserve(mesh.tris.size() + 1);
    out.face_indices.reserve(mesh.tris.size() * 3);
    out.face_normals.reserve(mesh.tris.size());
    out.face_offsets.push_back(0);

    for (const Triangle& tri : mesh.tris) {
        for (const Vec3& v : tri) {
            std::uint32_t index = 0;
            double best_d2 = 1e300;
            for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(out.vertices.size()); ++i) {
                double d2 = norm2(out.vertices[i] - v);
                if (d2 < best_d2) {
                    best_d2 = d2;
                    index = i;
                }
            }
            out.face_indices.push_back(index);
        }
        out.face_offsets.push_back(static_cast<std::uint32_t>(out.face_indices.size()));
        out.face_normals.push_back(normalize(cross(tri[1] - tri[0], tri[2] - tri[0])));
    }

    out.radius = mesh.radius;
    out.mean_edge = mesh.mean_edge;
    return out;
}

}  // namespace zdem
