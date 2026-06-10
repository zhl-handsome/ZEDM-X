#pragma once

#include <array>
#include <vector>

#include "core/mat3.hpp"
#include "core/vec3.hpp"

namespace zdem {

using Triangle = std::array<Vec3, 3>;

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<Triangle> tris;
    Vec3 center;
    double radius = 0.0;
    double mean_edge = 0.0;
    double bbox_diag = 0.0;
    double volume = 0.0;
    Mat3 inertia_unit;
};

}  // namespace zdem
