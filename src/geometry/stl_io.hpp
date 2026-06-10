#pragma once

#include <string>
#include <vector>

#include "geometry/mesh.hpp"

namespace zdem {

bool load_stl(const std::string& path, std::vector<Triangle>& tris);

}  // namespace zdem
