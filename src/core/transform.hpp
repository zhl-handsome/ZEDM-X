#pragma once

#include "core/quat.hpp"
#include "core/vec3.hpp"

namespace zdem {

struct Transform {
    Vec3 pos;
    Quat rot;
};

}  // namespace zdem
