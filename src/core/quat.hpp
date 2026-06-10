#pragma once

#include <cmath>

#include "core/vec3.hpp"

namespace zdem {

struct Quat {
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Quat quat_conj(const Quat& q) {
    return Quat{q.w, -q.x, -q.y, -q.z};
}

inline Quat quat_mul(const Quat& a, const Quat& b) {
    return Quat{
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
    };
}

inline Quat quat_normalize(const Quat& q) {
    double n = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    if (n < 1e-14) {
        return Quat{};
    }
    return Quat{q.w / n, q.x / n, q.y / n, q.z / n};
}

inline Vec3 quat_rotate(const Quat& q, const Vec3& v) {
    Quat p{0.0, v.x, v.y, v.z};
    Quat qn = quat_normalize(q);
    Quat r = quat_mul(quat_mul(qn, p), quat_conj(qn));
    return Vec3{r.x, r.y, r.z};
}

}  // namespace zdem
