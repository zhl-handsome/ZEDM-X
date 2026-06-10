#pragma once

#include <cmath>

#include "core/constants.hpp"

namespace zdem {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return Vec3{x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return Vec3{x - o.x, y - o.y, z - o.z}; }
    Vec3 operator-() const { return Vec3{-x, -y, -z}; }
    Vec3 operator*(double s) const { return Vec3{x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return Vec3{x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};

inline Vec3 operator*(double s, const Vec3& v) {
    return v * s;
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

inline double norm2(const Vec3& v) {
    return dot(v, v);
}

inline double norm(const Vec3& v) {
    return std::sqrt(norm2(v));
}

inline Vec3 normalize(const Vec3& v) {
    double n = norm(v);
    if (n < 1e-14) {
        return Vec3{0.0, 0.0, 0.0};
    }
    return v / n;
}

}  // namespace zdem
