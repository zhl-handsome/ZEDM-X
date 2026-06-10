#pragma once

#include <cmath>

#include "core/quat.hpp"
#include "core/vec3.hpp"

namespace zdem {

struct Mat3 {
    double m[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
};

inline Mat3 mat3_identity() {
    Mat3 r;
    r.m[0][0] = 1.0;
    r.m[1][1] = 1.0;
    r.m[2][2] = 1.0;
    return r;
}

inline Mat3 mat3_add(const Mat3& a, const Mat3& b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = a.m[i][j] + b.m[i][j];
        }
    }
    return r;
}

inline Mat3 mat3_sub(const Mat3& a, const Mat3& b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = a.m[i][j] - b.m[i][j];
        }
    }
    return r;
}

inline Mat3 mat3_scale(const Mat3& a, double s) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = a.m[i][j] * s;
        }
    }
    return r;
}

inline Mat3 mat3_outer(const Vec3& a, const Vec3& b) {
    Mat3 r;
    r.m[0][0] = a.x * b.x; r.m[0][1] = a.x * b.y; r.m[0][2] = a.x * b.z;
    r.m[1][0] = a.y * b.x; r.m[1][1] = a.y * b.y; r.m[1][2] = a.y * b.z;
    r.m[2][0] = a.z * b.x; r.m[2][1] = a.z * b.y; r.m[2][2] = a.z * b.z;
    return r;
}

inline Mat3 quat_to_mat3(const Quat& q) {
    Quat n = quat_normalize(q);
    double w = n.w, x = n.x, y = n.y, z = n.z;
    Mat3 R;
    R.m[0][0] = 1.0 - 2.0 * (y * y + z * z);
    R.m[0][1] = 2.0 * (x * y - z * w);
    R.m[0][2] = 2.0 * (x * z + y * w);
    R.m[1][0] = 2.0 * (x * y + z * w);
    R.m[1][1] = 1.0 - 2.0 * (x * x + z * z);
    R.m[1][2] = 2.0 * (y * z - x * w);
    R.m[2][0] = 2.0 * (x * z - y * w);
    R.m[2][1] = 2.0 * (y * z + x * w);
    R.m[2][2] = 1.0 - 2.0 * (x * x + y * y);
    return R;
}

inline Vec3 mat3_mul_vec3(const Mat3& m, const Vec3& v) {
    return Vec3{
        m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z,
        m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z,
        m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z,
    };
}

inline Mat3 mat3_mul(const Mat3& a, const Mat3& b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                r.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }
    return r;
}

inline Mat3 mat3_transpose(const Mat3& m) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = m.m[j][i];
        }
    }
    return r;
}

inline double mat3_det(const Mat3& a) {
    return a.m[0][0] * (a.m[1][1] * a.m[2][2] - a.m[1][2] * a.m[2][1])
         - a.m[0][1] * (a.m[1][0] * a.m[2][2] - a.m[1][2] * a.m[2][0])
         + a.m[0][2] * (a.m[1][0] * a.m[2][1] - a.m[1][1] * a.m[2][0]);
}

inline Mat3 mat3_inverse(const Mat3& a) {
    double det = mat3_det(a);
    if (std::abs(det) < 1e-18) {
        return mat3_identity();
    }
    double invdet = 1.0 / det;
    Mat3 r;
    r.m[0][0] =  (a.m[1][1] * a.m[2][2] - a.m[1][2] * a.m[2][1]) * invdet;
    r.m[0][1] = -(a.m[0][1] * a.m[2][2] - a.m[0][2] * a.m[2][1]) * invdet;
    r.m[0][2] =  (a.m[0][1] * a.m[1][2] - a.m[0][2] * a.m[1][1]) * invdet;
    r.m[1][0] = -(a.m[1][0] * a.m[2][2] - a.m[1][2] * a.m[2][0]) * invdet;
    r.m[1][1] =  (a.m[0][0] * a.m[2][2] - a.m[0][2] * a.m[2][0]) * invdet;
    r.m[1][2] = -(a.m[0][0] * a.m[1][2] - a.m[0][2] * a.m[1][0]) * invdet;
    r.m[2][0] =  (a.m[1][0] * a.m[2][1] - a.m[1][1] * a.m[2][0]) * invdet;
    r.m[2][1] = -(a.m[0][0] * a.m[2][1] - a.m[0][1] * a.m[2][0]) * invdet;
    r.m[2][2] =  (a.m[0][0] * a.m[1][1] - a.m[0][1] * a.m[1][0]) * invdet;
    return r;
}

}  // namespace zdem
