#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "core/mat3.hpp"
#include "core/quat.hpp"
#include "core/vec3.hpp"

namespace {

bool approx(double a, double b, double eps = 1e-12) {
    return std::abs(a - b) <= eps;
}

void require(bool ok, const std::string& msg) {
    if (!ok) {
        std::cerr << "FAIL: " << msg << "\n";
        std::exit(1);
    }
}

}  // namespace

int main() {
    using namespace zdem;

    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{-2.0, 0.5, 4.0};
    Vec3 c = a + b;
    require(approx(c.x, -1.0) && approx(c.y, 2.5) && approx(c.z, 7.0), "Vec3 addition");
    require(approx(dot(a, b), 11.0), "Vec3 dot product");

    Vec3 x{1.0, 0.0, 0.0};
    Vec3 y{0.0, 1.0, 0.0};
    Vec3 z = cross(x, y);
    require(approx(z.x, 0.0) && approx(z.y, 0.0) && approx(z.z, 1.0), "Vec3 cross product");
    require(approx(norm(Vec3{3.0, 4.0, 0.0}), 5.0), "Vec3 norm");

    Quat q = quat_normalize(Quat{2.0, 0.0, 0.0, 0.0});
    require(approx(q.w, 1.0) && approx(q.x, 0.0), "Quat normalize");
    Vec3 rotated = quat_rotate(Quat{1.0, 0.0, 0.0, 0.0}, Vec3{1.0, 2.0, 3.0});
    require(approx(rotated.x, 1.0) && approx(rotated.y, 2.0) && approx(rotated.z, 3.0), "Quat identity rotation");

    Mat3 I = mat3_identity();
    Vec3 iv = mat3_mul_vec3(I, Vec3{5.0, -1.0, 2.0});
    require(approx(iv.x, 5.0) && approx(iv.y, -1.0) && approx(iv.z, 2.0), "Mat3 identity multiply");
    require(approx(mat3_det(I), 1.0), "Mat3 determinant identity");

    std::cout << "test_math passed\n";
    return 0;
}
