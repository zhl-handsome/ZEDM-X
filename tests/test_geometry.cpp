#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "geometry/csr_mesh.hpp"
#include "geometry/mesh_build.hpp"

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

    std::vector<Triangle> tris;
    tris.push_back(Triangle{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}});

    Mesh mesh = build_mesh(tris, false);
    require(mesh.tris.size() == 1, "mesh triangle count");
    require(mesh.vertices.size() == 3, "mesh unique vertex count");
    require(mesh.mean_edge > 0.0, "mesh mean edge");

    Vec3 n = tri_normal(mesh.tris[0]);
    require(approx(n.x, 0.0) && approx(n.y, 0.0) && approx(n.z, 1.0), "triangle normal");

    CsrMesh csr = build_csr_mesh_from_triangles(mesh);
    require(csr.face_offsets.size() == 2, "csr offset count");
    require(csr.face_offsets[0] == 0 && csr.face_offsets[1] == 3, "csr offsets");
    require(csr.face_indices.size() == 3, "csr index count");
    require(csr.face_normals.size() == 1, "csr normal count");

    std::cout << "test_geometry passed\n";
    return 0;
}
