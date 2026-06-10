# M0 M1 C++ Core Baseline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Freeze the current single-thread C++ DEM prototype as a reproducible baseline, then extract the math and geometry foundation into tested C++17 modules without changing simulation behavior.

**Architecture:** Keep `src/main.cpp` as the simulation driver for this milestone, but move low-level math and geometry code into small `src/core` and `src/geometry` units. Add CTest-based unit tests and Python smoke tests so later broad-phase, OpenMP, CUDA, and MPI work has a stable correctness baseline.

**Tech Stack:** C++17, CMake 3.14+, CTest, Windows PowerShell, Anaconda Base Python at `D:/ProgramData/anaconda3/python.exe`.

---

## Current Context

- Repository root: `D:/Users/zhl/Desktop/Project/06-FCL-unsphere/ZDEM-X`
- Current executable target: `zdem_cpu`
- Current source file: `src/main.cpp`
- Current CMake file: `CMakeLists.txt`
- Existing configs: `config/example_sim.txt`, `config/wall_test.txt`
- Existing dirty files observed before this plan: `src/main.cpp`, `config/wall_test.txt`
- Do not revert or overwrite user changes in dirty files. Read diffs before editing.

## File Structure After M0 + M1

- Modify: `CMakeLists.txt`
  - Add `zdem_core` static library.
  - Add CTest tests.
  - Keep `zdem_cpu` executable.
- Modify: `src/main.cpp`
  - Remove duplicated low-level definitions after extraction.
  - Include extracted headers.
  - Keep simulation loop, config parsing, contact logic, VTK output in place for this milestone.
- Create: `src/core/constants.hpp`
  - Shared numeric constants.
- Create: `src/core/vec3.hpp`
  - `zdem::Vec3`, vector arithmetic, dot/cross/norm helpers.
- Create: `src/core/quat.hpp`
  - `zdem::Quat`, quaternion normalization, multiplication, conjugation, vector rotation.
- Create: `src/core/mat3.hpp`
  - `zdem::Mat3`, matrix arithmetic, determinant, inverse, quaternion-to-matrix conversion.
- Create: `src/core/transform.hpp`
  - `zdem::Transform`.
- Create: `src/geometry/mesh.hpp`
  - `zdem::Mesh`, `zdem::Triangle`.
- Create: `src/geometry/stl_io.hpp`
  - STL loader declaration.
- Create: `src/geometry/stl_io.cpp`
  - STL loader implementation moved from `src/main.cpp`.
- Create: `src/geometry/mesh_build.hpp`
  - Mesh construction and mass-property declarations.
- Create: `src/geometry/mesh_build.cpp`
  - `unique_vertices`, `compute_mass_properties`, `build_mesh`, `transform_tris`, `tri_normal`.
- Create: `src/geometry/csr_mesh.hpp`
  - CSR-style mesh data structure for future SoA/CUDA work.
- Create: `tests/test_math.cpp`
  - Unit tests for vector, quaternion, and matrix behavior.
- Create: `tests/test_geometry.cpp`
  - Unit tests for mesh construction and CSR mesh conversion.
- Create: `tests/smoke/run_smoke_tests.py`
  - Runs small configs and validates executable success/output.
- Create: `config/smoke_no_contact.txt`
  - Fast no-contact smoke config.
- Create: `config/smoke_wall_contact.txt`
  - Fast particle-wall smoke config.
- Create: `docs/m0_baseline.md`
  - Baseline commands, configs, and acceptance criteria.

---

## Task 1: Protect Current Workspace And Record Baseline

**Files:**
- Create: `docs/m0_baseline.md`

- [ ] **Step 1: Inspect dirty files before editing**

Run:

```powershell
git status --short
git diff -- src/main.cpp config/wall_test.txt
```

Expected:

```text
M config/wall_test.txt
M src/main.cpp
```

If the diff contains user edits unrelated to this milestone, keep them and edit around them.

- [ ] **Step 2: Build current executable**

Run:

```powershell
cmake -S . -B build
cmake --build build
```

Expected:

```text
zdem_cpu builds successfully
```

- [ ] **Step 3: Run the current example config**

Run:

```powershell
.\build\Debug\zdem_cpu.exe --config config/example_sim.txt
```

If the generator creates a different configuration, run the executable path printed by CMake under `build`.

Expected:

```text
The program prints step lines and ends with done.
```

- [ ] **Step 4: Create baseline document**

Create `docs/m0_baseline.md` with this content:

````markdown
# M0 Baseline

Date: 2026-06-10

## Purpose

This document records the reproducible baseline for the current C++ DEM prototype before the M1 math and geometry extraction.

## Commands

```powershell
cmake -S . -B build
cmake --build build
.\build\Debug\zdem_cpu.exe --config config/example_sim.txt
```

## Expected Behavior

- The executable builds as `zdem_cpu`.
- The example simulation prints progress lines containing `step=`.
- The run exits with code 0.
- The run prints `done` before exit.

## Known Current Limitations

- `src/main.cpp` contains math, geometry, contact detection, force calculation, integration, IO, and CLI code in one file.
- The baseline uses a brute-force particle-pair loop and triangle checks inside candidate pairs.
- Particle-wall contact and tangential history exist in the prototype but are not isolated behind stable interfaces.
- There is no test target in CMake before M1.
````

- [ ] **Step 5: Commit baseline doc only if user asks for commits**

Do not commit automatically. If the user asks for commits, run:

```powershell
git add docs/m0_baseline.md
git commit -m "docs: record M0 baseline"
```

Expected:

```text
Commit succeeds and includes only docs/m0_baseline.md.
```

---

## Task 2: Add Core Math Headers

**Files:**
- Create: `src/core/constants.hpp`
- Create: `src/core/vec3.hpp`
- Create: `src/core/quat.hpp`
- Create: `src/core/mat3.hpp`
- Create: `src/core/transform.hpp`

- [ ] **Step 1: Create constants header**

Create `src/core/constants.hpp`:

```cpp
#pragma once

namespace zdem {

constexpr double kEps = 1e-12;
constexpr double kTiny = 1e-30;
constexpr double kPi = 3.141592653589793238462643383279502884;

}  // namespace zdem
```

- [ ] **Step 2: Create vector header**

Create `src/core/vec3.hpp`:

```cpp
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
```

- [ ] **Step 3: Create quaternion header**

Create `src/core/quat.hpp`:

```cpp
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
```

- [ ] **Step 4: Create matrix header**

Create `src/core/mat3.hpp` by moving the current `Mat3` helpers from `src/main.cpp` into namespace `zdem`. Include these functions exactly: `mat3_identity`, `mat3_add`, `mat3_sub`, `mat3_scale`, `mat3_outer`, `quat_to_mat3`, `mat3_mul_vec3`, `mat3_mul`, `mat3_transpose`, `mat3_det`, `mat3_inverse`.

The header must start with:

```cpp
#pragma once

#include <cmath>

#include "core/quat.hpp"
#include "core/vec3.hpp"

namespace zdem {
```

The header must end with:

```cpp
}  // namespace zdem
```

- [ ] **Step 5: Create transform header**

Create `src/core/transform.hpp`:

```cpp
#pragma once

#include "core/quat.hpp"
#include "core/vec3.hpp"

namespace zdem {

struct Transform {
    Vec3 pos;
    Quat rot;
};

}  // namespace zdem
```

- [ ] **Step 6: Add math unit tests**

Create `tests/test_math.cpp`:

```cpp
#include <cmath>
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
```

- [ ] **Step 7: Wire math test in CMake**

Modify `CMakeLists.txt` to include:

```cmake
enable_testing()

add_executable(test_math
    tests/test_math.cpp
)
target_include_directories(test_math PRIVATE src)
add_test(NAME test_math COMMAND test_math)
```

- [ ] **Step 8: Run math test and verify it passes**

Run:

```powershell
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure -R test_math
```

Expected:

```text
100% tests passed
```

---

## Task 3: Extract Mesh, STL IO, And Transform Helpers

**Files:**
- Create: `src/geometry/mesh.hpp`
- Create: `src/geometry/stl_io.hpp`
- Create: `src/geometry/stl_io.cpp`
- Create: `src/geometry/mesh_build.hpp`
- Create: `src/geometry/mesh_build.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Create mesh header**

Create `src/geometry/mesh.hpp`:

```cpp
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
```

- [ ] **Step 2: Create STL IO declarations**

Create `src/geometry/stl_io.hpp`:

```cpp
#pragma once

#include <string>
#include <vector>

#include "geometry/mesh.hpp"

namespace zdem {

bool load_stl(const std::string& path, std::vector<Triangle>& tris);

}  // namespace zdem
```

- [ ] **Step 3: Create STL IO implementation**

Create `src/geometry/stl_io.cpp` by moving the current `load_stl` implementation from `src/main.cpp`. Change the triangle type from `std::array<Vec3, 3>` to `Triangle`, wrap the function in `namespace zdem`, and include:

```cpp
#include "geometry/stl_io.hpp"

#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
```

- [ ] **Step 4: Create mesh build declarations**

Create `src/geometry/mesh_build.hpp`:

```cpp
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
Vec3 tri_normal(const Triangle& tri);

}  // namespace zdem
```

- [ ] **Step 5: Create mesh build implementation**

Create `src/geometry/mesh_build.cpp` by moving these current functions from `src/main.cpp`: `unique_vertices`, `compute_mass_properties`, `build_mesh`, `transform_tris`, `tri_normal`.

Required includes:

```cpp
#include "geometry/mesh_build.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
```

Required namespace wrapper:

```cpp
namespace zdem {

// moved functions

}  // namespace zdem
```

- [ ] **Step 6: Add `zdem_core` library target**

Modify `CMakeLists.txt` so it contains:

```cmake
cmake_minimum_required(VERSION 3.14)
project(zdem_cpu LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_library(zdem_core
    src/geometry/stl_io.cpp
    src/geometry/mesh_build.cpp
)
target_include_directories(zdem_core PUBLIC src)

add_executable(zdem_cpu
    src/main.cpp
)
target_link_libraries(zdem_cpu PRIVATE zdem_core)

add_executable(test_math
    tests/test_math.cpp
)
target_link_libraries(test_math PRIVATE zdem_core)
add_test(NAME test_math COMMAND test_math)
```

- [ ] **Step 7: Update main includes**

At the top of `src/main.cpp`, add:

```cpp
#include "core/mat3.hpp"
#include "core/quat.hpp"
#include "core/transform.hpp"
#include "core/vec3.hpp"
#include "geometry/mesh.hpp"
#include "geometry/mesh_build.hpp"
#include "geometry/stl_io.hpp"
```

After includes, add:

```cpp
using namespace zdem;
```

- [ ] **Step 8: Remove duplicate moved definitions from main**

Remove from `src/main.cpp` only the definitions moved to headers or `.cpp` files:

```text
Vec3
dot
cross
norm
normalize
norm2
Quat
quat_conj
quat_mul
quat_normalize
quat_rotate
Mat3
mat3_identity
mat3_add
mat3_sub
mat3_scale
mat3_outer
quat_to_mat3
mat3_mul_vec3
mat3_mul
mat3_transpose
mat3_det
mat3_inverse
Transform
Mesh
unique_vertices
load_stl
compute_mass_properties
build_mesh
transform_tris
tri_normal
```

Do not remove `ParticleInit`, `Particle`, `WallInit`, `Wall`, geometry contact functions, GJK/EPA, config parsing, VTK output, or `main`.

- [ ] **Step 9: Build after extraction**

Run:

```powershell
cmake -S . -B build
cmake --build build
```

Expected:

```text
zdem_core, zdem_cpu, and test_math build successfully.
```

- [ ] **Step 10: Run the existing example after extraction**

Run:

```powershell
.\build\Debug\zdem_cpu.exe --config config/example_sim.txt
```

Expected:

```text
The run exits with code 0 and prints done.
```

---

## Task 4: Add CSR Mesh Structure For Future SoA/CUDA Work

**Files:**
- Create: `src/geometry/csr_mesh.hpp`
- Create: `tests/test_geometry.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Create CSR mesh header**

Create `src/geometry/csr_mesh.hpp`:

```cpp
#pragma once

#include <array>
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
```

- [ ] **Step 2: Create geometry tests**

Create `tests/test_geometry.cpp`:

```cpp
#include <cmath>
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
```

- [ ] **Step 3: Wire geometry test in CMake**

Add to `CMakeLists.txt`:

```cmake
add_executable(test_geometry
    tests/test_geometry.cpp
)
target_link_libraries(test_geometry PRIVATE zdem_core)
add_test(NAME test_geometry COMMAND test_geometry)
```

- [ ] **Step 4: Run geometry tests**

Run:

```powershell
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure -R "test_math|test_geometry"
```

Expected:

```text
100% tests passed
```

---

## Task 5: Add Fast Smoke Configs And Runner

**Files:**
- Create: `config/smoke_no_contact.txt`
- Create: `config/smoke_wall_contact.txt`
- Create: `tests/smoke/run_smoke_tests.py`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Create no-contact smoke config**

Create `config/smoke_no_contact.txt`:

```text
steps = 2
dt = 0.00001
split_contacts = 0
contact_debug = 0
gravity = 0 0 0
center_mesh = 1
vtk_prefix = smoke_no_contact
output_dir = build/test_output/smoke_no_contact
output_interval = 1

particle
stl = geometry/low-poly-banbana.stl
pos = -1 0 0
vel = 0 0 0
quat = 1 0 0 0
omega = 0 0 0
scale = 0.001
density = 2500
young = 1e7
poisson = 0.25
mu = 0.5
restitution = 0.5
end_particle

particle
stl = geometry/low-poly-banbana.stl
pos = 1 0 0
vel = 0 0 0
quat = 1 0 0 0
omega = 0 0 0
scale = 0.001
density = 2500
young = 1e7
poisson = 0.25
mu = 0.5
restitution = 0.5
end_particle
```

- [ ] **Step 2: Create wall smoke config**

Create `config/smoke_wall_contact.txt`:

```text
steps = 2
dt = 0.00001
split_contacts = 0
contact_debug = 0
gravity = 0 0 -9.81
center_mesh = 1
vtk_prefix = smoke_wall_contact
output_dir = build/test_output/smoke_wall_contact
output_interval = 1

wall
stl = geometry/wall_plane.stl
pos = 0 0 0
quat = 1 0 0 0
scale = 1.0
mu = 0.5
restitution = 0.1
end_wall

particle
stl = geometry/low-poly-banbana.stl
pos = 0 0 0.2
vel = 0 0 0
quat = 1 0 0 0
omega = 0 0 0
scale = 0.001
density = 2500
young = 1e7
poisson = 0.25
mu = 0.5
restitution = 0.1
end_particle
```

- [ ] **Step 3: Create smoke runner**

Create `tests/smoke/run_smoke_tests.py`:

```python
import argparse
import pathlib
import subprocess
import sys


def run_case(exe: pathlib.Path, root: pathlib.Path, config: str) -> None:
    cmd = [str(exe), "--config", config]
    completed = subprocess.run(
        cmd,
        cwd=str(root),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
    )
    if completed.returncode != 0:
        print(completed.stdout)
        raise SystemExit(f"{config} failed with exit code {completed.returncode}")
    if "done" not in completed.stdout:
        print(completed.stdout)
        raise SystemExit(f"{config} did not print done")
    if "step=" not in completed.stdout:
        print(completed.stdout)
        raise SystemExit(f"{config} did not print step progress")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True)
    parser.add_argument("--root", required=True)
    args = parser.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    root = pathlib.Path(args.root).resolve()
    cases = [
        "config/smoke_no_contact.txt",
        "config/smoke_wall_contact.txt",
    ]
    for case in cases:
        run_case(exe, root, case)
    print("smoke tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Wire smoke test in CMake**

Add to `CMakeLists.txt`:

```cmake
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(
        NAME smoke_configs
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tests/smoke/run_smoke_tests.py
                --exe $<TARGET_FILE:zdem_cpu>
                --root ${CMAKE_SOURCE_DIR}
    )
endif()
```

- [ ] **Step 5: Run smoke tests using Anaconda Python if CMake cannot find Python**

First run:

```powershell
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure -R smoke_configs
```

If CMake does not find Python, run directly:

```powershell
D:/ProgramData/anaconda3/python.exe tests/smoke/run_smoke_tests.py --exe build/Debug/zdem_cpu.exe --root .
```

Expected:

```text
smoke tests passed
```

---

## Task 6: Final Verification And Regression Check

**Files:**
- Modify only files created or intentionally refactored by Tasks 1-5.

- [ ] **Step 1: Run full configure/build/test**

Run:

```powershell
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected:

```text
100% tests passed
```

- [ ] **Step 2: Run original example manually**

Run:

```powershell
.\build\Debug\zdem_cpu.exe --config config/example_sim.txt
```

Expected:

```text
The executable prints step progress and done.
```

- [ ] **Step 3: Check changed files**

Run:

```powershell
git status --short
git diff --stat
```

Expected changed files:

```text
CMakeLists.txt
docs/m0_baseline.md
src/main.cpp
src/core/constants.hpp
src/core/vec3.hpp
src/core/quat.hpp
src/core/mat3.hpp
src/core/transform.hpp
src/geometry/mesh.hpp
src/geometry/stl_io.hpp
src/geometry/stl_io.cpp
src/geometry/mesh_build.hpp
src/geometry/mesh_build.cpp
src/geometry/csr_mesh.hpp
tests/test_math.cpp
tests/test_geometry.cpp
tests/smoke/run_smoke_tests.py
config/smoke_no_contact.txt
config/smoke_wall_contact.txt
```

Existing dirty files `src/main.cpp` and `config/wall_test.txt` may include user changes from before this milestone. Preserve them unless the user explicitly asks to discard them.

- [ ] **Step 4: Optional commit if user asks**

Run only after user approval:

```powershell
git add CMakeLists.txt docs/m0_baseline.md src/main.cpp src/core src/geometry tests config/smoke_no_contact.txt config/smoke_wall_contact.txt
git commit -m "refactor: extract C++ core math and geometry baseline"
```

Expected:

```text
Commit succeeds.
```

---

## Acceptance Criteria

- `zdem_cpu` still builds and runs `config/example_sim.txt`.
- `test_math` passes.
- `test_geometry` passes.
- `smoke_configs` passes or the direct Anaconda Python smoke command passes.
- `src/main.cpp` no longer owns `Vec3`, `Quat`, `Mat3`, `Transform`, `Mesh`, STL loading, or mesh build helpers.
- The extracted core APIs are in namespace `zdem`.
- No user edits are reverted.
- No OpenMP, CUDA, MPI, or algorithmic contact-model changes are introduced in M0 + M1.

## Self-Review Notes

- Scope covers M0 baseline freezing and M1 math/geometry extraction only.
- CPU broad-phase redesign, contact history redesign, OpenMP, CUDA, and MPI are intentionally outside this plan.
- Every new test has an exact command and expected result.
- The plan avoids changing physical behavior except for mechanical namespace/include extraction.
