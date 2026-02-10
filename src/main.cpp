#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return Vec3{x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return Vec3{x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return Vec3{x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return Vec3{x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};

static inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

static inline double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

static inline Vec3 normalize(const Vec3& v) {
    double n = norm(v);
    if (n < 1e-14) {
        return Vec3{0.0, 0.0, 0.0};
    }
    return v / n;
}

struct Quat {
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

static inline Quat quat_conj(const Quat& q) {
    return Quat{q.w, -q.x, -q.y, -q.z};
}

static inline Quat quat_mul(const Quat& a, const Quat& b) {
    return Quat{
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
    };
}

static inline Quat quat_normalize(const Quat& q) {
    double n = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    if (n < 1e-14) {
        return Quat{};
    }
    return Quat{q.w / n, q.x / n, q.y / n, q.z / n};
}

static inline Vec3 quat_rotate(const Quat& q, const Vec3& v) {
    Quat p{0.0, v.x, v.y, v.z};
    Quat qn = quat_normalize(q);
    Quat r = quat_mul(quat_mul(qn, p), quat_conj(qn));
    return Vec3{r.x, r.y, r.z};
}

struct Mat3 {
    double m[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
};

static inline Mat3 quat_to_mat3(const Quat& q) {
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

static inline Vec3 mat3_mul_vec3(const Mat3& m, const Vec3& v) {
    return Vec3{
        m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z,
        m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z,
        m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z,
    };
}

static inline Mat3 mat3_mul(const Mat3& a, const Mat3& b) {
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

static inline Mat3 mat3_transpose(const Mat3& m) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = m.m[j][i];
        }
    }
    return r;
}

struct Transform {
    Vec3 pos;
    Quat rot;
};

struct Mesh {
    std::vector<Vec3> vertices;
    Vec3 center;
    double radius = 0.0;
};

struct Particle {
    Transform tf;
    Vec3 vel;
    Vec3 omega;
    double mass = 1.0;
    double inv_mass = 1.0;
    Vec3 inertia_body;
    Vec3 inertia_body_inv;
    double radius = 0.0;
};

static std::vector<Vec3> unique_vertices(const std::vector<std::array<Vec3, 3>>& tris) {
    std::vector<Vec3> verts;
    const double tol = 1e-10;
    for (const auto& tri : tris) {
        for (const auto& v : tri) {
            bool found = false;
            for (const auto& u : verts) {
                Vec3 d = v - u;
                if (dot(d, d) <= tol) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                verts.push_back(v);
            }
        }
    }
    return verts;
}

static bool load_stl(const std::string& path, std::vector<std::array<Vec3, 3>>& tris) {
    tris.clear();
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }

    in.seekg(0, std::ios::end);
    std::streampos size = in.tellg();
    in.seekg(0, std::ios::beg);
    if (size < 84) {
        return false;
    }

    std::array<char, 80> header{};
    in.read(header.data(), header.size());
    uint32_t n_tri = 0;
    in.read(reinterpret_cast<char*>(&n_tri), sizeof(uint32_t));

    std::streampos expected = 84 + static_cast<std::streampos>(n_tri) * 50;
    if (expected == size) {
        tris.reserve(n_tri);
        for (uint32_t i = 0; i < n_tri; ++i) {
            float data[12];
            uint16_t attr = 0;
            in.read(reinterpret_cast<char*>(data), sizeof(float) * 12);
            in.read(reinterpret_cast<char*>(&attr), sizeof(uint16_t));
            std::array<Vec3, 3> tri{
                Vec3{data[3], data[4], data[5]},
                Vec3{data[6], data[7], data[8]},
                Vec3{data[9], data[10], data[11]}
            };
            tris.push_back(tri);
        }
        return true;
    }

    in.close();
    std::ifstream ia(path);
    if (!ia) {
        return false;
    }
    std::string line;
    std::vector<Vec3> vbuf;
    while (std::getline(ia, line)) {
        std::string s = line;
        for (auto& c : s) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        if (s.find("vertex") != std::string::npos) {
            std::istringstream iss(s);
            std::string tag;
            double x, y, z;
            if (iss >> tag >> x >> y >> z) {
                vbuf.emplace_back(x, y, z);
            }
        }
        if (s.find("endloop") != std::string::npos && vbuf.size() >= 3) {
            std::array<Vec3, 3> tri{vbuf[vbuf.size() - 3], vbuf[vbuf.size() - 2], vbuf[vbuf.size() - 1]};
            tris.push_back(tri);
        }
    }
    return !tris.empty();
}

static Mesh build_mesh(const std::vector<std::array<Vec3, 3>>& tris, bool center_mesh) {
    Mesh m;
    m.vertices = unique_vertices(tris);
    Vec3 center{0.0, 0.0, 0.0};
    for (const auto& v : m.vertices) {
        center += v;
    }
    if (!m.vertices.empty()) {
        center *= (1.0 / static_cast<double>(m.vertices.size()));
    }
    m.center = center;
    if (center_mesh) {
        for (auto& v : m.vertices) {
            v -= center;
        }
        m.center = Vec3{0.0, 0.0, 0.0};
    }
    double r2 = 0.0;
    for (const auto& v : m.vertices) {
        r2 = std::max(r2, dot(v, v));
    }
    m.radius = std::sqrt(r2);
    return m;
}

static Vec3 support_point(const Mesh& mesh, const Transform& tf, const Vec3& dir) {
    Quat qc = quat_conj(tf.rot);
    Vec3 local_dir = quat_rotate(qc, dir);
    double best = -std::numeric_limits<double>::infinity();
    Vec3 best_v;
    for (const auto& v : mesh.vertices) {
        double d = dot(v, local_dir);
        if (d > best) {
            best = d;
            best_v = v;
        }
    }
    return quat_rotate(tf.rot, best_v) + tf.pos;
}

static Vec3 support_minkowski(const Mesh& a, const Transform& ta,
                              const Mesh& b, const Transform& tb,
                              const Vec3& dir) {
    Vec3 p1 = support_point(a, ta, dir);
    Vec3 p2 = support_point(b, tb, dir * -1.0);
    return p1 - p2;
}

struct Simplex {
    std::vector<Vec3> pts;
};

static bool handle_simplex(Simplex& s, Vec3& dir) {
    if (s.pts.size() == 2) {
        Vec3 a = s.pts[1];
        Vec3 b = s.pts[0];
        Vec3 ab = b - a;
        Vec3 ao = Vec3{} - a;
        if (dot(ab, ao) > 0.0) {
            dir = cross(cross(ab, ao), ab);
        } else {
            s.pts = {a};
            dir = ao;
        }
        return false;
    }
    if (s.pts.size() == 3) {
        Vec3 a = s.pts[2];
        Vec3 b = s.pts[1];
        Vec3 c = s.pts[0];
        Vec3 ab = b - a;
        Vec3 ac = c - a;
        Vec3 ao = Vec3{} - a;
        Vec3 abc = cross(ab, ac);

        Vec3 ab_perp = cross(abc, ab);
        if (dot(ab_perp, ao) > 0.0) {
            s.pts = {b, a};
            dir = cross(cross(ab, ao), ab);
            return false;
        }
        Vec3 ac_perp = cross(ac, abc);
        if (dot(ac_perp, ao) > 0.0) {
            s.pts = {c, a};
            dir = cross(cross(ac, ao), ac);
            return false;
        }

        if (dot(abc, ao) > 0.0) {
            dir = abc;
        } else {
            s.pts = {b, c, a};
            dir = abc * -1.0;
        }
        return false;
    }
    if (s.pts.size() == 4) {
        Vec3 a = s.pts[3];
        Vec3 b = s.pts[2];
        Vec3 c = s.pts[1];
        Vec3 d = s.pts[0];
        Vec3 ao = Vec3{} - a;
        Vec3 ab = b - a;
        Vec3 ac = c - a;
        Vec3 ad = d - a;

        Vec3 abc = cross(ab, ac);
        Vec3 acd = cross(ac, ad);
        Vec3 adb = cross(ad, ab);

        if (dot(abc, ao) > 0.0) {
            s.pts = {c, b, a};
            dir = abc;
            return false;
        }
        if (dot(acd, ao) > 0.0) {
            s.pts = {d, c, a};
            dir = acd;
            return false;
        }
        if (dot(adb, ao) > 0.0) {
            s.pts = {b, d, a};
            dir = adb;
            return false;
        }
        return true;
    }
    return false;
}

static bool gjk_intersect(const Mesh& a, const Transform& ta,
                          const Mesh& b, const Transform& tb,
                          Simplex& simplex) {
    Vec3 dir = ta.pos - tb.pos;
    if (norm(dir) < 1e-12) {
        dir = Vec3{1.0, 0.0, 0.0};
    }
    simplex.pts.clear();
    simplex.pts.push_back(support_minkowski(a, ta, b, tb, dir));
    dir = simplex.pts.back() * -1.0;

    for (int iter = 0; iter < 32; ++iter) {
        Vec3 p = support_minkowski(a, ta, b, tb, dir);
        if (dot(p, dir) < 0.0) {
            return false;
        }
        simplex.pts.push_back(p);
        if (handle_simplex(simplex, dir)) {
            return true;
        }
    }
    return false;
}

struct EPAFace {
    int a = 0;
    int b = 0;
    int c = 0;
    Vec3 n;
    double d = 0.0;
};

static EPAFace make_face(const std::vector<Vec3>& verts, int a, int b, int c) {
    Vec3 v0 = verts[a];
    Vec3 v1 = verts[b];
    Vec3 v2 = verts[c];
    Vec3 n = normalize(cross(v1 - v0, v2 - v0));
    double d = dot(n, v0);
    if (d < 0.0) {
        std::swap(b, c);
        n = n * -1.0;
        d = -d;
    }
    EPAFace f;
    f.a = a; f.b = b; f.c = c; f.n = n; f.d = d;
    return f;
}

static bool epa_penetration(const Mesh& a, const Transform& ta,
                            const Mesh& b, const Transform& tb,
                            const Simplex& simplex,
                            Vec3& out_normal, double& out_depth) {
    if (simplex.pts.size() < 4) {
        return false;
    }
    std::vector<Vec3> verts = simplex.pts;
    std::vector<EPAFace> faces;
    faces.push_back(make_face(verts, 0, 1, 2));
    faces.push_back(make_face(verts, 0, 3, 1));
    faces.push_back(make_face(verts, 0, 2, 3));
    faces.push_back(make_face(verts, 1, 3, 2));

    const double tol = 1e-6;
    for (int iter = 0; iter < 32; ++iter) {
        int best = -1;
        double min_d = std::numeric_limits<double>::infinity();
        for (int i = 0; i < static_cast<int>(faces.size()); ++i) {
            if (faces[i].d < min_d) {
                min_d = faces[i].d;
                best = i;
            }
        }
        if (best < 0) {
            return false;
        }
        EPAFace f = faces[best];
        Vec3 p = support_minkowski(a, ta, b, tb, f.n);
        double d = dot(p, f.n);
        if (d - f.d < tol) {
            out_normal = f.n;
            out_depth = d;
            return true;
        }

        int new_index = static_cast<int>(verts.size());
        verts.push_back(p);

        std::vector<std::pair<int, int>> edges;
        std::vector<EPAFace> new_faces;
        for (const auto& face : faces) {
            Vec3 v0 = verts[face.a];
            if (dot(face.n, p - v0) > 0.0) {
                auto add_edge = [&](int a, int b) {
                    auto it = std::find(edges.begin(), edges.end(), std::make_pair(b, a));
                    if (it != edges.end()) {
                        edges.erase(it);
                    } else {
                        edges.emplace_back(a, b);
                    }
                };
                add_edge(face.a, face.b);
                add_edge(face.b, face.c);
                add_edge(face.c, face.a);
            } else {
                new_faces.push_back(face);
            }
        }

        for (const auto& e : edges) {
            new_faces.push_back(make_face(verts, e.first, e.second, new_index));
        }
        faces.swap(new_faces);
    }
    return false;
}

static Mat3 inertia_world(const Particle& p) {
    Mat3 R = quat_to_mat3(p.tf.rot);
    Mat3 Rt = mat3_transpose(R);
    Mat3 Ibody;
    Ibody.m[0][0] = p.inertia_body.x;
    Ibody.m[1][1] = p.inertia_body.y;
    Ibody.m[2][2] = p.inertia_body.z;
    return mat3_mul(mat3_mul(R, Ibody), Rt);
}

static Mat3 inertia_world_inv(const Particle& p) {
    Mat3 R = quat_to_mat3(p.tf.rot);
    Mat3 Rt = mat3_transpose(R);
    Mat3 Iinv;
    Iinv.m[0][0] = p.inertia_body_inv.x;
    Iinv.m[1][1] = p.inertia_body_inv.y;
    Iinv.m[2][2] = p.inertia_body_inv.z;
    return mat3_mul(mat3_mul(R, Iinv), Rt);
}

static void integrate_particle(Particle& p, const Vec3& force, const Vec3& torque, const Vec3& gravity, double dt) {
    Vec3 acc = force * p.inv_mass + gravity;
    p.vel += acc * dt;
    p.tf.pos += p.vel * dt;

    Mat3 Iinv = inertia_world_inv(p);
    Vec3 ang_acc = mat3_mul_vec3(Iinv, torque);
    p.omega += ang_acc * dt;

    Quat dq{0.0, p.omega.x, p.omega.y, p.omega.z};
    Quat qdot = quat_mul(dq, p.tf.rot);
    p.tf.rot.w += 0.5 * qdot.w * dt;
    p.tf.rot.x += 0.5 * qdot.x * dt;
    p.tf.rot.y += 0.5 * qdot.y * dt;
    p.tf.rot.z += 0.5 * qdot.z * dt;
    p.tf.rot = quat_normalize(p.tf.rot);
}

struct SimConfig {
    int n = 2;
    int steps = 1;
    double dt = 1e-4;
    double spacing = 2.5;
    double kn = 1e5;
    double cn = 0.0;
    double mass = 1.0;
    Vec3 gravity{0.0, -9.81, 0.0};
    bool center_mesh = true;
    std::string stl_path;
    std::string vtk_prefix = "particles";
    int output_interval = 1;
    Vec3 v0{0.0, 0.0, 0.0};
    Vec3 v1{0.0, 0.0, 0.0};
};

static bool parse_config_file(const std::string& path, SimConfig& cfg) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Failed to open config: " << path << "\n";
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::string cleaned;
        cleaned.reserve(line.size());
        for (char c : line) {
            if (c == '#') {
                break;
            }
            if (c == '=') {
                cleaned.push_back(' ');
            } else {
                cleaned.push_back(c);
            }
        }
        std::istringstream iss(cleaned);
        std::string key;
        if (!(iss >> key)) {
            continue;
        }
        if (key == "stl") {
            iss >> cfg.stl_path;
        } else if (key == "n") {
            iss >> cfg.n;
        } else if (key == "steps") {
            iss >> cfg.steps;
        } else if (key == "dt") {
            iss >> cfg.dt;
        } else if (key == "spacing") {
            iss >> cfg.spacing;
        } else if (key == "kn") {
            iss >> cfg.kn;
        } else if (key == "cn") {
            iss >> cfg.cn;
        } else if (key == "mass") {
            iss >> cfg.mass;
        } else if (key == "gravity") {
            iss >> cfg.gravity.x >> cfg.gravity.y >> cfg.gravity.z;
        } else if (key == "center_mesh") {
            int v = 1;
            iss >> v;
            cfg.center_mesh = (v != 0);
        } else if (key == "vtk_prefix") {
            iss >> cfg.vtk_prefix;
        } else if (key == "output_interval") {
            iss >> cfg.output_interval;
        } else if (key == "v0") {
            iss >> cfg.v0.x >> cfg.v0.y >> cfg.v0.z;
        } else if (key == "v1") {
            iss >> cfg.v1.x >> cfg.v1.y >> cfg.v1.z;
        }
    }
    if (cfg.output_interval < 1) {
        cfg.output_interval = 1;
    }
    return !cfg.stl_path.empty();
}

static bool parse_args(int argc, char** argv, std::string& config_path) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--config" && i + 1 < argc) {
            config_path = argv[++i];
        } else if (arg == "--help") {
            return false;
        } else {
            std::cerr << "Unknown arg: " << arg << "\n";
            return false;
        }
    }
    return !config_path.empty();
}

static void usage() {
    std::cout << "Usage: zdem_cpu --config path.txt\n";
}

static void write_vtk_particles(const std::string& path,
                                const std::vector<Particle>& particles,
                                const std::vector<Vec3>& forces,
                                const std::vector<Vec3>& torques,
                                const std::vector<int>& contact_counts) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to write VTK: " << path << "\n";
        return;
    }
    out << "# vtk DataFile Version 3.0\n";
    out << "zdem particles\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << particles.size() << " float\n";
    for (const auto& p : particles) {
        out << static_cast<float>(p.tf.pos.x) << " "
            << static_cast<float>(p.tf.pos.y) << " "
            << static_cast<float>(p.tf.pos.z) << "\n";
    }

    out << "POINT_DATA " << particles.size() << "\n";

    out << "SCALARS id int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        out << i << "\n";
    }

    out << "SCALARS mass float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        out << static_cast<float>(p.mass) << "\n";
    }

    out << "SCALARS radius float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        out << static_cast<float>(p.radius) << "\n";
    }

    out << "SCALARS contact_count int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int c : contact_counts) {
        out << c << "\n";
    }

    out << "VECTORS velocity float\n";
    for (const auto& p : particles) {
        out << static_cast<float>(p.vel.x) << " "
            << static_cast<float>(p.vel.y) << " "
            << static_cast<float>(p.vel.z) << "\n";
    }

    out << "VECTORS omega float\n";
    for (const auto& p : particles) {
        out << static_cast<float>(p.omega.x) << " "
            << static_cast<float>(p.omega.y) << " "
            << static_cast<float>(p.omega.z) << "\n";
    }

    out << "VECTORS force float\n";
    for (const auto& f : forces) {
        out << static_cast<float>(f.x) << " "
            << static_cast<float>(f.y) << " "
            << static_cast<float>(f.z) << "\n";
    }

    out << "VECTORS torque float\n";
    for (const auto& t : torques) {
        out << static_cast<float>(t.x) << " "
            << static_cast<float>(t.y) << " "
            << static_cast<float>(t.z) << "\n";
    }

    out << "FIELD FieldData 1\n";
    out << "orientation 4 " << particles.size() << " float\n";
    for (const auto& p : particles) {
        out << static_cast<float>(p.tf.rot.w) << " "
            << static_cast<float>(p.tf.rot.x) << " "
            << static_cast<float>(p.tf.rot.y) << " "
            << static_cast<float>(p.tf.rot.z) << "\n";
    }
}

int main(int argc, char** argv) {
    SimConfig cfg;
    std::string config_path;
    if (!parse_args(argc, argv, config_path)) {
        usage();
        return 1;
    }
    if (!parse_config_file(config_path, cfg)) {
        return 1;
    }

    std::vector<std::array<Vec3, 3>> tris;
    if (!load_stl(cfg.stl_path, tris)) {
        std::cerr << "Failed to load STL: " << cfg.stl_path << "\n";
        return 1;
    }

    Mesh mesh = build_mesh(tris, cfg.center_mesh);
    if (mesh.vertices.empty()) {
        std::cerr << "Mesh has no vertices.\n";
        return 1;
    }

    std::vector<Particle> particles(cfg.n);
    for (int i = 0; i < cfg.n; ++i) {
        Particle p;
        p.tf.pos = Vec3{cfg.spacing * i, 0.0, 0.0};
        p.tf.rot = Quat{};
        p.vel = Vec3{0.0, 0.0, 0.0};
        p.omega = Vec3{0.0, 0.0, 0.0};
        p.mass = cfg.mass;
        p.inv_mass = 1.0 / cfg.mass;
        double I = 0.4 * p.mass * mesh.radius * mesh.radius;
        p.inertia_body = Vec3{I, I, I};
        p.inertia_body_inv = Vec3{1.0 / I, 1.0 / I, 1.0 / I};
        p.radius = mesh.radius;
        particles[i] = p;
    }

    if (cfg.n >= 1) {
        particles[0].vel = cfg.v0;
    }
    if (cfg.n >= 2) {
        particles[1].vel = cfg.v1;
    }

    std::vector<Vec3> forces(cfg.n);
    std::vector<Vec3> torques(cfg.n);
    std::vector<int> contact_counts(cfg.n);

    for (int step = 0; step < cfg.steps; ++step) {
        std::fill(forces.begin(), forces.end(), Vec3{});
        std::fill(torques.begin(), torques.end(), Vec3{});
        std::fill(contact_counts.begin(), contact_counts.end(), 0);

        int contacts = 0;
        for (int i = 0; i < cfg.n; ++i) {
            for (int j = i + 1; j < cfg.n; ++j) {
                Particle& pa = particles[i];
                Particle& pb = particles[j];
                Vec3 dpos = pb.tf.pos - pa.tf.pos;
                double dist2 = dot(dpos, dpos);
                double rsum = pa.radius + pb.radius;
                if (dist2 > rsum * rsum) {
                    continue;
                }
                Simplex simplex;
                if (!gjk_intersect(mesh, pa.tf, mesh, pb.tf, simplex)) {
                    continue;
                }

                Vec3 normal{0.0, 1.0, 0.0};
                double depth = 0.0;
                if (!epa_penetration(mesh, pa.tf, mesh, pb.tf, simplex, normal, depth)) {
                    Vec3 fallback = normalize(dpos);
                    if (norm(fallback) > 0.0) {
                        normal = fallback;
                        depth = rsum - std::sqrt(dist2);
                    }
                }

                Vec3 rel_vel = pb.vel - pa.vel;
                double vn = dot(rel_vel, normal);
                double fn = cfg.kn * depth - cfg.cn * vn;
                if (fn < 0.0) {
                    fn = 0.0;
                }
                Vec3 f = normal * fn;
                forces[i] -= f;
                forces[j] += f;
                contact_counts[i] += 1;
                contact_counts[j] += 1;
                contacts++;
            }
        }

        for (int i = 0; i < cfg.n; ++i) {
            integrate_particle(particles[i], forces[i], torques[i], cfg.gravity, cfg.dt);
        }

        if (step % cfg.output_interval == 0) {
            std::ostringstream oss;
            oss << cfg.vtk_prefix << "_" << std::setw(6) << std::setfill('0') << step << ".vtk";
            write_vtk_particles(oss.str(), particles, forces, torques, contact_counts);
        }

        std::cout << "step=" << step << " contacts=" << contacts << "\n";
    }

    std::cout << "done\n";
    return 0;
}
