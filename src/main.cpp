#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <chrono>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

static inline double norm2(const Vec3& v) {
    return dot(v, v);
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

static inline Mat3 mat3_identity() {
    Mat3 r;
    r.m[0][0] = 1.0;
    r.m[1][1] = 1.0;
    r.m[2][2] = 1.0;
    return r;
}

static inline Mat3 mat3_add(const Mat3& a, const Mat3& b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = a.m[i][j] + b.m[i][j];
        }
    }
    return r;
}

static inline Mat3 mat3_sub(const Mat3& a, const Mat3& b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = a.m[i][j] - b.m[i][j];
        }
    }
    return r;
}

static inline Mat3 mat3_scale(const Mat3& a, double s) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r.m[i][j] = a.m[i][j] * s;
        }
    }
    return r;
}

static inline Mat3 mat3_outer(const Vec3& a, const Vec3& b) {
    Mat3 r;
    r.m[0][0] = a.x * b.x; r.m[0][1] = a.x * b.y; r.m[0][2] = a.x * b.z;
    r.m[1][0] = a.y * b.x; r.m[1][1] = a.y * b.y; r.m[1][2] = a.y * b.z;
    r.m[2][0] = a.z * b.x; r.m[2][1] = a.z * b.y; r.m[2][2] = a.z * b.z;
    return r;
}

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

static inline double mat3_det(const Mat3& a) {
    return a.m[0][0] * (a.m[1][1] * a.m[2][2] - a.m[1][2] * a.m[2][1])
         - a.m[0][1] * (a.m[1][0] * a.m[2][2] - a.m[1][2] * a.m[2][0])
         + a.m[0][2] * (a.m[1][0] * a.m[2][1] - a.m[1][1] * a.m[2][0]);
}

static inline Mat3 mat3_inverse(const Mat3& a) {
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

struct Transform {
    Vec3 pos;
    Quat rot;
};

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<Vec3, 3>> tris;
    Vec3 center;
    double radius = 0.0;
    double mean_edge = 0.0;
    double bbox_diag = 0.0;
    double volume = 0.0;
    Mat3 inertia_unit;
};

struct ParticleInit {
    std::string stl_path;
    Vec3 pos{0.0, 0.0, 0.0};
    Vec3 vel{0.0, 0.0, 0.0};
    Vec3 omega{0.0, 0.0, 0.0};
    Quat rot{1.0, 0.0, 0.0, 0.0};
    double scale = 0.0;
    double density = 1.0;
    double young = 1e7;
    double poisson = 0.25;
    double mu = 0.5;
    double restitution = 0.5;
};

struct Particle {
    Transform tf;
    Vec3 vel;
    Vec3 omega;
    Vec3 L;
    double mass = 1.0;
    double inv_mass = 1.0;
    Mat3 inertia_body;
    Mat3 inertia_body_inv;
    double radius = 0.0;
    double equiv_radius = 0.0;
    double young = 0.0;
    double poisson = 0.0;
    double mu = 0.0;
    double restitution = 0.0;
    int mesh_index = 0;
};

struct WallInit {
    std::string stl_path;
    Vec3 pos{0.0, 0.0, 0.0};
    Quat rot{1.0, 0.0, 0.0, 0.0};
    double scale = 1.0;
    double mu = 0.5;
    double restitution = 0.5;
};

struct Wall {
    Mesh mesh;
    Transform tf;
    double mu = 0.5;
    double restitution = 0.5;
    std::vector<Vec3> tri_normals;
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

static void compute_mass_properties(const std::vector<std::array<Vec3, 3>>& tris,
                                    double& out_volume,
                                    Vec3& out_centroid,
                                    Mat3& out_inertia) {
    double det_sum = 0.0;
    Vec3 csum{0.0, 0.0, 0.0};
    double Ixx = 0.0, Iyy = 0.0, Izz = 0.0;
    double Ixy = 0.0, Ixz = 0.0, Iyz = 0.0;

    for (const auto& tri : tris) {
        Vec3 a = tri[0];
        Vec3 b = tri[1];
        Vec3 c = tri[2];
        double det = dot(a, cross(b, c));
        det_sum += det;
        csum += (a + b + c) * det;

        double x1 = a.x, x2 = b.x, x3 = c.x;
        double y1 = a.y, y2 = b.y, y3 = c.y;
        double z1 = a.z, z2 = b.z, z3 = c.z;

        double f2x = x1 * x1 + x2 * x2 + x3 * x3 + x1 * x2 + x1 * x3 + x2 * x3;
        double f2y = y1 * y1 + y2 * y2 + y3 * y3 + y1 * y2 + y1 * y3 + y2 * y3;
        double f2z = z1 * z1 + z2 * z2 + z3 * z3 + z1 * z2 + z1 * z3 + z2 * z3;

        Ixx += det * (f2y + f2z);
        Iyy += det * (f2x + f2z);
        Izz += det * (f2x + f2y);

        double gxy = 2.0 * x1 * y1 + 2.0 * x2 * y2 + 2.0 * x3 * y3
                   + x1 * y2 + x2 * y1 + x1 * y3 + x3 * y1 + x2 * y3 + x3 * y2;
        double gxz = 2.0 * x1 * z1 + 2.0 * x2 * z2 + 2.0 * x3 * z3
                   + x1 * z2 + x2 * z1 + x1 * z3 + x3 * z1 + x2 * z3 + x3 * z2;
        double gyz = 2.0 * y1 * z1 + 2.0 * y2 * z2 + 2.0 * y3 * z3
                   + y1 * z2 + y2 * z1 + y1 * z3 + y3 * z1 + y2 * z3 + y3 * z2;
        Ixy += det * gxy;
        Ixz += det * gxz;
        Iyz += det * gyz;
    }

    if (std::abs(det_sum) < 1e-18) {
        out_volume = 0.0;
        out_centroid = Vec3{0.0, 0.0, 0.0};
        out_inertia = Mat3{};
        return;
    }

    if (det_sum < 0.0) {
        det_sum = -det_sum;
        csum *= -1.0;
        Ixx *= -1.0; Iyy *= -1.0; Izz *= -1.0;
        Ixy *= -1.0; Ixz *= -1.0; Iyz *= -1.0;
    }

    out_volume = det_sum / 6.0;
    out_centroid = csum / (4.0 * det_sum);

    Ixx /= 60.0;
    Iyy /= 60.0;
    Izz /= 60.0;
    Ixy /= -120.0;
    Ixz /= -120.0;
    Iyz /= -120.0;

    Mat3 I;
    I.m[0][0] = Ixx; I.m[1][1] = Iyy; I.m[2][2] = Izz;
    I.m[0][1] = Ixy; I.m[1][0] = Ixy;
    I.m[0][2] = Ixz; I.m[2][0] = Ixz;
    I.m[1][2] = Iyz; I.m[2][1] = Iyz;
    out_inertia = I;
}

static Mesh build_mesh(const std::vector<std::array<Vec3, 3>>& tris, bool center_mesh) {
    Mesh m;
    m.vertices = unique_vertices(tris);
    m.tris = tris;
    Vec3 centroid{0.0, 0.0, 0.0};
    Mat3 inertia_origin;
    compute_mass_properties(m.tris, m.volume, centroid, inertia_origin);
    m.center = centroid;
    Mat3 inertia_cm = inertia_origin;
    if (m.volume > 0.0) {
        Mat3 shift = mat3_sub(mat3_scale(mat3_identity(), dot(centroid, centroid)), mat3_outer(centroid, centroid));
        inertia_cm = mat3_sub(inertia_origin, mat3_scale(shift, m.volume));
    }
    if (center_mesh) {
        for (auto& v : m.vertices) {
            v -= centroid;
        }
        for (auto& tri : m.tris) {
            tri[0] -= centroid;
            tri[1] -= centroid;
            tri[2] -= centroid;
        }
        m.center = Vec3{0.0, 0.0, 0.0};
    }
    m.inertia_unit = inertia_cm;
    double r2 = 0.0;
    Vec3 mn{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Vec3 mx{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (const auto& v : m.vertices) {
        r2 = std::max(r2, dot(v, v));
        mn.x = std::min(mn.x, v.x);
        mn.y = std::min(mn.y, v.y);
        mn.z = std::min(mn.z, v.z);
        mx.x = std::max(mx.x, v.x);
        mx.y = std::max(mx.y, v.y);
        mx.z = std::max(mx.z, v.z);
    }
    m.radius = std::sqrt(r2);
    m.bbox_diag = norm(mx - mn);

    double edge_sum = 0.0;
    std::size_t edge_cnt = 0;
    for (const auto& tri : m.tris) {
        edge_sum += norm(tri[1] - tri[0]);
        edge_sum += norm(tri[2] - tri[1]);
        edge_sum += norm(tri[0] - tri[2]);
        edge_cnt += 3;
    }
    if (edge_cnt > 0) {
        m.mean_edge = edge_sum / static_cast<double>(edge_cnt);
    }
    return m;
}

static std::vector<std::array<Vec3, 3>> transform_tris(const Mesh& mesh, const Transform& tf) {
    std::vector<std::array<Vec3, 3>> out;
    out.reserve(mesh.tris.size());
    for (const auto& tri : mesh.tris) {
        std::array<Vec3, 3> t;
        t[0] = quat_rotate(tf.rot, tri[0]) + tf.pos;
        t[1] = quat_rotate(tf.rot, tri[1]) + tf.pos;
        t[2] = quat_rotate(tf.rot, tri[2]) + tf.pos;
        out.push_back(t);
    }
    return out;
}

static Vec3 tri_normal(const std::array<Vec3, 3>& tri) {
    Vec3 e1 = tri[1] - tri[0];
    Vec3 e2 = tri[2] - tri[0];
    return normalize(cross(e1, e2));
}

static bool plane_from_tri(const std::array<Vec3, 3>& tri, Vec3& n, double& d) {
    n = cross(tri[1] - tri[0], tri[2] - tri[0]);
    double ln = norm(n);
    if (ln < 1e-30) {
        return false;
    }
    n = n / ln;
    d = -dot(n, tri[0]);
    return true;
}

static bool point_in_tri(const Vec3& p, const std::array<Vec3, 3>& tri, const Vec3& n) {
    Vec3 a = tri[0];
    Vec3 b = tri[1];
    Vec3 c = tri[2];
    Vec3 ab = b - a;
    Vec3 bc = c - b;
    Vec3 ca = a - c;
    Vec3 ap = p - a;
    Vec3 bp = p - b;
    Vec3 cp = p - c;
    double c1 = dot(cross(ab, ap), n);
    double c2 = dot(cross(bc, bp), n);
    double c3 = dot(cross(ca, cp), n);
    return (c1 >= -1e-10) && (c2 >= -1e-10) && (c3 >= -1e-10);
}

static bool seg_plane_intersection(const Vec3& p0, const Vec3& p1, const Vec3& n, double d, Vec3& out) {
    Vec3 u = p1 - p0;
    double denom = dot(n, u);
    double num = -(dot(n, p0) + d);
    if (std::abs(denom) < 1e-30) {
        return false;
    }
    double t = num / denom;
    if (t < -1e-10 || t > 1.0 + 1e-10) {
        return false;
    }
    t = std::min(1.0, std::max(0.0, t));
    out = p0 + u * t;
    return true;
}

static std::vector<Vec3> unique_points(const std::vector<Vec3>& pts, double tol) {
    std::vector<Vec3> out;
    double tol2 = tol * tol;
    for (const auto& p : pts) {
        bool keep = true;
        for (const auto& q : out) {
            if (norm2(p - q) <= tol2) {
                keep = false;
                break;
            }
        }
        if (keep) {
            out.push_back(p);
        }
    }
    return out;
}

static bool tri_tri_intersection_segment(const std::array<Vec3, 3>& triA, const Vec3& nA, double dA,
                                         const std::array<Vec3, 3>& triB, const Vec3& nB, double dB,
                                         Vec3& out_p0, Vec3& out_p1) {
    Vec3 tau = cross(nA, nB);
    double lt = norm(tau);
    if (lt < 1e-20) {
        return false;
    }
    tau = tau / lt;

    std::vector<Vec3> candidates;
    for (int i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        Vec3 p;
        if (seg_plane_intersection(triA[i], triA[j], nB, dB, p) && point_in_tri(p, triB, nB)) {
            candidates.push_back(p);
        }
    }
    for (int i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        Vec3 p;
        if (seg_plane_intersection(triB[i], triB[j], nA, dA, p) && point_in_tri(p, triA, nA)) {
            candidates.push_back(p);
        }
    }

    candidates = unique_points(candidates, 1e-9);
    if (candidates.size() < 2) {
        return false;
    }
    double min_t = std::numeric_limits<double>::infinity();
    double max_t = -std::numeric_limits<double>::infinity();
    int i0 = -1, i1 = -1;
    for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
        double t = dot(candidates[i], tau);
        if (t < min_t) { min_t = t; i0 = i; }
        if (t > max_t) { max_t = t; i1 = i; }
    }
    if (i0 == i1) {
        return false;
    }
    Vec3 p0 = candidates[i0];
    Vec3 p1 = candidates[i1];
    if (norm2(p1 - p0) < 1e-18) {
        return false;
    }
    if (dot(p1 - p0, tau) < 0.0) {
        std::swap(p0, p1);
    }
    out_p0 = p0;
    out_p1 = p1;
    return true;
}

struct Int3 {
    int x = 0;
    int y = 0;
    int z = 0;
    bool operator==(const Int3& o) const { return x == o.x && y == o.y && z == o.z; }
};

struct Int3Hash {
    std::size_t operator()(const Int3& k) const {
        std::size_t h1 = std::hash<int>{}(k.x);
        std::size_t h2 = std::hash<int>{}(k.y);
        std::size_t h3 = std::hash<int>{}(k.z);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

static Int3 cell_index(const Vec3& p, double h) {
    return Int3{
        static_cast<int>(std::floor(p.x / h)),
        static_cast<int>(std::floor(p.y / h)),
        static_cast<int>(std::floor(p.z / h))
    };
}

struct UnionFind {
    std::vector<int> parent;
    std::vector<int> rank;

    explicit UnionFind(int n = 0) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }

    void unite(int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) return;
        if (rank[ra] < rank[rb]) parent[ra] = rb;
        else if (rank[ra] > rank[rb]) parent[rb] = ra;
        else { parent[rb] = ra; rank[ra] += 1; }
    }
};

static void snap_endpoints(const std::vector<std::pair<Vec3, Vec3>>& segments,
                           double tol,
                           std::vector<Vec3>& out_vertices,
                           std::vector<std::pair<int, int>>& out_edges) {
    double h = tol;
    int nseg = static_cast<int>(segments.size());
    std::vector<Vec3> pts(2 * nseg);
    for (int i = 0; i < nseg; ++i) {
        pts[2 * i] = segments[i].first;
        pts[2 * i + 1] = segments[i].second;
    }

    UnionFind uf(2 * nseg);
    std::unordered_map<Int3, std::vector<int>, Int3Hash> grid;

    for (int i = 0; i < 2 * nseg; ++i) {
        const Vec3& p = pts[i];
        Int3 c = cell_index(p, h);
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    Int3 cc{c.x + dx, c.y + dy, c.z + dz};
                    auto it = grid.find(cc);
                    if (it == grid.end()) continue;
                    for (int j : it->second) {
                        if (norm2(p - pts[j]) <= h * h) {
                            uf.unite(i, j);
                        }
                    }
                }
            }
        }
        grid[c].push_back(i);
    }

    std::unordered_map<int, int> root_to_vid;
    std::unordered_map<int, Vec3> sums;
    std::unordered_map<int, int> counts;
    out_vertices.clear();
    for (int i = 0; i < 2 * nseg; ++i) {
        int r = uf.find(i);
        if (root_to_vid.find(r) == root_to_vid.end()) {
            int vid = static_cast<int>(root_to_vid.size());
            root_to_vid[r] = vid;
            sums[r] = pts[i];
            counts[r] = 1;
        } else {
            sums[r] += pts[i];
            counts[r] += 1;
        }
    }
    out_vertices.resize(root_to_vid.size());
    for (const auto& kv : root_to_vid) {
        int r = kv.first;
        int vid = kv.second;
        out_vertices[vid] = sums[r] / static_cast<double>(counts[r]);
    }

    out_edges.clear();
    for (int s = 0; s < nseg; ++s) {
        int a = root_to_vid[uf.find(2 * s)];
        int b = root_to_vid[uf.find(2 * s + 1)];
        if (a != b) {
            out_edges.emplace_back(a, b);
        }
    }
}

static std::pair<double, double> point_segment_dist2_and_t(const Vec3& p, const Vec3& a, const Vec3& b) {
    Vec3 ab = b - a;
    double denom = dot(ab, ab);
    if (denom <= 1e-12) {
        return {norm2(p - a), 0.0};
    }
    double t = dot(p - a, ab) / denom;
    double t_clamped = std::min(1.0, std::max(0.0, t));
    Vec3 proj = a + ab * t_clamped;
    double d2 = norm2(p - proj);
    return {d2, t_clamped};
}

static double seg_seg_dist2(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2) {
    Vec3 d1 = q1 - p1;
    Vec3 d2 = q2 - p2;
    Vec3 r = p1 - p2;
    double a = dot(d1, d1);
    double e = dot(d2, d2);
    double f = dot(d2, r);

    if (a <= 1e-12 && e <= 1e-12) {
        return dot(p1 - p2, p1 - p2);
    }
    if (a <= 1e-12) {
        double t = std::min(1.0, std::max(0.0, f / e));
        Vec3 c = p2 + d2 * t;
        return dot(p1 - c, p1 - c);
    }
    if (e <= 1e-12) {
        double s = std::min(1.0, std::max(0.0, -dot(d1, r) / a));
        Vec3 c = p1 + d1 * s;
        return dot(c - p2, c - p2);
    }

    double b = dot(d1, d2);
    double c = dot(d1, r);
    double denom = a * e - b * b;
    double s = 0.0;
    if (std::abs(denom) > 0.0) {
        s = std::min(1.0, std::max(0.0, (b * f - c * e) / denom));
    }
    double t = (b * s + f) / e;
    if (t < 0.0) {
        t = 0.0;
        s = std::min(1.0, std::max(0.0, -c / a));
    } else if (t > 1.0) {
        t = 1.0;
        s = std::min(1.0, std::max(0.0, (b - c) / a));
    }
    Vec3 cp1 = p1 + d1 * s;
    Vec3 cp2 = p2 + d2 * t;
    Vec3 diff = cp1 - cp2;
    return dot(diff, diff);
}

static void label_segments_by_contact_greedy(const std::vector<std::pair<Vec3, Vec3>>& segments,
                                             double tol,
                                             std::vector<std::vector<int>>& comps) {
    double h = tol;
    int nseg = static_cast<int>(segments.size());
    comps.clear();
    if (nseg == 0) {
        return;
    }

    std::unordered_map<Int3, std::vector<int>, Int3Hash> grid;
    std::vector<std::pair<Vec3, Vec3>> seg_aabb;
    std::vector<std::array<Vec3, 3>> rep_points;
    seg_aabb.reserve(nseg);
    rep_points.reserve(nseg);

    for (int i = 0; i < nseg; ++i) {
        Vec3 p0 = segments[i].first;
        Vec3 p1 = segments[i].second;
        Vec3 mid = (p0 + p1) * 0.5;
        rep_points.push_back({p0, p1, mid});
        Vec3 mn{std::min(p0.x, p1.x) - h, std::min(p0.y, p1.y) - h, std::min(p0.z, p1.z) - h};
        Vec3 mx{std::max(p0.x, p1.x) + h, std::max(p0.y, p1.y) + h, std::max(p0.z, p1.z) + h};
        seg_aabb.push_back({mn, mx});
        for (const auto& rp : rep_points.back()) {
            Int3 c = cell_index(rp, h);
            grid[c].push_back(i);
        }
    }

    auto nearby_candidates = [&](int i) {
        std::vector<int> cand;
        std::unordered_map<int, bool> mark;
        for (const auto& rp : rep_points[i]) {
            Int3 c = cell_index(rp, h);
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        Int3 cc{c.x + dx, c.y + dy, c.z + dz};
                        auto it = grid.find(cc);
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            if (j == i) continue;
                            if (!mark[j]) {
                                mark[j] = true;
                                cand.push_back(j);
                            }
                        }
                    }
                }
            }
        }
        return cand;
    };

    auto aabb_overlap = [&](int i, int j) {
        const Vec3& mn1 = seg_aabb[i].first;
        const Vec3& mx1 = seg_aabb[i].second;
        const Vec3& mn2 = seg_aabb[j].first;
        const Vec3& mx2 = seg_aabb[j].second;
        return !(mx1.x < mn2.x || mn1.x > mx2.x ||
                 mx1.y < mn2.y || mn1.y > mx2.y ||
                 mx1.z < mn2.z || mn1.z > mx2.z);
    };

    double tol2 = h * h;
    std::vector<bool> visited(nseg, false);
    for (int seed = 0; seed < nseg; ++seed) {
        if (visited[seed]) continue;
        std::vector<int> queue;
        queue.push_back(seed);
        visited[seed] = true;
        std::vector<int> comp;
        comp.push_back(seed);
        while (!queue.empty()) {
            int i = queue.back();
            queue.pop_back();
            Vec3 p1 = segments[i].first;
            Vec3 q1 = segments[i].second;
            for (int j : nearby_candidates(i)) {
                if (visited[j]) continue;
                if (!aabb_overlap(i, j)) continue;
                Vec3 p2 = segments[j].first;
                Vec3 q2 = segments[j].second;
                if (seg_seg_dist2(p1, q1, p2, q2) <= tol2) {
                    visited[j] = true;
                    queue.push_back(j);
                    comp.push_back(j);
                }
            }
        }
        comps.push_back(comp);
    }
    std::sort(comps.begin(), comps.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
}

static void split_t_junctions(const std::vector<Vec3>& vertices,
                              const std::vector<std::pair<int, int>>& edges,
                              double tol,
                              std::vector<Vec3>& out_vertices,
                              std::vector<std::pair<int, int>>& out_edges) {
    double h = tol;
    double tol2 = h * h;
    std::unordered_map<Int3, std::vector<int>, Int3Hash> vgrid;
    for (int vid = 0; vid < static_cast<int>(vertices.size()); ++vid) {
        Int3 c = cell_index(vertices[vid], h);
        vgrid[c].push_back(vid);
    }

    out_vertices = vertices;
    out_edges.clear();
    for (const auto& e : edges) {
        int a = e.first;
        int b = e.second;
        Vec3 pa = vertices[a];
        Vec3 pb = vertices[b];
        Vec3 mn{std::min(pa.x, pb.x) - h, std::min(pa.y, pb.y) - h, std::min(pa.z, pb.z) - h};
        Vec3 mx{std::max(pa.x, pb.x) + h, std::max(pa.y, pb.y) + h, std::max(pa.z, pb.z) + h};
        Int3 cmin = cell_index(mn, h);
        Int3 cmax = cell_index(mx, h);

        std::vector<std::pair<double, int>> splits;
        for (int ix = cmin.x; ix <= cmax.x; ++ix) {
            for (int iy = cmin.y; iy <= cmax.y; ++iy) {
                for (int iz = cmin.z; iz <= cmax.z; ++iz) {
                    auto it = vgrid.find(Int3{ix, iy, iz});
                    if (it == vgrid.end()) continue;
                    for (int vid : it->second) {
                        if (vid == a || vid == b) continue;
                        auto dt = point_segment_dist2_and_t(vertices[vid], pa, pb);
                        if (dt.first <= tol2 && dt.second > 1e-6 && dt.second < 1.0 - 1e-6) {
                            splits.emplace_back(dt.second, vid);
                        }
                    }
                }
            }
        }
        if (splits.empty()) {
            out_edges.push_back(e);
            continue;
        }
        std::sort(splits.begin(), splits.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        std::vector<int> chain;
        chain.push_back(a);
        int last = -1;
        for (const auto& s : splits) {
            if (s.second == last) continue;
            chain.push_back(s.second);
            last = s.second;
        }
        chain.push_back(b);
        for (int i = 0; i + 1 < static_cast<int>(chain.size()); ++i) {
            if (chain[i] != chain[i + 1]) {
                out_edges.emplace_back(chain[i], chain[i + 1]);
            }
        }
    }
}

static std::vector<std::vector<int>> extract_loops(const std::vector<Vec3>& vertices,
                                                   const std::vector<std::pair<int, int>>& edges) {
    int nv = static_cast<int>(vertices.size());
    std::vector<std::vector<std::pair<int, int>>> adj(nv);
    for (int ei = 0; ei < static_cast<int>(edges.size()); ++ei) {
        int a = edges[ei].first;
        int b = edges[ei].second;
        adj[a].push_back({b, ei});
        adj[b].push_back({a, ei});
    }

    std::vector<bool> used(edges.size(), false);
    std::vector<int> start_vertices(nv);
    for (int i = 0; i < nv; ++i) start_vertices[i] = i;
    std::sort(start_vertices.begin(), start_vertices.end(), [&](int i, int j) {
        int di = static_cast<int>(adj[i].size());
        int dj = static_cast<int>(adj[j].size());
        if ((di == 2) != (dj == 2)) return di == 2;
        return di > dj;
    });

    std::vector<std::vector<int>> loops;
    for (int sv : start_vertices) {
        for (const auto& nb : adj[sv]) {
            int ei = nb.second;
            if (used[ei]) continue;
            std::vector<int> path;
            path.push_back(sv);
            int prev = sv;
            int cur = nb.first;
            used[ei] = true;
            path.push_back(cur);

            int steps = 0;
            while (steps < 100000) {
                steps += 1;
                if (cur == sv) {
                    if (path.size() > 1) {
                        path.pop_back();
                        loops.push_back(path);
                    }
                    break;
                }
                std::vector<std::pair<int, int>> cand;
                for (const auto& e : adj[cur]) {
                    if (!used[e.second]) cand.push_back(e);
                }
                if (cand.empty()) {
                    break;
                }
                std::pair<int, int> pick = cand[0];
                for (const auto& c : cand) {
                    if (c.first != prev) { pick = c; break; }
                }
                used[pick.second] = true;
                prev = cur;
                cur = pick.first;
                path.push_back(cur);
            }
        }
    }
    return loops;
}

static std::pair<Vec3, Vec3> accumulate_Sn_Gn_from_polyline(const std::vector<Vec3>& loop_pts) {
    Vec3 Sn{0.0, 0.0, 0.0};
    Vec3 Gn{0.0, 0.0, 0.0};
    int m = static_cast<int>(loop_pts.size());
    for (int i = 0; i < m; ++i) {
        const Vec3& p0 = loop_pts[i];
        const Vec3& p1 = loop_pts[(i + 1) % m];
        Sn += cross(p0, p1) * 0.5;
        Vec3 dx = p1 - p0;
        double s = dot(p0, p1) + (1.0 / 3.0) * dot(dx, dx);
        Gn += dx * (-(1.0 / 3.0) * s);
    }
    return {Sn, Gn};
}

static bool contact_point_xc0(const Vec3& Sn, const Vec3& Gn, Vec3& xc0, Vec3& nA, double& area) {
    area = norm(Sn);
    if (area < 1e-14) {
        xc0 = Vec3{0.0, 0.0, 0.0};
        nA = Vec3{0.0, 0.0, 0.0};
        return false;
    }
    nA = Sn * (-1.0 / area);
    xc0 = cross(nA, Gn) / area;
    return true;
}

// ============== Sutherland-Hodgman Polygon Clipping for Wall Contact ==============

// 2D point for clipping algorithm
struct Vec2 {
    double x = 0.0;
    double y = 0.0;
    Vec2() = default;
    Vec2(double x_, double y_) : x(x_), y(y_) {}
    Vec2 operator+(const Vec2& o) const { return Vec2{x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2& o) const { return Vec2{x - o.x, y - o.y}; }
    Vec2 operator*(double s) const { return Vec2{x * s, y * s}; }
};

static inline double cross2d(const Vec2& a, const Vec2& b) {
    return a.x * b.y - a.y * b.x;
}

// Check if point p is inside the clipping edge (on the left side of edge from cp1 to cp2)
static inline bool is_inside_clip(const Vec2& p, const Vec2& cp1, const Vec2& cp2) {
    return cross2d(cp2 - cp1, p - cp1) >= 0.0;
}

// Compute intersection of line segment (s, e) with clipping edge (cp1, cp2)
static inline Vec2 compute_intersection_clip(const Vec2& s, const Vec2& e,
                                              const Vec2& cp1, const Vec2& cp2) {
    Vec2 dc = cp1 - cp2;
    Vec2 dp = s - e;
    double n1 = cross2d(cp1, cp2);
    double n2 = cross2d(s, e);
    double n3 = cross2d(dc, dp);
    if (std::abs(n3) < 1e-30) {
        return s;
    }
    double t = (n1 * dp.x - n2 * dc.x) / n3;
    double u = (n1 * dp.y - n2 * dc.y) / n3;
    return Vec2{t, u};
}

// Sutherland-Hodgman polygon clipping algorithm
static std::vector<Vec2> sutherland_hodgman_clip(const std::vector<Vec2>& subject,
                                                  const std::vector<Vec2>& clip) {
    std::vector<Vec2> output = subject;
    if (clip.size() < 3 || subject.size() < 3) {
        return output;
    }

    Vec2 cp1 = clip.back();
    for (const auto& cp2 : clip) {
        if (output.empty()) break;

        std::vector<Vec2> input = output;
        output.clear();

        Vec2 s = input.back();
        for (const auto& e : input) {
            if (is_inside_clip(e, cp1, cp2)) {
                if (!is_inside_clip(s, cp1, cp2)) {
                    output.push_back(compute_intersection_clip(s, e, cp1, cp2));
                }
                output.push_back(e);
            } else if (is_inside_clip(s, cp1, cp2)) {
                output.push_back(compute_intersection_clip(s, e, cp1, cp2));
            }
            s = e;
        }
        cp1 = cp2;
    }
    return output;
}

// Compute 2D polygon area using shoelace formula
static double polygon_area_2d(const std::vector<Vec2>& poly) {
    if (poly.size() < 3) return 0.0;
    double area = 0.0;
    int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += cross2d(poly[i], poly[j]);
    }
    return std::abs(area) * 0.5;
}

// Compute 2D polygon centroid
static Vec2 polygon_centroid_2d(const std::vector<Vec2>& poly) {
    if (poly.empty()) return Vec2{0.0, 0.0};
    Vec2 c{0.0, 0.0};
    for (const auto& p : poly) {
        c.x += p.x;
        c.y += p.y;
    }
    double n = static_cast<double>(poly.size());
    return Vec2{c.x / n, c.y / n};
}

// Build orthonormal basis from a normal vector
static void build_basis(const Vec3& n, Vec3& u, Vec3& v) {
    // Choose a vector not parallel to n
    if (std::abs(n.x) < 0.9) {
        u = normalize(cross(n, Vec3{1.0, 0.0, 0.0}));
    } else {
        u = normalize(cross(n, Vec3{0.0, 1.0, 0.0}));
    }
    v = cross(n, u);
}

// Project 3D point to 2D using basis vectors
static Vec2 project_to_2d(const Vec3& p, const Vec3& origin, const Vec3& u, const Vec3& v) {
    Vec3 d = p - origin;
    return Vec2{dot(d, u), dot(d, v)};
}

// Unproject 2D point back to 3D
static Vec3 unproject_to_3d(const Vec2& p2d, const Vec3& origin, const Vec3& u, const Vec3& v) {
    return origin + u * p2d.x + v * p2d.y;
}

// Get particle cross-section polygon at a plane (defined by wall triangle)
// Returns vertices of the polygon formed by intersecting particle mesh with the plane
// Also returns the maximum penetration depth (negative distance behind plane)
static std::vector<Vec3> get_particle_plane_section(const Mesh& mesh, const Transform& tf,
                                                     const Vec3& plane_n, double plane_d,
                                                     double tol, double& out_penetration) {
    std::vector<Vec3> section_pts;
    std::vector<std::array<Vec3, 3>> tris = transform_tris(mesh, tf);

    out_penetration = 0.0;
    double min_dist = std::numeric_limits<double>::infinity();

    for (const auto& tri : tris) {
        // Check each edge for intersection with plane
        for (int i = 0; i < 3; ++i) {
            int j = (i + 1) % 3;
            Vec3 p0 = tri[i];
            Vec3 p1 = tri[j];

            double d0 = dot(plane_n, p0) + plane_d;
            double d1 = dot(plane_n, p1) + plane_d;

            // Track minimum distance (penetration if negative)
            if (d0 < min_dist) min_dist = d0;
            if (d1 < min_dist) min_dist = d1;

            // Check if edge actually crosses the plane (strict condition)
            if (d0 * d1 < -tol * tol) {
                // Edge crosses the plane - compute intersection point
                double t = d0 / (d0 - d1);
                t = std::min(1.0, std::max(0.0, t));
                Vec3 ip = p0 + (p1 - p0) * t;
                section_pts.push_back(ip);
            }
        }
    }

    // Penetration depth is negative of minimum distance
    out_penetration = -min_dist;

    // Remove duplicate points
    return unique_points(section_pts, tol);
}

// Compute wall-particle contact using Sutherland-Hodgman clipping
// Returns: overlap area, contact point (in 3D), contact normal (wall normal), penetration depth
static bool compute_wall_contact(const Mesh& particle_mesh, const Transform& particle_tf,
                                  const std::array<Vec3, 3>& wall_tri, const Vec3& wall_normal,
                                  double particle_radius, double tol,
                                  double& out_area, Vec3& out_contact_point, Vec3& out_normal,
                                  double& out_penetration) {
    // Compute wall plane
    Vec3 plane_n = wall_normal;
    double plane_d = -dot(plane_n, wall_tri[0]);

    // Get particle cross-section at wall plane and penetration depth
    double penetration;
    std::vector<Vec3> section_3d = get_particle_plane_section(particle_mesh, particle_tf,
                                                               plane_n, plane_d, tol, penetration);

    // Only consider contact if there is actual penetration
    if (penetration <= tol) {
        return false;
    }

    if (section_3d.size() < 3) {
        return false;
    }

    // Build 2D coordinate system on wall plane
    Vec3 basis_u, basis_v;
    build_basis(plane_n, basis_u, basis_v);
    Vec3 origin = wall_tri[0];

    // Project particle section to 2D
    std::vector<Vec2> subject_2d;
    for (const auto& p : section_3d) {
        subject_2d.push_back(project_to_2d(p, origin, basis_u, basis_v));
    }

    // Project wall triangle to 2D
    std::vector<Vec2> clip_2d;
    for (const auto& p : wall_tri) {
        clip_2d.push_back(project_to_2d(p, origin, basis_u, basis_v));
    }

    // Clip particle section by wall triangle
    std::vector<Vec2> clipped = sutherland_hodgman_clip(subject_2d, clip_2d);
    if (clipped.size() < 3) {
        return false;
    }

    // Compute overlap area
    double area = polygon_area_2d(clipped);
    if (area < tol * tol * 0.01) {
        return false;
    }

    // Compute centroid
    Vec2 centroid_2d = polygon_centroid_2d(clipped);
    Vec3 contact_point = unproject_to_3d(centroid_2d, origin, basis_u, basis_v);

    out_area = area;
    out_contact_point = contact_point;
    out_normal = plane_n;
    out_penetration = penetration;
    return true;
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
    return mat3_mul(mat3_mul(R, p.inertia_body), Rt);
}


static Mat3 inertia_world_inv(const Particle& p) {
    Mat3 R = quat_to_mat3(p.tf.rot);
    Mat3 Rt = mat3_transpose(R);
    return mat3_mul(mat3_mul(R, p.inertia_body_inv), Rt);
}

static void integrate_particle(Particle& p, const Vec3& force, const Vec3& torque, const Vec3& gravity, double dt) {
    Vec3 acc = force * p.inv_mass + gravity;
    p.vel += acc * dt;
    p.tf.pos += p.vel * dt;

    Mat3 Iinv = inertia_world_inv(p);
    p.L += torque * dt;
    p.omega = mat3_mul_vec3(Iinv, p.L);

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
    Vec3 gravity{0.0, -9.81, 0.0};
    bool split_contacts = true;
    bool center_mesh = true;
    bool contact_debug = false;
    std::string stl_path;
    std::string vtk_prefix = "particles";
    int output_interval = 1;
    std::string output_dir = "output";
    Vec3 v0{0.0, 0.0, 0.0};
    Vec3 v1{0.0, 0.0, 0.0};
    std::vector<ParticleInit> particle_inits;
    std::vector<WallInit> wall_inits;
};

static bool parse_config_file(const std::string& path, SimConfig& cfg) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Failed to open config: " << path << "\n";
        return false;
    }
    std::string line;
    bool in_particle = false;
    bool in_wall = false;
    ParticleInit current_particle;
    WallInit current_wall;
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
        if (key == "particle") {
            if (in_particle) {
                cfg.particle_inits.push_back(current_particle);
            }
            if (in_wall) {
                cfg.wall_inits.push_back(current_wall);
                in_wall = false;
            }
            current_particle = ParticleInit{};
            in_particle = true;
            continue;
        }
        if (key == "end_particle" || key == "particle_end") {
            if (in_particle) {
                cfg.particle_inits.push_back(current_particle);
                in_particle = false;
            }
            continue;
        }
        if (key == "wall") {
            if (in_particle) {
                cfg.particle_inits.push_back(current_particle);
                in_particle = false;
            }
            if (in_wall) {
                cfg.wall_inits.push_back(current_wall);
            }
            current_wall = WallInit{};
            in_wall = true;
            continue;
        }
        if (key == "end_wall" || key == "wall_end") {
            if (in_wall) {
                cfg.wall_inits.push_back(current_wall);
                in_wall = false;
            }
            continue;
        }

        if (in_particle) {
            if (key == "stl") {
                iss >> current_particle.stl_path;
            } else if (key == "pos") {
                iss >> current_particle.pos.x >> current_particle.pos.y >> current_particle.pos.z;
            } else if (key == "vel") {
                iss >> current_particle.vel.x >> current_particle.vel.y >> current_particle.vel.z;
            } else if (key == "omega") {
                iss >> current_particle.omega.x >> current_particle.omega.y >> current_particle.omega.z;
            } else if (key == "quat") {
                iss >> current_particle.rot.w >> current_particle.rot.x >> current_particle.rot.y >> current_particle.rot.z;
            } else if (key == "scale") {
                iss >> current_particle.scale;
            } else if (key == "density") {
                iss >> current_particle.density;
            } else if (key == "young") {
                iss >> current_particle.young;
            } else if (key == "poisson") {
                iss >> current_particle.poisson;
            } else if (key == "mu") {
                iss >> current_particle.mu;
            } else if (key == "restitution") {
                iss >> current_particle.restitution;
            }
        } else if (in_wall) {
            if (key == "stl") {
                iss >> current_wall.stl_path;
            } else if (key == "pos") {
                iss >> current_wall.pos.x >> current_wall.pos.y >> current_wall.pos.z;
            } else if (key == "quat") {
                iss >> current_wall.rot.w >> current_wall.rot.x >> current_wall.rot.y >> current_wall.rot.z;
            } else if (key == "scale") {
                iss >> current_wall.scale;
            } else if (key == "mu") {
                iss >> current_wall.mu;
            } else if (key == "restitution") {
                iss >> current_wall.restitution;
            }
        } else if (key == "stl") {
            iss >> cfg.stl_path;
        } else if (key == "n") {
            iss >> cfg.n;
        } else if (key == "steps") {
            iss >> cfg.steps;
        } else if (key == "dt") {
            iss >> cfg.dt;
        } else if (key == "spacing") {
            iss >> cfg.spacing;
        } else if (key == "split_contacts") {
            int v = 1;
            iss >> v;
            cfg.split_contacts = (v != 0);
        } else if (key == "contact_debug") {
            int v = 0;
            iss >> v;
            cfg.contact_debug = (v != 0);
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
        } else if (key == "output_dir") {
            iss >> cfg.output_dir;
        } else if (key == "v0") {
            iss >> cfg.v0.x >> cfg.v0.y >> cfg.v0.z;
        } else if (key == "v1") {
            iss >> cfg.v1.x >> cfg.v1.y >> cfg.v1.z;
        }
    }
    if (in_particle) {
        cfg.particle_inits.push_back(current_particle);
    }
    if (in_wall) {
        cfg.wall_inits.push_back(current_wall);
    }
    if (cfg.output_interval < 1) {
        cfg.output_interval = 1;
    }
    if (!cfg.stl_path.empty()) {
        return true;
    }
    return !cfg.particle_inits.empty();
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
                                const std::vector<Mesh>& meshes,
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
    std::size_t total_tris = 0;
    for (const auto& p : particles) {
        total_tris += meshes[p.mesh_index].tris.size();
    }
    const std::size_t total_points = total_tris * 3;
    out << "POINTS " << total_points << " float\n";
    for (const auto& p : particles) {
        const auto& m = meshes[p.mesh_index];
        for (const auto& tri : m.tris) {
            Vec3 v0 = quat_rotate(p.tf.rot, tri[0]) + p.tf.pos;
            Vec3 v1 = quat_rotate(p.tf.rot, tri[1]) + p.tf.pos;
            Vec3 v2 = quat_rotate(p.tf.rot, tri[2]) + p.tf.pos;
            out << static_cast<float>(v0.x) << " " << static_cast<float>(v0.y) << " " << static_cast<float>(v0.z) << "\n";
            out << static_cast<float>(v1.x) << " " << static_cast<float>(v1.y) << " " << static_cast<float>(v1.z) << "\n";
            out << static_cast<float>(v2.x) << " " << static_cast<float>(v2.y) << " " << static_cast<float>(v2.z) << "\n";
        }
    }

    out << "POLYGONS " << total_tris << " " << total_tris * 4 << "\n";
    for (std::size_t i = 0; i < total_tris; ++i) {
        std::size_t base = i * 3;
        out << "3 " << base << " " << base + 1 << " " << base + 2 << "\n";
    }

    out << "CELL_DATA " << total_tris << "\n";

    out << "SCALARS id int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        std::size_t tris_per_particle = meshes[particles[i].mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << i << "\n";
        }
    }

    out << "SCALARS mass float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.mass) << "\n";
        }
    }

    out << "SCALARS radius float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.radius) << "\n";
        }
    }

    out << "SCALARS contact_count int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < contact_counts.size(); ++i) {
        std::size_t tris_per_particle = meshes[particles[i].mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << contact_counts[i] << "\n";
        }
    }

    out << "VECTORS velocity float\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.vel.x) << " "
                << static_cast<float>(p.vel.y) << " "
                << static_cast<float>(p.vel.z) << "\n";
        }
    }

    out << "VECTORS omega float\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.omega.x) << " "
                << static_cast<float>(p.omega.y) << " "
                << static_cast<float>(p.omega.z) << "\n";
        }
    }

    out << "VECTORS force float\n";
    for (std::size_t i = 0; i < particles.size(); ++i) {
        std::size_t tris_per_particle = meshes[particles[i].mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            const Vec3& f = forces[i];
            out << static_cast<float>(f.x) << " "
                << static_cast<float>(f.y) << " "
                << static_cast<float>(f.z) << "\n";
        }
    }

    out << "VECTORS torque float\n";
    for (std::size_t i = 0; i < particles.size(); ++i) {
        std::size_t tris_per_particle = meshes[particles[i].mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            const Vec3& tt = torques[i];
            out << static_cast<float>(tt.x) << " "
                << static_cast<float>(tt.y) << " "
                << static_cast<float>(tt.z) << "\n";
        }
    }

    out << "FIELD FieldData 1\n";
    out << "orientation 4 " << total_tris << " float\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.tf.rot.w) << " "
                << static_cast<float>(p.tf.rot.x) << " "
                << static_cast<float>(p.tf.rot.y) << " "
                << static_cast<float>(p.tf.rot.z) << "\n";
        }
    }
}

static void write_particle_state_txt(const std::string& path,
                                     const std::vector<Particle>& particles,
                                     const std::vector<Mesh>& meshes,
                                     const std::vector<ParticleInit>& inits) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to write particle state: " << path << "\n";
        return;
    }
    out.setf(std::ios::fixed);
    out << std::setprecision(9);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const Particle& p = particles[i];
        const ParticleInit& init = inits[i];
        out << "particle\n";
        out << "stl = " << (init.stl_path.empty() ? "" : init.stl_path) << "\n";
        out << "pos = " << p.tf.pos.x << " " << p.tf.pos.y << " " << p.tf.pos.z << "\n";
        out << "vel = " << p.vel.x << " " << p.vel.y << " " << p.vel.z << "\n";
        out << "quat = " << p.tf.rot.w << " " << p.tf.rot.x << " " << p.tf.rot.y << " " << p.tf.rot.z << "\n";
        out << "omega = " << p.omega.x << " " << p.omega.y << " " << p.omega.z << "\n";
        out << "scale = " << init.scale << "\n";
        out << "density = " << init.density << "\n";
        out << "young = " << init.young << "\n";
        out << "poisson = " << init.poisson << "\n";
        out << "mu = " << init.mu << "\n";
        out << "restitution = " << init.restitution << "\n";
        out << "end_particle\n\n";
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

    std::unordered_map<std::string, std::vector<std::array<Vec3, 3>>> tri_cache;
    std::unordered_map<std::string, int> mesh_cache;
    std::vector<Mesh> meshes;

    auto load_mesh_for = [&](const std::string& path, double scale) -> int {
        std::string key = path + "|" + std::to_string(scale);
        auto it_cache = mesh_cache.find(key);
        if (it_cache != mesh_cache.end()) {
            return it_cache->second;
        }

        auto it = tri_cache.find(path);
        if (it == tri_cache.end()) {
            std::vector<std::array<Vec3, 3>> tris;
            if (!load_stl(path, tris)) {
                return -1;
            }
            tri_cache[path] = tris;
            it = tri_cache.find(path);
        }
        std::vector<std::array<Vec3, 3>> tris = it->second;
        for (auto& tri : tris) {
            tri[0] *= scale;
            tri[1] *= scale;
            tri[2] *= scale;
        }
        Mesh m = build_mesh(tris, cfg.center_mesh);
        meshes.push_back(m);
        int idx = static_cast<int>(meshes.size() - 1);
        mesh_cache[key] = idx;
        return idx;
    };

    std::vector<ParticleInit> inits;
    if (!cfg.particle_inits.empty()) {
        inits = cfg.particle_inits;
    } else {
        inits.resize(cfg.n);
        for (int i = 0; i < cfg.n; ++i) {
            inits[i].stl_path = cfg.stl_path;
            inits[i].pos = Vec3{cfg.spacing * i, 0.0, 0.0};
            inits[i].vel = (i == 0) ? cfg.v0 : (i == 1 ? cfg.v1 : Vec3{0.0, 0.0, 0.0});
            inits[i].rot = Quat{};
            inits[i].omega = Vec3{0.0, 0.0, 0.0};
            inits[i].scale = 1.0;
            inits[i].density = 1.0;
            inits[i].young = 1e7;
            inits[i].poisson = 0.25;
            inits[i].mu = 0.5;
            inits[i].restitution = 0.5;
        }
    }

    std::vector<Particle> particles(inits.size());
    for (std::size_t i = 0; i < inits.size(); ++i) {
        const auto& init = inits[i];
        std::string stl_path = init.stl_path.empty() ? cfg.stl_path : init.stl_path;
        if (stl_path.empty()) {
            std::cerr << "Missing stl for particle " << i << "\n";
            return 1;
        }
        double scale = init.scale > 0.0 ? init.scale : 1.0;
        int mesh_index = load_mesh_for(stl_path, scale);
        if (mesh_index < 0) {
            std::cerr << "Failed to load STL: " << stl_path << "\n";
            return 1;
        }
        const Mesh& mesh = meshes[mesh_index];
        if (mesh.vertices.empty()) {
            std::cerr << "Mesh has no vertices.\n";
            return 1;
        }
        Particle p;
        p.tf.pos = init.pos;
        p.tf.rot = init.rot;
        p.vel = init.vel;
        p.omega = init.omega;
        if (mesh.volume > 0.0) {
            double density = init.density > 0.0 ? init.density : 1.0;
            p.mass = density * mesh.volume;
            p.inv_mass = 1.0 / p.mass;
            p.inertia_body = mat3_scale(mesh.inertia_unit, density);
            p.inertia_body_inv = mat3_inverse(p.inertia_body);
        } else {
            double density = init.density > 0.0 ? init.density : 1.0;
            p.mass = density;
            p.inv_mass = 1.0 / p.mass;
            double I = 0.4 * p.mass * mesh.radius * mesh.radius;
            Mat3 Ibody = mat3_identity();
            Ibody.m[0][0] = I;
            Ibody.m[1][1] = I;
            Ibody.m[2][2] = I;
            p.inertia_body = Ibody;
            p.inertia_body_inv = mat3_inverse(Ibody);
        }
        p.radius = mesh.radius;
        p.equiv_radius = mesh.volume > 0.0 ? std::cbrt(3.0 * mesh.volume / (4.0 * 3.141592653589793)) : mesh.radius;
        p.young = init.young;
        p.poisson = init.poisson;
        p.mu = init.mu;
        p.restitution = init.restitution;
        p.mesh_index = mesh_index;
        p.L = mat3_mul_vec3(inertia_world(p), p.omega);
        particles[i] = p;
    }

    // Initialize walls
    std::vector<Wall> walls;
    for (const auto& wi : cfg.wall_inits) {
        if (wi.stl_path.empty()) {
            std::cerr << "Wall missing stl path.\n";
            continue;
        }
        std::vector<std::array<Vec3, 3>> wall_tris;
        if (!load_stl(wi.stl_path, wall_tris)) {
            std::cerr << "Failed to load wall STL: " << wi.stl_path << "\n";
            continue;
        }
        double scale = wi.scale > 0.0 ? wi.scale : 1.0;
        for (auto& tri : wall_tris) {
            tri[0] = tri[0] * scale;
            tri[1] = tri[1] * scale;
            tri[2] = tri[2] * scale;
        }
        Wall wall;
        wall.mesh = build_mesh(wall_tris, cfg.center_mesh);
        wall.tf.pos = wi.pos;
        wall.tf.rot = wi.rot;
        wall.mu = wi.mu;
        wall.restitution = wi.restitution;

        // Precompute wall triangle normals in world space
        std::vector<std::array<Vec3, 3>> world_tris = transform_tris(wall.mesh, wall.tf);
        wall.tri_normals.resize(world_tris.size());
        for (std::size_t t = 0; t < world_tris.size(); ++t) {
            wall.tri_normals[t] = tri_normal(world_tris[t]);
        }
        walls.push_back(wall);
        std::cout << "Loaded wall: " << wi.stl_path << " (tris=" << wall.mesh.tris.size() << ")\n";
    }

    // Tangential displacement for wall contacts: key = (particle_idx << 16) | wall_idx
    std::unordered_map<long long, Vec3> wall_tangential_disp;

    std::vector<Vec3> forces(particles.size());
    std::vector<Vec3> torques(particles.size());
    std::vector<int> contact_counts(particles.size());
    std::unordered_map<long long, Vec3> tangential_disp;

    std::filesystem::create_directories(cfg.output_dir);

    for (int step = 0; step < cfg.steps; ++step) {
        auto step_t0 = std::chrono::steady_clock::now();
        std::fill(forces.begin(), forces.end(), Vec3{});
        std::fill(torques.begin(), torques.end(), Vec3{});
        std::fill(contact_counts.begin(), contact_counts.end(), 0);

        int contacts = 0;
        std::unordered_set<long long> active_pairs;
        std::size_t total_segments = 0;
        std::size_t total_loops = 0;
        std::size_t total_comps = 0;
        std::size_t tri_checks = 0;
        auto t_tri = std::chrono::steady_clock::duration::zero();
        auto t_rebuild = std::chrono::steady_clock::duration::zero();
        auto t_output = std::chrono::steady_clock::duration::zero();
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(particles.size()); ++j) {
                Particle& pa = particles[i];
                Particle& pb = particles[j];
                long long pair_key = (static_cast<long long>(i) << 32) | static_cast<unsigned long long>(j);
                Vec3 dpos = pb.tf.pos - pa.tf.pos;
                double dist2 = dot(dpos, dpos);
                double rsum = pa.radius + pb.radius;
                if (dist2 > rsum * rsum) {
                    continue;
                }

                auto t0 = std::chrono::steady_clock::now();
                const Mesh& meshA = meshes[pa.mesh_index];
                const Mesh& meshB = meshes[pb.mesh_index];
                std::vector<std::array<Vec3, 3>> trisA = transform_tris(meshA, pa.tf);
                std::vector<std::array<Vec3, 3>> trisB = transform_tris(meshB, pb.tf);

                std::vector<Vec3> nA_all(trisA.size());
                std::vector<Vec3> nB_all(trisB.size());
                std::vector<std::pair<Vec3, double>> planesA(trisA.size());
                std::vector<std::pair<Vec3, double>> planesB(trisB.size());

                for (std::size_t ta = 0; ta < trisA.size(); ++ta) {
                    Vec3 n; double d;
                    if (!plane_from_tri(trisA[ta], n, d)) {
                        nA_all[ta] = Vec3{0.0, 0.0, 0.0};
                        planesA[ta] = {Vec3{0.0, 0.0, 0.0}, 0.0};
                    } else {
                        nA_all[ta] = n;
                        planesA[ta] = {n, d};
                    }
                }
                for (std::size_t tb = 0; tb < trisB.size(); ++tb) {
                    Vec3 n; double d;
                    if (!plane_from_tri(trisB[tb], n, d)) {
                        nB_all[tb] = Vec3{0.0, 0.0, 0.0};
                        planesB[tb] = {Vec3{0.0, 0.0, 0.0}, 0.0};
                    } else {
                        nB_all[tb] = n;
                        planesB[tb] = {n, d};
                    }
                }

                std::vector<Vec3> aabbA_min(trisA.size());
                std::vector<Vec3> aabbA_max(trisA.size());
                std::vector<Vec3> aabbB_min(trisB.size());
                std::vector<Vec3> aabbB_max(trisB.size());
                for (std::size_t ta = 0; ta < trisA.size(); ++ta) {
                    Vec3 mn{std::min({trisA[ta][0].x, trisA[ta][1].x, trisA[ta][2].x}),
                            std::min({trisA[ta][0].y, trisA[ta][1].y, trisA[ta][2].y}),
                            std::min({trisA[ta][0].z, trisA[ta][1].z, trisA[ta][2].z})};
                    Vec3 mx{std::max({trisA[ta][0].x, trisA[ta][1].x, trisA[ta][2].x}),
                            std::max({trisA[ta][0].y, trisA[ta][1].y, trisA[ta][2].y}),
                            std::max({trisA[ta][0].z, trisA[ta][1].z, trisA[ta][2].z})};
                    aabbA_min[ta] = mn;
                    aabbA_max[ta] = mx;
                }
                for (std::size_t tb = 0; tb < trisB.size(); ++tb) {
                    Vec3 mn{std::min({trisB[tb][0].x, trisB[tb][1].x, trisB[tb][2].x}),
                            std::min({trisB[tb][0].y, trisB[tb][1].y, trisB[tb][2].y}),
                            std::min({trisB[tb][0].z, trisB[tb][1].z, trisB[tb][2].z})};
                    Vec3 mx{std::max({trisB[tb][0].x, trisB[tb][1].x, trisB[tb][2].x}),
                            std::max({trisB[tb][0].y, trisB[tb][1].y, trisB[tb][2].y}),
                            std::max({trisB[tb][0].z, trisB[tb][1].z, trisB[tb][2].z})};
                    aabbB_min[tb] = mn;
                    aabbB_max[tb] = mx;
                }

                double tol = meshA.mean_edge > 0.0 ? (meshA.mean_edge * 0.1) : (meshA.bbox_diag * 1e-2);
                if (tol <= 0.0) {
                    tol = 1e-6;
                }

                std::vector<std::pair<Vec3, Vec3>> segments;
                for (std::size_t ta = 0; ta < trisA.size(); ++ta) {
                    if (norm2(nA_all[ta]) < 1e-30) {
                        continue;
                    }
                    for (std::size_t tb = 0; tb < trisB.size(); ++tb) {
                        tri_checks++;
                        if (aabbA_max[ta].x < aabbB_min[tb].x || aabbA_min[ta].x > aabbB_max[tb].x ||
                            aabbA_max[ta].y < aabbB_min[tb].y || aabbA_min[ta].y > aabbB_max[tb].y ||
                            aabbA_max[ta].z < aabbB_min[tb].z || aabbA_min[ta].z > aabbB_max[tb].z) {
                            continue;
                        }
                        if (norm2(nB_all[tb]) < 1e-30) {
                            continue;
                        }
                        Vec3 p0, p1;
                        if (tri_tri_intersection_segment(trisA[ta], planesA[ta].first, planesA[ta].second,
                                                         trisB[tb], planesB[tb].first, planesB[tb].second,
                                                         p0, p1)) {
                            segments.emplace_back(p0, p1);
                        }
                    }
                }
                t_tri += std::chrono::steady_clock::now() - t0;

                if (segments.empty()) {
                    continue;
                }

                t0 = std::chrono::steady_clock::now();
                std::vector<std::vector<int>> comps;
                if (cfg.split_contacts) {
                    label_segments_by_contact_greedy(segments, tol, comps);
                } else {
                    comps.resize(1);
                    comps[0].reserve(segments.size());
                    for (int s = 0; s < static_cast<int>(segments.size()); ++s) {
                        comps[0].push_back(s);
                    }
                }
                total_segments += segments.size();
                total_comps += comps.size();

                for (const auto& comp : comps) {
                    std::vector<std::pair<Vec3, Vec3>> segs;
                    segs.reserve(comp.size());
                    for (int idx : comp) {
                        segs.push_back(segments[idx]);
                    }

                    std::vector<Vec3> V;
                    std::vector<std::pair<int, int>> edges;
                    snap_endpoints(segs, tol, V, edges);

                    std::vector<Vec3> V2;
                    std::vector<std::pair<int, int>> edges2;
                    split_t_junctions(V, edges, tol, V2, edges2);

                    std::vector<std::vector<int>> loops_vids = extract_loops(V2, edges2);
                    if (loops_vids.empty()) {
                        continue;
                    }
                    total_loops += loops_vids.size();

                    for (const auto& lv : loops_vids) {
                        if (lv.size() < 3) {
                            continue;
                        }
                        std::vector<Vec3> loop_pts;
                        loop_pts.reserve(lv.size());
                        for (int vid : lv) {
                            loop_pts.push_back(V2[vid]);
                        }
                        Vec3 Sn, Gn;
                        auto acc = accumulate_Sn_Gn_from_polyline(loop_pts);
                        Sn = acc.first;
                        Gn = acc.second;
                        Vec3 xc0, nA;
                        double area = 0.0;
                        if (!contact_point_xc0(Sn, Gn, xc0, nA, area)) {
                            continue;
                        }
                        double area_eps = meshA.mean_edge > 0.0 ? (meshA.mean_edge * meshA.mean_edge * 1e-4) : 1e-12;
                        if (area < area_eps) {
                            continue;
                        }
                        Vec3 dirAB = pb.tf.pos - pa.tf.pos;
                        if (dot(nA, dirAB) < 0.0) {
                            // std::reverse(loop_pts.begin(), loop_pts.end());
                            // acc = accumulate_Sn_Gn_from_polyline(loop_pts);
                            // Sn = acc.first;
                            // Gn = acc.second;
                            // if (!contact_point_xc0(Sn, Gn, xc0, nA, area)) {
                            //     continue;
                            // }
                            nA = nA * -1.0;
                        }
                        Vec3 rA = xc0 - pa.tf.pos;
                        Vec3 rB = xc0 - pb.tf.pos;
                        Vec3 vA = pa.vel + cross(pa.omega, rA);
                        Vec3 vB = pb.vel + cross(pb.omega, rB);
                        Vec3 vrel = vB - vA;
                        double vn = dot(vrel, nA);

                        double R1 = std::max(pa.equiv_radius, 1e-12);
                        double R2 = std::max(pb.equiv_radius, 1e-12);
                        double E1 = std::max(pa.young, 1e-12);
                        double E2 = std::max(pb.young, 1e-12);
                        double nu1 = pa.poisson;
                        double nu2 = pb.poisson;
                        double A_n = area;

                        // Normal stiffness (Feng 2021 energy-conserving model)
                        // kn = (E1/R1 + E2/R2) * An
                        double kn = (E1 / R1 + E2 / R2) * A_n;
                        // Tangential stiffness
                        // kt = (E1/(2*(1+nu1)*R1) + E2/(2*(1+nu2)*R2)) * An
                        double kt = (E1 / (2.0 * (1.0 + nu1) * R1) + E2 / (2.0 * (1.0 + nu2) * R2)) * A_n;
                        double m_eff = 1.0 / (pa.inv_mass + pb.inv_mass);
                        double e = std::min(pa.restitution, pb.restitution);
                        e = std::min(0.9999, std::max(1e-6, e));
                        double loge = std::log(e);
                        double cn = 2.0 * loge * std::sqrt(kn * m_eff) / std::sqrt(loge * loge + 3.141592653589793 * 3.141592653589793);
                        double ct = 2.0 * loge * std::sqrt(kt * m_eff) / std::sqrt(loge * loge + 3.141592653589793 * 3.141592653589793);
                        cn = std::abs(cn);
                        ct = std::abs(ct);

                        // Estimate overlap volume from contact loop geometry
                        // Characteristic length from contact area: Lc = sqrt(Area)
                        // Overlap volume estimate: deltaV ≈ Area^(3/2) (assuming quasi-circular contact)
                        double deltaV = std::pow(A_n, 1.5);

                        // Normal force: Fn = kn * deltaV^(1/3)
                        // deltaV^(1/3) ≈ Area^(1/2) acts as equivalent penetration depth
                        double fn = kn * std::cbrt(deltaV);
                         /*if (vn < 0.0) {
                            fn += -cn * vn;
                        }*/
                        if (fn < 0.0) fn = 0.0;
                        Vec3 fn_vec = nA * fn;

                        Vec3 vt = vrel - nA * vn;
                        Vec3 ds = tangential_disp[pair_key];
                        ds += vt * cfg.dt;
                        ds -= nA * dot(ds, nA);
                        Vec3 ft = ds * (-kt) + vt * (-ct);
                        double mu = std::min(pa.mu, pb.mu);
                        double ft_norm = norm(ft);
                        double ft_max = mu * fn;
                        if (ft_norm > ft_max && ft_norm > 1e-14) {
                            ft = ft * (ft_max / ft_norm);
                            if (kt > 1e-12) {
                                ds = ft * (-1.0 / kt);
                            }
                        }
                        tangential_disp[pair_key] = ds;

                        Vec3 f = fn_vec + ft;
                        forces[i] -= f;
                        forces[j] += f;
                        torques[i] += cross(rA, f * -1.0);
                        torques[j] += cross(rB, f);
                        active_pairs.insert(pair_key);
                        contact_counts[i] += 1;
                        contact_counts[j] += 1;
                        if (cfg.contact_debug && contacts == 0) {
                            Vec3 tauA = cross(rA, f);
                            Vec3 tauB = cross(rB, f * -1.0);
                            std::cout << "  contact_debug: i=" << i << " j=" << j
                                      << " xc0=" << xc0.x << "," << xc0.y << "," << xc0.z
                                      << " n=" << nA.x << "," << nA.y << "," << nA.z
                                      << " F=" << f.x << "," << f.y << "," << f.z
                                      << " rA=" << rA.x << "," << rA.y << "," << rA.z
                                      << " tauA=" << tauA.x << "," << tauA.y << "," << tauA.z
                                      << " rB=" << rB.x << "," << rB.y << "," << rB.z
                                      << " tauB=" << tauB.x << "," << tauB.y << "," << tauB.z
                                      << "\n";
                        }
                        contacts++;
                    }
                }
                t_rebuild += std::chrono::steady_clock::now() - t0;
            }
        }

        // ============== Particle-Wall Contact Detection ==============
        std::unordered_set<long long> active_wall_pairs;
        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            Particle& p = particles[i];
            const Mesh& pmesh = meshes[p.mesh_index];

            for (int w = 0; w < static_cast<int>(walls.size()); ++w) {
                const Wall& wall = walls[w];

                // Broad phase: check if particle bounding sphere overlaps with wall bounding sphere
                Vec3 dpos = p.tf.pos - wall.tf.pos;
                double dist2 = dot(dpos, dpos);
                double rsum = p.radius + wall.mesh.radius;
                if (dist2 > rsum * rsum) {
                    continue;
                }

                // Get transformed particle triangles
                std::vector<std::array<Vec3, 3>> particle_tris = transform_tris(pmesh, p.tf);
                std::vector<std::array<Vec3, 3>> wall_tris = transform_tris(wall.mesh, wall.tf);

                double tol = pmesh.mean_edge > 0.0 ? (pmesh.mean_edge * 0.1) : (pmesh.bbox_diag * 1e-2);
                if (tol <= 0.0) tol = 1e-6;

                // For each wall triangle, check for contact
                for (std::size_t t = 0; t < wall_tris.size(); ++t) {
                    const auto& wtri = wall_tris[t];
                    Vec3 wn = wall.tri_normals[t];

                    // Quick AABB check between particle and wall triangle
                    Vec3 wmn{std::min({wtri[0].x, wtri[1].x, wtri[2].x}),
                             std::min({wtri[0].y, wtri[1].y, wtri[2].y}),
                             std::min({wtri[0].z, wtri[1].z, wtri[2].z})};
                    Vec3 wmx{std::max({wtri[0].x, wtri[1].x, wtri[2].x}),
                             std::max({wtri[0].y, wtri[1].y, wtri[2].y}),
                             std::max({wtri[0].z, wtri[1].z, wtri[2].z})};

                    // Particle AABB
                    Vec3 pmn = p.tf.pos - Vec3{p.radius, p.radius, p.radius};
                    Vec3 pmx = p.tf.pos + Vec3{p.radius, p.radius, p.radius};

                    if (pmx.x < wmn.x || pmn.x > wmx.x ||
                        pmx.y < wmn.y || pmn.y > wmx.y ||
                        pmx.z < wmn.z || pmn.z > wmx.z) {
                        continue;
                    }

                    // Compute wall contact using Sutherland-Hodgman clipping
                    double area, penetration;
                    Vec3 contact_point, contact_normal;
                    if (!compute_wall_contact(pmesh, p.tf, wtri, wn, p.radius, tol,
                                              area, contact_point, contact_normal, penetration)) {
                        continue;
                    }

                    // Ensure contact normal points from wall to particle (outward from wall)
                    Vec3 to_particle = p.tf.pos - contact_point;
                    if (dot(contact_normal, to_particle) < 0.0) {
                        contact_normal = contact_normal * -1.0;
                    }

                    // Compute contact force
                    Vec3 rP = contact_point - p.tf.pos;
                    Vec3 vP = p.vel + cross(p.omega, rP);
                    Vec3 vrel = vP;  // Wall velocity is zero
                    double vn = dot(vrel, contact_normal);

                    // Wall contact stiffness (wall is rigid, so only particle properties matter)
                    double Rp = std::max(p.equiv_radius, 1e-12);
                    double Ep = std::max(p.young, 1e-12);
                    double nu_p = p.poisson;
                    double A_n = area;

                    // kn = (Ep/Rp) * An for wall contact
                    double kn = (Ep / Rp) * A_n;
                    // kt = (Ep/(2*(1+nu_p)*Rp)) * An
                    double kt = (Ep / (2.0 * (1.0 + nu_p) * Rp)) * A_n;

                    // Effective mass for wall contact (wall mass is infinite)
                    double m_eff = p.mass;

                    double e = p.restitution;
                    e = std::min(0.9999, std::max(1e-6, e));
                    double loge = std::log(e);
                    double cn = 2.0 * loge * std::sqrt(kn * m_eff) / std::sqrt(loge * loge + 3.141592653589793 * 3.141592653589793);
                    double ct = 2.0 * loge * std::sqrt(kt * m_eff) / std::sqrt(loge * loge + 3.141592653589793 * 3.141592653589793);
                    cn = std::abs(cn);
                    ct = std::abs(ct);

                    // Normal force using actual penetration depth
                    // fn = kn * penetration (Hookean spring model)
                    double fn = kn * penetration;
                    // Add damping force when approaching
                    if (vn < 0.0) {
                        fn += -cn * vn;
                    }
                    if (fn < 0.0) fn = 0.0;
                    Vec3 fn_vec = contact_normal * fn;

                    // Tangential force
                    Vec3 vt = vrel - contact_normal * vn;
                    long long wall_pair_key = (static_cast<long long>(i) << 16) | static_cast<long long>(w);
                    Vec3 ds = wall_tangential_disp[wall_pair_key];
                    ds += vt * cfg.dt;
                    ds -= contact_normal * dot(ds, contact_normal);
                    Vec3 ft = ds * (-kt) + vt * (-ct);

                    double mu = std::min(p.mu, wall.mu);
                    double ft_norm = norm(ft);
                    double ft_max = mu * fn;
                    if (ft_norm > ft_max && ft_norm > 1e-14) {
                        ft = ft * (ft_max / ft_norm);
                        if (kt > 1e-12) {
                            ds = ft * (-1.0 / kt);
                        }
                    }
                    wall_tangential_disp[wall_pair_key] = ds;

                    // Apply force to particle (wall doesn't move)
                    Vec3 f = fn_vec + ft;
                    forces[i] += f;
                    torques[i] += cross(rP, f);

                    active_wall_pairs.insert(wall_pair_key);
                    contact_counts[i] += 1;
                    contacts++;

                    if (cfg.contact_debug && contacts < 3) {
                        std::cout << "  wall_contact: particle=" << i << " wall=" << w
                                  << " tri=" << t
                                  << " area=" << area
                                  << " penetration=" << penetration
                                  << " xc=" << contact_point.x << "," << contact_point.y << "," << contact_point.z
                                  << " n=" << contact_normal.x << "," << contact_normal.y << "," << contact_normal.z
                                  << " F=" << f.x << "," << f.y << "," << f.z
                                  << "\n";
                    }
                }
            }
        }

        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            integrate_particle(particles[i], forces[i], torques[i], cfg.gravity, cfg.dt);
        }

        for (auto it = tangential_disp.begin(); it != tangential_disp.end(); ) {
            if (active_pairs.find(it->first) == active_pairs.end()) {
                it = tangential_disp.erase(it);
            } else {
                ++it;
            }
        }

        // Clean up inactive wall contact tangential displacements
        for (auto it = wall_tangential_disp.begin(); it != wall_tangential_disp.end(); ) {
            if (active_wall_pairs.find(it->first) == active_wall_pairs.end()) {
                it = wall_tangential_disp.erase(it);
            } else {
                ++it;
            }
        }

        if (step % cfg.output_interval == 0) {
            std::ostringstream oss;
            oss << cfg.vtk_prefix << "_" << std::setw(6) << std::setfill('0') << step << ".vtk";
            std::filesystem::path out_path = std::filesystem::path(cfg.output_dir) / oss.str();
            auto t_out0 = std::chrono::steady_clock::now();
            write_vtk_particles(out_path.string(), meshes, particles, forces, torques, contact_counts);

            // Output particle state to txt
            std::ostringstream oss_txt;
            oss_txt << cfg.vtk_prefix << "_" << std::setw(6) << std::setfill('0') << step << ".txt";
            std::filesystem::path txt_path = std::filesystem::path(cfg.output_dir) / oss_txt.str();
            write_particle_state_txt(txt_path.string(), particles, meshes, inits);

            t_output += std::chrono::steady_clock::now() - t_out0;

            auto step_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - step_t0).count();
            auto tri_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_tri).count();
            auto rebuild_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_rebuild).count();
            auto out_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_output).count();
            std::cout << "step=" << step
                  << " contacts=" << contacts
                  << " step_ms=" << step_ms
                  << " tri_ms=" << tri_ms
                  << " rebuild_ms=" << rebuild_ms
                  << " output_ms=" << out_ms
                  << " tri_checks=" << tri_checks
                  << " segments=" << total_segments
                  << " comps=" << total_comps
                  << " loops=" << total_loops
                  << "\n";
        }
    }

    std::cout << "done\n";
    return 0;
}
