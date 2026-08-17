#include "geometry/mesh_build.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace zdem {

std::vector<Vec3> unique_vertices(const std::vector<Triangle>& tris) {
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

void compute_mass_properties(const std::vector<Triangle>& tris,
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

Mesh build_mesh(const std::vector<Triangle>& tris, bool center_mesh) {
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

std::vector<Triangle> transform_tris(const Mesh& mesh, const Transform& tf) {
    std::vector<Triangle> out;
    out.reserve(mesh.tris.size());
    for (const auto& tri : mesh.tris) {
        Triangle t;
        t[0] = quat_rotate(tf.rot, tri[0]) + tf.pos;
        t[1] = quat_rotate(tf.rot, tri[1]) + tf.pos;
        t[2] = quat_rotate(tf.rot, tri[2]) + tf.pos;
        out.push_back(t);
    }
    return out;
}

Vec3 tri_normal(const Triangle& tri) {
    Vec3 e1 = tri[1] - tri[0];
    Vec3 e2 = tri[2] - tri[0];
    return normalize(cross(e1, e2));
}

bool point_inside_mesh(const std::vector<Triangle>& tris, const Vec3& p) {
    // Ray parity: shoot +x and count surface crossings. A watertight mesh
    // gives odd = inside. Degenerate hits (ray ~parallel to a triangle) are
    // counted conservatively; the STLs here are not adversarial.
    int crossings = 0;
    for (const auto& tri : tris) {
        const Vec3& a = tri[0];
        const Vec3& b = tri[1];
        const Vec3& c = tri[2];
        // Bounding-box shortcut on y/z before the full test.
        double ymin = std::min({a.y, b.y, c.y});
        double ymax = std::max({a.y, b.y, c.y});
        double zmin = std::min({a.z, b.z, c.z});
        double zmax = std::max({a.z, b.z, c.z});
        if (p.y < ymin || p.y > ymax || p.z < zmin || p.z > zmax) continue;
        // Solve p + t*(1,0,0) = a + u*(b-a) + v*(c-a) for the y,z rows:
        //   e1.y*u + e2.y*v = d.y ;  e1.z*u + e2.z*v = d.z   (d = p - a)
        Vec3 e1 = b - a;
        Vec3 e2 = c - a;
        double det = e1.y * e2.z - e1.z * e2.y;
        if (std::abs(det) < 1e-14) continue;  // ray parallel to triangle plane
        Vec3 d = p - a;
        double u = (d.y * e2.z - d.z * e2.y) / det;
        double v = (e1.y * d.z - e1.z * d.y) / det;
        if (u < -1e-12 || v < -1e-12 || u + v > 1.0 + 1e-12) continue;
        // Intersection x coordinate: x = a.x + u*e1.x + v*e2.x; need t = x - p.x > 0.
        double xhit = a.x + u * e1.x + v * e2.x;
        if (xhit > p.x) crossings++;
    }
    return (crossings % 2) == 1;
}

double point_mesh_distance(const std::vector<Triangle>& tris, const Vec3& p,
                           Vec3& out_closest_point, Vec3& out_closest_normal) {
    double best = std::numeric_limits<double>::max();
    Vec3 best_q{0.0, 0.0, 0.0};
    Vec3 best_n{0.0, 0.0, 0.0};
    for (const auto& tri : tris) {
        // Closest point on triangle (Ericson, Real-Time Collision Detection 5.1.5).
        Vec3 ab = tri[1] - tri[0];
        Vec3 ac = tri[2] - tri[0];
        Vec3 ap = p - tri[0];
        double d1 = dot(ab, ap);
        double d2 = dot(ac, ap);
        Vec3 q;
        if (d1 <= 0.0 && d2 <= 0.0) {
            q = tri[0];                                        // vertex A region
        } else {
            Vec3 bp = p - tri[1];
            double d3 = dot(ab, bp);
            double d4 = dot(ac, bp);
            if (d3 >= 0.0 && d4 <= d3) {
                q = tri[1];                                    // vertex B region
            } else {
                Vec3 cp = p - tri[2];
                double d5 = dot(ab, cp);
                double d6 = dot(ac, cp);
                if (d6 >= 0.0 && d5 <= d6) {
                    q = tri[2];                                // vertex C region
                } else if (d1 * d4 - d3 * d2 <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
                    q = tri[0] + ab * (d1 / (d1 - d3));        // edge AB
                } else if (d5 * d2 - d6 * d1 <= 0.0 && d6 >= 0.0 && d5 <= d3) {
                    q = tri[2] + ac * ((d6 - d5) / (d3 - d5)); // edge AC
                } else {
                    double va = d3 * d6 - d5 * d4;
                    double vb = d5 * d2 - d1 * d6;
                    double vc = d1 * d4 - d3 * d2;
                    double denom = 1.0 / (va + vb + vc);
                    double v = vb * denom;
                    double w = vc * denom;
                    q = tri[0] + ab * v + ac * w;              // face region
                }
            }
        }
        double dist2 = norm2(p - q);
        if (dist2 < best) {
            best = dist2;
            best_q = q;
            best_n = tri_normal(tri);
        }
    }
    out_closest_point = best_q;
    out_closest_normal = best_n;
    return (best == std::numeric_limits<double>::max()) ? 0.0 : std::sqrt(best);
}

}  // namespace zdem
