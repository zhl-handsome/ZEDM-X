#include "physics/wall_contact.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "geometry/mesh_build.hpp"
#include "core/quat.hpp"

using namespace zdem;

// Verbatim copy of the CPU point_in_tri (src/main.cpp): the wall footprint
// test. Duplicated here because main.cpp's copy is still used by the
// tri-tri intersection pipeline that remains there.
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

int wall_contact_particle(const Particle& p, const Mesh& pmesh,
                          const std::vector<Wall>& walls,
                          double tangential_damping,
                          Vec3& f, Vec3& t,
                          bool contact_debug,
                          int dbg_particle,
                          int dbg_contacts) {
    int n_contact_groups = 0;
    for (int w = 0; w < static_cast<int>(walls.size()); ++w) {
        const Wall& wall = walls[w];

        // Broad phase: check if particle bounding sphere overlaps with wall bounding sphere
        Vec3 dpos = p.tf.pos - wall.tf.pos;
        double dist2 = dot(dpos, dpos);
        double rsum = p.radius + wall.mesh.radius;
        if (dist2 > rsum * rsum) {
            continue;
        }

        std::vector<std::array<Vec3, 3>> wall_tris = transform_tris(wall.mesh, wall.tf);

        double tol = pmesh.mean_edge > 0.0 ? (pmesh.mean_edge * 0.1) : (pmesh.bbox_diag * 1e-2);
        if (tol <= 0.0) tol = 1e-6;

        // ============== Wall contact: per-vertex penalty (v8) ==============
        // Same contact element as particle-particle: every vertex that
        // crossed a wall plane is one penalty spring acting exactly at
        // that vertex. Coplanar wall triangles passing the AABB test
        // are grouped so a vertex in the shared footprint is charged
        // once (a split plane wall is TWO coplanar triangles -- without
        // grouping every contact force would be doubled). The old
        // clip-loop pipeline drove the force by the deepest vertex but
        // applied it at a depth-weighted average point -- different
        // material points, the spin-mismatch injection fixed in the
        // pair contact; after the integrator stopped dissipating
        // numerically it left the e=0.3 wall case oscillating forever.
        struct WallPlaneGroup {
            Vec3 n; double d;
            std::vector<std::array<Vec3, 3>> footprints;
        };
        std::vector<WallPlaneGroup> groups;
        for (std::size_t t_idx = 0; t_idx < wall_tris.size(); ++t_idx) {
            const auto& wtri = wall_tris[t_idx];
            Vec3 wmn{std::min({wtri[0].x, wtri[1].x, wtri[2].x}),
                     std::min({wtri[0].y, wtri[1].y, wtri[2].y}),
                     std::min({wtri[0].z, wtri[1].z, wtri[2].z})};
            Vec3 wmx{std::max({wtri[0].x, wtri[1].x, wtri[2].x}),
                     std::max({wtri[0].y, wtri[1].y, wtri[2].y}),
                     std::max({wtri[0].z, wtri[1].z, wtri[2].z})};
            Vec3 pmn = p.tf.pos - Vec3{p.radius, p.radius, p.radius};
            Vec3 pmx = p.tf.pos + Vec3{p.radius, p.radius, p.radius};
            if (pmx.x < wmn.x || pmn.x > wmx.x ||
                pmx.y < wmn.y || pmn.y > wmx.y ||
                pmx.z < wmn.z || pmn.z > wmx.z) {
                continue;
            }
            Vec3 pn = tri_normal(wtri);
            if (dot(pn, p.tf.pos - wtri[0]) < 0.0) {
                pn = pn * -1.0;  // push toward the particle side
            }
            double pd = -dot(pn, wtri[0]);
            const double pq = 1e6;
            bool merged = false;
            for (auto& g : groups) {
                if (std::llround(g.n.x * pq) == std::llround(pn.x * pq) &&
                    std::llround(g.n.y * pq) == std::llround(pn.y * pq) &&
                    std::llround(g.n.z * pq) == std::llround(pn.z * pq) &&
                    std::llround(g.d * pq) == std::llround(pd * pq)) {
                    g.footprints.push_back(wtri);
                    merged = true;
                    break;
                }
            }
            if (!merged) {
                WallPlaneGroup g;
                g.n = pn; g.d = pd;
                g.footprints.push_back(wtri);
                groups.push_back(g);
            }
        }

        double Rp = std::max(p.equiv_radius, 1e-12);
        double Ep = std::max(p.young, 1e-12);
        double nu_p = p.poisson;
        double Estar = Ep / (2.0 * (1.0 - nu_p * nu_p));
        double k_hw = (4.0 / 3.0) * Estar * std::sqrt(Rp);
        double kt_w = (Ep * Rp) / (2.0 * (1.0 + nu_p));
        double e_w = std::min(0.9999, std::max(1e-6, p.restitution));
        double loge_w = std::log(e_w);
        double pi2 = 3.141592653589793 * 3.141592653589793;
        double mu_w = std::min(p.mu, wall.mu);

        for (const auto& g : groups) {
            // Count crossings first so the damping uses each vertex share
            // of the particle mass.
            int n_cross = 0;
            for (const auto& v : pmesh.vertices) {
                Vec3 wv = quat_rotate(p.tf.rot, v) + p.tf.pos;
                double s = dot(g.n, wv) + g.d;
                if (s > -1e-12) continue;
                Vec3 p_proj = wv - g.n * s;  // drop onto the wall plane
                for (const auto& fp : g.footprints) {
                    if (point_in_tri(p_proj, fp, g.n)) { ++n_cross; break; }
                }
            }
            if (n_cross == 0) continue;
            double m_eff_v = p.mass / static_cast<double>(n_cross);
            double ct_v = tangential_damping * std::sqrt(kt_w * m_eff_v);
            bool any_force = false;
            for (const auto& v : pmesh.vertices) {
                Vec3 wv = quat_rotate(p.tf.rot, v) + p.tf.pos;
                double s = dot(g.n, wv) + g.d;
                if (s > -1e-12) continue;
                Vec3 p_proj = wv - g.n * s;
                bool inside = false;
                for (const auto& fp : g.footprints) {
                    if (point_in_tri(p_proj, fp, g.n)) { inside = true; break; }
                }
                if (!inside) continue;
                double d_v = -s;
                Vec3 r_v = wv - p.tf.pos;
                Vec3 v_v = p.vel + cross(p.omega, r_v);
                double vn_v = dot(v_v, g.n);  // <0 = penetrating deeper
                double kn_eff_v = 1.5 * k_hw * std::sqrt(std::max(d_v, 1e-12));
                double cn_v = std::abs(2.0 * loge_w * std::sqrt(kn_eff_v * m_eff_v) /
                    std::sqrt(loge_w * loge_w + pi2));
                double fn_v = k_hw * d_v * std::sqrt(d_v);
                if (vn_v < 0.0) fn_v += -cn_v * vn_v;
                if (fn_v < 0.0) fn_v = 0.0;
                Vec3 vt_v = v_v - g.n * vn_v;
                Vec3 ft_v = vt_v * (-ct_v);
                double ftv_norm = norm(ft_v);
                double ftv_max = mu_w * fn_v;
                if (ftv_norm > ftv_max && ftv_norm > 1e-14) {
                    ft_v = ft_v * (ftv_max / ftv_norm);
                }
                Vec3 f_v = g.n * fn_v + ft_v;
                f += f_v;
                t += cross(r_v, f_v);
                any_force = true;
            }
            if (any_force) {
                n_contact_groups++;
                if (contact_debug && dbg_contacts < 3) {
                    std::cout << "  wall_contact: particle=" << dbg_particle << " wall=" << w
                              << " n_cross=" << n_cross
                              << " F=" << f.x << "," << f.y << "," << f.z
                              << "\n";
                }
            }
        }
    }
    return n_contact_groups;
}
