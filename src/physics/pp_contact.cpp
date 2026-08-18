#include "physics/pp_contact.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "core/quat.hpp"
#include "geometry/mesh_build.hpp"

using namespace zdem;

namespace {

// Per-pair penalty constants (Hertz rate, tangential rate/damping, friction
// bound) shared by every per-vertex contact element of one pair. Computed
// once per pair in the exact expression order of the original main.cpp block.
struct PPConsts {
    double k_hr;
    double kt_r;
    double m_eff_v;
    double loge_r;
    double pi2;
    double ct_v;
    double mu_r;
};

}  // namespace

// One penalty contact element: a contained vertex (wv) and its closest point
// on the other surface (cp) at depth d. The spring force acts on exactly
// those two material points; f_own/t_own accumulate the owner side,
// f_oth/t_oth the other side (same per-vertex order as the original
// forces[i_own]/forces[i_other] interleaving).
static void apply_penalty(const Vec3& wv, const Vec3& cp, double d,
                          const Particle& p_own, const Particle& p_other,
                          const PPConsts& c,
                          Vec3& f_own, Vec3& t_own, Vec3& f_oth, Vec3& t_oth) {
    Vec3 u = (wv - cp) * (1.0 / std::max(d, 1e-15));  // surface -> interior
    Vec3 n_push = u * -1.0;                            // separation direction
    Vec3 v_own = p_own.vel + cross(p_own.omega, wv - p_own.tf.pos);
    Vec3 v_oth = p_other.vel + cross(p_other.omega, cp - p_other.tf.pos);
    Vec3 vrel = v_own - v_oth;
    double vn = dot(vrel, n_push);  // <0 = penetrating deeper
    double kn_eff_v = 1.5 * c.k_hr * std::sqrt(std::max(d, 1e-12));
    double cn_v = std::abs(2.0 * c.loge_r * std::sqrt(kn_eff_v * c.m_eff_v) /
        std::sqrt(c.loge_r * c.loge_r + c.pi2));
    double fn_v = c.k_hr * d * std::sqrt(d);
    if (vn < 0.0) fn_v += -cn_v * vn;
    if (fn_v < 0.0) fn_v = 0.0;
    Vec3 fn_vec = n_push * fn_v;
    f_own += fn_vec;
    f_oth -= fn_vec;
    t_own += cross(wv - p_own.tf.pos, fn_vec);
    t_oth += cross(cp - p_other.tf.pos, fn_vec * -1.0);
    // Tangential friction at the same pair of points
    Vec3 vt = vrel - n_push * vn;
    Vec3 ft_v = vt * (-c.ct_v);
    double ftv_norm = norm(ft_v);
    double ftv_max = c.mu_r * fn_v;
    if (ftv_norm > ftv_max && ftv_norm > 1e-14) {
        ft_v = ft_v * (ftv_max / ftv_norm);
    }
    f_own += ft_v;
    f_oth -= ft_v;
    t_own += cross(wv - p_own.tf.pos, ft_v);
    t_oth += cross(cp - p_other.tf.pos, ft_v * -1.0);
}

int pp_contact_pair(const Particle& pa, const Particle& pb,
                    const Mesh& ma, const Mesh& mb,
                    Vec3& f_i, Vec3& t_i, Vec3& f_j, Vec3& t_j,
                    double tangential_damping) {
    std::vector<std::array<Vec3, 3>> trisA = transform_tris(ma, pa.tf);
    std::vector<std::array<Vec3, 3>> trisB = transform_tris(mb, pb.tf);

    // ============== Containment scan ==============
    // Vertices that have crossed INTO the other mesh: the earliest
    // contact signal for concave shapes. The segment pipeline stays
    // blind while a vertex slides through the other body's bay -- the
    // first loop used to appear with the vertex already 30 mm deep,
    // and loading the full Hertz force on that overlap teleported
    // (2/5)*kh*pen^2.5 ~ 57 J of potential energy into a 0.1 J impact.
    // Each hit records (vertex, closest surface point, depth): that
    // pair of material points is one penalty contact element.
    std::vector<Vec3> incA, incB;          // contained vertices
    std::vector<Vec3> incA_cp, incB_cp;    // their closest points on the other surface
    std::vector<double> incA_d, incB_d;    // depths |vertex - closest|
    {
        Vec3 bmn = trisB[0][0], bmx = trisB[0][0];
        for (const auto& t : trisB) {
            for (const auto& c : t) {
                bmn.x = std::min(bmn.x, c.x); bmn.y = std::min(bmn.y, c.y); bmn.z = std::min(bmn.z, c.z);
                bmx.x = std::max(bmx.x, c.x); bmx.y = std::max(bmx.y, c.y); bmx.z = std::max(bmx.z, c.z);
            }
        }
        for (const auto& v : ma.vertices) {
            Vec3 wv = quat_rotate(pa.tf.rot, v) + pa.tf.pos;
            if (wv.x < bmn.x || wv.x > bmx.x || wv.y < bmn.y || wv.y > bmx.y || wv.z < bmn.z || wv.z > bmx.z) continue;
            if (!point_inside_mesh(trisB, wv)) continue;
            Vec3 cp, nf; double d = point_mesh_distance(trisB, wv, cp, nf);
            if (d <= 1e-12) continue;
            incA.push_back(wv); incA_cp.push_back(cp); incA_d.push_back(d);
        }
        Vec3 amn = trisA[0][0], amx = trisA[0][0];
        for (const auto& t : trisA) {
            for (const auto& c : t) {
                amn.x = std::min(amn.x, c.x); amn.y = std::min(amn.y, c.y); amn.z = std::min(amn.z, c.z);
                amx.x = std::max(amx.x, c.x); amx.y = std::max(amx.y, c.y); amx.z = std::max(amx.z, c.z);
            }
        }
        for (const auto& v : mb.vertices) {
            Vec3 wv = quat_rotate(pb.tf.rot, v) + pb.tf.pos;
            if (wv.x < amn.x || wv.x > amx.x || wv.y < amn.y || wv.y > amx.y || wv.z < amn.z || wv.z > amx.z) continue;
            if (!point_inside_mesh(trisA, wv)) continue;
            Vec3 cp, nf; double d = point_mesh_distance(trisA, wv, cp, nf);
            if (d <= 1e-12) continue;
            incB.push_back(wv); incB_cp.push_back(cp); incB_d.push_back(d);
        }
    }

    f_i = Vec3{}; t_i = Vec3{}; f_j = Vec3{}; t_j = Vec3{};
    int n_inc = static_cast<int>(incA.size() + incB.size());
    if (!incA.empty() || !incB.empty()) {
        // ============== Per-vertex penalty contact ==============
        // Each contained vertex and its closest point on the other
        // surface form ONE contact element: the spring length is the
        // depth d, and the force acts on exactly the two material
        // points whose separation defines d. The spring's rate of
        // change therefore equals the relative velocity of the very
        // points the force acts on -- dU/dt = F.v holds by
        // construction. Single-point forms (loop contact point,
        // weighted average, deepest vertex) all mismatched the
        // driving quantity's time derivative against the force's
        // application point under spin and injected up to kJ.
        // The push direction is -(vertex - closest)/d: winding-free
        // (this STL's face windings are NOT consistent, 19/102 faces
        // point inward), always from the surface point away from the
        // interior point, i.e. the separation direction.
        double Rp1 = std::max(pa.equiv_radius, 1e-12);
        double Rp2 = std::max(pb.equiv_radius, 1e-12);
        double E1r = std::max(pa.young, 1e-12);
        double E2r = std::max(pb.young, 1e-12);
        double Estar_r = 1.0 / ((1.0 - pa.poisson * pa.poisson) / E1r + (1.0 - pb.poisson * pb.poisson) / E2r);
        double Rstar_r = (Rp1 * Rp2) / std::max(Rp1 + Rp2, 1e-12);
        double k_hr = (4.0 / 3.0) * Estar_r * std::sqrt(Rstar_r);
        double kt_r = 0.5 * ((E1r * Rp1) / (2.0 * (1.0 + pa.poisson)) + (E2r * Rp2) / (2.0 * (1.0 + pb.poisson)));
        double m_eff_r = 1.0 / (pa.inv_mass + pb.inv_mass);
        double m_eff_v = m_eff_r / std::max(n_inc, 1);  // each vertex carries its mass share
        double e_r = std::min(pa.restitution, pb.restitution);
        e_r = std::min(0.9999, std::max(1e-6, e_r));
        double loge_r = std::log(e_r);
        double pi2 = 3.141592653589793 * 3.141592653589793;
        double ct_v = tangential_damping * std::sqrt(kt_r * m_eff_v);
        double mu_r = std::min(pa.mu, pb.mu);
        PPConsts c{k_hr, kt_r, m_eff_v, loge_r, pi2, ct_v, mu_r};
        for (std::size_t k = 0; k < incA.size(); ++k) {
            apply_penalty(incA[k], incA_cp[k], incA_d[k], pa, pb, c, f_i, t_i, f_j, t_j);
        }
        for (std::size_t k = 0; k < incB.size(); ++k) {
            apply_penalty(incB[k], incB_cp[k], incB_d[k], pb, pa, c, f_j, t_j, f_i, t_i);
        }
    }
    return n_inc;
}
