// src/gpu/kernels/wall_contact.cuh -- per-vertex wall penalty contact (GPU Task 5)
//
// Port of the CPU wall contact (src/main.cpp "Particle-Wall Contact
// Detection", per-vertex penalty v8) expression for expression:
//   phase 1  count crossing vertices (n_cross) so the Tsuji/tangential
//            damping uses each vertex's share of the particle mass
//   phase 2  per crossing vertex: Hertz spring k_hw*d^1.5 + Tsuji damping,
//            tangential -ct_v*vt_v clamped to mu_w*fn_v, accumulate
//            force += fn_vec + ft and torque += r_v x (fn_vec + ft)
// Contact bookkeeping: one count per (particle, group) block that produced a
// contact -- CPU increments contact_counts[i] and the step `contacts` counter
// once per group with any_force; n_cross > 0 and any_force are the same set
// (phase 2 recomputes the identical deterministic test of phase 1, and every
// crossing vertex accumulates force, fn is only clamped to 0 never removed).
//
// Sign contract (differs from CPU only in generality): the CPU re-orients the
// group normal toward the particle every step, so its s = dot(n,wv)+d is
// always <= 0 for a crossing vertex (n points at the particle). Device wall
// groups were oriented toward the FIRST particle once, at upload (Task 3).
// This kernel derives the particle's side from its center of mass,
// side = sign(dot(n,pos)+d), and works the CPU's rule in that side's frame:
// a vertex penetrates iff s*side <= -1e-12 (beyond the plane, seen from the
// particle side), depth d_v = -s*side, push direction n_eff = n*side. That is
// algebraically the CPU's per-particle orientation for a particle on EITHER
// side of the plane, and numerically identical for the single-side scenarios
// of every test config (particle stays on the upload side: side=+1, the
// crossing test is the CPU's `s > -1e-12 -> skip` verbatim). A plain |s| rule
// would misfire: any vertex on the particle's own side whose projection lands
// in a footprint would count as penetrating at arbitrary distance. The
// footprint test uses the stored group normal; its outcome depends only on
// the winding consistency between footprint and that normal, not on side.
//
// Broad phase: the CPU skips a (particle, wall) pair when bounding spheres
// miss and drops wall triangles whose AABB misses the particle box. The
// kernel skips both: groups are few (one block each is cheap) and a vertex
// whose projection is inside a footprint is necessarily within the particle
// bounding sphere of the wall plane, so in the static-wall scenarios the AABB
// filter never removes a footprint hit; point_in_tri is the geometric filter.
#pragma once
#include "gpu/kernels/integrate.cuh"  // real, r_sqrt/r_log/r_fabs, quat_rotate_device

// CPU point_in_tri (src/main.cpp): the three edge cross products dotted with
// the plane normal must all be >= -1e-10 (inclusive, CPU tolerance literal).
__device__ inline bool point_in_tri_device(const real p[3], const real tri[9],
                                           const real n[3]) {
    real ab[3] = {tri[3] - tri[0], tri[4] - tri[1], tri[5] - tri[2]};
    real bc[3] = {tri[6] - tri[3], tri[7] - tri[4], tri[8] - tri[5]};
    real ca[3] = {tri[0] - tri[6], tri[1] - tri[7], tri[2] - tri[8]};
    real ap[3] = {p[0] - tri[0], p[1] - tri[1], p[2] - tri[2]};
    real bp[3] = {p[0] - tri[3], p[1] - tri[4], p[2] - tri[5]};
    real cp[3] = {p[0] - tri[6], p[1] - tri[7], p[2] - tri[8]};
    real c1 = (ab[1]*ap[2] - ab[2]*ap[1]) * n[0]
            + (ab[2]*ap[0] - ab[0]*ap[2]) * n[1]
            + (ab[0]*ap[1] - ab[1]*ap[0]) * n[2];
    real c2 = (bc[1]*bp[2] - bc[2]*bp[1]) * n[0]
            + (bc[2]*bp[0] - bc[0]*bp[2]) * n[1]
            + (bc[0]*bp[1] - bc[1]*bp[0]) * n[2];
    real c3 = (ca[1]*cp[2] - ca[2]*cp[1]) * n[0]
            + (ca[2]*cp[0] - ca[0]*cp[2]) * n[1]
            + (ca[0]*cp[1] - ca[1]*cp[0]) * n[2];
    return (c1 >= real(-1e-10)) && (c2 >= real(-1e-10)) && (c3 >= real(-1e-10));
}

// One block per (particle, wall-group): blockIdx.x = particle, blockIdx.y =
// group; the block's threads stride the particle's mesh vertices.
__global__ void wall_contact_kernel(
    const real* pos, const real* quat, const real* vel, const real* omega,
    const real* mass, const real* equiv_radius, const real* young,
    const real* poisson, const real* mu, const real* restitution,
    const int* mesh_index, const real* m_verts, const int* m_voffset,
    const real* gn, const real* gd, const real* fp, const int* fp_start,
    const real* gmu,
    const real tangential_damping,
    real* force, real* torque, int* contact_count, int* contacts,
    int n_particles, int n_groups) {
    const int i = blockIdx.x;
    const int g = blockIdx.y;
    if (i >= n_particles || g >= n_groups) return;

    const int mi = mesh_index[i];
    const int v_begin = m_voffset[mi];
    const int v_end = m_voffset[mi + 1];

    const real n[3] = {gn[3*g], gn[3*g + 1], gn[3*g + 2]};
    const real d = gd[g];
    const int fp_begin = fp_start[g];
    const int fp_end = fp_start[g + 1];

    const real q[4] = {quat[4*i], quat[4*i + 1], quat[4*i + 2], quat[4*i + 3]};
    const real ppos[3] = {pos[3*i], pos[3*i + 1], pos[3*i + 2]};

    // Particle side of the group plane (see sign contract above): all s are
    // worked in this side's frame, s_rel = s*side, so the crossing rule and
    // the force expressions mirror the CPU line for line.
    const real s_p = n[0]*ppos[0] + n[1]*ppos[1] + n[2]*ppos[2] + d;
    const real side = (s_p >= real(0)) ? real(1) : real(-1);

    // ---- phase 1: count vertices that crossed the plane inside a footprint
    // (CPU crossing rule: process iff s > -1e-12 -> skip; here the same rule
    // on s_rel). The count uses a shared counter + two barriers
    // (__syncthreadsCount is not declared by the CUDA 12.8 headers; this is
    // the same semantics).
    __shared__ int sh_count;
    if (threadIdx.x == 0) sh_count = 0;
    __syncthreads();
    int crossing = 0;
    for (int v = v_begin + threadIdx.x; v < v_end; v += blockDim.x) {
        real wv[3];
        quat_rotate_device(q, m_verts + 3*v, wv);
        wv[0] += ppos[0]; wv[1] += ppos[1]; wv[2] += ppos[2];
        real s = n[0]*wv[0] + n[1]*wv[1] + n[2]*wv[2] + d;
        if (s * side > real(-1e-12)) continue;
        real p_proj[3] = {wv[0] - n[0]*s, wv[1] - n[1]*s, wv[2] - n[2]*s};
        for (int f = fp_begin; f < fp_end; ++f) {
            if (point_in_tri_device(p_proj, fp + 9*f, n)) { crossing = 1; break; }
        }
    }
    if (crossing) atomicAdd(&sh_count, 1);
    __syncthreads();
    const int n_cross = sh_count;
    if (n_cross == 0) return;  // uniform over the block

    // ---- constants (CPU lines: Rp/Ep/Estar/k_hw/kt_w/e_w/mu_w block).
    // Recomputed per thread: a handful of flops, cheaper than a shared-memory
    // staging + extra barrier, and every thread gets the identical value.
    const real Rp = fmax(equiv_radius[i], real(1e-12));
    const real Ep = fmax(young[i], real(1e-12));
    const real nu = poisson[i];
    const real Estar = Ep / (real(2) * (real(1) - nu*nu));
    const real k_hw = (real(4) / real(3)) * Estar * r_sqrt(Rp);
    const real kt_w = (Ep * Rp) / (real(2) * (real(1) + nu));
    const real e_w = fmin(real(0.9999), fmax(real(1e-6), restitution[i]));
    const real loge = r_log(e_w);
    const real pi2 = real(3.141592653589793) * real(3.141592653589793);
    const real mu_w = fmin(mu[i], gmu[g]);
    const real m_eff_v = mass[i] / real(n_cross);   // vertex mass share
    const real ct_v = tangential_damping * r_sqrt(kt_w * m_eff_v);

    // ---- phase 2: penalty force per crossing vertex (same test, recomputed)
    const real pvel[3] = {vel[3*i], vel[3*i + 1], vel[3*i + 2]};
    const real pom[3] = {omega[3*i], omega[3*i + 1], omega[3*i + 2]};
    for (int v = v_begin + threadIdx.x; v < v_end; v += blockDim.x) {
        real wv[3];
        quat_rotate_device(q, m_verts + 3*v, wv);
        wv[0] += ppos[0]; wv[1] += ppos[1]; wv[2] += ppos[2];
        real s = n[0]*wv[0] + n[1]*wv[1] + n[2]*wv[2] + d;
        if (s * side > real(-1e-12)) continue;
        real p_proj[3] = {wv[0] - n[0]*s, wv[1] - n[1]*s, wv[2] - n[2]*s};
        bool inside = false;
        for (int f = fp_begin; f < fp_end; ++f) {
            if (point_in_tri_device(p_proj, fp + 9*f, n)) { inside = true; break; }
        }
        if (!inside) continue;

        const real n_eff[3] = {n[0]*side, n[1]*side, n[2]*side};  // push toward the particle side
        const real d_v = -s * side;                               // penetration depth (>= 1e-12)
        const real r_v[3] = {wv[0] - ppos[0], wv[1] - ppos[1], wv[2] - ppos[2]};
        // v_v = vel + omega x r_v
        const real v_v[3] = {pvel[0] + pom[1]*r_v[2] - pom[2]*r_v[1],
                             pvel[1] + pom[2]*r_v[0] - pom[0]*r_v[2],
                             pvel[2] + pom[0]*r_v[1] - pom[1]*r_v[0]};
        const real vn_v = n_eff[0]*v_v[0] + n_eff[1]*v_v[1] + n_eff[2]*v_v[2];
        const real kn_eff_v = real(1.5) * k_hw * r_sqrt(fmax(d_v, real(1e-12)));
        const real cn_v = r_fabs(real(2) * loge * r_sqrt(kn_eff_v * m_eff_v) /
                                 r_sqrt(loge*loge + pi2));
        real fn_v = k_hw * d_v * r_sqrt(d_v);
        if (vn_v < real(0)) fn_v += -cn_v * vn_v;   // damping only while closing
        if (fn_v < real(0)) fn_v = real(0);
        // tangential: -ct_v * vt on the contact plane, Coulomb clamp mu_w*fn_v
        const real vt[3] = {v_v[0] - n_eff[0]*vn_v,
                            v_v[1] - n_eff[1]*vn_v,
                            v_v[2] - n_eff[2]*vn_v};
        real ft_v[3] = {-ct_v * vt[0], -ct_v * vt[1], -ct_v * vt[2]};
        const real ftv_norm = r_sqrt(ft_v[0]*ft_v[0] + ft_v[1]*ft_v[1] + ft_v[2]*ft_v[2]);
        const real ftv_max = mu_w * fn_v;
        if (ftv_norm > ftv_max && ftv_norm > real(1e-14)) {
            const real scale = ftv_max / ftv_norm;
            ft_v[0] *= scale; ft_v[1] *= scale; ft_v[2] *= scale;
        }
        const real f_v[3] = {n_eff[0]*fn_v + ft_v[0],
                             n_eff[1]*fn_v + ft_v[1],
                             n_eff[2]*fn_v + ft_v[2]};
        const real tq[3] = {r_v[1]*f_v[2] - r_v[2]*f_v[1],
                            r_v[2]*f_v[0] - r_v[0]*f_v[2],
                            r_v[0]*f_v[1] - r_v[1]*f_v[0]};
        for (int k = 0; k < 3; ++k) {
            atomicAdd(&force[3*i + k], f_v[k]);
            atomicAdd(&torque[3*i + k], tq[k]);
        }
    }

    // one contact per (particle, group) with force: same counter semantics
    // as the CPU (contact_counts[i]++ / step contacts++ per active group)
    if (threadIdx.x == 0) {
        atomicAdd(&contact_count[i], 1);
        atomicAdd(contacts, 1);
    }
}
