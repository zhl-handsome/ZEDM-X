// src/gpu/kernels/pp_contact.cuh -- particle-particle containment penalty
// contact (GPU Task 7).
//
// Port of the CPU particle-particle path (src/main.cpp, "Containment scan"
// block + the per-vertex penalty apply, through `continue`), expression for
// expression:
//   1. one block (128 threads) per broadphase candidate pair (i, j)
//   2. world AABB of both meshes (threads stride vertices, block reduce)
//   3. containment scan: every vertex of A inside B's world AABB gets a ray
//      parity test (point_inside_mesh) and, if inside, an Ericson
//      closest-point depth (point_mesh_distance); symmetric for B in A
//   4. no contained vertices on either side -> no contact, block returns
//   5. per contained vertex: Hertz penalty spring + Tsuji damping along the
//      winding-free separation direction -(wv-cp)/d, distributed tangential
//      damping-friction clamped to mu*fn, at exactly the two material points
//      (wv on the owner, cp on the other) that define the spring length
//   6. contact_count[i]++ / [j]++ and the step contacts counter, once per
//      pair -- same bookkeeping as the CPU
//
// Determinism notes (why FP64 parity holds at the tolerance level, and where
// bit-exactness stops):
//   - the containment decisions (ray parity = integer crossing count,
//     order-independent; closest point = strict-< min tracking over the
//     triangle loop, same serial order 0..F-1 as the CPU per vertex) are
//     computed ENTIRELY by one thread per vertex, so every (inside?, cp, d)
//     decision is bit-identical to the CPU in FP64
//   - the mesh world AABB min/max reduction across threads is exact (no
//     rounding), any reduction order matches the CPU loop
//   - the force ACCUMULATION order (atomicAdd across contained vertices and
//     across pairs/wall blocks) is not the CPU's sequential order: sums can
//     differ in the last ulp. With ~10 contained vertices the drift over
//     20k steps stays far inside the parity tolerances (pos 1e-6, energy
//     1e-9); bit-exactness is NOT claimed here (it IS for the single-source
//     wall and free-spin paths).
//
// Shared hit list: kPPMaxHits contained vertices per side. The banana mesh
// has 53 vertices, so 64 is unreachable overflow; for a hypothetical mesh
// with more contained vertices than the cap the DEEPEST ones are kept
// (arrival-order slots, overflow path evicts the shallowest stored hit under
// a block-local spin lock -- greedy eviction keeps the kPPMaxHits deepest
// set) and n_inc counts the stored hits only, a documented divergence from
// the CPU (which processes all).
#pragma once
#include "gpu/kernels/integrate.cuh"  // real, r_sqrt/r_log/r_fabs, quat_rotate_device

namespace {
constexpr int kPPThreads = 128;  // threads per pair block
constexpr int kPPMaxHits = 64;   // contained-vertex slots per side (banana: 53 verts)
}  // namespace

// Largest finite real (initial value for min/max reductions; CPU uses
// std::numeric_limits<double>::max() for the distance sentinel).
__device__ inline real pp_huge() {
    if constexpr (sizeof(real) == 4) return real(3.402823466e+38f);
    else return real(1.7976931348623157e+308);
}

// World-space corner c of triangle t: quat_rotate(q, verts[vidx]) + pos --
// the exact transform_tris expression (quat_rotate normalizes internally on
// both sides).
__device__ inline void pp_world_corner(const real q[4], const real p[3],
                                       const real* verts, const int* tris,
                                       int t, int c, real out[3]) {
    const int vidx = tris[3 * t + c];
    real w[3];
    quat_rotate_device(q, verts + 3 * (size_t)vidx, w);
    out[0] = w[0] + p[0];
    out[1] = w[1] + p[1];
    out[2] = w[2] + p[2];
}

// CPU point_inside_mesh (src/geometry/mesh_build.cpp): +x ray parity.
// Serial triangle loop in the same 0..F-1 order; the crossing count is an
// integer, so the result never depends on iteration order.
__device__ inline bool point_inside_mesh_device(
    const real q[4], const real p[3], const real* verts, const int* tris,
    int t_begin, int t_end, const real pt[3]) {
    int crossings = 0;
    for (int t = t_begin; t < t_end; ++t) {
        real a[3], b[3], c[3];
        pp_world_corner(q, p, verts, tris, t, 0, a);
        pp_world_corner(q, p, verts, tris, t, 1, b);
        pp_world_corner(q, p, verts, tris, t, 2, c);
        // Bounding-box shortcut on y/z before the full test.
        const real ymin = fmin(fmin(a[1], b[1]), c[1]);
        const real ymax = fmax(fmax(a[1], b[1]), c[1]);
        const real zmin = fmin(fmin(a[2], b[2]), c[2]);
        const real zmax = fmax(fmax(a[2], b[2]), c[2]);
        if (pt[1] < ymin || pt[1] > ymax || pt[2] < zmin || pt[2] > zmax) continue;
        // Solve pt + t*(1,0,0) = a + u*(b-a) + v*(c-a) for the y,z rows.
        const real e1y = b[1] - a[1], e1z = b[2] - a[2];
        const real e2y = c[1] - a[1], e2z = c[2] - a[2];
        const real det = e1y * e2z - e1z * e2y;
        if (r_fabs(det) < real(1e-14)) continue;  // ray parallel to triangle plane
        const real dy = pt[1] - a[1], dz = pt[2] - a[2];
        const real u = (dy * e2z - dz * e2y) / det;
        const real v = (e1y * dz - e1z * dy) / det;
        if (u < real(-1e-12) || v < real(-1e-12) || u + v > real(1) + real(1e-12)) continue;
        const real xhit = a[0] + u * (b[0] - a[0]) + v * (c[0] - a[0]);
        if (xhit > pt[0]) crossings++;
    }
    return (crossings % 2) == 1;
}

// CPU point_mesh_distance (Ericson, Real-Time Collision Detection 5.1.5) over
// ALL triangles, strict-< min tracking in the same serial order (the CPU's
// unused out_closest_normal is dropped here). Returns the distance, writes
// the closest surface point to out_cp.
__device__ inline real point_mesh_distance_device(
    const real q[4], const real p[3], const real* verts, const int* tris,
    int t_begin, int t_end, const real pt[3], real out_cp[3]) {
    real best2 = pp_huge();
    real bq[3] = {real(0), real(0), real(0)};
    for (int t = t_begin; t < t_end; ++t) {
        real a[3], b[3], c[3];
        pp_world_corner(q, p, verts, tris, t, 0, a);
        pp_world_corner(q, p, verts, tris, t, 1, b);
        pp_world_corner(q, p, verts, tris, t, 2, c);
        real ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        real ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        real ap[3] = {pt[0] - a[0], pt[1] - a[1], pt[2] - a[2]};
        const real d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
        const real d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];
        real cpt[3];  // closest point on this triangle (NOT the quaternion q!)
        if (d1 <= real(0) && d2 <= real(0)) {
            cpt[0] = a[0]; cpt[1] = a[1]; cpt[2] = a[2];          // vertex A region
        } else {
            real bp[3] = {pt[0] - b[0], pt[1] - b[1], pt[2] - b[2]};
            const real d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
            const real d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];
            if (d3 >= real(0) && d4 <= d3) {
                cpt[0] = b[0]; cpt[1] = b[1]; cpt[2] = b[2];      // vertex B region
            } else {
                real cp[3] = {pt[0] - c[0], pt[1] - c[1], pt[2] - c[2]};
                const real d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
                const real d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];
                if (d6 >= real(0) && d5 <= d6) {
                    cpt[0] = c[0]; cpt[1] = c[1]; cpt[2] = c[2];  // vertex C region
                } else if (d1*d4 - d3*d2 <= real(0) && d1 >= real(0) && d3 <= real(0)) {
                    const real s = d1 / (d1 - d3);                 // edge AB
                    cpt[0] = a[0] + ab[0]*s; cpt[1] = a[1] + ab[1]*s; cpt[2] = a[2] + ab[2]*s;
                } else if (d5*d2 - d6*d1 <= real(0) && d6 >= real(0) && d5 <= d3) {
                    const real s = (d6 - d5) / (d3 - d5);          // edge AC
                    cpt[0] = c[0] + ac[0]*s; cpt[1] = c[1] + ac[1]*s; cpt[2] = c[2] + ac[2]*s;
                } else {
                    const real va = d3 * d6 - d5 * d4;             // face region
                    const real vb = d5 * d2 - d1 * d6;
                    const real vc = d1 * d4 - d3 * d2;
                    const real denom = real(1) / (va + vb + vc);
                    const real v = vb * denom;
                    const real w = vc * denom;
                    cpt[0] = a[0] + ab[0]*v + ac[0]*w;
                    cpt[1] = a[1] + ab[1]*v + ac[1]*w;
                    cpt[2] = a[2] + ab[2]*v + ac[2]*w;
                }
            }
        }
        const real dx = pt[0] - cpt[0], dy = pt[1] - cpt[1], dz = pt[2] - cpt[2];
        const real dist2 = dx*dx + dy*dy + dz*dz;
        if (dist2 < best2) {
            best2 = dist2;
            bq[0] = cpt[0]; bq[1] = cpt[1]; bq[2] = cpt[2];
        }
    }
    out_cp[0] = bq[0]; out_cp[1] = bq[1]; out_cp[2] = bq[2];
    return r_sqrt(best2);
}

// World AABB of one mesh: threads stride body vertices (quat_rotate + pos,
// the exact transform_tris corner expression; every registry vertex is
// referenced by a triangle, so this equals the CPU's per-triangle-corner
// loop), per-thread local min/max staged to shared, merged by thread 0.
// min/max are exact -- no reduction-order rounding. Contains two barriers;
// out_min/out_max are shared arrays written by thread 0 only.
__device__ inline void pp_block_world_aabb(
    const real q[4], const real p[3], const real* verts, int vb, int ve,
    real* sh_min, real* sh_max, real out_min[3], real out_max[3]) {
    real lmn[3] = {pp_huge(), pp_huge(), pp_huge()};
    real lmx[3] = {-pp_huge(), -pp_huge(), -pp_huge()};
    for (int v = vb + (int)threadIdx.x; v < ve; v += (int)blockDim.x) {
        real w[3];
        quat_rotate_device(q, verts + 3 * (size_t)v, w);
        for (int k = 0; k < 3; ++k) {
            const real x = w[k] + p[k];
            lmn[k] = fmin(lmn[k], x);
            lmx[k] = fmax(lmx[k], x);
        }
    }
    for (int k = 0; k < 3; ++k) {
        sh_min[3 * threadIdx.x + k] = lmn[k];
        sh_max[3 * threadIdx.x + k] = lmx[k];
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        for (int k = 0; k < 3; ++k) {
            real mn = pp_huge(), mx = -pp_huge();
            for (int t = 0; t < (int)blockDim.x; ++t) {
                mn = fmin(mn, sh_min[3 * t + k]);
                mx = fmax(mx, sh_max[3 * t + k]);
            }
            out_min[k] = mn;
            out_max[k] = mx;
        }
    }
    __syncthreads();
}

// Append one contained vertex (wv, cp, d) to a side's shared hit list.
// Slots are arrival order (atomicAdd); beyond the cap the DEEPEST set is
// kept by evicting the shallowest stored hit under a block-local spin lock
// (greedy eviction: the final stored set is the kPPMaxHits deepest hits).
__device__ inline void pp_store_hit(
    int* s_cnt, int* s_lock, real* hwv, real* hcp, real* hd,
    const real wv[3], const real cp[3], real d) {
    const int slot = atomicAdd(s_cnt, 1);
    if (slot < kPPMaxHits) {
        hwv[3 * slot + 0] = wv[0]; hwv[3 * slot + 1] = wv[1]; hwv[3 * slot + 2] = wv[2];
        hcp[3 * slot + 0] = cp[0]; hcp[3 * slot + 1] = cp[1]; hcp[3 * slot + 2] = cp[2];
        hd[slot] = d;
        return;
    }
    while (atomicExch(s_lock, 1) == 1) { /* spin */ }
    int shallow = 0;
    real dmin = hd[0];
    for (int k = 1; k < kPPMaxHits; ++k) {
        if (hd[k] < dmin) { dmin = hd[k]; shallow = k; }
    }
    if (d > dmin) {  // deeper than the shallowest stored hit -> evict it
        hwv[3 * shallow + 0] = wv[0]; hwv[3 * shallow + 1] = wv[1]; hwv[3 * shallow + 2] = wv[2];
        hcp[3 * shallow + 0] = cp[0]; hcp[3 * shallow + 1] = cp[1]; hcp[3 * shallow + 2] = cp[2];
        hd[shallow] = d;
    }
    atomicExch(s_lock, 0);
}

// One block of kPPThreads per candidate pair. blockDim.x MUST be
// kPPThreads: the shared staging below is sized by kPPThreads while the
// launch site (dem_gpu.cu) passes its own kThreads -- two constants that
// only coincide by convention. The guard at the top of the kernel turns a
// mismatched launch into a no-op instead of out-of-bounds shared writes
// (s_red is indexed 3*threadIdx.x for every launched thread).
__global__ void pp_contact_kernel(
    /* particle state (const) */ const real* pos, const real* quat,
    const real* vel, const real* omega, const real* inv_mass,
    const real* radius, const real* equiv_radius, const real* young,
    const real* poisson, const real* mu, const real* restitution,
    const int* mesh_index,
    /* mesh registry */ const real* m_verts, const int* m_tris,
    const int* m_voffset, const int* m_toffset,
    /* pairs */ const int* d_pairs, int n_pairs,
    /* params */ real tangential_damping,
    /* out */ real* force, real* torque, int* contact_count, int* contacts) {
    if (blockDim.x != kPPThreads) return;  // shared staging sized kPPThreads
    (void)radius;  // broad phase already applied the bounding-sphere test
    const int p = blockIdx.x;
    if (p >= n_pairs) return;
    const int i = d_pairs[2 * p + 0];   // i < j (broadphase enumeration)
    const int j = d_pairs[2 * p + 1];
    const int tid = (int)threadIdx.x;

    // ---- pair state + mesh descriptors (every thread loads its own copy)
    const real qi[4] = {quat[4*i], quat[4*i+1], quat[4*i+2], quat[4*i+3]};
    const real qj[4] = {quat[4*j], quat[4*j+1], quat[4*j+2], quat[4*j+3]};
    const real pi[3] = {pos[3*i], pos[3*i+1], pos[3*i+2]};
    const real pj[3] = {pos[3*j], pos[3*j+1], pos[3*j+2]};
    const int mi = mesh_index[i], mj = mesh_index[j];
    const int ivb = m_voffset[mi], ive = m_voffset[mi + 1];
    const int jvb = m_voffset[mj], jve = m_voffset[mj + 1];
    const int itb = m_toffset[mi], ite = m_toffset[mi + 1];
    const int jtb = m_toffset[mj], jte = m_toffset[mj + 1];

    __shared__ int s_cntA, s_cntB, s_lockA, s_lockB;
    __shared__ real s_aabbA_min[3], s_aabbA_max[3];
    __shared__ real s_aabbB_min[3], s_aabbB_max[3];
    __shared__ real s_red[2 * kPPThreads * 3];  // AABB staging: min half + max half
    __shared__ real s_hitA_wv[3 * kPPMaxHits], s_hitA_cp[3 * kPPMaxHits], s_hitA_d[kPPMaxHits];
    __shared__ real s_hitB_wv[3 * kPPMaxHits], s_hitB_cp[3 * kPPMaxHits], s_hitB_d[kPPMaxHits];
    if (tid == 0) { s_cntA = 0; s_cntB = 0; s_lockA = 0; s_lockB = 0; }
    __syncthreads();

    // ---- phase 1: both meshes' world AABBs (CPU: bmn/bmx over tris corners)
    pp_block_world_aabb(qi, pi, m_verts, ivb, ive,
                        s_red, s_red + kPPThreads * 3, s_aabbA_min, s_aabbA_max);
    pp_block_world_aabb(qj, pj, m_verts, jvb, jve,
                        s_red, s_red + kPPThreads * 3, s_aabbB_min, s_aabbB_max);

    // ---- phase 2: containment scan. Thread-per-vertex (strided); one
    // thread does the whole serial triangle work per vertex, in the CPU's
    // order, so every (inside?, cp, d) triple is bit-identical to the CPU.
    for (int v = ivb + tid; v < ive; v += kPPThreads) {
        real wv[3];
        quat_rotate_device(qi, m_verts + 3 * (size_t)v, wv);
        wv[0] += pi[0]; wv[1] += pi[1]; wv[2] += pi[2];
        if (wv[0] < s_aabbB_min[0] || wv[0] > s_aabbB_max[0] ||
            wv[1] < s_aabbB_min[1] || wv[1] > s_aabbB_max[1] ||
            wv[2] < s_aabbB_min[2] || wv[2] > s_aabbB_max[2]) continue;
        if (!point_inside_mesh_device(qj, pj, m_verts, m_tris, jtb, jte, wv)) continue;
        real cp[3];
        const real d = point_mesh_distance_device(qj, pj, m_verts, m_tris, jtb, jte, wv, cp);
        if (d <= real(1e-12)) continue;
        pp_store_hit(&s_cntA, &s_lockA, s_hitA_wv, s_hitA_cp, s_hitA_d, wv, cp, d);
    }
    for (int v = jvb + tid; v < jve; v += kPPThreads) {
        real wv[3];
        quat_rotate_device(qj, m_verts + 3 * (size_t)v, wv);
        wv[0] += pj[0]; wv[1] += pj[1]; wv[2] += pj[2];
        if (wv[0] < s_aabbA_min[0] || wv[0] > s_aabbA_max[0] ||
            wv[1] < s_aabbA_min[1] || wv[1] > s_aabbA_max[1] ||
            wv[2] < s_aabbA_min[2] || wv[2] > s_aabbA_max[2]) continue;
        if (!point_inside_mesh_device(qi, pi, m_verts, m_tris, itb, ite, wv)) continue;
        real cp[3];
        const real d = point_mesh_distance_device(qi, pi, m_verts, m_tris, itb, ite, wv, cp);
        if (d <= real(1e-12)) continue;
        pp_store_hit(&s_cntB, &s_lockB, s_hitB_wv, s_hitB_cp, s_hitB_d, wv, cp, d);
    }
    __syncthreads();

    const int cntA = s_cntA < kPPMaxHits ? s_cntA : kPPMaxHits;
    const int cntB = s_cntB < kPPMaxHits ? s_cntB : kPPMaxHits;
    const int n_inc = cntA + cntB;   // total contained vertices both sides
    if (n_inc == 0) return;          // uniform over the block: no contact

    // ---- phase 3: pair constants (CPU Rp/Estar/Rstar/k_hr/kt_r/m_eff_v
    // block; recomputed per thread -- identical math on identical inputs).
    const real Rp1 = fmax(equiv_radius[i], real(1e-12));
    const real Rp2 = fmax(equiv_radius[j], real(1e-12));
    const real E1r = fmax(young[i], real(1e-12));
    const real E2r = fmax(young[j], real(1e-12));
    const real nu1 = poisson[i], nu2 = poisson[j];
    const real Estar = real(1) / ((real(1) - nu1 * nu1) / E1r + (real(1) - nu2 * nu2) / E2r);
    const real Rstar = (Rp1 * Rp2) / fmax(Rp1 + Rp2, real(1e-12));
    const real k_hr = (real(4) / real(3)) * Estar * r_sqrt(Rstar);
    const real kt_r = real(0.5) * ((E1r * Rp1) / (real(2) * (real(1) + nu1)) +
                                   (E2r * Rp2) / (real(2) * (real(1) + nu2)));
    const real m_eff = real(1) / (inv_mass[i] + inv_mass[j]);
    const real m_eff_v = m_eff / (real)n_inc;  // each vertex carries its mass share
    real e_r = fmin(restitution[i], restitution[j]);
    e_r = fmin(real(0.9999), fmax(real(1e-6), e_r));
    const real loge = r_log(e_r);
    const real pi2 = real(3.141592653589793) * real(3.141592653589793);
    const real ct_v = tangential_damping * r_sqrt(kt_r * m_eff_v);
    const real mu_r = fmin(mu[i], mu[j]);

    // ---- phase 4: penalty per contained vertex (CPU apply_penalty lambda,
    // verbatim). Threads stride the shared hit list; A-side hits act on
    // (owner=i at wv, other=j at cp), B-side symmetrically.
    for (int k = tid; k < n_inc; k += kPPThreads) {
        const bool sideA = k < cntA;
        const int idx = sideA ? k : k - cntA;
        const real* hwv = sideA ? s_hitA_wv : s_hitB_wv;
        const real* hcp = sideA ? s_hitA_cp : s_hitB_cp;
        const real* hd = sideA ? s_hitA_d : s_hitB_d;
        const int own = sideA ? i : j;
        const int oth = sideA ? j : i;

        const real wv[3] = {hwv[3*idx], hwv[3*idx+1], hwv[3*idx+2]};
        const real cp[3] = {hcp[3*idx], hcp[3*idx+1], hcp[3*idx+2]};
        const real d = hd[idx];

        const real* vel_o = vel + 3 * (size_t)own;
        const real* om_o = omega + 3 * (size_t)own;
        const real* pos_o = pos + 3 * (size_t)own;
        const real* vel_t = vel + 3 * (size_t)oth;
        const real* om_t = omega + 3 * (size_t)oth;
        const real* pos_t = pos + 3 * (size_t)oth;

        // push direction: surface point -> interior point, negated
        // (winding-free separation direction; this STL has 19/102 faces in)
        const real idmax = real(1) / fmax(d, real(1e-15));
        const real n_push[3] = {-(wv[0] - cp[0]) * idmax,
                                -(wv[1] - cp[1]) * idmax,
                                -(wv[2] - cp[2]) * idmax};
        const real r_own[3] = {wv[0] - pos_o[0], wv[1] - pos_o[1], wv[2] - pos_o[2]};
        const real r_oth[3] = {cp[0] - pos_t[0], cp[1] - pos_t[1], cp[2] - pos_t[2]};
        // v_own = vel + omega x r_own ; v_oth likewise at the cp material point
        const real v_own[3] = {vel_o[0] + om_o[1]*r_own[2] - om_o[2]*r_own[1],
                               vel_o[1] + om_o[2]*r_own[0] - om_o[0]*r_own[2],
                               vel_o[2] + om_o[0]*r_own[1] - om_o[1]*r_own[0]};
        const real v_oth[3] = {vel_t[0] + om_t[1]*r_oth[2] - om_t[2]*r_oth[1],
                               vel_t[1] + om_t[2]*r_oth[0] - om_t[0]*r_oth[2],
                               vel_t[2] + om_t[0]*r_oth[1] - om_t[1]*r_oth[0]};
        const real vrel[3] = {v_own[0] - v_oth[0], v_own[1] - v_oth[1], v_own[2] - v_oth[2]};
        const real vn = vrel[0]*n_push[0] + vrel[1]*n_push[1] + vrel[2]*n_push[2];  // <0 = closing

        // Hertz spring + Tsuji damping (restitution-closed) on this vertex
        const real kn_eff_v = real(1.5) * k_hr * r_sqrt(fmax(d, real(1e-12)));
        const real cn_v = r_fabs(real(2) * loge * r_sqrt(kn_eff_v * m_eff_v) /
                                 r_sqrt(loge * loge + pi2));
        real fn_v = k_hr * d * r_sqrt(d);
        if (vn < real(0)) fn_v += -cn_v * vn;   // damping only while closing
        if (fn_v < real(0)) fn_v = real(0);
        const real fn_vec[3] = {n_push[0]*fn_v, n_push[1]*fn_v, n_push[2]*fn_v};
        for (int c = 0; c < 3; ++c) {
            atomicAdd(&force[3*own + c], fn_vec[c]);
            atomicAdd(&force[3*oth + c], -fn_vec[c]);
        }
        // torques[i_own] += (wv - pos_own) x fn_vec ;
        // torques[i_other] += (cp - pos_other) x (-fn_vec)
        atomicAdd(&torque[3*own + 0], r_own[1]*fn_vec[2] - r_own[2]*fn_vec[1]);
        atomicAdd(&torque[3*own + 1], r_own[2]*fn_vec[0] - r_own[0]*fn_vec[2]);
        atomicAdd(&torque[3*own + 2], r_own[0]*fn_vec[1] - r_own[1]*fn_vec[0]);
        atomicAdd(&torque[3*oth + 0], -(r_oth[1]*fn_vec[2] - r_oth[2]*fn_vec[1]));
        atomicAdd(&torque[3*oth + 1], -(r_oth[2]*fn_vec[0] - r_oth[0]*fn_vec[2]));
        atomicAdd(&torque[3*oth + 2], -(r_oth[0]*fn_vec[1] - r_oth[1]*fn_vec[0]));

        // tangential friction at the same pair of points, Coulomb-clamped
        const real vt[3] = {vrel[0] - n_push[0]*vn,
                            vrel[1] - n_push[1]*vn,
                            vrel[2] - n_push[2]*vn};
        real ft_v[3] = {vt[0] * (-ct_v), vt[1] * (-ct_v), vt[2] * (-ct_v)};
        const real ftv_norm = r_sqrt(ft_v[0]*ft_v[0] + ft_v[1]*ft_v[1] + ft_v[2]*ft_v[2]);
        const real ftv_max = mu_r * fn_v;
        if (ftv_norm > ftv_max && ftv_norm > real(1e-14)) {
            const real scale = ftv_max / ftv_norm;
            ft_v[0] *= scale; ft_v[1] *= scale; ft_v[2] *= scale;
        }
        for (int c = 0; c < 3; ++c) {
            atomicAdd(&force[3*own + c], ft_v[c]);
            atomicAdd(&force[3*oth + c], -ft_v[c]);
        }
        atomicAdd(&torque[3*own + 0], r_own[1]*ft_v[2] - r_own[2]*ft_v[1]);
        atomicAdd(&torque[3*own + 1], r_own[2]*ft_v[0] - r_own[0]*ft_v[2]);
        atomicAdd(&torque[3*own + 2], r_own[0]*ft_v[1] - r_own[1]*ft_v[0]);
        atomicAdd(&torque[3*oth + 0], -(r_oth[1]*ft_v[2] - r_oth[2]*ft_v[1]));
        atomicAdd(&torque[3*oth + 1], -(r_oth[2]*ft_v[0] - r_oth[0]*ft_v[2]));
        atomicAdd(&torque[3*oth + 2], -(r_oth[0]*ft_v[1] - r_oth[1]*ft_v[0]));
    }

    // ---- phase 5: contact bookkeeping once per pair (CPU contact_counts
    // +=1 for BOTH i and j, step contacts counter +=1)
    if (tid == 0) {
        atomicAdd(&contact_count[i], 1);
        atomicAdd(&contact_count[j], 1);
        atomicAdd(contacts, 1);
    }
}
