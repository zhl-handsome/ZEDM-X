// src/gpu/kernels/integrate.cuh -- force clear + time integration (GPU Task 4)
//
// integrate_kernel replicates the CPU integrator (src/main.cpp
// integrate_particle / quat_step_world / inertia_world_inv) expression for
// expression so FP64 runs stay bit-comparable with zdem_cpu:
//   - semi-implicit Euler translation: vel += (F*inv_m + g)*dt; pos += vel*dt
//   - L += torque*dt
//   - midpoint attitude: w0 = Iinv_world(q_old)*L -> half-step exponential
//     -> wh = Iinv_world(q_half)*L -> full-step exponential from q_old with wh
//     -> omega = Iinv_world(q_NEW)*L  (world Iinv evaluated at THREE
//        attitudes: old, half, new)
// The world inverse inertia is (R * Ibody^-1) * R^T with R built from the
// NORMALIZED quaternion (CPU quat_to_mat3 normalizes first) and the same
// two-step product association as CPU mat3_mul(mat3_mul(R, Iinv), Rt).
#pragma once
#include "gpu/real.hpp"

// --------------------------------------------------------------------------
// precision-switching math wrappers: float -> sinf/cosf/sqrtf, double stays
// with the double-precision CUDA math library.
// --------------------------------------------------------------------------
__device__ inline real r_sin(real x) {
    if constexpr (sizeof(real) == 4) return sinf(x);
    else return sin(x);
}

__device__ inline real r_cos(real x) {
    if constexpr (sizeof(real) == 4) return cosf(x);
    else return cos(x);
}

__device__ inline real r_sqrt(real x) {
    if constexpr (sizeof(real) == 4) return sqrtf(x);
    else return sqrt(x);
}

__device__ inline real r_log(real x) {
    if constexpr (sizeof(real) == 4) return logf(x);
    else return log(x);
}

__device__ inline real r_fabs(real x) {
    if constexpr (sizeof(real) == 4) return fabsf(x);
    else return fabs(x);
}

// --------------------------------------------------------------------------
// quaternion helpers (same conventions as src/core/quat.hpp)
// --------------------------------------------------------------------------
__device__ inline void quat_mul_device(const real a[4], const real b[4], real out[4]) {
    out[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    out[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
    out[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
    out[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
}

// Rotate v by q (q = {w,x,y,z}): v + 2*qv x (qv x v + w*v). q is normalized
// first so config-supplied non-unit quaternions match CPU quat_rotate (which
// normalizes internally). A degenerate norm == 0 returns v unchanged.
__device__ inline void quat_rotate_device(const real q[4], const real v[3], real out[3]) {
    real n = r_sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    if (n == real(0)) {
        out[0] = v[0]; out[1] = v[1]; out[2] = v[2];
        return;
    }
    real qw = q[0] / n;
    real qv[3] = {q[1] / n, q[2] / n, q[3] / n};
    real c[3] = {qv[1]*v[2] - qv[2]*v[1],
                 qv[2]*v[0] - qv[0]*v[2],
                 qv[0]*v[1] - qv[1]*v[0]};
    real cc[3] = {qv[1]*c[2] - qv[2]*c[1],
                  qv[2]*c[0] - qv[0]*c[2],
                  qv[0]*c[1] - qv[1]*c[0]};
    out[0] = v[0] + real(2)*(qw*c[0] + cc[0]);
    out[1] = v[1] + real(2)*(qw*c[1] + cc[1]);
    out[2] = v[2] + real(2)*(qw*c[2] + cc[2]);
}

// Exact exponential rotation increment for a world-frame angular velocity:
// dq = {cos(th/2), w_hat*sin(th/2)}, out = dq (x) q. Norm-preserving per step.
// Guard on wm == 0 (NOT th < 1e-300: in float the literal underflows to 0 so
// the guard would never fire and sin(0)/0 would produce NaN).
__device__ inline void quat_step_world_device(const real q[4], const real w[3],
                                              real dt, real out[4]) {
    real wm = r_sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
    real th = wm * dt;
    real dq[4];
    if (wm == real(0)) {
        dq[0] = real(1); dq[1] = real(0); dq[2] = real(0); dq[3] = real(0);
    } else {
        real half_th = real(0.5) * th;
        real s = r_sin(half_th) / wm;
        dq[0] = r_cos(half_th);
        dq[1] = w[0]*s; dq[2] = w[1]*s; dq[3] = w[2]*s;
    }
    quat_mul_device(dq, q, out);
}

// --------------------------------------------------------------------------
// world-frame inverse inertia
// --------------------------------------------------------------------------
// Row-major R from quaternion, exactly the CPU quat_to_mat3 expansion with
// the quaternion NORMALIZED first (CPU normalizes inside quat_to_mat3; a
// degenerate near-zero quaternion yields the identity rotation).
__device__ inline void quat_to_mat3_device(const real q[4], real R[9]) {
    real n = r_sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    real w, x, y, z;
    if (n < real(1e-14)) {
        w = real(1); x = real(0); y = real(0); z = real(0);
    } else {
        w = q[0]/n; x = q[1]/n; y = q[2]/n; z = q[3]/n;
    }
    R[0] = real(1) - real(2)*(y*y + z*z);
    R[1] = real(2)*(x*y - z*w);
    R[2] = real(2)*(x*z + y*w);
    R[3] = real(2)*(x*y + z*w);
    R[4] = real(1) - real(2)*(x*x + z*z);
    R[5] = real(2)*(y*z - x*w);
    R[6] = real(2)*(x*z - y*w);
    R[7] = real(2)*(y*z + x*w);
    R[8] = real(1) - real(2)*(x*x + y*y);
}

// W = (R * Ibody^-1) * R^T, ibinv row-major [9], same two-step association as
// CPU mat3_mul(mat3_mul(R, inertia_body_inv), transpose(R)).
__device__ inline void inertia_world_inv_device(const real q[4], const real ibinv[9],
                                                real W[9]) {
    real R[9];
    quat_to_mat3_device(q, R);
    real T[9];  // T = R * Ibody^-1
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            real s = real(0);
            for (int k = 0; k < 3; ++k) s += R[3*a + k] * ibinv[3*k + b];
            T[3*a + b] = s;
        }
    }
    // W = T * R^T  (R^T[k][b] = R[b][k] in row-major storage)
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            real s = real(0);
            for (int k = 0; k < 3; ++k) s += T[3*a + k] * R[3*b + k];
            W[3*a + b] = s;
        }
    }
}

__device__ inline void mat3_mul_vec3_device(const real m[9], const real v[3], real out[3]) {
    out[0] = m[0]*v[0] + m[1]*v[1] + m[2]*v[2];
    out[1] = m[3]*v[0] + m[4]*v[1] + m[5]*v[2];
    out[2] = m[6]*v[0] + m[7]*v[1] + m[8]*v[2];
}

// --------------------------------------------------------------------------
// kernels
// --------------------------------------------------------------------------
// NOTE: no `inline` on the kernels -- nvcc ignores it for __global__ anyway.
// contacts[0] is the per-step contact counter (wall groups here, pp pairs in
// Task 7): zeroed here, incremented once per contacting block by the contact
// kernels, read back by the driver for the step log. Thread 0 of block 0
// zeroes it; the counter is only touched by contact kernels AFTER this
// kernel completes (stream ordered), so no race.
__global__ void clear_forces_kernel(real* force, real* torque, int* cc,
                                    int* contacts, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i == 0) contacts[0] = 0;
    if (i >= n) return;
    for (int k = 0; k < 3; ++k) {
        force[3*i + k] = real(0);
        torque[3*i + k] = real(0);
    }
    cc[i] = 0;
}

__global__ void integrate_kernel(real* pos, real* vel, real* omega, real* quat, real* L,
                                        const real* force, const real* torque,
                                        const real* inv_mass, const real* inertia_body_inv,
                                        const real* gravity, real dt, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    // semi-implicit Euler translation
    for (int k = 0; k < 3; ++k) {
        real acc = force[3*i + k] * inv_mass[i] + gravity[k];
        vel[3*i + k] += acc * dt;
        pos[3*i + k] += vel[3*i + k] * dt;
    }
    for (int k = 0; k < 3; ++k) L[3*i + k] += torque[3*i + k] * dt;

    // midpoint attitude (world Iinv evaluated at old / half / new attitude)
    real q0[4] = {quat[4*i], quat[4*i + 1], quat[4*i + 2], quat[4*i + 3]};
    real Li[3] = {L[3*i], L[3*i + 1], L[3*i + 2]};
    const real* ibinv = inertia_body_inv + 9*i;

    real W[9], w0[3], w1[3], wh[3], qh[4], q1[4];
    inertia_world_inv_device(q0, ibinv, W);
    mat3_mul_vec3_device(W, Li, w0);
    quat_step_world_device(q0, w0, real(0.5)*dt, qh);       // half-step attitude
    inertia_world_inv_device(qh, ibinv, W);
    mat3_mul_vec3_device(W, Li, wh);                        // omega at half step
    quat_step_world_device(q0, wh, dt, q1);                 // full step from q0
    for (int k = 0; k < 4; ++k) quat[4*i + k] = q1[k];
    inertia_world_inv_device(q1, ibinv, W);
    mat3_mul_vec3_device(W, Li, w1);                        // omega at NEW attitude
    for (int k = 0; k < 3; ++k) omega[3*i + k] = w1[k];
}
