#include "physics/integrate.hpp"

#include <cmath>

#include "core/mat3.hpp"
#include "core/quat.hpp"

using namespace zdem;

// Moved verbatim from src/main.cpp (2026-08 MPI extraction); the helpers below
// moved with it because they are used only by integrate_particle.

static Mat3 inertia_world_inv(const Particle& p) {
    Mat3 R = quat_to_mat3(p.tf.rot);
    Mat3 Rt = mat3_transpose(R);
    return mat3_mul(mat3_mul(R, p.inertia_body_inv), Rt);
}

// Exact exponential rotation increment for a world-frame angular velocity:
// dq = {cos(th/2), w_hat*sin(th/2)}, q' = dq (x) q. Unlike the linearized
// q += 0.5*qdot*dt this preserves the quaternion norm exactly per step.
static Quat quat_step_world(const Quat& q, const Vec3& w, double dt) {
    double wmag = std::sqrt(w.x * w.x + w.y * w.y + w.z * w.z);
    double th = wmag * dt;
    Quat dq;
    if (th < 1e-300) {
        dq.w = 1.0; dq.x = 0.0; dq.y = 0.0; dq.z = 0.0;
    } else {
        double s = std::sin(0.5 * th) / wmag;
        dq.w = std::cos(0.5 * th);
        dq.x = w.x * s; dq.y = w.y * s; dq.z = w.z * s;
    }
    return quat_mul(dq, q);
}

void integrate_particle(Particle& p, const Vec3& force, const Vec3& torque, const Vec3& gravity, double dt) {
    Vec3 acc = force * p.inv_mass + gravity;
    p.vel += acc * dt;
    p.tf.pos += p.vel * dt;

    p.L += torque * dt;

    // Midpoint attitude integration with exact exponential rotation steps.
    // The old explicit form (omega from the PRE-rotation inertia, quaternion
    // advanced linearly) injects rotational energy for an asymmetric body:
    // a torque-free spinner at 400 rad/s grew to 2875 rad/s (isolated test).
    // Here the half-step attitude is used to re-evaluate omega before the
    // full step, and each rotation increment is a norm-preserving exponential.
    Quat q0 = p.tf.rot;
    Vec3 w0 = mat3_mul_vec3(inertia_world_inv(p), p.L);
    p.tf.rot = quat_step_world(q0, w0, 0.5 * dt);
    Vec3 wh = mat3_mul_vec3(inertia_world_inv(p), p.L);
    p.tf.rot = quat_step_world(q0, wh, dt);
    p.omega = mat3_mul_vec3(inertia_world_inv(p), p.L);
}
