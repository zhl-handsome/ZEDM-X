// src/gpu/dem_gpu.cu -- GPU DEM driver (v8 per-vertex penalty contact).
// Task 4 scope: parse config -> host build -> device upload -> main loop
// (clear_forces -> [contact kernels in Tasks 5/7] -> integrate) with
// periodic VTK output through the shared host writer.
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include "gpu/cuda_check.hpp"
#include "gpu/dem_state.cuh"
#include "gpu/kernels/integrate.cuh"
#include "gpu/kernels/wall_contact.cuh"
#include "gpu/real.hpp"
#include "host/config_io.hpp"
#include "host/vtk_io.hpp"

namespace {

constexpr int kThreads = 128;

__host__ inline const char* precision_name() {
    return sizeof(real) == 8 ? "double" : "float";
}

void usage() {
    std::printf("Usage: zdem_gpu --config path.txt [--check-sums]\n");
}

// Per-output-interval snapshot: download the full frame, refresh the host
// particle registry copy with device state, rotate mesh vertices to world in
// write_vtk_particles (same CELL_DATA fields and file naming as zdem_cpu).
void output_vtk(const GpuSim& sim, const SimConfig& cfg, int step) {
    std::vector<double> pos, vel, omega, quat, force, torque;
    std::vector<int> cc;
    sim.download_frame(pos, vel, omega, quat, force, torque, cc);

    const int n = sim.P.n;
    std::vector<Particle> parts = sim.host_particles;
    std::vector<Vec3> forces(n), torques(n);
    for (int i = 0; i < n; ++i) {
        Particle& p = parts[i];
        p.tf.pos = Vec3{pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]};
        p.vel = Vec3{vel[3 * i], vel[3 * i + 1], vel[3 * i + 2]};
        p.omega = Vec3{omega[3 * i], omega[3 * i + 1], omega[3 * i + 2]};
        p.tf.rot = Quat{quat[4 * i], quat[4 * i + 1], quat[4 * i + 2], quat[4 * i + 3]};
        forces[i] = Vec3{force[3 * i], force[3 * i + 1], force[3 * i + 2]};
        torques[i] = Vec3{torque[3 * i], torque[3 * i + 1], torque[3 * i + 2]};
    }

    std::ostringstream oss;
    oss << cfg.vtk_prefix << "_" << std::setw(6) << std::setfill('0') << step << ".vtk";
    std::filesystem::path out_path = std::filesystem::path(cfg.output_dir) / oss.str();
    write_vtk_particles(out_path.string(), sim.host_meshes, parts, forces, torques, cc);
}

}  // namespace

int main(int argc, char** argv) {
    std::string config_path;
    bool check_sums = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--config" && i + 1 < argc) {
            config_path = argv[++i];
        } else if (arg == "--check-sums") {
            check_sums = true;
        } else if (arg == "--help") {
            usage();
            return 0;
        } else {
            std::fprintf(stderr, "Unknown arg: %s\n", arg.c_str());
            usage();
            return 1;
        }
    }
    if (config_path.empty()) {
        usage();
        return 1;
    }

    SimConfig cfg;
    if (!parse_config_file(config_path, cfg)) {
        return 1;
    }

    std::printf("zdem_gpu (CUDA, precision=%s)\n", precision_name());

    GpuSim sim;
    sim.upload(cfg);

    if (check_sums) {
        // Device roundtrip: copy the raw real buffers back and sum as double.
        const std::size_t n = static_cast<std::size_t>(sim.P.n);
        double dev_pos = 0.0;
        double dev_mass = 0.0;
        if (n > 0) {
            std::vector<real> rpos(3 * n);
            CUDA_CHECK(cudaMemcpy(rpos.data(), sim.P.pos,
                                  3 * n * sizeof(real), cudaMemcpyDeviceToHost));
            for (real v : rpos) dev_pos += static_cast<double>(v);
            std::vector<real> rmass(n);
            CUDA_CHECK(cudaMemcpy(rmass.data(), sim.P.mass,
                                  n * sizeof(real), cudaMemcpyDeviceToHost));
            for (real v : rmass) dev_mass += static_cast<double>(v);
        }
        std::printf("checksum pos=%.17g mass=%.17g\n", dev_pos, dev_mass);
        std::printf("hostsum pos=%.17g mass=%.17g\n",
                    sim.host_pos_sum, sim.host_mass_sum);
        // Exercise the generic download path as well (pos SoA -> double).
        std::vector<double> pos, vel, omega;
        sim.download_soa(pos, vel, omega);
        double dl_pos = 0.0;
        for (double v : pos) dl_pos += v;
        std::printf("download_soa n=%zu pos=%.17g\n", pos.size() / 3, dl_pos);
    }

    std::filesystem::create_directories(cfg.output_dir);
    const int n = sim.P.n;
    const int blocks = (n + kThreads - 1) / kThreads;
    std::printf("simulating: n=%d steps=%d dt=%g output_interval=%d\n",
                n, cfg.steps, cfg.dt, cfg.output_interval);

    for (int step = 0; step < cfg.steps; ++step) {
        auto step_t0 = std::chrono::steady_clock::now();
        if (n > 0) {
            clear_forces_kernel<<<blocks, kThreads>>>(
                sim.P.force, sim.P.torque, sim.P.contact_count, sim.P.contacts, n);
            CUDA_CHECK(cudaGetLastError());
            // Wall contact: one block per (particle, wall-group). Same counter
            // semantics as the CPU step log (wall groups with contact; the
            // particle-particle kernel joins the same accumulator in Task 7).
            if (sim.W.n_groups > 0) {
                dim3 wc_grid(n, sim.W.n_groups);
                wall_contact_kernel<<<wc_grid, kThreads>>>(
                    sim.P.pos, sim.P.quat, sim.P.vel, sim.P.omega,
                    sim.P.mass, sim.P.equiv_radius, sim.P.young,
                    sim.P.poisson, sim.P.mu, sim.P.restitution,
                    sim.P.mesh_index, sim.M.d_verts, sim.M.d_voffset,
                    sim.W.d_gn, sim.W.d_gd, sim.W.d_fp, sim.W.d_fp_start,
                    sim.W.d_mu, sim.tangential_damping,
                    sim.P.force, sim.P.torque, sim.P.contact_count, sim.P.contacts,
                    n, sim.W.n_groups);
                CUDA_CHECK(cudaGetLastError());
            }
            integrate_kernel<<<blocks, kThreads>>>(
                sim.P.pos, sim.P.vel, sim.P.omega, sim.P.quat, sim.P.L,
                sim.P.force, sim.P.torque, sim.P.inv_mass, sim.P.inertia_body_inv,
                sim.d_gravity, sim.dt, n);
            CUDA_CHECK(cudaGetLastError());
        }
        if (step % cfg.output_interval == 0) {
            CUDA_CHECK(cudaDeviceSynchronize());
            int contacts = 0;
            CUDA_CHECK(cudaMemcpy(&contacts, sim.P.contacts, sizeof(int),
                                  cudaMemcpyDeviceToHost));
            output_vtk(sim, cfg, step);
            double step_ms = std::chrono::duration<double, std::milli>(
                                 std::chrono::steady_clock::now() - step_t0)
                                 .count();
            std::printf("step=%d contacts=%d step_ms=%.1f\n", step, contacts, step_ms);
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    std::printf("done\n");

    sim.free_all();
    return 0;
}
