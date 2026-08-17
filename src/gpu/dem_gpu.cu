// src/gpu/dem_gpu.cu -- GPU DEM driver (v8 per-vertex penalty contact).
// Task 3 scope: parse config -> host build -> device upload -> optional
// --check-sums verification (device roundtrip sums vs host staged sums).
// Contact kernels and time integration arrive in later tasks.
#include <cstdio>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include "gpu/cuda_check.hpp"
#include "gpu/dem_state.cuh"
#include "gpu/real.hpp"
#include "host/config_io.hpp"

namespace {

__host__ inline const char* precision_name() {
    return sizeof(real) == 8 ? "double" : "float";
}

void usage() {
    std::printf("Usage: zdem_gpu --config path.txt [--check-sums]\n");
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

    sim.free_all();
    return 0;
}
