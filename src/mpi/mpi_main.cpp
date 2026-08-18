// src/mpi/mpi_main.cpp -- MPI DEM driver (v8, Newton off, static bricks).
// Skeleton scope: init -> config -> build -> per-rank filter -> report.
// The time loop lands in a later task; cfg.steps is not consumed yet.
#include <cstdio>
#include <string>
#include <vector>

#include "mpi.h"

#include "host/config_io.hpp"
#include "host/sim_build.hpp"
#include "mpi/decomp.hpp"

static bool parse_args(int argc, char** argv, std::string& config_path) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--config" && i + 1 < argc) {
            config_path = argv[++i];
        } else {
            return false;
        }
    }
    return !config_path.empty();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int wrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    std::string config_path;
    if (!parse_args(argc, argv, config_path)) {
        if (wrank == 0) std::fprintf(stderr, "Usage: zdem_mpi --config path.txt\n");
        MPI_Finalize();
        return 1;
    }
    SimConfig cfg;
    if (!parse_config_file(config_path, cfg)) {
        MPI_Finalize();
        return 1;
    }

    SimBuild sim;
    std::string err;
    if (!build_sim(cfg, sim, err)) {
        if (wrank == 0) std::fprintf(stderr, "build_sim failed: %s\n", err.c_str());
        MPI_Finalize();
        return 1;
    }

    Decomp d = make_decomp(cfg, sim);

    // Initial distribution: every rank builds the full config, then keeps the
    // particles whose center falls in its brick (gid = config index, carried
    // alongside the particle).
    std::vector<Particle> local;
    std::vector<int> gids;
    for (int i = 0; i < (int)sim.particles.size(); ++i) {
        if (in_sub(d, sim.particles[i].tf.pos)) {
            local.push_back(sim.particles[i]);
            gids.push_back(i);
        }
    }

    // Conservation check: every particle is owned by exactly one rank.
    long long loc = (long long)local.size();
    long long tot = 0;
    MPI_Allreduce(&loc, &tot, 1, MPI_LONG_LONG, MPI_SUM, d.cart);
    if (tot != (long long)sim.particles.size()) {
        if (wrank == 0) {
            std::fprintf(stderr, "initial distribution lost/duplicated particles: sum(nlocal)=%lld != N=%lld\n",
                         tot, (long long)sim.particles.size());
        }
        MPI_Abort(d.cart, 1);
    }

    if (d.rank == 0) {
        std::printf("zdem_mpi: nprocs=%d dims=%dx%dx%d box=[%.3f %.3f]-[%.3f %.3f]-[%.3f %.3f] ghost=%.3f N=%lld\n",
                    d.nprocs, d.dims[0], d.dims[1], d.dims[2],
                    d.box_lo[0], d.box_hi[0], d.box_lo[1], d.box_hi[1], d.box_lo[2], d.box_hi[2],
                    d.ghost_depth, tot);
    }
    std::printf("  rank %d: nlocal=%zu gids=[%d..%d]\n",
                d.rank, local.size(),
                local.empty() ? -1 : gids.front(),
                local.empty() ? -1 : gids.back());
    std::fflush(stdout);

    MPI_Finalize();
    return 0;
}
