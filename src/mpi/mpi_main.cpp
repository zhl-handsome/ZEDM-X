// src/mpi/mpi_main.cpp -- MPI DEM driver (v8, Newton off, static bricks).
// Full time loop: local-local pp pairs (same i<j order and double-sided
// whole-block adds as the CPU driver), ghost pp pairs (own half only,
// Newton off; halo is empty until Task 5), wall contact for local particles
// against the replicated wall set, local integration, then gather-based VTK
// output on rank 0. Per-step semantics (wall contact_count accumulation,
// per-step contacts counter, output condition/frame set/log line) are
// copied from src/main.cpp.
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "mpi.h"

#include "host/config_io.hpp"
#include "host/sim_build.hpp"
#include "mpi/decomp.hpp"
#include "mpi/gather.hpp"
#include "mpi/ghost.hpp"
#include "physics/integrate.hpp"
#include "physics/pp_contact.hpp"
#include "physics/wall_contact.hpp"

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

    // Only rank 0 writes frames (gather output), so only it needs the dir.
    if (d.rank == 0) {
        std::filesystem::create_directories(cfg.output_dir);
    }

    GhostLayer ghost;
    exchange_ghosts(d, local, gids, ghost);   // n=1: empty; Task 5 fills this

    std::vector<Vec3> forces, torques;
    std::vector<int> cc;
    for (int step = 0; step < cfg.steps; ++step) {
        forces.assign(local.size(), Vec3{});
        torques.assign(local.size(), Vec3{});
        cc.assign(local.size(), 0);
        int total_contacts = 0;   // per-step counter, same as the CPU driver

        // ---- particle-particle (Newton off) ----
        for (int li = 0; li < (int)local.size(); ++li) {
            // Local pairs: same i<j order, bounding-sphere prefilter and
            // double-sided whole-block adds as the CPU loop; at n=1 this
            // covers every pair exactly like the CPU.
            for (int lj = li + 1; lj < (int)local.size(); ++lj) {
                Vec3 dpos = local[lj].tf.pos - local[li].tf.pos;
                double rsum = local[li].radius + local[lj].radius;
                if (dot(dpos, dpos) > rsum * rsum) {
                    continue;
                }
                Vec3 f_i, t_i, f_j, t_j;
                if (pp_contact_pair(local[li], local[lj],
                                    sim.meshes[local[li].mesh_index],
                                    sim.meshes[local[lj].mesh_index],
                                    f_i, t_i, f_j, t_j,
                                    cfg.tangential_damping) > 0) {
                    forces[li] += f_i; torques[li] += t_i;
                    forces[lj] += f_j; torques[lj] += t_j;
                    cc[li]++; cc[lj]++;
                    total_contacts++;
                }
            }
            // Ghost pairs: this rank owns li, so keep only the own half
            // (f_i/t_i); the ghost's owner adds the mirrored other half.
            for (int gj = 0; gj < (int)ghost.particles.size(); ++gj) {
                if (ghost.gids[gj] == gids[li]) continue;   // defensive; migration keeps gids unique
                const Particle& g = ghost.particles[gj];
                Vec3 dv = g.tf.pos - local[li].tf.pos;
                double rs = local[li].radius + g.radius;
                if (dot(dv, dv) > rs * rs) {
                    continue;
                }
                Vec3 f_i, t_i, f_j, t_j;
                if (pp_contact_pair(local[li], g,
                                    sim.meshes[local[li].mesh_index],
                                    sim.meshes[g.mesh_index],
                                    f_i, t_i, f_j, t_j,
                                    cfg.tangential_damping) > 0) {
                    forces[li] += f_i; torques[li] += t_i;
                    cc[li]++;
                    if (gids[li] < ghost.gids[gj]) total_contacts++;   // count each cross-rank pair once (mirror CPU i<j)
                }
            }
        }

        // ---- walls: local particles x all walls (walls are replicated) ----
        // CPU semantics (main.cpp): wall_contact_particle accumulates the
        // per-vertex forces into forces[i]/torques[i] in place and returns
        // the contacted wall-group count; contact_counts[i] and the per-step
        // contacts counter both add the raw group count (not clamped to 1).
        for (int li = 0; li < (int)local.size(); ++li) {
            int n_wall_contacts = wall_contact_particle(local[li],
                                                        sim.meshes[local[li].mesh_index],
                                                        sim.walls,
                                                        cfg.tangential_damping,
                                                        forces[li], torques[li]);
            cc[li] += n_wall_contacts;
            total_contacts += n_wall_contacts;
        }

        // ---- integrate locals ----
        for (int li = 0; li < (int)local.size(); ++li) {
            integrate_particle(local[li], forces[li], torques[li], cfg.gravity, cfg.dt);
        }

        // ---- output: same condition, frame set and log fields as the CPU ----
        if (step % cfg.output_interval == 0) {
            gather_write_frame(d, cfg, sim, local, gids, cc, forces, torques, step);
            int glob_contacts = reduce_add(d, total_contacts);
            if (d.rank == 0) {
                std::printf("step=%d contacts=%d\n", step, glob_contacts);
                std::fflush(stdout);
            }
        }

        // ---- migration (Task 6) and ghost rebuild (Task 5) hook here ----
    }

    if (d.rank == 0) {
        std::printf("done\n");
    }

    MPI_Finalize();
    return 0;
}
