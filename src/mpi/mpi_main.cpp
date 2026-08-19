// src/mpi/mpi_main.cpp -- MPI DEM driver (v8, Newton off, static bricks).
// Full time loop: local-local pp pairs (same i<j order and double-sided
// whole-block adds as the CPU driver), ghost pp pairs (own half only,
// Newton off; the halo is rebuilt by the chained three-axis exchange after
// every step), wall contact for local particles against the replicated
// wall set, local integration, gather-based VTK output on rank 0, then
// centroid migration (particles that left their brick change rank) and the
// halo rebuild so the next step sees fresh ownership.
// Per-step semantics (wall contact_count accumulation, per-step contacts
// counter, output condition/frame set/log line) are copied from
// src/main.cpp.
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "mpi.h"

#include "host/config_io.hpp"
#include "host/sim_build.hpp"
#include "geometry/mesh_build.hpp"
#include "mpi/decomp.hpp"
#include "mpi/gather.hpp"
#include "mpi/ghost.hpp"
#include "mpi/migrate.hpp"
#include "physics/broadphase.hpp"
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
    // This is also the standing invariant for migration: this phase never
    // inserts or deletes particles, so the global N stays fixed for the
    // whole run.
    long long tot = assert_global_count(d, (long long)local.size(),
                                        (long long)sim.particles.size());

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
    exchange_ghosts(d, local, gids, ghost);   // step 0's halo (rebuilt at the end of every step)

    std::vector<Vec3> forces, torques;
    std::vector<int> cc;
    // Per-rank step timing, reset at every output step: phys_ns covers the
    // force-computation + integrate sections (world-triangle cache,
    // broadphase, pp, walls, integrate), comm_ns the migration + ghost
    // exchange sections. Reported per output step as the per-rank MAX,
    // averaged over the interval, so slow-rank skew (hybrid cores, load
    // imbalance) shows up in the log.
    std::chrono::steady_clock::duration phys_ns{}, comm_ns{};
    int steps_since_out = 0;
    for (int step = 0; step < cfg.steps; ++step) {
        forces.assign(local.size(), Vec3{});
        torques.assign(local.size(), Vec3{});
        cc.assign(local.size(), 0);
        int total_contacts = 0;   // per-step counter, same as the CPU driver
        const auto phys_t0 = std::chrono::steady_clock::now();

        // World triangles cached once per step (local + ghost copies):
        // transform_tris is pure and tf is frozen within the step, so this
        // is bit-identical to the per-pair recomputation inside
        // pp_contact_pair -- it just stops re-transforming a mesh for every
        // pair it appears in. Rebuilt after the end-of-step migration/
        // exchange so indices always match the current local/ghost arrays.
        std::vector<std::vector<std::array<Vec3, 3>>> wtris_local(local.size());
        for (int li = 0; li < (int)local.size(); ++li) {
            wtris_local[li] = transform_tris(sim.meshes[local[li].mesh_index], local[li].tf);
        }
        std::vector<std::vector<std::array<Vec3, 3>>> wtris_ghost(ghost.particles.size());
        for (int gj = 0; gj < (int)ghost.particles.size(); ++gj) {
            wtris_ghost[gj] = transform_tris(sim.meshes[ghost.particles[gj].mesh_index],
                                             ghost.particles[gj].tf);
        }

        // ---- particle-particle (Newton off) ----
        // Broadphase: one shared spatial-hash call over the combined
        // [local; ghost] array replaces the former O(N^2) local double loop
        // plus full ghost scan. The sorted (i, j) pair list reproduces the
        // old per-li nesting exactly: for each local i, the local partners
        // (j < n_local, ascending, j > i) come first, then the ghosts
        // (j >= n_local, ascending) -- lexicographic order on the combined
        // array IS that order. Ghost-ghost pairs (first index >= n_local)
        // are dropped inside broadphase_pairs; the drivers never computed
        // them. The prechecks below are bit-identical no-op copies of the
        // test already applied inside broadphase_pairs (kept conservative
        // per the task brief).
        std::vector<Particle> all;
        all.reserve(local.size() + ghost.particles.size());
        all.insert(all.end(), local.begin(), local.end());
        all.insert(all.end(), ghost.particles.begin(), ghost.particles.end());
        const int n_local = (int)local.size();
        std::vector<std::pair<int, int>> pp_pairs;
        broadphase_pairs(all, n_local, pp_pairs);
        for (const auto& [li, lj] : pp_pairs) {
            if (lj < n_local) {
                // Local pair: same i<j order, bounding-sphere prefilter and
                // double-sided whole-block adds as the CPU loop; at n=1 this
                // covers every pair exactly like the CPU.
                Vec3 dpos = local[lj].tf.pos - local[li].tf.pos;
                double rsum = local[li].radius + local[lj].radius;
                if (dot(dpos, dpos) > rsum * rsum) {
                    continue;
                }
                Vec3 f_i, t_i, f_j, t_j;
                if (pp_contact_pair(local[li], local[lj],
                                    sim.meshes[local[li].mesh_index],
                                    sim.meshes[local[lj].mesh_index],
                                    wtris_local[li], wtris_local[lj],
                                    f_i, t_i, f_j, t_j,
                                    cfg.tangential_damping) > 0) {
                    forces[li] += f_i; torques[li] += t_i;
                    forces[lj] += f_j; torques[lj] += t_j;
                    cc[li]++; cc[lj]++;
                    total_contacts++;
                }
            } else {
                // Ghost pair: this rank owns li, so keep only the own half
                // (f_i/t_i); the ghost's owner adds the mirrored other half.
                const int gj = lj - n_local;
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
                                    wtris_local[li], wtris_ghost[gj],
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
        phys_ns += std::chrono::steady_clock::now() - phys_t0;

        // ---- output: same condition, frame set and log fields as the CPU ----
        if (step % cfg.output_interval == 0) {
            gather_write_frame(d, cfg, sim, local, gids, cc, forces, torques, step);
            int glob_contacts = reduce_add(d, total_contacts);
            // Per-rank timing maxima: ONE Allreduce carries both
            // accumulators (packed as a 2-element array). Collective
            // hygiene: every rank calls the reduce unconditionally, only
            // the printf is rank-gated (this codebase once deadlocked on
            // a reduce inside a rank-0 branch).
            long long tim_in[2] = {
                std::chrono::duration_cast<std::chrono::nanoseconds>(phys_ns).count(),
                std::chrono::duration_cast<std::chrono::nanoseconds>(comm_ns).count()};
            long long tim_max[2] = {0, 0};
            MPI_Allreduce(tim_in, tim_max, 2, MPI_LONG_LONG, MPI_MAX, d.cart);
            phys_ns = {};   // reset right after the Allreduce
            comm_ns = {};
            // Interval-average per-step milliseconds (integer division);
            // step 0 has zero elapsed steps, so it reports zeros.
            long long ms_max = 0, comm_ms_max = 0;
            if (steps_since_out > 0) {
                ms_max = tim_max[0] / steps_since_out / 1000000LL;
                comm_ms_max = tim_max[1] / steps_since_out / 1000000LL;
            }
            if (d.rank == 0) {
                std::printf("step=%d contacts=%d ms_max=%lld comm_ms_max=%lld\n",
                            step, glob_contacts, ms_max, comm_ms_max);
                std::fflush(stdout);
            }
            // Per-output-step ownership distribution (design-spec load
            // diagnostics): every rank contributes its size, rank 0 prints
            // rank_n=a,b,c. Migration shows up here as counts shifting
            // between ranks after a boundary crossing.
            std::vector<int> nlocal_hint(1, (int)local.size());
            log_rank_n(d, nlocal_hint);
            steps_since_out = 0;
        }

        // ---- migration + halo rebuild: particles whose centroid left its
        // brick move to the owner rank, then the halo is rebuilt so the
        // next step's force computation sees fresh ownership and ghost
        // positions. Output already happened above (per the plan's Global
        // Constraints: force/cc arrays are indexed by the pre-migration
        // local order). ----
        // The halo rebuilt here is consumed only by the NEXT step's force
        // pass; after the last step nothing reads it, so skip the final
        // exchange (one all-rank communication saved at run end).
        const auto comm_t0 = std::chrono::steady_clock::now();
        migrate_particles(d, local, gids);
        if (step + 1 < cfg.steps) {
            exchange_ghosts(d, local, gids, ghost);
        }
        comm_ns += std::chrono::steady_clock::now() - comm_t0;
        ++steps_since_out;
    }

    if (d.rank == 0) {
        std::printf("done\n");
    }

    MPI_Finalize();
    return 0;
}
