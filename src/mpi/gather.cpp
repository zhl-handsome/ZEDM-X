#include "mpi/gather.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

void gather_write_frame(const Decomp& d, const SimConfig& cfg, const SimBuild& sim,
                        const std::vector<Particle>& local, const std::vector<int>& gids,
                        const std::vector<int>& cc, const std::vector<Vec3>& forces,
                        const std::vector<Vec3>& torques, int step) {
    // Pack this rank's particles into flat snapshot records.
    std::vector<FrameSnapshot> snaps(local.size());
    for (std::size_t i = 0; i < local.size(); ++i) {
        snaps[i].p = local[i];
        snaps[i].gid = gids[i];
        snaps[i].cc = cc[i];
        snaps[i].force = forces[i];
        snaps[i].torque = torques[i];
    }

    // counts/displs in MPI_BYTE units, sized from each rank's record count
    // (Allgather of the per-rank element counts).
    const int nlocal = static_cast<int>(snaps.size());
    const int snap_bytes = static_cast<int>(sizeof(FrameSnapshot));
    std::vector<int> elems(d.nprocs, 0);
    MPI_Allgather(&nlocal, 1, MPI_INT, elems.data(), 1, MPI_INT, d.cart);
    std::vector<int> bcounts(d.nprocs, 0);
    std::vector<int> bdispls(d.nprocs, 0);
    int total_elems = 0;
    for (int r = 0; r < d.nprocs; ++r) {
        bcounts[r] = elems[r] * snap_bytes;
        bdispls[r] = total_elems * snap_bytes;
        total_elems += elems[r];
    }

    std::vector<FrameSnapshot> all;
    if (d.rank == 0) {
        all.resize(total_elems);
    }
    MPI_Gatherv(snaps.data(), nlocal * snap_bytes, MPI_BYTE,
                all.data(), bcounts.data(), bdispls.data(), MPI_BYTE, 0, d.cart);
    if (d.rank != 0) {
        return;
    }

    // gid order == config order == the CPU particle order, so the assembled
    // frame (mesh layout, id field, per-cell arrays) is identical.
    std::sort(all.begin(), all.end(),
              [](const FrameSnapshot& a, const FrameSnapshot& b) { return a.gid < b.gid; });

    std::vector<Particle> particles;
    std::vector<Vec3> f;
    std::vector<Vec3> t;
    std::vector<int> counts;
    particles.reserve(all.size());
    f.reserve(all.size());
    t.reserve(all.size());
    counts.reserve(all.size());
    for (const FrameSnapshot& s : all) {
        particles.push_back(s.p);
        f.push_back(s.force);
        t.push_back(s.torque);
        counts.push_back(s.cc);
    }

    // Same construction as the CPU driver: vtk_prefix_%06d.vtk under
    // cfg.output_dir (setw(6) + setfill('0') == "%06d").
    char stem[32];
    std::snprintf(stem, sizeof(stem), "%06d", step);
    std::filesystem::path out_path =
        std::filesystem::path(cfg.output_dir) / (cfg.vtk_prefix + "_" + stem + ".vtk");
    write_vtk_particles(out_path.string(), sim.meshes, particles, f, t, counts);
}

int reduce_add(const Decomp& d, int local_val) {
    int total = 0;
    MPI_Reduce(&local_val, &total, 1, MPI_INT, MPI_SUM, 0, d.cart);
    return d.rank == 0 ? total : 0;
}

void log_rank_n(const Decomp& d, const std::vector<int>& local_sizes_hint) {
    const int n_send = static_cast<int>(local_sizes_hint.size());
    std::vector<int> counts(d.nprocs, 0);
    MPI_Allgather(&n_send, 1, MPI_INT, counts.data(), 1, MPI_INT, d.cart);
    std::vector<int> displs(d.nprocs, 0);
    int total = 0;
    for (int r = 0; r < d.nprocs; ++r) {
        displs[r] = total;
        total += counts[r];
    }
    std::vector<int> all(static_cast<std::size_t>(total > 0 ? total : 1), 0);
    MPI_Allgatherv(local_sizes_hint.data(), n_send, MPI_INT,
                   all.data(), counts.data(), displs.data(), MPI_INT, d.cart);
    if (d.rank != 0) {
        return;
    }
    std::printf("rank_n=");
    for (int r = 0; r < d.nprocs; ++r) {
        for (int k = 0; k < counts[r]; ++k) {
            std::printf("%s%d", (r == 0 && k == 0) ? "" : ",", all[displs[r] + k]);
        }
    }
    std::printf("\n");
}
