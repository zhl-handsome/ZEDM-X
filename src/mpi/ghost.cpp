#include "mpi/ghost.hpp"

#include <algorithm>
#include <cstdio>
#include <vector>

namespace {

// One ghost particle on the wire. POD, shipped as raw MPI_BYTES -- the same
// trick gather.cpp uses for FrameSnapshot.
struct PackedPart {
    Particle p;
    int gid;
};

inline double axis_val(const Vec3& pos, int axis) {
    return axis == 0 ? pos.x : (axis == 1 ? pos.y : pos.z);
}

// Swap a variable-length PackedPart list with one face-neighbor pair:
// first a 1-int count Sendrecv, then the fixed-size records. Two-phase
// counting avoids MPI_Probe entirely and is portable across MPI stacks.
//
// dir=+1 sends our to_hi window to the +axis neighbor and receives its
// to_hi window back (particles just below our sub_lo face); dir=-1 mirrors
// that across the lo face. MPI_PROC_NULL at non-periodic box edges turns
// each Sendrecv into a local no-op (recv count stays 0).
//
// Tags: every message moving in the +axis direction is sent by a dir=+1
// call and received by a dir=+1 call on the far side (Cart_shift's
// source/dest are symmetric across the shared face), so a tag derived from
// (axis, dir) matches on both ends. Count and data get separate tags.
std::vector<PackedPart> exchange_face(MPI_Comm cart, int axis, int dir,
                                      const std::vector<PackedPart>& sendbuf) {
    int src = MPI_PROC_NULL, dst = MPI_PROC_NULL;
    MPI_Cart_shift(cart, axis, dir, &src, &dst);

    const int tag_count = 4 * axis + (dir > 0 ? 0 : 2);
    const int tag_data = tag_count + 1;

    const int n_send = static_cast<int>(sendbuf.size());
    int n_recv = 0;
    MPI_Sendrecv(&n_send, 1, MPI_INT, dst, tag_count,
                 &n_recv, 1, MPI_INT, src, tag_count, cart, MPI_STATUS_IGNORE);

    std::vector<PackedPart> recvbuf(n_recv);
    const int unit = static_cast<int>(sizeof(PackedPart));
    MPI_Sendrecv(n_send > 0 ? sendbuf.data() : nullptr,
                 n_send * unit, MPI_BYTE, dst, tag_data,
                 n_recv > 0 ? recvbuf.data() : nullptr,
                 n_recv * unit, MPI_BYTE, src, tag_data,
                 cart, MPI_STATUS_IGNORE);
    return recvbuf;
}

}  // namespace

void exchange_ghosts(const Decomp& d, const std::vector<Particle>& local,
                     const std::vector<int>& gids, GhostLayer& out) {
    // Chained three-axis halo exchange (LAMMPS comm_brick style). The
    // working set W starts as the locally owned particles; ghosts received
    // on axis k join W before axis k+1 is packed, so a diagonal ghost
    // arrives in two hops (axis i, then axis j > i).
    //
    // A particle reaches this rank through at most one path: hops along the
    // axes on which owner and receiver differ, in increasing axis order --
    // any second path would need a second hop on some already-finished
    // axis, which the per-axis ordering forbids. make_decomp also
    // guarantees brick thickness >= ghost_depth, so only face-adjacent
    // bricks can ever exchange. The duplicate check at the end guards both
    // invariants at runtime.
    std::vector<PackedPart> W;
    W.reserve(local.size());
    for (std::size_t i = 0; i < local.size(); ++i) {
        W.push_back(PackedPart{local[i], gids[i]});
    }
    const std::size_t n_local = W.size();

    for (int axis = 0; axis < 3; ++axis) {
        // Pack both face windows from the same W snapshot. Strict bounds: a
        // particle exactly ghost_depth inside a face cannot overlap anything
        // across it (contact range r_i + r_j <= 2*max_radius == ghost_depth).
        std::vector<PackedPart> to_hi, to_lo;
        for (const PackedPart& w : W) {
            const double v = axis_val(w.p.tf.pos, axis);
            if (v > d.sub_hi[axis] - d.ghost_depth) to_hi.push_back(w);
            if (v < d.sub_lo[axis] + d.ghost_depth) to_lo.push_back(w);
        }
        std::vector<PackedPart> got_lo = exchange_face(d.cart, axis, +1, to_hi);
        std::vector<PackedPart> got_hi = exchange_face(d.cart, axis, -1, to_lo);
        W.insert(W.end(), got_lo.begin(), got_lo.end());
        W.insert(W.end(), got_hi.begin(), got_hi.end());
    }

    // Everything past the local prefix of W is a received ghost.
    out.particles.clear();
    out.gids.clear();
    out.particles.reserve(W.size() - n_local);
    out.gids.reserve(W.size() - n_local);
    for (std::size_t i = n_local; i < W.size(); ++i) {
        out.particles.push_back(W[i].p);
        out.gids.push_back(W[i].gid);
    }

    // Development guard, kept in Release: sorting one small gid copy is
    // cheap next to the pp pass. A duplicate gid (two ghosts, or a ghost
    // colliding with a local) means the chain shipped the same particle
    // along two paths -- the Newton-off force halves would be applied
    // twice, so die loudly instead of corrupting the dynamics.
    {
        std::vector<int> check;
        check.reserve(W.size());
        for (const PackedPart& w : W) check.push_back(w.gid);
        std::sort(check.begin(), check.end());
        std::vector<int>::const_iterator dup =
            std::adjacent_find(check.begin(), check.end());
        if (dup != check.end()) {
            std::fprintf(stderr,
                         "zdem_mpi fatal: rank %d got duplicate particle gid %d "
                         "after ghost exchange (brick/ghost invariant broken)\n",
                         d.rank, *dup);
            std::fflush(stderr);
            MPI_Abort(d.cart, 1);
        }
    }
}
