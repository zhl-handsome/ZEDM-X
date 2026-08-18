#pragma once

#include <mpi.h>

#include <vector>

#include "host/sim_build.hpp"

// Shared particle-shipping primitives for the two exchange families:
// the ghost halo rebuild (ghost.cpp) and centroid migration (migrate.cpp).
// Both ship variable-length particle lists between face neighbors, so both
// use the same count-then-data Sendrecv pattern; only the tag scheme
// differs (each caller derives its own tags and passes the count tag in).

// One particle on the wire. POD, shipped as raw MPI_BYTES -- the same
// trick gather.cpp uses for FrameSnapshot (homogeneous-cluster assumption:
// same executable, same ABI).
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
// dir=+1 sends sendbuf to the +axis neighbor and receives that neighbor's
// outbound list back (what crossed its hi face into our brick); dir=-1
// mirrors that across the lo face. MPI_PROC_NULL at non-periodic box edges
// turns each Sendrecv into a local no-op (recv count stays 0) -- callers
// that must not silently lose particles (migration) detect out-of-box
// positions before calling.
//
// Tags: tag_count is the caller-supplied base tag, tag_data = tag_count+1.
// Every message moving in the +axis direction is sent by a dir=+1 call and
// received by a dir=+1 call on the far side (Cart_shift's source/dest are
// symmetric across the shared face), so a tag derived from (axis, dir)
// matches on both ends. Count and data get separate tags.
inline std::vector<PackedPart> sendrecv_particles(MPI_Comm cart, int axis, int dir,
                                                  int tag_count,
                                                  const std::vector<PackedPart>& sendbuf) {
    int src = MPI_PROC_NULL, dst = MPI_PROC_NULL;
    MPI_Cart_shift(cart, axis, dir, &src, &dst);

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
