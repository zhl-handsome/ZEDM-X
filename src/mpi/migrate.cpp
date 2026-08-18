#include "mpi/migrate.hpp"

#include <cstdarg>
#include <cstdio>
#include <vector>

#include "mpi/comm_util.hpp"

namespace {

// Migration tags: 32+4*axis+{0..3}, disjoint from the ghost exchange's
// 0..11 range. (Each Sendrecv completes before either side issues its next
// call, so overlap is impossible anyway; the disjoint range just keeps the
// two exchange families distinguishable in trace logs.)
const int kMigrateTagBase = 32;

// The detecting rank logs (the offending particle lives on ONE rank, so
// unlike decomp.cpp's identical-input die() only that rank has the facts;
// a rank-0-only print would silently lose the message), then the whole job
// dies. Mirrors decomp.cpp's die() otherwise: migration failure is global
// and unrecoverable.
[[noreturn]] void die(const Decomp& d, const char* fmt, ...) {
    std::fprintf(stderr, "zdem_mpi fatal (rank %d): ", d.rank);
    va_list ap;
    va_start(ap, fmt);
    std::vfprintf(stderr, fmt, ap);
    va_end(ap);
    std::fprintf(stderr, "\n");
    std::fflush(stderr);
    MPI_Abort(d.cart, 1);
}

// Same left-closed right-open convention as the brick indexing: the box is
// [box_lo, box_hi)^3, so box_hi itself counts as outside (its owner by
// floor would be the last brick, but positions only ever land there by
// rounding at the closed end).
bool out_of_box(const Decomp& d, const Vec3& pos) {
    return pos.x < d.box_lo[0] || pos.x >= d.box_hi[0] ||
           pos.y < d.box_lo[1] || pos.y >= d.box_hi[1] ||
           pos.z < d.box_lo[2] || pos.z >= d.box_hi[2];
}

void check_in_box(const Decomp& d, const std::vector<PackedPart>& W) {
    for (const PackedPart& w : W) {
        if (out_of_box(d, w.p.tf.pos)) {
            die(d,
                "particle gid %d left the global box: pos=(%.6f %.6f %.6f) "
                "box=[%.3f %.3f]-[%.3f %.3f]-[%.3f %.3f]",
                w.gid, w.p.tf.pos.x, w.p.tf.pos.y, w.p.tf.pos.z,
                d.box_lo[0], d.box_hi[0], d.box_lo[1], d.box_hi[1],
                d.box_lo[2], d.box_hi[2]);
        }
    }
}

}  // namespace

void migrate_particles(const Decomp& d, std::vector<Particle>& local, std::vector<int>& gids) {
    // Working set W: particles currently on this rank. Round 0 starts from
    // the owned set; later rounds re-test whatever is still displaced
    // (multi-brick movers -- one round moves at most one brick per axis).
    std::vector<PackedPart> W;
    W.reserve(local.size());
    for (std::size_t i = 0; i < local.size(); ++i) {
        W.push_back(PackedPart{local[i], gids[i]});
    }

    for (int round = 0; round < 8; ++round) {
        // Out-of-box must abort BEFORE any send: at a box edge the face
        // neighbor is MPI_PROC_NULL and a send there silently drops the
        // particle. Positions do not change inside migrate_particles, so
        // this also covers everything shipped in this round.
        check_in_box(d, W);

        if (round > 0) {
            // Continue only while some rank still holds a displaced
            // particle. The Allreduce keeps the round count identical on
            // every rank -- a rank that quit early would miss particles
            // shipped to it by a slower chain.
            int local_out = 0;
            for (const PackedPart& w : W) {
                if (!in_sub(d, w.p.tf.pos)) {
                    local_out = 1;
                    break;
                }
            }
            int any_out = 0;
            MPI_Allreduce(&local_out, &any_out, 1, MPI_INT, MPI_MAX, d.cart);
            if (any_out == 0) {
                break;
            }
        }

        // Chained three-axis handoff, same exchange pattern as
        // exchange_ghosts: a particle received on axis k is re-tested on
        // axis k+1, so a corner crossing moves two bricks in one round.
        // Sent particles leave W immediately (no duplicate copies exist at
        // any instant).
        for (int axis = 0; axis < 3; ++axis) {
            std::vector<PackedPart> keep, to_hi, to_lo;
            keep.reserve(W.size());
            for (const PackedPart& w : W) {
                const double v = axis_val(w.p.tf.pos, axis);
                if (v >= d.sub_hi[axis]) {
                    to_hi.push_back(w);   // crossed the +axis face
                } else if (v < d.sub_lo[axis]) {
                    to_lo.push_back(w);   // crossed the -axis face
                } else {
                    keep.push_back(w);
                }
            }
            std::vector<PackedPart> got_lo = sendrecv_particles(
                d.cart, axis, +1, kMigrateTagBase + 4 * axis + 0, to_hi);
            std::vector<PackedPart> got_hi = sendrecv_particles(
                d.cart, axis, -1, kMigrateTagBase + 4 * axis + 2, to_lo);
            W.swap(keep);
            W.insert(W.end(), got_lo.begin(), got_lo.end());
            W.insert(W.end(), got_hi.begin(), got_hi.end());
        }
    }

    // Post-loop validation. The box check is defensive: every shipped
    // particle was in-box at this round's top check on its sender (positions
    // never change inside migrate_particles), so this should be
    // unreachable -- kept because it is one cheap scan guarding the most
    // dangerous failure (a particle dropped at a PROC_NULL edge). Then the
    // round cap: 8 rounds still left something displaced.
    check_in_box(d, W);
    for (const PackedPart& w : W) {
        if (!in_sub(d, w.p.tf.pos)) {
            die(d,
                "particle gid %d still outside its owner brick after 8 migration "
                "rounds: pos=(%.6f %.6f %.6f) sub=[%.3f %.3f]-[%.3f %.3f]-[%.3f %.3f] "
                "(dt too large for the brick size?)",
                w.gid, w.p.tf.pos.x, w.p.tf.pos.y, w.p.tf.pos.z,
                d.sub_lo[0], d.sub_hi[0], d.sub_lo[1], d.sub_hi[1],
                d.sub_lo[2], d.sub_hi[2]);
        }
    }

    // Publish the settled set back to the driver. Ownership moves but the
    // global set never changes in this phase (no insert/delete), so the
    // caller's expected total stays valid across the whole run.
    local.clear();
    gids.clear();
    local.reserve(W.size());
    gids.reserve(W.size());
    for (const PackedPart& w : W) {
        local.push_back(w.p);
        gids.push_back(w.gid);
    }
}

long long assert_global_count(const Decomp& d, long long nlocal, long long expected_total) {
    long long tot = 0;
    MPI_Allreduce(&nlocal, &tot, 1, MPI_LONG_LONG, MPI_SUM, d.cart);
    if (tot != expected_total) {
        if (d.rank == 0) {
            std::fprintf(stderr,
                         "zdem_mpi fatal: global particle count broken: "
                         "sum(nlocal)=%lld != N=%lld (particle lost or duplicated)\n",
                         tot, expected_total);
            std::fflush(stderr);
        }
        MPI_Abort(d.cart, 1);
    }
    return tot;
}
