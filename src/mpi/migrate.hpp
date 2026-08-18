#pragma once

#include <vector>

#include "host/sim_build.hpp"
#include "mpi/decomp.hpp"

// Centroid migration: after every step, particles whose centroid left this
// rank's brick change ownership (LAMMPS comm_brick style, chained over the
// three axes like exchange_ghosts). Per-axis decision uses the same floor
// convention as in_sub/owner_rank (bricks are left-closed right-open):
//   pos[axis] >= sub_hi[axis] -> +axis face neighbor
//   pos[axis] <  sub_lo[axis] -> -axis face neighbor
// Only the displaced locals are shipped (via mpi/comm_util.hpp's
// count-then-data Sendrecv, disjoint tag range from the ghost exchange);
// received particles join the working set and are re-tested on the next
// axis, so a corner crossing settles within one round. A particle that
// crossed several bricks in one step needs another full round; rounds
// repeat up to 8 (the round count is agreed collectively -- an Allreduce
// keeps every rank in the loop while ANY rank still holds a displaced
// particle, because that particle may be on its way here).
//
// Fatal paths (the rank holding the offending particle logs, then
// MPI_Abort -- the message would be lost if only rank 0 could print, since
// the particle lives on exactly one rank): a particle outside the global
// box [box_lo, box_hi) -- checked BEFORE any send, since at a box edge the
// face neighbor is MPI_PROC_NULL and the send would silently drop it -- or
// a particle still outside its owner brick after 8 rounds.
//
// gid uniqueness survives migration structurally: a sent particle leaves
// the sender's working set immediately (at most one send per axis per
// round), so no rank ever holds two copies. assert_global_count guards the
// global total as a standing invariant (this phase never inserts or
// deletes particles).
void migrate_particles(const Decomp& d, std::vector<Particle>& local, std::vector<int>& gids);

// Allreduce SUM of nlocal over d.cart; a total != expected_total means a
// particle was lost or duplicated -- rank 0 logs, MPI_Abort. Returns the
// global count so callers can reuse it (e.g. the startup log line).
// NOTE: the task sketch had assert_global_count(d, expected_total); the
// local size has to enter the reduction somehow, so it is an explicit
// parameter.
long long assert_global_count(const Decomp& d, long long nlocal, long long expected_total);
