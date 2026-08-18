#include "mpi/ghost.hpp"

// Task 4 placeholder: single-rank (V0) runs own every particle, so there is
// no halo to exchange. Clearing `out` keeps the interface contract ("after
// the call, out holds the current halo") without any MPI calls. Task 5
// replaces this with the neighbor isend/irecv exchange.
void exchange_ghosts(const Decomp& d, const std::vector<Particle>& local,
                     const std::vector<int>& gids, GhostLayer& out) {
    (void)d;
    (void)local;
    (void)gids;
    out.particles.clear();
    out.gids.clear();
}
