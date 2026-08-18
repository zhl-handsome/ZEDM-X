#pragma once

#include <vector>

#include "host/vtk_io.hpp"
#include "mpi/decomp.hpp"

// Halo particles owned by neighboring bricks, consumed by the Newton-off
// pp pass (each side keeps only its own force half).
struct GhostLayer {
    std::vector<Particle> particles;
    std::vector<int> gids;
};

// Fill `out` with the current halo of the neighboring bricks: every
// non-locally-owned particle within ghost_depth of this rank's brick
// boundary. Chained three-axis exchange (LAMMPS comm_brick style): ghosts
// received on axis k join the send set for axis k+1, so diagonal ghosts
// arrive in two hops. Ownership does not move (migration is a later task),
// so `local`/`gids` stay this rank's initial brick contents.
void exchange_ghosts(const Decomp& d, const std::vector<Particle>& local,
                     const std::vector<int>& gids, GhostLayer& out);
