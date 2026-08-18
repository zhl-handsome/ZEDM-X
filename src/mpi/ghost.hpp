#pragma once

#include <vector>

#include "host/vtk_io.hpp"
#include "mpi/decomp.hpp"

// Halo particles owned by neighboring bricks, consumed by the Newton-off
// pp pass (each side keeps only its own force half). Task 4 placeholder:
// single-rank (V0) runs own every particle, so exchange_ghosts() only
// clears `out`; Task 5 implements the real neighbor exchange.
struct GhostLayer {
    std::vector<Particle> particles;
    std::vector<int> gids;
};

// Fill `out` with the current halo of the neighboring bricks (particles
// within ghost_depth of this rank's brick boundary). Task 4: no-op that
// clears `out` (n=1 has no neighbors); no MPI calls are made.
void exchange_ghosts(const Decomp& d, const std::vector<Particle>& local,
                     const std::vector<int>& gids, GhostLayer& out);
