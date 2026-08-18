#pragma once

#include <vector>

#include "core/vec3.hpp"
#include "geometry/mesh.hpp"
#include "host/sim_build.hpp"
#include "host/vtk_io.hpp"

using namespace zdem;

// One particle x all walls: v8 grouping / penetration / penalty (extracted
// from main.cpp, structure unchanged). f/t are output accumulators: the CPU
// caller passes forces[i]/torques[i] directly and this function accumulates
// per vertex in the original order (identical FP order to the pre-extraction
// forces[i] += ... sequence).
// Returns the number of wall groups this particle contacted (contact_counts
// semantics, same as CPU). The three trailing parameters are optional
// diagnostics (the contact_debug stdout line; dbg_contacts is the step-wide
// contacts snapshot at call time, used only for the "contacts < 3" print gate).
int wall_contact_particle(const Particle& p, const Mesh& pmesh,
                          const std::vector<Wall>& walls,
                          double tangential_damping,
                          Vec3& f, Vec3& t,
                          bool contact_debug = false,
                          int dbg_particle = 0,
                          int dbg_contacts = 0);
