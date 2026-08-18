#pragma once

#include "core/vec3.hpp"
#include "host/vtk_io.hpp"

using namespace zdem;

void integrate_particle(Particle& p, const Vec3& force, const Vec3& torque,
                        const Vec3& gravity, double dt);
