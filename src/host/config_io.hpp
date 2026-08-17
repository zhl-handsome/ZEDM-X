#pragma once

#include <string>
#include <vector>

#include "core/quat.hpp"
#include "core/vec3.hpp"

using namespace zdem;

struct ParticleInit {
    std::string stl_path;
    Vec3 pos{0.0, 0.0, 0.0};
    Vec3 vel{0.0, 0.0, 0.0};
    Vec3 omega{0.0, 0.0, 0.0};
    Quat rot{1.0, 0.0, 0.0, 0.0};
    double scale = 0.0;
    double density = 1.0;
    double young = 1e7;
    double poisson = 0.25;
    double mu = 0.5;
    double restitution = 0.5;
};

struct WallInit {
    std::string stl_path;
    Vec3 pos{0.0, 0.0, 0.0};
    Quat rot{1.0, 0.0, 0.0, 0.0};
    double scale = 1.0;
    double mu = 0.5;
    double restitution = 0.5;
};

struct SimConfig {
    int n = 2;
    int steps = 1;
    double dt = 1e-4;
    double spacing = 2.5;
    Vec3 gravity{0.0, -9.81, 0.0};
    bool split_contacts = true;
    bool center_mesh = true;
    bool contact_debug = false;
    double tangential_damping = 1.0;  // 切向阻尼系数: ct = tangential_damping * sqrt(kt*m_eff)
    std::string stl_path;
    std::string vtk_prefix = "particles";
    int output_interval = 1;
    std::string output_dir = "output";
    Vec3 v0{0.0, 0.0, 0.0};
    Vec3 v1{0.0, 0.0, 0.0};
    std::vector<ParticleInit> particle_inits;
    std::vector<WallInit> wall_inits;
};

bool parse_config_file(const std::string& path, SimConfig& cfg);
