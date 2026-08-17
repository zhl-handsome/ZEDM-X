#include "host/config_io.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

bool parse_config_file(const std::string& path, SimConfig& cfg) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Failed to open config: " << path << "\n";
        return false;
    }
    std::string line;
    bool in_particle = false;
    bool in_wall = false;
    ParticleInit current_particle;
    WallInit current_wall;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::string cleaned;
        cleaned.reserve(line.size());
        for (char c : line) {
            if (c == '#') {
                break;
            }
            if (c == '=') {
                cleaned.push_back(' ');
            } else {
                cleaned.push_back(c);
            }
        }
        std::istringstream iss(cleaned);
        std::string key;
        if (!(iss >> key)) {
            continue;
        }
        if (key == "particle") {
            if (in_particle) {
                cfg.particle_inits.push_back(current_particle);
            }
            if (in_wall) {
                cfg.wall_inits.push_back(current_wall);
                in_wall = false;
            }
            current_particle = ParticleInit{};
            in_particle = true;
            continue;
        }
        if (key == "end_particle" || key == "particle_end") {
            if (in_particle) {
                cfg.particle_inits.push_back(current_particle);
                in_particle = false;
            }
            continue;
        }
        if (key == "wall") {
            if (in_particle) {
                cfg.particle_inits.push_back(current_particle);
                in_particle = false;
            }
            if (in_wall) {
                cfg.wall_inits.push_back(current_wall);
            }
            current_wall = WallInit{};
            in_wall = true;
            continue;
        }
        if (key == "end_wall" || key == "wall_end") {
            if (in_wall) {
                cfg.wall_inits.push_back(current_wall);
                in_wall = false;
            }
            continue;
        }

        if (in_particle) {
            if (key == "stl") {
                iss >> current_particle.stl_path;
            } else if (key == "pos") {
                iss >> current_particle.pos.x >> current_particle.pos.y >> current_particle.pos.z;
            } else if (key == "vel") {
                iss >> current_particle.vel.x >> current_particle.vel.y >> current_particle.vel.z;
            } else if (key == "omega") {
                iss >> current_particle.omega.x >> current_particle.omega.y >> current_particle.omega.z;
            } else if (key == "quat") {
                iss >> current_particle.rot.w >> current_particle.rot.x >> current_particle.rot.y >> current_particle.rot.z;
            } else if (key == "scale") {
                iss >> current_particle.scale;
            } else if (key == "density") {
                iss >> current_particle.density;
            } else if (key == "young") {
                iss >> current_particle.young;
            } else if (key == "poisson") {
                iss >> current_particle.poisson;
            } else if (key == "mu") {
                iss >> current_particle.mu;
            } else if (key == "restitution") {
                iss >> current_particle.restitution;
            }
        } else if (in_wall) {
            if (key == "stl") {
                iss >> current_wall.stl_path;
            } else if (key == "pos") {
                iss >> current_wall.pos.x >> current_wall.pos.y >> current_wall.pos.z;
            } else if (key == "quat") {
                iss >> current_wall.rot.w >> current_wall.rot.x >> current_wall.rot.y >> current_wall.rot.z;
            } else if (key == "scale") {
                iss >> current_wall.scale;
            } else if (key == "mu") {
                iss >> current_wall.mu;
            } else if (key == "restitution") {
                iss >> current_wall.restitution;
            }
        } else if (key == "stl") {
            iss >> cfg.stl_path;
        } else if (key == "n") {
            iss >> cfg.n;
        } else if (key == "steps") {
            iss >> cfg.steps;
        } else if (key == "dt") {
            iss >> cfg.dt;
        } else if (key == "spacing") {
            iss >> cfg.spacing;
        } else if (key == "split_contacts") {
            int v = 1;
            iss >> v;
            cfg.split_contacts = (v != 0);
        } else if (key == "contact_debug") {
            int v = 0;
            iss >> v;
            cfg.contact_debug = (v != 0);
        } else if (key == "gravity") {
            iss >> cfg.gravity.x >> cfg.gravity.y >> cfg.gravity.z;
        } else if (key == "center_mesh") {
            int v = 1;
            iss >> v;
            cfg.center_mesh = (v != 0);
        } else if (key == "vtk_prefix") {
            iss >> cfg.vtk_prefix;
        } else if (key == "output_interval") {
            iss >> cfg.output_interval;
        } else if (key == "output_dir") {
            iss >> cfg.output_dir;
        } else if (key == "v0") {
            iss >> cfg.v0.x >> cfg.v0.y >> cfg.v0.z;
        } else if (key == "v1") {
            iss >> cfg.v1.x >> cfg.v1.y >> cfg.v1.z;
        } else if (key == "tangential_damping") {
            iss >> cfg.tangential_damping;
        }
    }
    if (in_particle) {
        cfg.particle_inits.push_back(current_particle);
    }
    if (in_wall) {
        cfg.wall_inits.push_back(current_wall);
    }
    if (cfg.output_interval < 1) {
        cfg.output_interval = 1;
    }
    if (!cfg.stl_path.empty()) {
        return true;
    }
    return !cfg.particle_inits.empty();
}
