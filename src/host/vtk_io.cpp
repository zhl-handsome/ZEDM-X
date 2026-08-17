#include "host/vtk_io.hpp"

#include <fstream>
#include <iostream>

#include "core/quat.hpp"

void write_vtk_particles(const std::string& path,
                         const std::vector<Mesh>& meshes,
                         const std::vector<Particle>& particles,
                         const std::vector<Vec3>& forces,
                         const std::vector<Vec3>& torques,
                         const std::vector<int>& contact_counts) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to write VTK: " << path << "\n";
        return;
    }
    out << "# vtk DataFile Version 3.0\n";
    out << "zdem particles\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    std::size_t total_tris = 0;
    for (const auto& p : particles) {
        total_tris += meshes[p.mesh_index].tris.size();
    }
    const std::size_t total_points = total_tris * 3;
    out << "POINTS " << total_points << " float\n";
    for (const auto& p : particles) {
        const auto& m = meshes[p.mesh_index];
        for (const auto& tri : m.tris) {
            Vec3 v0 = quat_rotate(p.tf.rot, tri[0]) + p.tf.pos;
            Vec3 v1 = quat_rotate(p.tf.rot, tri[1]) + p.tf.pos;
            Vec3 v2 = quat_rotate(p.tf.rot, tri[2]) + p.tf.pos;
            out << static_cast<float>(v0.x) << " " << static_cast<float>(v0.y) << " " << static_cast<float>(v0.z) << "\n";
            out << static_cast<float>(v1.x) << " " << static_cast<float>(v1.y) << " " << static_cast<float>(v1.z) << "\n";
            out << static_cast<float>(v2.x) << " " << static_cast<float>(v2.y) << " " << static_cast<float>(v2.z) << "\n";
        }
    }

    out << "POLYGONS " << total_tris << " " << total_tris * 4 << "\n";
    for (std::size_t i = 0; i < total_tris; ++i) {
        std::size_t base = i * 3;
        out << "3 " << base << " " << base + 1 << " " << base + 2 << "\n";
    }

    out << "CELL_DATA " << total_tris << "\n";

    out << "SCALARS id int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        std::size_t tris_per_particle = meshes[particles[i].mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << i << "\n";
        }
    }

    out << "SCALARS mass float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.mass) << "\n";
        }
    }

    out << "SCALARS radius float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.radius) << "\n";
        }
    }

    out << "SCALARS contact_count int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < contact_counts.size(); ++i) {
        std::size_t tris_per_particle = meshes[particles[i].mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << contact_counts[i] << "\n";
        }
    }

    out << "VECTORS velocity float\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.vel.x) << " "
                << static_cast<float>(p.vel.y) << " "
                << static_cast<float>(p.vel.z) << "\n";
        }
    }

    out << "VECTORS omega float\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.omega.x) << " "
                << static_cast<float>(p.omega.y) << " "
                << static_cast<float>(p.omega.z) << "\n";
        }
    }

    out << "VECTORS force float\n";
    for (std::size_t i = 0; i < particles.size(); ++i) {
        std::size_t tris_per_particle = meshes[particles[i].mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            const Vec3& f = forces[i];
            out << static_cast<float>(f.x) << " "
                << static_cast<float>(f.y) << " "
                << static_cast<float>(f.z) << "\n";
        }
    }

    out << "VECTORS torque float\n";
    for (std::size_t i = 0; i < particles.size(); ++i) {
        std::size_t tris_per_particle = meshes[particles[i].mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            const Vec3& tt = torques[i];
            out << static_cast<float>(tt.x) << " "
                << static_cast<float>(tt.y) << " "
                << static_cast<float>(tt.z) << "\n";
        }
    }

    out << "FIELD FieldData 1\n";
    out << "orientation 4 " << total_tris << " float\n";
    for (const auto& p : particles) {
        std::size_t tris_per_particle = meshes[p.mesh_index].tris.size();
        for (std::size_t t = 0; t < tris_per_particle; ++t) {
            out << static_cast<float>(p.tf.rot.w) << " "
                << static_cast<float>(p.tf.rot.x) << " "
                << static_cast<float>(p.tf.rot.y) << " "
                << static_cast<float>(p.tf.rot.z) << "\n";
        }
    }
}
