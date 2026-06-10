#include "geometry/stl_io.hpp"

#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>

namespace zdem {

bool load_stl(const std::string& path, std::vector<Triangle>& tris) {
    tris.clear();
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }

    in.seekg(0, std::ios::end);
    std::streampos size = in.tellg();
    in.seekg(0, std::ios::beg);
    if (size < 84) {
        return false;
    }

    std::array<char, 80> header{};
    in.read(header.data(), header.size());
    uint32_t n_tri = 0;
    in.read(reinterpret_cast<char*>(&n_tri), sizeof(uint32_t));

    std::streampos expected = 84 + static_cast<std::streampos>(n_tri) * 50;
    if (expected == size) {
        tris.reserve(n_tri);
        for (uint32_t i = 0; i < n_tri; ++i) {
            float data[12];
            uint16_t attr = 0;
            in.read(reinterpret_cast<char*>(data), sizeof(float) * 12);
            in.read(reinterpret_cast<char*>(&attr), sizeof(uint16_t));
            Triangle tri{
                Vec3{data[3], data[4], data[5]},
                Vec3{data[6], data[7], data[8]},
                Vec3{data[9], data[10], data[11]}
            };
            tris.push_back(tri);
        }
        return true;
    }

    in.close();
    std::ifstream ia(path);
    if (!ia) {
        return false;
    }
    std::string line;
    std::vector<Vec3> vbuf;
    while (std::getline(ia, line)) {
        std::string s = line;
        for (auto& c : s) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        if (s.find("vertex") != std::string::npos) {
            std::istringstream iss(s);
            std::string tag;
            double x, y, z;
            if (iss >> tag >> x >> y >> z) {
                vbuf.emplace_back(x, y, z);
            }
        }
        if (s.find("endloop") != std::string::npos && vbuf.size() >= 3) {
            Triangle tri{vbuf[vbuf.size() - 3], vbuf[vbuf.size() - 2], vbuf[vbuf.size() - 1]};
            tris.push_back(tri);
        }
    }
    return !tris.empty();
}

}  // namespace zdem
