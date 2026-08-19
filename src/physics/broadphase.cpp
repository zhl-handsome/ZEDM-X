// src/physics/broadphase.cpp -- uniform-grid spatial hash for bounding-sphere
// candidate pairs (O(N^2) -> O(N) for sparse scenes). Shared by the CPU and
// MPI drivers; see broadphase.hpp for the pair-set / order contract.
//
// Per-step allocations (hash map, cell vectors, out) match the drivers'
// existing per-step world-tris caching style; keep it simple until a profile
// says otherwise.
#include "physics/broadphase.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <unordered_map>

namespace {

struct CellKey {
    long long c[3];
    bool operator==(const CellKey& o) const {
        return c[0] == o.c[0] && c[1] == o.c[1] && c[2] == o.c[2];
    }
};

struct CellKeyHash {
    std::size_t operator()(const CellKey& k) const {
        // 64-bit FNV-1a over the three cell coords. std::size_t is 64-bit on
        // the supported platforms, and unordered_map re-checks equality with
        // CellKey::operator==, so hash collisions cannot alter the pair set.
        unsigned long long h = 1469598103934665603ULL;
        for (int a = 0; a < 3; ++a) {
            h ^= (unsigned long long)k.c[a];
            h *= 1099511628211ULL;
        }
        return (std::size_t)h;
    }
};

inline double coord_of(const Particle& p, int axis) {
    return axis == 0 ? p.tf.pos.x : (axis == 1 ? p.tf.pos.y : p.tf.pos.z);
}

}  // namespace

void broadphase_pairs(const std::vector<Particle>& ps, int n_local,
                      std::vector<std::pair<int, int>>& out) {
    out.clear();
    const int n = (int)ps.size();
    if (n < 2) {
        return;
    }

    // Cell edge. Coverage needs cell >= r_i + r_j for every pair, i.e.
    // cell >= 2*max_radius. The (1 + 2^-30) margin makes that argument
    // airtight in floating point: cell indices are
    // floor((pos - lo) / cell), and each computed quotient carries at most
    // ~2 ulps of rounding (subtraction + division) versus the exact one, so
    // two centers less than 2*max_radius apart get indices differing by at
    // most 1 per axis as long as the domain spans well under ~1e9 cells
    // (this repo: 256 particles over ~14 cells). With a bare 2*max_radius
    // cell there is a sub-ulp-wide window where the quotients of a touching
    // pair could straddle two integer boundaries each and the
    // 27-neighborhood scan would miss it.
    double max_radius = 0.0;
    for (int i = 0; i < n; ++i) {
        max_radius = std::max(max_radius, ps[i].radius);
    }
    const double cell = 2.0 * max_radius * (1.0 + 1.0 / 1073741824.0);
    if (!(cell > 0.0)) {
        // Degenerate: all radii zero. Mesh bounding-sphere radii are always
        // positive in practice; an all-zero set has no bounding-sphere pair.
        return;
    }

    // Grid reference corner = min position over ps. Indexing against min
    // keeps every (coord - lo) >= 0 and std::floor behaves plainly; casting
    // negative values to long long would truncate toward zero instead.
    double lo[3];
    for (int a = 0; a < 3; ++a) {
        lo[a] = coord_of(ps[0], a);
    }
    for (int i = 1; i < n; ++i) {
        for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], coord_of(ps[i], a));
        }
    }

    // Bucket every particle into its cell; cell_of caches each particle's
    // own key so the neighbor scan below does not recompute it.
    std::unordered_map<CellKey, std::vector<int>, CellKeyHash> grid;
    grid.reserve((std::size_t)n * 2);
    std::vector<CellKey> cell_of((std::size_t)n);
    for (int i = 0; i < n; ++i) {
        CellKey k;
        for (int a = 0; a < 3; ++a) {
            k.c[a] = (long long)std::floor((coord_of(ps[i], a) - lo[a]) / cell);
        }
        grid[k].push_back(i);
        cell_of[(std::size_t)i] = k;
    }

    // 27-neighborhood scan. Each overlapping pair is emitted exactly once:
    // only while scanning its smaller index (j > i check) and each particle
    // occupies exactly one cell. The distance test below is bit-identical to
    // the drivers' former in-loop precheck -- same operands, same expression
    // order -- so the collected set equals the O(N^2) i<j double loop's set
    // exactly (the precheck kept in the MPI consumer is a no-op copy of it).
    for (int i = 0; i < n_local; ++i) {
        const CellKey& ki = cell_of[(std::size_t)i];
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    CellKey k;
                    k.c[0] = ki.c[0] + dx;
                    k.c[1] = ki.c[1] + dy;
                    k.c[2] = ki.c[2] + dz;
                    auto it = grid.find(k);
                    if (it == grid.end()) {
                        continue;
                    }
                    for (int j : it->second) {
                        if (j <= i) {
                            continue;
                        }
                        Vec3 dpos = ps[(std::size_t)j].tf.pos - ps[(std::size_t)i].tf.pos;
                        double rsum = ps[(std::size_t)i].radius + ps[(std::size_t)j].radius;
                        if (dot(dpos, dpos) > rsum * rsum) {
                            continue;
                        }
                        out.emplace_back(i, j);
                    }
                }
            }
        }
    }

    // Lexicographic (i, j) order == the i<j double-loop enumeration order.
    std::sort(out.begin(), out.end());
}
