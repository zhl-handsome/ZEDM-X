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

}  // namespace

void broadphase_pairs(const std::vector<Particle>& local,
                      const std::vector<Particle>* ghosts,
                      std::vector<std::pair<int, int>>& out) {
    out.clear();
    const int n_local = (int)local.size();
    const int n = n_local + (ghosts ? (int)ghosts->size() : 0);
    if (n < 2) {
        return;
    }
    // Combined-index view: j < n_local addresses local[j], else the ghost
    // (*ghosts)[j - n_local]. Iterating 0..n below therefore enumerates the
    // SAME combined sequence the drivers' former [local; ghosts] copy had.
    const auto p_at = [&](int j) -> const Particle& {
        return j < n_local ? local[(std::size_t)j]
                           : (*ghosts)[(std::size_t)(j - n_local)];
    };

    // Cell edge. Coverage needs cell >= r_i + r_j for every pair, i.e.
    // cell >= 2*max_radius. The (1 + 2^-30) margin makes that argument
    // airtight in floating point: cell indices are
    // floor((pos - lo) / cell), and each computed quotient q carries at
    // most ~2 ulps of rounding (subtraction + division) versus the exact
    // one -- an absolute error of about |q|*2^-51 per quotient -- so two
    // centers less than 2*max_radius apart get indices differing by at
    // most 1 per axis as long as the domain stays under ~2e6 cells per
    // axis (|q| < 2^21 keeps 2*|q|*2^-52 below the 2^-30 slack; this
    // repo: 256 particles over ~14 cells). Beyond that bound pairs could
    // be silently missed. With a bare 2*max_radius cell there is also a
    // sub-ulp-wide window where the quotients of a touching pair could
    // straddle two integer boundaries each and the 27-neighborhood scan
    // would miss it.
    double max_radius = 0.0;
    for (int i = 0; i < n; ++i) {
        max_radius = std::max(max_radius, p_at(i).radius);
    }
    const double cell = 2.0 * max_radius * (1.0 + 1.0 / 1073741824.0);
    if (!(cell > 0.0)) {
        // Degenerate: all radii zero. Mesh bounding-sphere radii are always
        // positive in practice; an all-zero set has no bounding-sphere pair.
        return;
    }

    // Grid reference corner = min position over the combined set. Indexing
    // against min keeps every (coord - lo) >= 0 and std::floor behaves
    // plainly; casting negative values to long long would truncate toward
    // zero instead.
    double lo[3];
    for (int a = 0; a < 3; ++a) {
        lo[a] = a == 0 ? p_at(0).tf.pos.x : (a == 1 ? p_at(0).tf.pos.y
                                                   : p_at(0).tf.pos.z);
    }
    for (int i = 1; i < n; ++i) {
        const Vec3& p = p_at(i).tf.pos;
        lo[0] = std::min(lo[0], p.x);
        lo[1] = std::min(lo[1], p.y);
        lo[2] = std::min(lo[2], p.z);
    }

    // Bucket every particle into its cell; cell_of caches each particle's
    // own key so the neighbor scan below does not recompute it. Local
    // particles are inserted first (combined indices 0..n_local), then the
    // ghosts -- the same sequence the former combined array carried.
    std::unordered_map<CellKey, std::vector<int>, CellKeyHash> grid;
    grid.reserve((std::size_t)n * 2);
    std::vector<CellKey> cell_of((std::size_t)n);
    for (int i = 0; i < n; ++i) {
        const Vec3& pos = p_at(i).tf.pos;
        CellKey k;
        k.c[0] = (long long)std::floor((pos.x - lo[0]) / cell);
        k.c[1] = (long long)std::floor((pos.y - lo[1]) / cell);
        k.c[2] = (long long)std::floor((pos.z - lo[2]) / cell);
        grid[k].push_back(i);
        cell_of[(std::size_t)i] = k;
    }

    // 27-neighborhood scan. Each overlapping pair is emitted exactly once:
    // only while scanning its smaller index (j > i check) and each particle
    // occupies exactly one cell. The distance test below is bit-identical to
    // the drivers' former in-loop precheck -- same operands, same expression
    // order -- so the collected set equals the O(N^2) i<j double loop's set
    // exactly (the precheck kept in the MPI consumer is a no-op copy of it).
    // The loop covers only local i, so ghost-ghost pairs (both indices >=
    // n_local) are structurally absent, as in the drivers.
    for (int i = 0; i < n_local; ++i) {
        const CellKey& ki = cell_of[(std::size_t)i];
        const Vec3& pi = p_at(i).tf.pos;
        const double ri = p_at(i).radius;
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
                        Vec3 dpos = p_at(j).tf.pos - pi;
                        double rsum = ri + p_at(j).radius;
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
