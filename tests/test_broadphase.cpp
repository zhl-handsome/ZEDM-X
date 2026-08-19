// tests/test_broadphase.cpp -- order-preservation unit test for the shared
// spatial-hash broadphase (physics/broadphase). Locks the contract the CPU
// and MPI drivers rely on for byte-level output stability:
//
//   1. The emitted pair LIST (set AND lexicographic order) equals an
//      O(N^2) i<j double loop over the combined [local; ghosts] sequence
//      using the drivers' verbatim bounding-sphere predicate.
//   2. Ghost semantics: first index always local, ghost partners carry the
//      combined index (n_local + gj), ghost-ghost pairs never appear.
//   3. Degenerates: N < 2 -> empty, all-zero radii -> empty (documented
//      behavior; positions must be distinct there -- duplicates would be
//      kept by the reference via dot == 0), ghosts = nullptr (CPU form).
//   4. Touching boundary: pairs at exactly dist == rsum are KEPT (the
//      predicate is a non-strict > test); constructed with dyadic
//      coordinates so the comparison is exact in floating point.
//
// Layouts are generated from a fixed-seed LCG (no random device): the test
// is deterministic run-to-run and machine-to-machine. Each layout is built
// exactly once in main() and reused, so the LCG consumption order is fixed.
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "host/vtk_io.hpp"
#include "physics/broadphase.hpp"

namespace {

int g_failures = 0;

void require(bool ok, const std::string& msg) {
    if (!ok) {
        std::cerr << "FAIL: " << msg << "\n";
        ++g_failures;
    }
}

// Fixed-seed LCG -> doubles in [0, 1). 53-bit resolution, no <random>
// device: identical sequences on every platform and every run.
unsigned long long g_lcg = 0x9E3779B97F4A7C15ULL;
double urand() {
    g_lcg = g_lcg * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((g_lcg >> 11) & 0x1FFFFFFFFFFFFFULL) * (1.0 / 9007199254740992.0);
}

Particle make_p(const Vec3& pos, double radius) {
    Particle p;   // mass/inertia/etc. are irrelevant to the broadphase
    p.tf.pos = pos;
    p.radius = radius;
    return p;
}

// O(N^2) reference: the drivers' original i<j double loop over the COMBINED
// array (local first, then ghosts). i only runs over the local range; the
// distance predicate is verbatim from the drivers (non-strict >).
std::vector<std::pair<int, int>> reference_pairs(const std::vector<Particle>& combined,
                                                 int n_local) {
    std::vector<std::pair<int, int>> ref;
    const int n = (int)combined.size();
    for (int i = 0; i < n_local; ++i) {
        for (int j = i + 1; j < n; ++j) {
            Vec3 dpos = combined[j].tf.pos - combined[i].tf.pos;
            double rsum = combined[i].radius + combined[j].radius;
            if (dot(dpos, dpos) > rsum * rsum) {
                continue;
            }
            ref.emplace_back(i, j);
        }
    }
    return ref;
}

// Run the broadphase on (local, ghosts) and compare against the reference
// ELEMENT-WISE: same pair count, same (i, j) at every position -- locks
// both the pair set and the lexicographic order. Also asserts the ghost
// semantics on every emitted pair (first local, no ghost-ghost).
void check_matches_reference(const std::vector<Particle>& local,
                             const std::vector<Particle>* ghosts,
                             const std::string& label) {
    std::vector<Particle> combined(local);
    if (ghosts) {
        combined.insert(combined.end(), ghosts->begin(), ghosts->end());
    }
    const int n_local = (int)local.size();
    const std::vector<std::pair<int, int>> ref = reference_pairs(combined, n_local);

    std::vector<std::pair<int, int>> got;
    broadphase_pairs(local, ghosts, got);

    require(got.size() == ref.size(),
            label + ": pair count " + std::to_string(got.size()) + " != reference " +
                std::to_string(ref.size()));
    const std::size_t common = got.size() < ref.size() ? got.size() : ref.size();
    for (std::size_t k = 0; k < common; ++k) {
        require(got[k] == ref[k],
                label + ": pair[" + std::to_string(k) + "] got (" +
                    std::to_string(got[k].first) + "," + std::to_string(got[k].second) +
                    ") != reference (" + std::to_string(ref[k].first) + "," +
                    std::to_string(ref[k].second) + ")");
        require(got[k].first < n_local,
                label + ": pair[" + std::to_string(k) + "] first index " +
                    std::to_string(got[k].first) + " >= n_local " +
                    std::to_string(n_local));
        require(!(got[k].first >= n_local && got[k].second >= n_local),
                label + ": ghost-ghost pair emitted at [" + std::to_string(k) + "]");
    }
}

// -------- layout generators (deterministic via the fixed-seed LCG) --------

// Two points at EXACTLY touching distance: radii 0.5 at x = -0.5 / +0.5 ->
// dist^2 == 1.0 == rsum^2, all dyadic so the comparison is exact.
std::vector<Particle> layout_touching_pair() {
    return {make_p(Vec3{-0.5, 0.0, 0.0}, 0.5), make_p(Vec3{0.5, 0.0, 0.0}, 0.5)};
}

// Small clustered set: mixed radii, negative coordinates, two coincident
// positions (dot == 0 -> the duplicate pair must be kept).
std::vector<Particle> layout_clustered5() {
    std::vector<Particle> ps;
    ps.push_back(make_p(Vec3{-1.0, -2.0, -0.5}, 0.4));
    ps.push_back(make_p(Vec3{-0.8, -1.9, -0.4}, 0.25));
    ps.push_back(make_p(Vec3{-0.5, -2.3, -0.2}, 0.6));
    ps.push_back(make_p(Vec3{-1.0, -2.0, -0.5}, 0.3));   // duplicate position
    ps.push_back(make_p(Vec3{+2.0, +2.0, +2.0}, 0.35));  // far outlier
    return ps;
}

// 40 particles on a jittered 5x4x2 lattice, single radius 0.35 (rsum 0.7):
// pitch 1.0 with +-0.3 jitter puts adjacent sites at distances 0.4..1.6,
// a genuine hit/miss mix that also straddles the 0.7-wide cell boundaries.
// Centered so roughly half the coordinates are negative.
std::vector<Particle> layout_uniform40() {
    std::vector<Particle> ps;
    for (int i = 0; i < 40; ++i) {
        const int ix = i % 5;
        const int iy = (i / 5) % 4;
        const int iz = i / 20;
        double x = (ix - 2.0) + (urand() - 0.5) * 0.6;
        double y = (iy - 1.5) + (urand() - 0.5) * 0.6;
        double z = (iz - 0.5) + (urand() - 0.5) * 0.6;
        ps.push_back(make_p(Vec3{x, y, z}, 0.35));
    }
    return ps;
}

// 40 particles in three separated blobs with mixed radii 0.1..0.6:
// dense intra-blob pairs, no inter-blob pairs.
std::vector<Particle> layout_blobs40() {
    std::vector<Particle> ps;
    const Vec3 centers[3] = {Vec3{-5.0, -5.0, 0.0}, Vec3{5.0, -5.0, 0.0},
                             Vec3{0.0, 6.0, 0.0}};
    for (int i = 0; i < 40; ++i) {
        const Vec3& c = centers[i % 3];
        double x = c.x + (urand() - 0.5) * 0.8;
        double y = c.y + (urand() - 0.5) * 0.8;
        double z = c.z + (urand() - 0.5) * 0.8;
        double r = 0.1 + 0.5 * urand();   // mixed radii
        ps.push_back(make_p(Vec3{x, y, z}, r));
    }
    return ps;
}

// Distinct positions, all radii zero -> documented-degenerate empty output.
std::vector<Particle> layout_zero_radii() {
    std::vector<Particle> ps;
    for (int i = 0; i < 5; ++i) {
        ps.push_back(make_p(Vec3{0.3 * i, 0.0, 0.0}, 0.0));
    }
    return ps;
}

}  // namespace

int main() {
    // Build every layout exactly once (fixes the LCG consumption order).
    const std::vector<Particle> touch2 = layout_touching_pair();
    const std::vector<Particle> cluster5 = layout_clustered5();
    const std::vector<Particle> uniform40 = layout_uniform40();
    const std::vector<Particle> blobs40 = layout_blobs40();
    const std::vector<Particle> zero_r = layout_zero_radii();

    // ---- 1 + 3: order preservation across layouts; CPU (ghosts=nullptr)
    // form, plus the N=0 / N=1 degenerates as layouts.
    {
        std::vector<std::pair<int, int>> got;

        std::vector<Particle> none;
        check_matches_reference(none, nullptr, "cpu/empty0");
        std::vector<Particle> single = {make_p(Vec3{0.0, 0.0, 0.0}, 0.5)};
        check_matches_reference(single, nullptr, "cpu/single1");
        check_matches_reference(touch2, nullptr, "cpu/touch2");
        check_matches_reference(cluster5, nullptr, "cpu/cluster5");
        check_matches_reference(uniform40, nullptr, "cpu/uniform40");
        check_matches_reference(blobs40, nullptr, "cpu/blobs40");

        // N < 2 and all-zero radii must yield empty output exactly.
        broadphase_pairs(none, nullptr, got);
        require(got.empty(), "deg: N=0 not empty");
        broadphase_pairs(single, nullptr, got);
        require(got.empty(), "deg: N=1 not empty");
        std::vector<Particle> one_ghost = {make_p(Vec3{1.0, 2.0, 3.0}, 0.4)};
        broadphase_pairs(none, &one_ghost, got);
        require(got.empty(), "deg: empty local + 1 ghost not empty");
        broadphase_pairs(zero_r, nullptr, got);
        require(got.empty(), "deg: all-zero radii not empty");
        // Ghost variant of the documented degenerate (ghost radii zero too).
        broadphase_pairs(zero_r, &zero_r, got);
        require(got.empty(), "deg: all-zero radii with ghosts not empty");
    }

    // ---- 4: touching boundary -- dist == rsum exactly must be KEPT
    // (non-strict > in the predicate); dyadic coords keep the test exact.
    {
        std::vector<std::pair<int, int>> got;
        broadphase_pairs(touch2, nullptr, got);
        require(got.size() == 1 && got[0].first == 0 && got[0].second == 1,
                "touch: exact-distance pair dropped (got " +
                    std::to_string(got.size()) + " pairs)");
    }

    // ---- 2: ghost semantics -- combined indices, local-first, no gg pairs.
    // Split each layout into local + ghosts at several cut points and let
    // the reference combined loop arbitrate the pair set/order; covers
    // local-local, local-ghost and (never) ghost-ghost combinations.
    {
        struct Case {
            const char* name;
            const std::vector<Particle>* ps;
            int n_local;
        };
        const Case cases[] = {
            {"touch2/1+1", &touch2, 1},
            {"cluster5/2+3", &cluster5, 2},
            {"cluster5/5+0", &cluster5, 5},   // empty halo (non-null pointer)
            {"cluster5/0+5", &cluster5, 0},   // empty local, ghosts only
            {"uniform40/25+15", &uniform40, 25},
            {"blobs40/13+27", &blobs40, 13},
        };
        for (const Case& c : cases) {
            std::vector<Particle> local(c.ps->begin(), c.ps->begin() + c.n_local);
            std::vector<Particle> ghosts(c.ps->begin() + c.n_local, c.ps->end());
            check_matches_reference(local, &ghosts, std::string("ghost/") + c.name);
        }
    }

    if (g_failures != 0) {
        std::cerr << "test_broadphase: " << g_failures << " failure(s)\n";
        return 1;
    }
    std::cout << "test_broadphase passed\n";
    return 0;
}
