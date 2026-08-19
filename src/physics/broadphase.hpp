#pragma once

#include <utility>
#include <vector>

#include "host/vtk_io.hpp"

// Bounding-sphere overlap candidate pairs over [local; optional ghosts] (the
// MPI driver passes its halo as ghosts; the CPU driver passes nullptr). No
// combined array is built: particle access resolves the combined index j as
// j < local.size() ? local[j] : (*ghosts)[j - local.size()]. Cell =
// 2*max(radius over both arrays) so any overlapping pair sits in adjacent
// cells; 27-neighborhood scan with j>i dedup, then LEXICOGRAPHIC SORT so the
// output order equals the O(N^2) i<j double loop over the combined array
// EXACTLY -- downstream force accumulation order and bit-level output are
// unchanged. Cell insertion visits local then ghosts; the outer emission
// loop runs only local i (first index always < local.size()), so ghost-ghost
// pairs are structurally never generated (matches the drivers).
void broadphase_pairs(const std::vector<Particle>& local,
                      const std::vector<Particle>* ghosts,
                      std::vector<std::pair<int, int>>& out);
