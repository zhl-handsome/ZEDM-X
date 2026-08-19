#pragma once

#include <utility>
#include <vector>

#include "host/vtk_io.hpp"

// Bounding-sphere overlap candidate pairs over ps[0..n) (combined local+ghost
// array in MPI). Cell = 2*max(radius) so any overlapping pair sits in adjacent
// cells; 27-neighborhood scan with j>i dedup, then LEXICOGRAPHIC SORT so the
// output order equals the O(N^2) i<j double loop EXACTLY -- downstream force
// accumulation order and bit-level output are unchanged. Pairs with first
// index >= n_local are dropped (ghost-ghost never computed, matches drivers).
void broadphase_pairs(const std::vector<Particle>& ps, int n_local,
                      std::vector<std::pair<int, int>>& out);
