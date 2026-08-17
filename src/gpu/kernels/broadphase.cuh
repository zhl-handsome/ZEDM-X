// src/gpu/kernels/broadphase.cuh -- spatial-hash broad phase (GPU Task 6)
//
// Produces the candidate particle-pair list (bounding-sphere overlap) that
// Task 7's particle-particle kernel will consume. Detection only: no forces,
// no state mutation -- safe to run anywhere between clear_forces and
// integrate; the driver keeps it right before integrate (Task 7's slot).
//
// Pipeline (per step, all buffers persistent -- allocated once by
// alloc_broadphase(), no per-step cudaMalloc):
//   1. compute_cell_ids_kernel : cell = floor((pos - origin)/cell) per
//      particle -> 64-bit linear cell id
//   2. cub::DeviceRadixSort::SortPairs(keys = cell id, values = particle
//      index) -> particles grouped by cell, sorted within
//   3. candidate_pairs_kernel, COUNT pass (pairs == nullptr): thread i
//      re-derives its own cell, scans the 27 neighbor cells (27 is
//      exhaustive: cell edge = 2*max(radius), so a candidate pair
//      dist < ri+rj <= cell has per-axis cell index difference <= 1), finds
//      each neighbor cell's particle run by binary search (lower_bound) on
//      the sorted keys, counts j > i with dist2 < (ri+rj)^2 -> counts[i]
//   4. cub::DeviceScan::ExclusiveSum over [counts, 0] (n+1 items, counts[n]
//      pinned 0 at alloc) -> offsets; offsets[n] = total pair count
//   5. candidate_pairs_kernel, WRITE pass (counts == nullptr): the same
//      kernel/enumeration (identical code path -> identical order), writes
//      pair (i, j) at row offsets[i] + running local index into BP.d_pairs
//
// Hash: cell coords are relative to the UPLOAD-TIME origin, so falling
// particles produce negative coords. Linear 64-bit id
//     ((cx+K)<<42) | ((cy+K)<<21) | (cz+K),  K = 2^20
// Each axis field is clamped to [0, 2^21-1] (NOT masked): positions more
// than 2^20 cells (~200 km at cell = 0.2 m) outside the origin grid merge
// into the edge cells. Merging distinct cells into one run cannot drop a
// pair: particle A still scans B's true cell coords, whose clamped id
// resolves to the merged run containing B; extra run members are removed by
// the distance test. Radix sort on 64-bit keys is fine (top bit always 0).
//
// Capacity: BP.d_pairs holds 64*n pairs (worst-case cap). Dense same-size
// packs touch ~12 neighbors per particle, so 64 rows per particle covers
// every realistic overlap; a pathological cluster beyond that clamps the
// stored list and sets the sticky BP.d_overflow flag (logical count in
// offsets[n] stays exact). The driver warns on overflow when it reads the
// pair count back (dump/step log); Task 7 kernels must check capacity too.
//
// Host-side consumers: readback_pairs() copies the compact list to the host
// (--dump-pairs debug flag, step log pair count). Device consumers (Task 7):
// pairs live in BP.d_pairs, the per-step count in BP.d_offsets[n].
//
// CPU parity note: the candidate test is the CPU broad phase verbatim
// (src/main.cpp bounding-sphere skip): strict dist2 < (ri+rj)^2 on the
// OUTSIDE bounding-sphere radius (p.radius, not equiv_radius).
#pragma once
#include <cstdio>
#include <cuda_runtime.h>
#include <vector>

#include <cub/cub.cuh>  // DeviceRadixSort, DeviceScan (ships with CUDA 12.8)

#include "gpu/cuda_check.hpp"
#include "gpu/dem_state.cuh"           // GpuSim, BroadPhaseWorkspace
#include "gpu/kernels/integrate.cuh"   // real, r_floor

namespace {

constexpr int kBPThreads = 128;

// Cell coordinate of one scalar: floor((x - origin)/cell). floor (not
// truncation) so negative coords round the right way; the SAME expression
// runs in every kernel that derives cells from positions, so a value on a
// cell boundary always resolves to the same cell regardless of which kernel
// evaluates it.
__device__ inline long long bp_cell_index(real x, real o, real cell) {
    return static_cast<long long>(r_floor((x - o) / cell));
}

// 64-bit linear id from three (possibly negative) cell coords; see header
// comment for the offset/clamp rationale.
__device__ inline unsigned long long bp_linear_id(long long cx, long long cy,
                                                  long long cz) {
    const long long K = 1ll << 20;
    const long long M = (1ll << 21) - 1;
    long long fx = cx + K; if (fx < 0) fx = 0; if (fx > M) fx = M;
    long long fy = cy + K; if (fy < 0) fy = 0; if (fy > M) fy = M;
    long long fz = cz + K; if (fz < 0) fz = 0; if (fz > M) fz = M;
    return (static_cast<unsigned long long>(fx) << 42) |
           (static_cast<unsigned long long>(fy) << 21) |
           static_cast<unsigned long long>(fz);
}

// First index in the sorted key array with keys[k] >= id (std::lower_bound
// equivalent); if id is present this is the first element of its run.
__device__ inline int bp_lower_bound(const unsigned long long* keys, int n,
                                     unsigned long long id) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (keys[mid] < id) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

}  // namespace

// Step 1: one thread per particle -> unsorted cell ids.
__global__ void compute_cell_ids_kernel(const real* pos, int n, real ox,
                                        real oy, real oz, real cell,
                                        unsigned long long* keys) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    keys[i] = bp_linear_id(bp_cell_index(pos[3 * i + 0], ox, cell),
                           bp_cell_index(pos[3 * i + 1], oy, cell),
                           bp_cell_index(pos[3 * i + 2], oz, cell));
}

// COUNT pass (pairs == nullptr -> stores counts[i]) / WRITE pass
// (counts == nullptr -> stores pairs rows at offsets[i] + local running
// index). Both passes MUST enumerate in the same order -- they run the
// identical code path over the same sorted arrays, so the write slot
// offsets[i] + cnt addresses exactly the pair counted in the count pass.
__global__ void candidate_pairs_kernel(
    const real* pos, const real* radius,
    const unsigned long long* keys_sorted, const int* idx_sorted, int n,
    real ox, real oy, real oz, real cell,
    int* counts, const int* offsets, int* pairs, int capacity, int* overflow) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    const long long ci = bp_cell_index(pos[3 * i + 0], ox, cell);
    const long long cj = bp_cell_index(pos[3 * i + 1], oy, cell);
    const long long ck = bp_cell_index(pos[3 * i + 2], oz, cell);
    const real xi = pos[3 * i + 0];
    const real yi = pos[3 * i + 1];
    const real zi = pos[3 * i + 2];
    const real ri = radius[i];

    int cnt = 0;
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                const unsigned long long id =
                    bp_linear_id(ci + dx, cj + dy, ck + dz);
                int k = bp_lower_bound(keys_sorted, n, id);
                for (; k < n && keys_sorted[k] == id; ++k) {
                    const int j = idx_sorted[k];
                    if (j <= i) continue;  // each unordered pair once, from min(i,j)
                    const real ddx = pos[3 * j + 0] - xi;
                    const real ddy = pos[3 * j + 1] - yi;
                    const real ddz = pos[3 * j + 2] - zi;
                    const real d2 = ddx * ddx + ddy * ddy + ddz * ddz;
                    const real rsum = ri + radius[j];
                    if (d2 < rsum * rsum) {
                        if (pairs != nullptr) {  // write pass only
                            const int slot = offsets[i] + cnt;
                            if (slot < capacity) {
                                pairs[2 * slot + 0] = i;
                                pairs[2 * slot + 1] = j;
                            } else {
                                atomicExch(overflow, 1);
                            }
                        }
                        ++cnt;
                    }
                }
            }
        }
    }
    if (counts != nullptr) counts[i] = cnt;  // count pass only
}

// Host-side view of one step's result (device consumers read BP.d_pairs /
// BP.d_offsets[n] directly; this struct is the --dump-pairs / step-log path).
struct BroadPhaseResult {
    int* d_pairs = nullptr;   // device pair buffer (same as sim.BP.d_pairs)
    int n_pairs = 0;          // logical total this step (can exceed stored)
    int stored = 0;           // rows actually present (== n_pairs unless clamped)
    int overflow = 0;         // sticky: capacity exceeded at some step
    std::vector<int> pairs;   // host copy, 2*stored ints, row = pair index
};

// One-time allocation after GpuSim::upload (called from dem_gpu.cu so CUB
// stays out of dem_state.cu). No-op for n <= 0 or a degenerate cell size.
inline void alloc_broadphase(GpuSim& sim) {
    BroadPhaseWorkspace& ws = sim.BP;
    const int n = sim.P.n;
    if (n <= 0 || sim.bp_cell <= real(0)) return;
    const std::size_t nu = static_cast<std::size_t>(n);

    ws.capacity = 64 * n;
    CUDA_CHECK(cudaMalloc(&ws.d_keys_in, nu * sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&ws.d_keys_sorted, nu * sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&ws.d_idx_in, nu * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ws.d_idx_sorted, nu * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ws.d_counts, (nu + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ws.d_offsets, (nu + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ws.d_pairs,
                          2 * static_cast<std::size_t>(ws.capacity) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ws.d_overflow, sizeof(int)));

    // Sort values = identity permutation, constant for the whole run.
    std::vector<int> identity(n);
    for (int i = 0; i < n; ++i) identity[i] = i;
    CUDA_CHECK(cudaMemcpy(ws.d_idx_in, identity.data(), nu * sizeof(int),
                          cudaMemcpyHostToDevice));
    // counts[n] stays 0 forever: the count pass writes [0..n-1] only, and the
    // trailing zero makes ExclusiveSum over n+1 items put the total at
    // offsets[n] without a second array.
    CUDA_CHECK(cudaMemset(ws.d_counts, 0, (nu + 1) * sizeof(int)));
    CUDA_CHECK(cudaMemset(ws.d_overflow, 0, sizeof(int)));

    // CUB two-call sizing: query both consumers, allocate the max once.
    std::size_t sort_bytes = 0, scan_bytes = 0;
    CUDA_CHECK(cub::DeviceRadixSort::SortPairs(
        nullptr, sort_bytes, ws.d_keys_in, ws.d_keys_sorted,
        ws.d_idx_in, ws.d_idx_sorted, n));
    CUDA_CHECK(cub::DeviceScan::ExclusiveSum(
        nullptr, scan_bytes, ws.d_counts, ws.d_offsets, n + 1));
    ws.temp_bytes = sort_bytes > scan_bytes ? sort_bytes : scan_bytes;
    CUDA_CHECK(cudaMalloc(&ws.d_temp, ws.temp_bytes));
}

// Per-step pipeline (steps 1-5 above). Launches are stream-ordered on the
// default stream; callers that need the result synchronize (readback_pairs
// does; Task 7 kernels just continue on the same stream).
inline void run_broadphase(GpuSim& sim) {
    const int n = sim.P.n;
    if (n <= 0 || sim.bp_cell <= real(0) || sim.BP.d_temp == nullptr) return;
    const int blocks = (n + kBPThreads - 1) / kBPThreads;
    BroadPhaseWorkspace& ws = sim.BP;

    compute_cell_ids_kernel<<<blocks, kBPThreads>>>(
        sim.P.pos, n, sim.bp_origin[0], sim.bp_origin[1], sim.bp_origin[2],
        sim.bp_cell, ws.d_keys_in);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cub::DeviceRadixSort::SortPairs(
        ws.d_temp, ws.temp_bytes, ws.d_keys_in, ws.d_keys_sorted,
        ws.d_idx_in, ws.d_idx_sorted, n));

    // count pass: counts filled, pairs == nullptr disables the write branch
    candidate_pairs_kernel<<<blocks, kBPThreads>>>(
        sim.P.pos, sim.P.radius, ws.d_keys_sorted, ws.d_idx_sorted, n,
        sim.bp_origin[0], sim.bp_origin[1], sim.bp_origin[2], sim.bp_cell,
        ws.d_counts, nullptr, nullptr, 0, nullptr);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cub::DeviceScan::ExclusiveSum(
        ws.d_temp, ws.temp_bytes, ws.d_counts, ws.d_offsets, n + 1));

    // write pass: counts == nullptr disables the store, offsets drives slots
    candidate_pairs_kernel<<<blocks, kBPThreads>>>(
        sim.P.pos, sim.P.radius, ws.d_keys_sorted, ws.d_idx_sorted, n,
        sim.bp_origin[0], sim.bp_origin[1], sim.bp_origin[2], sim.bp_cell,
        nullptr, ws.d_offsets, ws.d_pairs, ws.capacity, ws.d_overflow);
    CUDA_CHECK(cudaGetLastError());
}

// Host readback of the current step's pair list (synchronizes). Clamped to
// capacity; out.overflow flags the clamp (sticky since alloc).
inline void readback_pairs(const GpuSim& sim, BroadPhaseResult& out) {
    out.d_pairs = sim.BP.d_pairs;
    out.n_pairs = 0;
    out.stored = 0;
    out.overflow = 0;
    out.pairs.clear();
    const int n = sim.P.n;
    if (n <= 0 || sim.BP.d_offsets == nullptr) return;
    CUDA_CHECK(cudaDeviceSynchronize());
    int total = 0, of = 0;
    CUDA_CHECK(cudaMemcpy(&total, sim.BP.d_offsets + n, sizeof(int),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&of, sim.BP.d_overflow, sizeof(int),
                          cudaMemcpyDeviceToHost));
    if (total < 0) total = 0;
    out.n_pairs = total;
    out.overflow = of;
    const int stored = total < sim.BP.capacity ? total : sim.BP.capacity;
    out.stored = stored;
    if (stored > 0) {
        out.pairs.resize(2 * static_cast<std::size_t>(stored));
        CUDA_CHECK(cudaMemcpy(out.pairs.data(), sim.BP.d_pairs,
                              2 * static_cast<std::size_t>(stored) * sizeof(int),
                              cudaMemcpyDeviceToHost));
    }
}

// Cheap count-only readback for the step log (call after a sync, e.g. at
// output steps where the driver already synchronized).
inline int peek_pair_count(const GpuSim& sim) {
    if (sim.P.n <= 0 || sim.BP.d_offsets == nullptr) return 0;
    int total = 0;
    CUDA_CHECK(cudaMemcpy(&total, sim.BP.d_offsets + sim.P.n, sizeof(int),
                          cudaMemcpyDeviceToHost));
    return total < 0 ? 0 : total;
}
