#include "mpi/decomp.hpp"

#include <cstdarg>
#include <cmath>
#include <cstdio>

namespace {

// Rank 0 logs, then the whole job dies. Every rank calls make_decomp() with
// identical inputs, so any decomposition failure is global and unrecoverable.
[[noreturn]] void die(const char* fmt, ...) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::fprintf(stderr, "zdem_mpi fatal: ");
        va_list ap;
        va_start(ap, fmt);
        std::vfprintf(stderr, fmt, ap);
        va_end(ap);
        std::fprintf(stderr, "\n");
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
}

static int smallest_prime_factor(int n) {
    for (int f = 2; (long long)f * f <= (long long)n; ++f) {
        if (n % f == 0) return f;
    }
    return n;
}

// Single source of truth for brick indexing. Owner decisions and (in later
// tasks) ghost/migrate neighbor tests must share this exact floor expression
// so a particle sitting on a brick boundary is owned by exactly one rank
// (left-closed, right-open per axis). Positions outside the global box clamp
// into the first/last brick; callers that care detect out-of-box positions
// separately before asking for the owner.
static int cell_index(const Decomp& d, double v, int axis) {
    double len = d.box_hi[axis] - d.box_lo[axis];
    double t = std::floor((v - d.box_lo[axis]) / len * (double)d.dims[axis]);
    int i = (int)t;
    if (i < 0) i = 0;
    if (i > d.dims[axis] - 1) i = d.dims[axis] - 1;
    return i;
}

}  // namespace

Decomp make_decomp(const SimConfig& cfg, const SimBuild& sim) {
    Decomp d;
    MPI_Comm_rank(MPI_COMM_WORLD, &d.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d.nprocs);
    d.ghost_depth = 2.0 * sim.max_radius;

    // Global box: explicit mpi_box wins, else initial particle bbox padded by
    // margin on every side (default 5*max_radius, cfg.mpi_margin overrides).
    if (cfg.has_mpi_box) {
        for (int a = 0; a < 3; ++a) {
            d.box_lo[a] = cfg.mpi_box[2 * a];
            d.box_hi[a] = cfg.mpi_box[2 * a + 1];
        }
    } else {
        if (sim.particles.empty()) {
            die("make_decomp: no particles to size the global box, set mpi_box explicitly");
        }
        double margin = cfg.mpi_margin >= 0.0 ? cfg.mpi_margin : 5.0 * sim.max_radius;
        for (int a = 0; a < 3; ++a) {
            d.box_lo[a] = 1e300;
            d.box_hi[a] = -1e300;
        }
        for (const Particle& p : sim.particles) {
            const double v[3] = {p.tf.pos.x, p.tf.pos.y, p.tf.pos.z};
            for (int a = 0; a < 3; ++a) {
                if (v[a] < d.box_lo[a]) d.box_lo[a] = v[a];
                if (v[a] > d.box_hi[a]) d.box_hi[a] = v[a];
            }
        }
        for (int a = 0; a < 3; ++a) {
            d.box_lo[a] -= margin;
            d.box_hi[a] += margin;
        }
    }

    double len[3];
    for (int a = 0; a < 3; ++a) {
        len[a] = d.box_hi[a] - d.box_lo[a];
        if (!(len[a] > 0.0)) {
            die("make_decomp: degenerate global box along axis %d (length %g)", a, len[a]);
        }
    }

    // Greedy dims split: start at {1,1,1}, repeatedly multiply the remaining
    // process count (by its smallest prime factor) onto the axis with the
    // largest dims[a]*len[a] so brick footprints stay balanced. cap[a] bounds
    // how far an axis may be split; a brick-thickness violation halves the
    // offending axis and the split retries. Terminating because every retry
    // strictly lowers some cap, and no eligible axis left is fatal.
    int cap[3] = {d.nprocs, d.nprocs, d.nprocs};
    for (;;) {
        int dims[3] = {1, 1, 1};
        int remaining = d.nprocs;
        while (remaining > 1) {
            int f = smallest_prime_factor(remaining);
            int pick = -1;
            double best = -1.0;
            for (int a = 0; a < 3; ++a) {
                if ((long long)dims[a] * f > (long long)cap[a]) continue;
                double v = (double)dims[a] * len[a];
                if (v > best) {
                    best = v;
                    pick = a;
                }
            }
            if (pick < 0) {
                die("make_decomp: cannot split %d ranks over the box (ghost depth %.3f); box too small or too many ranks",
                    d.nprocs, d.ghost_depth);
            }
            dims[pick] *= f;
            remaining /= f;
        }
        bool ok = true;
        bool fatal = false;
        for (int a = 0; a < 3; ++a) {
            if (len[a] / (double)dims[a] < d.ghost_depth) {
                ok = false;
                if (dims[a] == 1) {
                    fatal = true;
                } else {
                    cap[a] = dims[a] / 2;
                }
            }
        }
        if (ok) {
            d.dims[0] = dims[0];
            d.dims[1] = dims[1];
            d.dims[2] = dims[2];
            break;
        }
        if (fatal) {
            die("make_decomp: brick thickness < ghost depth (%.3f) with dims==1 along some axis; box too small or too many ranks",
                d.ghost_depth);
        }
        // Retry with the violating axes capped at half their previous dims.
    }

    int periods[3] = {0, 0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 3, d.dims, periods, 0, &d.cart);
    MPI_Cart_coords(d.cart, d.rank, 3, d.coords);
    for (int a = 0; a < 3; ++a) {
        double thick = len[a] / (double)d.dims[a];
        d.sub_lo[a] = d.box_lo[a] + (double)d.coords[a] * thick;
        d.sub_hi[a] = d.sub_lo[a] + thick;
    }
    return d;
}

int owner_rank(const Decomp& d, const Vec3& pos) {
    int c[3] = {cell_index(d, pos.x, 0), cell_index(d, pos.y, 1), cell_index(d, pos.z, 2)};
    int r = 0;
    MPI_Cart_rank(d.cart, c, &r);
    return r;
}

bool in_sub(const Decomp& d, const Vec3& pos) {
    const double v[3] = {pos.x, pos.y, pos.z};
    for (int a = 0; a < 3; ++a) {
        if (cell_index(d, v[a], a) != d.coords[a]) return false;
    }
    return true;
}
