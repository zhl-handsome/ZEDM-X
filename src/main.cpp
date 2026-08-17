#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "core/mat3.hpp"
#include "core/quat.hpp"
#include "core/transform.hpp"
#include "core/vec3.hpp"
#include "geometry/mesh.hpp"
#include "geometry/mesh_build.hpp"
#include "geometry/stl_io.hpp"

using namespace zdem;

struct ParticleInit {
    std::string stl_path;
    Vec3 pos{0.0, 0.0, 0.0};
    Vec3 vel{0.0, 0.0, 0.0};
    Vec3 omega{0.0, 0.0, 0.0};
    Quat rot{1.0, 0.0, 0.0, 0.0};
    double scale = 0.0;
    double density = 1.0;
    double young = 1e7;
    double poisson = 0.25;
    double mu = 0.5;
    double restitution = 0.5;
};

struct Particle {
    Transform tf;
    Vec3 vel;
    Vec3 omega;
    Vec3 L;
    double mass = 1.0;
    double inv_mass = 1.0;
    Mat3 inertia_body;
    Mat3 inertia_body_inv;
    double radius = 0.0;
    double equiv_radius = 0.0;
    double young = 0.0;
    double poisson = 0.0;
    double mu = 0.0;
    double restitution = 0.0;
    int mesh_index = 0;
};

struct WallInit {
    std::string stl_path;
    Vec3 pos{0.0, 0.0, 0.0};
    Quat rot{1.0, 0.0, 0.0, 0.0};
    double scale = 1.0;
    double mu = 0.5;
    double restitution = 0.5;
};

struct Wall {
    Mesh mesh;
    Transform tf;
    double mu = 0.5;
    double restitution = 0.5;
    std::vector<Vec3> tri_normals;
};


static bool plane_from_tri(const std::array<Vec3, 3>& tri, Vec3& n, double& d) {
    n = cross(tri[1] - tri[0], tri[2] - tri[0]);
    double ln = norm(n);
    if (ln < 1e-30) {
        return false;
    }
    n = n / ln;
    d = -dot(n, tri[0]);
    return true;
}

static bool point_in_tri(const Vec3& p, const std::array<Vec3, 3>& tri, const Vec3& n) {
    Vec3 a = tri[0];
    Vec3 b = tri[1];
    Vec3 c = tri[2];
    Vec3 ab = b - a;
    Vec3 bc = c - b;
    Vec3 ca = a - c;
    Vec3 ap = p - a;
    Vec3 bp = p - b;
    Vec3 cp = p - c;
    double c1 = dot(cross(ab, ap), n);
    double c2 = dot(cross(bc, bp), n);
    double c3 = dot(cross(ca, cp), n);
    return (c1 >= -1e-10) && (c2 >= -1e-10) && (c3 >= -1e-10);
}

static bool seg_plane_intersection(const Vec3& p0, const Vec3& p1, const Vec3& n, double d, Vec3& out) {
    Vec3 u = p1 - p0;
    double denom = dot(n, u);
    double num = -(dot(n, p0) + d);
    if (std::abs(denom) < 1e-30) {
        return false;
    }
    double t = num / denom;
    if (t < -1e-10 || t > 1.0 + 1e-10) {
        return false;
    }
    t = std::min(1.0, std::max(0.0, t));
    out = p0 + u * t;
    return true;
}

static std::vector<Vec3> unique_points(const std::vector<Vec3>& pts, double tol) {
    std::vector<Vec3> out;
    double tol2 = tol * tol;
    for (const auto& p : pts) {
        bool keep = true;
        for (const auto& q : out) {
            if (norm2(p - q) <= tol2) {
                keep = false;
                break;
            }
        }
        if (keep) {
            out.push_back(p);
        }
    }
    return out;
}

static bool tri_tri_intersection_segment(const std::array<Vec3, 3>& triA, const Vec3& nA, double dA,
                                         const std::array<Vec3, 3>& triB, const Vec3& nB, double dB,
                                         Vec3& out_p0, Vec3& out_p1) {
    Vec3 tau = cross(nA, nB);
    double lt = norm(tau);
    if (lt < 1e-20) {
        return false;
    }
    tau = tau / lt;

    std::vector<Vec3> candidates;
    for (int i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        Vec3 p;
        if (seg_plane_intersection(triA[i], triA[j], nB, dB, p) && point_in_tri(p, triB, nB)) {
            candidates.push_back(p);
        }
    }
    for (int i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        Vec3 p;
        if (seg_plane_intersection(triB[i], triB[j], nA, dA, p) && point_in_tri(p, triA, nA)) {
            candidates.push_back(p);
        }
    }

    candidates = unique_points(candidates, 1e-9);
    if (candidates.size() < 2) {
        return false;
    }
    double min_t = std::numeric_limits<double>::infinity();
    double max_t = -std::numeric_limits<double>::infinity();
    int i0 = -1, i1 = -1;
    for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
        double t = dot(candidates[i], tau);
        if (t < min_t) { min_t = t; i0 = i; }
        if (t > max_t) { max_t = t; i1 = i; }
    }
    if (i0 == i1) {
        return false;
    }
    Vec3 p0 = candidates[i0];
    Vec3 p1 = candidates[i1];
    if (norm2(p1 - p0) < 1e-18) {
        return false;
    }
    if (dot(p1 - p0, tau) < 0.0) {
        std::swap(p0, p1);
    }
    out_p0 = p0;
    out_p1 = p1;
    return true;
}

struct Int3 {
    int x = 0;
    int y = 0;
    int z = 0;
    bool operator==(const Int3& o) const { return x == o.x && y == o.y && z == o.z; }
};

struct Int3Hash {
    std::size_t operator()(const Int3& k) const {
        std::size_t h1 = std::hash<int>{}(k.x);
        std::size_t h2 = std::hash<int>{}(k.y);
        std::size_t h3 = std::hash<int>{}(k.z);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

static Int3 cell_index(const Vec3& p, double h) {
    return Int3{
        static_cast<int>(std::floor(p.x / h)),
        static_cast<int>(std::floor(p.y / h)),
        static_cast<int>(std::floor(p.z / h))
    };
}

struct UnionFind {
    std::vector<int> parent;
    std::vector<int> rank;

    explicit UnionFind(int n = 0) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }

    void unite(int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) return;
        if (rank[ra] < rank[rb]) parent[ra] = rb;
        else if (rank[ra] > rank[rb]) parent[rb] = ra;
        else { parent[rb] = ra; rank[ra] += 1; }
    }
};

static void snap_endpoints(const std::vector<std::pair<Vec3, Vec3>>& segments,
                           double tol,
                           std::vector<Vec3>& out_vertices,
                           std::vector<std::pair<int, int>>& out_edges) {
    double h = tol;
    int nseg = static_cast<int>(segments.size());
    std::vector<Vec3> pts(2 * nseg);
    for (int i = 0; i < nseg; ++i) {
        pts[2 * i] = segments[i].first;
        pts[2 * i + 1] = segments[i].second;
    }

    UnionFind uf(2 * nseg);
    std::unordered_map<Int3, std::vector<int>, Int3Hash> grid;

    for (int i = 0; i < 2 * nseg; ++i) {
        const Vec3& p = pts[i];
        Int3 c = cell_index(p, h);
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    Int3 cc{c.x + dx, c.y + dy, c.z + dz};
                    auto it = grid.find(cc);
                    if (it == grid.end()) continue;
                    for (int j : it->second) {
                        if (norm2(p - pts[j]) <= h * h) {
                            uf.unite(i, j);
                        }
                    }
                }
            }
        }
        grid[c].push_back(i);
    }

    std::unordered_map<int, int> root_to_vid;
    std::unordered_map<int, Vec3> sums;
    std::unordered_map<int, int> counts;
    out_vertices.clear();
    for (int i = 0; i < 2 * nseg; ++i) {
        int r = uf.find(i);
        if (root_to_vid.find(r) == root_to_vid.end()) {
            int vid = static_cast<int>(root_to_vid.size());
            root_to_vid[r] = vid;
            sums[r] = pts[i];
            counts[r] = 1;
        } else {
            sums[r] += pts[i];
            counts[r] += 1;
        }
    }
    out_vertices.resize(root_to_vid.size());
    for (const auto& kv : root_to_vid) {
        int r = kv.first;
        int vid = kv.second;
        out_vertices[vid] = sums[r] / static_cast<double>(counts[r]);
    }

    out_edges.clear();
    for (int s = 0; s < nseg; ++s) {
        int a = root_to_vid[uf.find(2 * s)];
        int b = root_to_vid[uf.find(2 * s + 1)];
        if (a != b) {
            out_edges.emplace_back(a, b);
        }
    }
}

static std::pair<double, double> point_segment_dist2_and_t(const Vec3& p, const Vec3& a, const Vec3& b) {
    Vec3 ab = b - a;
    double denom = dot(ab, ab);
    if (denom <= 1e-12) {
        return {norm2(p - a), 0.0};
    }
    double t = dot(p - a, ab) / denom;
    double t_clamped = std::min(1.0, std::max(0.0, t));
    Vec3 proj = a + ab * t_clamped;
    double d2 = norm2(p - proj);
    return {d2, t_clamped};
}

static double seg_seg_dist2(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2) {
    Vec3 d1 = q1 - p1;
    Vec3 d2 = q2 - p2;
    Vec3 r = p1 - p2;
    double a = dot(d1, d1);
    double e = dot(d2, d2);
    double f = dot(d2, r);

    if (a <= 1e-12 && e <= 1e-12) {
        return dot(p1 - p2, p1 - p2);
    }
    if (a <= 1e-12) {
        double t = std::min(1.0, std::max(0.0, f / e));
        Vec3 c = p2 + d2 * t;
        return dot(p1 - c, p1 - c);
    }
    if (e <= 1e-12) {
        double s = std::min(1.0, std::max(0.0, -dot(d1, r) / a));
        Vec3 c = p1 + d1 * s;
        return dot(c - p2, c - p2);
    }

    double b = dot(d1, d2);
    double c = dot(d1, r);
    double denom = a * e - b * b;
    double s = 0.0;
    if (std::abs(denom) > 0.0) {
        s = std::min(1.0, std::max(0.0, (b * f - c * e) / denom));
    }
    double t = (b * s + f) / e;
    if (t < 0.0) {
        t = 0.0;
        s = std::min(1.0, std::max(0.0, -c / a));
    } else if (t > 1.0) {
        t = 1.0;
        s = std::min(1.0, std::max(0.0, (b - c) / a));
    }
    Vec3 cp1 = p1 + d1 * s;
    Vec3 cp2 = p2 + d2 * t;
    Vec3 diff = cp1 - cp2;
    return dot(diff, diff);
}

static void label_segments_by_contact_greedy(const std::vector<std::pair<Vec3, Vec3>>& segments,
                                             double tol,
                                             std::vector<std::vector<int>>& comps) {
    double h = tol;
    int nseg = static_cast<int>(segments.size());
    comps.clear();
    if (nseg == 0) {
        return;
    }

    std::unordered_map<Int3, std::vector<int>, Int3Hash> grid;
    std::vector<std::pair<Vec3, Vec3>> seg_aabb;
    std::vector<std::array<Vec3, 3>> rep_points;
    seg_aabb.reserve(nseg);
    rep_points.reserve(nseg);

    for (int i = 0; i < nseg; ++i) {
        Vec3 p0 = segments[i].first;
        Vec3 p1 = segments[i].second;
        Vec3 mid = (p0 + p1) * 0.5;
        rep_points.push_back({p0, p1, mid});
        Vec3 mn{std::min(p0.x, p1.x) - h, std::min(p0.y, p1.y) - h, std::min(p0.z, p1.z) - h};
        Vec3 mx{std::max(p0.x, p1.x) + h, std::max(p0.y, p1.y) + h, std::max(p0.z, p1.z) + h};
        seg_aabb.push_back({mn, mx});
        for (const auto& rp : rep_points.back()) {
            Int3 c = cell_index(rp, h);
            grid[c].push_back(i);
        }
    }

    auto nearby_candidates = [&](int i) {
        std::vector<int> cand;
        std::unordered_map<int, bool> mark;
        for (const auto& rp : rep_points[i]) {
            Int3 c = cell_index(rp, h);
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        Int3 cc{c.x + dx, c.y + dy, c.z + dz};
                        auto it = grid.find(cc);
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            if (j == i) continue;
                            if (!mark[j]) {
                                mark[j] = true;
                                cand.push_back(j);
                            }
                        }
                    }
                }
            }
        }
        return cand;
    };

    auto aabb_overlap = [&](int i, int j) {
        const Vec3& mn1 = seg_aabb[i].first;
        const Vec3& mx1 = seg_aabb[i].second;
        const Vec3& mn2 = seg_aabb[j].first;
        const Vec3& mx2 = seg_aabb[j].second;
        return !(mx1.x < mn2.x || mn1.x > mx2.x ||
                 mx1.y < mn2.y || mn1.y > mx2.y ||
                 mx1.z < mn2.z || mn1.z > mx2.z);
    };

    double tol2 = h * h;
    std::vector<bool> visited(nseg, false);
    for (int seed = 0; seed < nseg; ++seed) {
        if (visited[seed]) continue;
        std::vector<int> queue;
        queue.push_back(seed);
        visited[seed] = true;
        std::vector<int> comp;
        comp.push_back(seed);
        while (!queue.empty()) {
            int i = queue.back();
            queue.pop_back();
            Vec3 p1 = segments[i].first;
            Vec3 q1 = segments[i].second;
            for (int j : nearby_candidates(i)) {
                if (visited[j]) continue;
                if (!aabb_overlap(i, j)) continue;
                Vec3 p2 = segments[j].first;
                Vec3 q2 = segments[j].second;
                if (seg_seg_dist2(p1, q1, p2, q2) <= tol2) {
                    visited[j] = true;
                    queue.push_back(j);
                    comp.push_back(j);
                }
            }
        }
        comps.push_back(comp);
    }
    std::sort(comps.begin(), comps.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
}

static void split_t_junctions(const std::vector<Vec3>& vertices,
                              const std::vector<std::pair<int, int>>& edges,
                              double tol,
                              std::vector<Vec3>& out_vertices,
                              std::vector<std::pair<int, int>>& out_edges) {
    double h = tol;
    double tol2 = h * h;
    std::unordered_map<Int3, std::vector<int>, Int3Hash> vgrid;
    for (int vid = 0; vid < static_cast<int>(vertices.size()); ++vid) {
        Int3 c = cell_index(vertices[vid], h);
        vgrid[c].push_back(vid);
    }

    out_vertices = vertices;
    out_edges.clear();
    for (const auto& e : edges) {
        int a = e.first;
        int b = e.second;
        Vec3 pa = vertices[a];
        Vec3 pb = vertices[b];
        Vec3 mn{std::min(pa.x, pb.x) - h, std::min(pa.y, pb.y) - h, std::min(pa.z, pb.z) - h};
        Vec3 mx{std::max(pa.x, pb.x) + h, std::max(pa.y, pb.y) + h, std::max(pa.z, pb.z) + h};
        Int3 cmin = cell_index(mn, h);
        Int3 cmax = cell_index(mx, h);

        std::vector<std::pair<double, int>> splits;
        for (int ix = cmin.x; ix <= cmax.x; ++ix) {
            for (int iy = cmin.y; iy <= cmax.y; ++iy) {
                for (int iz = cmin.z; iz <= cmax.z; ++iz) {
                    auto it = vgrid.find(Int3{ix, iy, iz});
                    if (it == vgrid.end()) continue;
                    for (int vid : it->second) {
                        if (vid == a || vid == b) continue;
                        auto dt = point_segment_dist2_and_t(vertices[vid], pa, pb);
                        if (dt.first <= tol2 && dt.second > 1e-6 && dt.second < 1.0 - 1e-6) {
                            splits.emplace_back(dt.second, vid);
                        }
                    }
                }
            }
        }
        if (splits.empty()) {
            out_edges.push_back(e);
            continue;
        }
        std::sort(splits.begin(), splits.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        std::vector<int> chain;
        chain.push_back(a);
        int last = -1;
        for (const auto& s : splits) {
            if (s.second == last) continue;
            chain.push_back(s.second);
            last = s.second;
        }
        chain.push_back(b);
        for (int i = 0; i + 1 < static_cast<int>(chain.size()); ++i) {
            if (chain[i] != chain[i + 1]) {
                out_edges.emplace_back(chain[i], chain[i + 1]);
            }
        }
    }
}

static std::vector<std::vector<int>> extract_loops(const std::vector<Vec3>& vertices,
                                                   const std::vector<std::pair<int, int>>& edges) {
    int nv = static_cast<int>(vertices.size());
    std::vector<std::vector<std::pair<int, int>>> adj(nv);
    for (int ei = 0; ei < static_cast<int>(edges.size()); ++ei) {
        int a = edges[ei].first;
        int b = edges[ei].second;
        adj[a].push_back({b, ei});
        adj[b].push_back({a, ei});
    }

    std::vector<bool> used(edges.size(), false);
    std::vector<int> start_vertices(nv);
    for (int i = 0; i < nv; ++i) start_vertices[i] = i;
    std::sort(start_vertices.begin(), start_vertices.end(), [&](int i, int j) {
        int di = static_cast<int>(adj[i].size());
        int dj = static_cast<int>(adj[j].size());
        if ((di == 2) != (dj == 2)) return di == 2;
        return di > dj;
    });

    std::vector<std::vector<int>> loops;
    for (int sv : start_vertices) {
        for (const auto& nb : adj[sv]) {
            int ei = nb.second;
            if (used[ei]) continue;
            std::vector<int> path;
            path.push_back(sv);
            int prev = sv;
            int cur = nb.first;
            used[ei] = true;
            path.push_back(cur);

            int steps = 0;
            while (steps < 100000) {
                steps += 1;
                if (cur == sv) {
                    if (path.size() > 1) {
                        path.pop_back();
                        loops.push_back(path);
                    }
                    break;
                }
                std::vector<std::pair<int, int>> cand;
                for (const auto& e : adj[cur]) {
                    if (!used[e.second]) cand.push_back(e);
                }
                if (cand.empty()) {
                    break;
                }
                std::pair<int, int> pick = cand[0];
                for (const auto& c : cand) {
                    if (c.first != prev) { pick = c; break; }
                }
                used[pick.second] = true;
                prev = cur;
                cur = pick.first;
                path.push_back(cur);
            }
        }
    }
    return loops;
}

static std::pair<Vec3, Vec3> accumulate_Sn_Gn_from_polyline(const std::vector<Vec3>& loop_pts) {
    Vec3 Sn{0.0, 0.0, 0.0};
    Vec3 Gn{0.0, 0.0, 0.0};
    int m = static_cast<int>(loop_pts.size());
    for (int i = 0; i < m; ++i) {
        const Vec3& p0 = loop_pts[i];
        const Vec3& p1 = loop_pts[(i + 1) % m];
        Sn += cross(p0, p1) * 0.5;
        Vec3 dx = p1 - p0;
        double s = dot(p0, p1) + (1.0 / 3.0) * dot(dx, dx);
        Gn += dx * (-(1.0 / 3.0) * s);
    }
    return {Sn, Gn};
}

static bool contact_point_xc0(const Vec3& Sn, const Vec3& Gn, Vec3& xc0, Vec3& nA, double& area) {
    area = norm(Sn);
    if (area < 1e-14) {
        xc0 = Vec3{0.0, 0.0, 0.0};
        nA = Vec3{0.0, 0.0, 0.0};
        return false;
    }
    nA = Sn * (-1.0 / area);
    xc0 = cross(nA, Gn) / area;
    return true;
}

// ============== Sutherland-Hodgman Polygon Clipping for Wall Contact ==============

// 2D point for clipping algorithm
struct Vec2 {
    double x = 0.0;
    double y = 0.0;
    Vec2() = default;
    Vec2(double x_, double y_) : x(x_), y(y_) {}
    Vec2 operator+(const Vec2& o) const { return Vec2{x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2& o) const { return Vec2{x - o.x, y - o.y}; }
    Vec2 operator*(double s) const { return Vec2{x * s, y * s}; }
};

static inline double cross2d(const Vec2& a, const Vec2& b) {
    return a.x * b.y - a.y * b.x;
}

// Check if point p is inside the clipping edge (on the left side of edge from cp1 to cp2)
static inline bool is_inside_clip(const Vec2& p, const Vec2& cp1, const Vec2& cp2) {
    return cross2d(cp2 - cp1, p - cp1) >= 0.0;
}

// Compute intersection of line segment (s, e) with clipping edge (cp1, cp2)
static inline Vec2 compute_intersection_clip(const Vec2& s, const Vec2& e,
                                              const Vec2& cp1, const Vec2& cp2) {
    Vec2 dc = cp1 - cp2;
    Vec2 dp = s - e;
    double n1 = cross2d(cp1, cp2);
    double n2 = cross2d(s, e);
    double n3 = cross2d(dc, dp);
    if (std::abs(n3) < 1e-30) {
        return s;
    }
    double t = (n1 * dp.x - n2 * dc.x) / n3;
    double u = (n1 * dp.y - n2 * dc.y) / n3;
    return Vec2{t, u};
}

// Sutherland-Hodgman polygon clipping algorithm
static std::vector<Vec2> sutherland_hodgman_clip(const std::vector<Vec2>& subject,
                                                  const std::vector<Vec2>& clip) {
    std::vector<Vec2> output = subject;
    if (clip.size() < 3 || subject.size() < 3) {
        return output;
    }

    Vec2 cp1 = clip.back();
    for (const auto& cp2 : clip) {
        if (output.empty()) break;

        std::vector<Vec2> input = output;
        output.clear();

        Vec2 s = input.back();
        for (const auto& e : input) {
            if (is_inside_clip(e, cp1, cp2)) {
                if (!is_inside_clip(s, cp1, cp2)) {
                    output.push_back(compute_intersection_clip(s, e, cp1, cp2));
                }
                output.push_back(e);
            } else if (is_inside_clip(s, cp1, cp2)) {
                output.push_back(compute_intersection_clip(s, e, cp1, cp2));
            }
            s = e;
        }
        cp1 = cp2;
    }
    return output;
}

// Compute 2D polygon area using shoelace formula
static double polygon_area_2d(const std::vector<Vec2>& poly) {
    if (poly.size() < 3) return 0.0;
    double area = 0.0;
    int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += cross2d(poly[i], poly[j]);
    }
    return std::abs(area) * 0.5;
}

// Compute 2D polygon centroid
static Vec2 polygon_centroid_2d(const std::vector<Vec2>& poly) {
    if (poly.empty()) return Vec2{0.0, 0.0};
    Vec2 c{0.0, 0.0};
    for (const auto& p : poly) {
        c.x += p.x;
        c.y += p.y;
    }
    double n = static_cast<double>(poly.size());
    return Vec2{c.x / n, c.y / n};
}

// Build orthonormal basis from a normal vector
static void build_basis(const Vec3& n, Vec3& u, Vec3& v) {
    // Choose a vector not parallel to n
    if (std::abs(n.x) < 0.9) {
        u = normalize(cross(n, Vec3{1.0, 0.0, 0.0}));
    } else {
        u = normalize(cross(n, Vec3{0.0, 1.0, 0.0}));
    }
    v = cross(n, u);
}

// Project 3D point to 2D using basis vectors
static Vec2 project_to_2d(const Vec3& p, const Vec3& origin, const Vec3& u, const Vec3& v) {
    Vec3 d = p - origin;
    return Vec2{dot(d, u), dot(d, v)};
}

// Unproject 2D point back to 3D
static Vec3 unproject_to_3d(const Vec2& p2d, const Vec3& origin, const Vec3& u, const Vec3& v) {
    return origin + u * p2d.x + v * p2d.y;
}

// Point-in-triangle (2D, inclusive of the boundary) -- footprint test for
// the finite wall triangles in the per-vertex wall penalty contact.
static bool point_in_tri_2d(const Vec2& p, const Vec2& a, const Vec2& b, const Vec2& c) {
    auto sign = [](const Vec2& p1, const Vec2& p2, const Vec2& p3) {
        return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
    };
    double d1 = sign(p, a, b), d2 = sign(p, b, c), d3 = sign(p, c, a);
    bool has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    bool has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
    return !(has_neg && has_pos);
}

static std::vector<std::vector<Vec3>> get_particle_plane_section_loops(const Mesh& mesh,
                                                                       const Transform& tf,
                                                                       const Vec3& plane_n,
                                                                       double plane_d,
                                                                       double tol) {
    std::vector<std::pair<Vec3, Vec3>> segments;
    std::vector<std::array<Vec3, 3>> tris = transform_tris(mesh, tf);

    for (const auto& tri : tris) {
        std::vector<Vec3> pts;
        for (int i = 0; i < 3; ++i) {
            int j = (i + 1) % 3;
            Vec3 p0 = tri[i];
            Vec3 p1 = tri[j];
            double d0 = dot(plane_n, p0) + plane_d;
            double d1 = dot(plane_n, p1) + plane_d;
            bool on0 = std::abs(d0) <= tol;
            bool on1 = std::abs(d1) <= tol;

            if (on0 && on1) {
                if (norm2(p1 - p0) > tol * tol) {
                    segments.emplace_back(p0, p1);
                }
                continue;
            }
            if (on0) {
                pts.push_back(p0);
            } else if (on1) {
                pts.push_back(p1);
            } else if (d0 * d1 < 0.0) {
                double t = d0 / (d0 - d1);
                t = std::min(1.0, std::max(0.0, t));
                pts.push_back(p0 + (p1 - p0) * t);
            }
        }

        pts = unique_points(pts, tol);
        if (pts.size() == 2 && norm2(pts[1] - pts[0]) > tol * tol) {
            segments.emplace_back(pts[0], pts[1]);
        } else if (pts.size() > 2) {
            for (std::size_t i = 0; i < pts.size(); ++i) {
                for (std::size_t j = i + 1; j < pts.size(); ++j) {
                    if (norm2(pts[j] - pts[i]) > tol * tol) {
                        segments.emplace_back(pts[i], pts[j]);
                    }
                }
            }
        }
    }

    if (segments.empty()) {
        return {};
    }

    std::vector<Vec3> V;
    std::vector<std::pair<int, int>> edges;
    snap_endpoints(segments, tol, V, edges);

    std::vector<Vec3> V2;
    std::vector<std::pair<int, int>> edges2;
    split_t_junctions(V, edges, tol, V2, edges2);

    std::vector<std::vector<int>> loop_vids = extract_loops(V2, edges2);
    std::vector<std::vector<Vec3>> loops;
    loops.reserve(loop_vids.size());
    for (const auto& lv : loop_vids) {
        if (lv.size() < 3) {
            continue;
        }
        std::vector<Vec3> loop;
        loop.reserve(lv.size());
        for (int vid : lv) {
            loop.push_back(V2[vid]);
        }
        loops.push_back(loop);
    }
    return loops;
}


static Vec3 support_point(const Mesh& mesh, const Transform& tf, const Vec3& dir) {
    Quat qc = quat_conj(tf.rot);
    Vec3 local_dir = quat_rotate(qc, dir);
    double best = -std::numeric_limits<double>::infinity();
    Vec3 best_v;
    for (const auto& v : mesh.vertices) {
        double d = dot(v, local_dir);
        if (d > best) {
            best = d;
            best_v = v;
        }
    }
    return quat_rotate(tf.rot, best_v) + tf.pos;
}

static Vec3 support_minkowski(const Mesh& a, const Transform& ta,
                              const Mesh& b, const Transform& tb,
                              const Vec3& dir) {
    Vec3 p1 = support_point(a, ta, dir);
    Vec3 p2 = support_point(b, tb, dir * -1.0);
    return p1 - p2;
}

struct Simplex {
    std::vector<Vec3> pts;
};

static bool handle_simplex(Simplex& s, Vec3& dir) {
    if (s.pts.size() == 2) {
        Vec3 a = s.pts[1];
        Vec3 b = s.pts[0];
        Vec3 ab = b - a;
        Vec3 ao = Vec3{} - a;
        if (dot(ab, ao) > 0.0) {
            dir = cross(cross(ab, ao), ab);
        } else {
            s.pts = {a};
            dir = ao;
        }
        return false;
    }
    if (s.pts.size() == 3) {
        Vec3 a = s.pts[2];
        Vec3 b = s.pts[1];
        Vec3 c = s.pts[0];
        Vec3 ab = b - a;
        Vec3 ac = c - a;
        Vec3 ao = Vec3{} - a;
        Vec3 abc = cross(ab, ac);

        Vec3 ab_perp = cross(abc, ab);
        if (dot(ab_perp, ao) > 0.0) {
            s.pts = {b, a};
            dir = cross(cross(ab, ao), ab);
            return false;
        }
        Vec3 ac_perp = cross(ac, abc);
        if (dot(ac_perp, ao) > 0.0) {
            s.pts = {c, a};
            dir = cross(cross(ac, ao), ac);
            return false;
        }

        if (dot(abc, ao) > 0.0) {
            dir = abc;
        } else {
            s.pts = {b, c, a};
            dir = abc * -1.0;
        }
        return false;
    }
    if (s.pts.size() == 4) {
        Vec3 a = s.pts[3];
        Vec3 b = s.pts[2];
        Vec3 c = s.pts[1];
        Vec3 d = s.pts[0];
        Vec3 ao = Vec3{} - a;
        Vec3 ab = b - a;
        Vec3 ac = c - a;
        Vec3 ad = d - a;

        Vec3 abc = cross(ab, ac);
        Vec3 acd = cross(ac, ad);
        Vec3 adb = cross(ad, ab);

        if (dot(abc, ao) > 0.0) {
            s.pts = {c, b, a};
            dir = abc;
            return false;
        }
        if (dot(acd, ao) > 0.0) {
            s.pts = {d, c, a};
            dir = acd;
            return false;
        }
        if (dot(adb, ao) > 0.0) {
            s.pts = {b, d, a};
            dir = adb;
            return false;
        }
        return true;
    }
    return false;
}

static bool gjk_intersect(const Mesh& a, const Transform& ta,
                          const Mesh& b, const Transform& tb,
                          Simplex& simplex) {
    Vec3 dir = ta.pos - tb.pos;
    if (norm(dir) < 1e-12) {
        dir = Vec3{1.0, 0.0, 0.0};
    }
    simplex.pts.clear();
    simplex.pts.push_back(support_minkowski(a, ta, b, tb, dir));
    dir = simplex.pts.back() * -1.0;

    for (int iter = 0; iter < 32; ++iter) {
        Vec3 p = support_minkowski(a, ta, b, tb, dir);
        if (dot(p, dir) < 0.0) {
            return false;
        }
        simplex.pts.push_back(p);
        if (handle_simplex(simplex, dir)) {
            return true;
        }
    }
    return false;
}

struct EPAFace {
    int a = 0;
    int b = 0;
    int c = 0;
    Vec3 n;
    double d = 0.0;
};

static EPAFace make_face(const std::vector<Vec3>& verts, int a, int b, int c) {
    Vec3 v0 = verts[a];
    Vec3 v1 = verts[b];
    Vec3 v2 = verts[c];
    Vec3 n = normalize(cross(v1 - v0, v2 - v0));
    double d = dot(n, v0);
    if (d < 0.0) {
        std::swap(b, c);
        n = n * -1.0;
        d = -d;
    }
    EPAFace f;
    f.a = a; f.b = b; f.c = c; f.n = n; f.d = d;
    return f;
}

static bool epa_penetration(const Mesh& a, const Transform& ta,
                            const Mesh& b, const Transform& tb,
                            const Simplex& simplex,
                            Vec3& out_normal, double& out_depth) {
    if (simplex.pts.size() < 4) {
        return false;
    }
    std::vector<Vec3> verts = simplex.pts;
    std::vector<EPAFace> faces;
    faces.push_back(make_face(verts, 0, 1, 2));
    faces.push_back(make_face(verts, 0, 3, 1));
    faces.push_back(make_face(verts, 0, 2, 3));
    faces.push_back(make_face(verts, 1, 3, 2));

    const double tol = 1e-6;
    for (int iter = 0; iter < 32; ++iter) {
        int best = -1;
        double min_d = std::numeric_limits<double>::infinity();
        for (int i = 0; i < static_cast<int>(faces.size()); ++i) {
            if (faces[i].d < min_d) {
                min_d = faces[i].d;
                best = i;
            }
        }
        if (best < 0) {
            return false;
        }
        EPAFace f = faces[best];
        Vec3 p = support_minkowski(a, ta, b, tb, f.n);
        double d = dot(p, f.n);
        if (d - f.d < tol) {
            out_normal = f.n;
            out_depth = d;
            return true;
        }

        int new_index = static_cast<int>(verts.size());
        verts.push_back(p);

        std::vector<std::pair<int, int>> edges;
        std::vector<EPAFace> new_faces;
        for (const auto& face : faces) {
            Vec3 v0 = verts[face.a];
            if (dot(face.n, p - v0) > 0.0) {
                auto add_edge = [&](int a, int b) {
                    auto it = std::find(edges.begin(), edges.end(), std::make_pair(b, a));
                    if (it != edges.end()) {
                        edges.erase(it);
                    } else {
                        edges.emplace_back(a, b);
                    }
                };
                add_edge(face.a, face.b);
                add_edge(face.b, face.c);
                add_edge(face.c, face.a);
            } else {
                new_faces.push_back(face);
            }
        }

        for (const auto& e : edges) {
            new_faces.push_back(make_face(verts, e.first, e.second, new_index));
        }
        faces.swap(new_faces);
    }
    return false;
}

static Mat3 inertia_world(const Particle& p) {
    Mat3 R = quat_to_mat3(p.tf.rot);
    Mat3 Rt = mat3_transpose(R);
    return mat3_mul(mat3_mul(R, p.inertia_body), Rt);
}


static Mat3 inertia_world_inv(const Particle& p) {
    Mat3 R = quat_to_mat3(p.tf.rot);
    Mat3 Rt = mat3_transpose(R);
    return mat3_mul(mat3_mul(R, p.inertia_body_inv), Rt);
}

// Exact exponential rotation increment for a world-frame angular velocity:
// dq = {cos(th/2), w_hat*sin(th/2)}, q' = dq (x) q. Unlike the linearized
// q += 0.5*qdot*dt this preserves the quaternion norm exactly per step.
static Quat quat_step_world(const Quat& q, const Vec3& w, double dt) {
    double wmag = std::sqrt(w.x * w.x + w.y * w.y + w.z * w.z);
    double th = wmag * dt;
    Quat dq;
    if (th < 1e-300) {
        dq.w = 1.0; dq.x = 0.0; dq.y = 0.0; dq.z = 0.0;
    } else {
        double s = std::sin(0.5 * th) / wmag;
        dq.w = std::cos(0.5 * th);
        dq.x = w.x * s; dq.y = w.y * s; dq.z = w.z * s;
    }
    return quat_mul(dq, q);
}

static void integrate_particle(Particle& p, const Vec3& force, const Vec3& torque, const Vec3& gravity, double dt) {
    Vec3 acc = force * p.inv_mass + gravity;
    p.vel += acc * dt;
    p.tf.pos += p.vel * dt;

    p.L += torque * dt;

    // Midpoint attitude integration with exact exponential rotation steps.
    // The old explicit form (omega from the PRE-rotation inertia, quaternion
    // advanced linearly) injects rotational energy for an asymmetric body:
    // a torque-free spinner at 400 rad/s grew to 2875 rad/s (isolated test).
    // Here the half-step attitude is used to re-evaluate omega before the
    // full step, and each rotation increment is a norm-preserving exponential.
    Quat q0 = p.tf.rot;
    Vec3 w0 = mat3_mul_vec3(inertia_world_inv(p), p.L);
    p.tf.rot = quat_step_world(q0, w0, 0.5 * dt);
    Vec3 wh = mat3_mul_vec3(inertia_world_inv(p), p.L);
    p.tf.rot = quat_step_world(q0, wh, dt);
    p.omega = mat3_mul_vec3(inertia_world_inv(p), p.L);
}

struct SimConfig {
    int n = 2;
    int steps = 1;
    double dt = 1e-4;
    double spacing = 2.5;
    Vec3 gravity{0.0, -9.81, 0.0};
    bool split_contacts = true;
    bool center_mesh = true;
    bool contact_debug = false;
    double tangential_damping = 1.0;  // 切向阻尼系数: ct = tangential_damping * sqrt(kt*m_eff)
    std::string stl_path;
    std::string vtk_prefix = "particles";
    int output_interval = 1;
    std::string output_dir = "output";
    Vec3 v0{0.0, 0.0, 0.0};
    Vec3 v1{0.0, 0.0, 0.0};
    std::vector<ParticleInit> particle_inits;
    std::vector<WallInit> wall_inits;
};

static bool parse_config_file(const std::string& path, SimConfig& cfg) {
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

static bool parse_args(int argc, char** argv, std::string& config_path) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--config" && i + 1 < argc) {
            config_path = argv[++i];
        } else if (arg == "--help") {
            return false;
        } else {
            std::cerr << "Unknown arg: " << arg << "\n";
            return false;
        }
    }
    return !config_path.empty();
}

static void usage() {
    std::cout << "Usage: zdem_cpu --config path.txt\n";
}

static void write_vtk_particles(const std::string& path,
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

static void write_particle_state_txt(const std::string& path,
                                     const std::vector<Particle>& particles,
                                     const std::vector<Mesh>& meshes,
                                     const std::vector<ParticleInit>& inits) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to write particle state: " << path << "\n";
        return;
    }
    out.setf(std::ios::fixed);
    out << std::setprecision(9);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const Particle& p = particles[i];
        const ParticleInit& init = inits[i];
        out << "particle\n";
        out << "stl = " << (init.stl_path.empty() ? "" : init.stl_path) << "\n";
        out << "pos = " << p.tf.pos.x << " " << p.tf.pos.y << " " << p.tf.pos.z << "\n";
        out << "vel = " << p.vel.x << " " << p.vel.y << " " << p.vel.z << "\n";
        out << "quat = " << p.tf.rot.w << " " << p.tf.rot.x << " " << p.tf.rot.y << " " << p.tf.rot.z << "\n";
        out << "omega = " << p.omega.x << " " << p.omega.y << " " << p.omega.z << "\n";
        out << "scale = " << init.scale << "\n";
        out << "density = " << init.density << "\n";
        out << "young = " << init.young << "\n";
        out << "poisson = " << init.poisson << "\n";
        out << "mu = " << init.mu << "\n";
        out << "restitution = " << init.restitution << "\n";
        out << "end_particle\n\n";
    }
}

int main(int argc, char** argv) {
    SimConfig cfg;
    std::string config_path;
    if (!parse_args(argc, argv, config_path)) {
        usage();
        return 1;
    }
    if (!parse_config_file(config_path, cfg)) {
        return 1;
    }

    std::unordered_map<std::string, std::vector<std::array<Vec3, 3>>> tri_cache;
    std::unordered_map<std::string, int> mesh_cache;
    std::vector<Mesh> meshes;

    auto load_mesh_for = [&](const std::string& path, double scale) -> int {
        std::string key = path + "|" + std::to_string(scale);
        auto it_cache = mesh_cache.find(key);
        if (it_cache != mesh_cache.end()) {
            return it_cache->second;
        }

        auto it = tri_cache.find(path);
        if (it == tri_cache.end()) {
            std::vector<std::array<Vec3, 3>> tris;
            if (!load_stl(path, tris)) {
                return -1;
            }
            tri_cache[path] = tris;
            it = tri_cache.find(path);
        }
        std::vector<std::array<Vec3, 3>> tris = it->second;
        for (auto& tri : tris) {
            tri[0] *= scale;
            tri[1] *= scale;
            tri[2] *= scale;
        }
        Mesh m = build_mesh(tris, cfg.center_mesh);
        meshes.push_back(m);
        int idx = static_cast<int>(meshes.size() - 1);
        mesh_cache[key] = idx;
        return idx;
    };

    std::vector<ParticleInit> inits;
    if (!cfg.particle_inits.empty()) {
        inits = cfg.particle_inits;
    } else {
        inits.resize(cfg.n);
        for (int i = 0; i < cfg.n; ++i) {
            inits[i].stl_path = cfg.stl_path;
            inits[i].pos = Vec3{cfg.spacing * i, 0.0, 0.0};
            inits[i].vel = (i == 0) ? cfg.v0 : (i == 1 ? cfg.v1 : Vec3{0.0, 0.0, 0.0});
            inits[i].rot = Quat{};
            inits[i].omega = Vec3{0.0, 0.0, 0.0};
            inits[i].scale = 1.0;
            inits[i].density = 1.0;
            inits[i].young = 1e7;
            inits[i].poisson = 0.25;
            inits[i].mu = 0.5;
            inits[i].restitution = 0.5;
        }
    }

    std::vector<Particle> particles(inits.size());
    for (std::size_t i = 0; i < inits.size(); ++i) {
        const auto& init = inits[i];
        std::string stl_path = init.stl_path.empty() ? cfg.stl_path : init.stl_path;
        if (stl_path.empty()) {
            std::cerr << "Missing stl for particle " << i << "\n";
            return 1;
        }
        double scale = init.scale > 0.0 ? init.scale : 1.0;
        int mesh_index = load_mesh_for(stl_path, scale);
        if (mesh_index < 0) {
            std::cerr << "Failed to load STL: " << stl_path << "\n";
            return 1;
        }
        const Mesh& mesh = meshes[mesh_index];
        if (mesh.vertices.empty()) {
            std::cerr << "Mesh has no vertices.\n";
            return 1;
        }
        Particle p;
        p.tf.pos = init.pos;
        p.tf.rot = init.rot;
        p.vel = init.vel;
        p.omega = init.omega;
        if (mesh.volume > 0.0) {
            double density = init.density > 0.0 ? init.density : 1.0;
            p.mass = density * mesh.volume;
            p.inv_mass = 1.0 / p.mass;
            p.inertia_body = mat3_scale(mesh.inertia_unit, density);
            p.inertia_body_inv = mat3_inverse(p.inertia_body);
        } else {
            double density = init.density > 0.0 ? init.density : 1.0;
            p.mass = density;
            p.inv_mass = 1.0 / p.mass;
            double I = 0.4 * p.mass * mesh.radius * mesh.radius;
            Mat3 Ibody = mat3_identity();
            Ibody.m[0][0] = I;
            Ibody.m[1][1] = I;
            Ibody.m[2][2] = I;
            p.inertia_body = Ibody;
            p.inertia_body_inv = mat3_inverse(Ibody);
        }
        p.radius = mesh.radius;
        p.equiv_radius = mesh.volume > 0.0 ? std::cbrt(3.0 * mesh.volume / (4.0 * 3.141592653589793)) : mesh.radius;
        p.young = init.young;
        p.poisson = init.poisson;
        p.mu = init.mu;
        p.restitution = init.restitution;
        p.mesh_index = mesh_index;
        p.L = mat3_mul_vec3(inertia_world(p), p.omega);
        particles[i] = p;
    }

    // Initialize walls
    std::vector<Wall> walls;
    for (const auto& wi : cfg.wall_inits) {
        if (wi.stl_path.empty()) {
            std::cerr << "Wall missing stl path.\n";
            continue;
        }
        std::vector<std::array<Vec3, 3>> wall_tris;
        if (!load_stl(wi.stl_path, wall_tris)) {
            std::cerr << "Failed to load wall STL: " << wi.stl_path << "\n";
            continue;
        }
        double scale = wi.scale > 0.0 ? wi.scale : 1.0;
        for (auto& tri : wall_tris) {
            tri[0] = tri[0] * scale;
            tri[1] = tri[1] * scale;
            tri[2] = tri[2] * scale;
        }
        Wall wall;
        wall.mesh = build_mesh(wall_tris, cfg.center_mesh);
        wall.tf.pos = wi.pos;
        wall.tf.rot = wi.rot;
        wall.mu = wi.mu;
        wall.restitution = wi.restitution;

        // Precompute wall triangle normals in world space
        std::vector<std::array<Vec3, 3>> world_tris = transform_tris(wall.mesh, wall.tf);
        wall.tri_normals.resize(world_tris.size());
        for (std::size_t t = 0; t < world_tris.size(); ++t) {
            wall.tri_normals[t] = tri_normal(world_tris[t]);
        }
        walls.push_back(wall);
        std::cout << "Loaded wall: " << wi.stl_path << " (tris=" << wall.mesh.tris.size() << ")\n";
    }

    std::vector<Vec3> forces(particles.size());
    std::vector<Vec3> torques(particles.size());
    std::vector<int> contact_counts(particles.size());

    std::filesystem::create_directories(cfg.output_dir);

    for (int step = 0; step < cfg.steps; ++step) {
        auto step_t0 = std::chrono::steady_clock::now();
        std::fill(forces.begin(), forces.end(), Vec3{});
        std::fill(torques.begin(), torques.end(), Vec3{});
        std::fill(contact_counts.begin(), contact_counts.end(), 0);

        int contacts = 0;
        std::size_t total_segments = 0;
        std::size_t total_loops = 0;
        std::size_t total_comps = 0;
        std::size_t tri_checks = 0;
        auto t_tri = std::chrono::steady_clock::duration::zero();
        auto t_rebuild = std::chrono::steady_clock::duration::zero();
        auto t_output = std::chrono::steady_clock::duration::zero();
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(particles.size()); ++j) {
                Particle& pa = particles[i];
                Particle& pb = particles[j];
                long long pair_key = (static_cast<long long>(i) << 32) | static_cast<unsigned long long>(j);
                Vec3 dpos = pb.tf.pos - pa.tf.pos;
                double dist2 = dot(dpos, dpos);
                double rsum = pa.radius + pb.radius;
                if (dist2 > rsum * rsum) {
                    continue;
                }

                auto t0 = std::chrono::steady_clock::now();
                const Mesh& meshA = meshes[pa.mesh_index];
                const Mesh& meshB = meshes[pb.mesh_index];
                std::vector<std::array<Vec3, 3>> trisA = transform_tris(meshA, pa.tf);
                std::vector<std::array<Vec3, 3>> trisB = transform_tris(meshB, pb.tf);

                std::vector<Vec3> nA_all(trisA.size());
                std::vector<Vec3> nB_all(trisB.size());
                std::vector<std::pair<Vec3, double>> planesA(trisA.size());
                std::vector<std::pair<Vec3, double>> planesB(trisB.size());

                for (std::size_t ta = 0; ta < trisA.size(); ++ta) {
                    Vec3 n; double d;
                    if (!plane_from_tri(trisA[ta], n, d)) {
                        nA_all[ta] = Vec3{0.0, 0.0, 0.0};
                        planesA[ta] = {Vec3{0.0, 0.0, 0.0}, 0.0};
                    } else {
                        nA_all[ta] = n;
                        planesA[ta] = {n, d};
                    }
                }
                for (std::size_t tb = 0; tb < trisB.size(); ++tb) {
                    Vec3 n; double d;
                    if (!plane_from_tri(trisB[tb], n, d)) {
                        nB_all[tb] = Vec3{0.0, 0.0, 0.0};
                        planesB[tb] = {Vec3{0.0, 0.0, 0.0}, 0.0};
                    } else {
                        nB_all[tb] = n;
                        planesB[tb] = {n, d};
                    }
                }

                std::vector<Vec3> aabbA_min(trisA.size());
                std::vector<Vec3> aabbA_max(trisA.size());
                std::vector<Vec3> aabbB_min(trisB.size());
                std::vector<Vec3> aabbB_max(trisB.size());
                for (std::size_t ta = 0; ta < trisA.size(); ++ta) {
                    Vec3 mn{std::min({trisA[ta][0].x, trisA[ta][1].x, trisA[ta][2].x}),
                            std::min({trisA[ta][0].y, trisA[ta][1].y, trisA[ta][2].y}),
                            std::min({trisA[ta][0].z, trisA[ta][1].z, trisA[ta][2].z})};
                    Vec3 mx{std::max({trisA[ta][0].x, trisA[ta][1].x, trisA[ta][2].x}),
                            std::max({trisA[ta][0].y, trisA[ta][1].y, trisA[ta][2].y}),
                            std::max({trisA[ta][0].z, trisA[ta][1].z, trisA[ta][2].z})};
                    aabbA_min[ta] = mn;
                    aabbA_max[ta] = mx;
                }
                for (std::size_t tb = 0; tb < trisB.size(); ++tb) {
                    Vec3 mn{std::min({trisB[tb][0].x, trisB[tb][1].x, trisB[tb][2].x}),
                            std::min({trisB[tb][0].y, trisB[tb][1].y, trisB[tb][2].y}),
                            std::min({trisB[tb][0].z, trisB[tb][1].z, trisB[tb][2].z})};
                    Vec3 mx{std::max({trisB[tb][0].x, trisB[tb][1].x, trisB[tb][2].x}),
                            std::max({trisB[tb][0].y, trisB[tb][1].y, trisB[tb][2].y}),
                            std::max({trisB[tb][0].z, trisB[tb][1].z, trisB[tb][2].z})};
                    aabbB_min[tb] = mn;
                    aabbB_max[tb] = mx;
                }

                double tol = meshA.mean_edge > 0.0 ? (meshA.mean_edge * 0.1) : (meshA.bbox_diag * 1e-2);
                if (tol <= 0.0) {
                    tol = 1e-6;
                }

                std::vector<std::pair<Vec3, Vec3>> segments;
                for (std::size_t ta = 0; ta < trisA.size(); ++ta) {
                    if (norm2(nA_all[ta]) < 1e-30) {
                        continue;
                    }
                    for (std::size_t tb = 0; tb < trisB.size(); ++tb) {
                        tri_checks++;
                        if (aabbA_max[ta].x < aabbB_min[tb].x || aabbA_min[ta].x > aabbB_max[tb].x ||
                            aabbA_max[ta].y < aabbB_min[tb].y || aabbA_min[ta].y > aabbB_max[tb].y ||
                            aabbA_max[ta].z < aabbB_min[tb].z || aabbA_min[ta].z > aabbB_max[tb].z) {
                            continue;
                        }
                        if (norm2(nB_all[tb]) < 1e-30) {
                            continue;
                        }
                        Vec3 p0, p1;
                        if (tri_tri_intersection_segment(trisA[ta], planesA[ta].first, planesA[ta].second,
                                                         trisB[tb], planesB[tb].first, planesB[tb].second,
                                                         p0, p1)) {
                            segments.emplace_back(p0, p1);
                        }
                    }
                }
                t_tri += std::chrono::steady_clock::now() - t0;

                // ============== Containment scan ==============
                // Vertices that have crossed INTO the other mesh: the earliest
                // contact signal for concave shapes. The segment pipeline stays
                // blind while a vertex slides through the other body's bay -- the
                // first loop used to appear with the vertex already 30 mm deep,
                // and loading the full Hertz force on that overlap teleported
                // (2/5)*kh*pen^2.5 ~ 57 J of potential energy into a 0.1 J impact.
                // Each hit records (vertex, closest surface point, depth): that
                // pair of material points is one penalty contact element.
                std::vector<Vec3> incA, incB;          // contained vertices
                std::vector<Vec3> incA_cp, incB_cp;    // their closest points on the other surface
                std::vector<double> incA_d, incB_d;    // depths |vertex - closest|
                {
                    Vec3 bmn = trisB[0][0], bmx = trisB[0][0];
                    for (const auto& t : trisB) {
                        for (const auto& c : t) {
                            bmn.x = std::min(bmn.x, c.x); bmn.y = std::min(bmn.y, c.y); bmn.z = std::min(bmn.z, c.z);
                            bmx.x = std::max(bmx.x, c.x); bmx.y = std::max(bmx.y, c.y); bmx.z = std::max(bmx.z, c.z);
                        }
                    }
                    for (const auto& v : meshA.vertices) {
                        Vec3 wv = quat_rotate(pa.tf.rot, v) + pa.tf.pos;
                        if (wv.x < bmn.x || wv.x > bmx.x || wv.y < bmn.y || wv.y > bmx.y || wv.z < bmn.z || wv.z > bmx.z) continue;
                        if (!point_inside_mesh(trisB, wv)) continue;
                        Vec3 cp, nf; double d = point_mesh_distance(trisB, wv, cp, nf);
                        if (d <= 1e-12) continue;
                        incA.push_back(wv); incA_cp.push_back(cp); incA_d.push_back(d);
                    }
                    Vec3 amn = trisA[0][0], amx = trisA[0][0];
                    for (const auto& t : trisA) {
                        for (const auto& c : t) {
                            amn.x = std::min(amn.x, c.x); amn.y = std::min(amn.y, c.y); amn.z = std::min(amn.z, c.z);
                            amx.x = std::max(amx.x, c.x); amx.y = std::max(amx.y, c.y); amx.z = std::max(amx.z, c.z);
                        }
                    }
                    for (const auto& v : meshB.vertices) {
                        Vec3 wv = quat_rotate(pb.tf.rot, v) + pb.tf.pos;
                        if (wv.x < amn.x || wv.x > amx.x || wv.y < amn.y || wv.y > amx.y || wv.z < amn.z || wv.z > amx.z) continue;
                        if (!point_inside_mesh(trisA, wv)) continue;
                        Vec3 cp, nf; double d = point_mesh_distance(trisA, wv, cp, nf);
                        if (d <= 1e-12) continue;
                        incB.push_back(wv); incB_cp.push_back(cp); incB_d.push_back(d);
                    }
                }

                if (!incA.empty() || !incB.empty()) {
                    // ============== Per-vertex penalty contact ==============
                    // Each contained vertex and its closest point on the other
                    // surface form ONE contact element: the spring length is the
                    // depth d, and the force acts on exactly the two material
                    // points whose separation defines d. The spring's rate of
                    // change therefore equals the relative velocity of the very
                    // points the force acts on -- dU/dt = F.v holds by
                    // construction. Single-point forms (loop contact point,
                    // weighted average, deepest vertex) all mismatched the
                    // driving quantity's time derivative against the force's
                    // application point under spin and injected up to kJ.
                    // The push direction is -(vertex - closest)/d: winding-free
                    // (this STL's face windings are NOT consistent, 19/102 faces
                    // point inward), always from the surface point away from the
                    // interior point, i.e. the separation direction.
                    double Rp1 = std::max(pa.equiv_radius, 1e-12);
                    double Rp2 = std::max(pb.equiv_radius, 1e-12);
                    double E1r = std::max(pa.young, 1e-12);
                    double E2r = std::max(pb.young, 1e-12);
                    double Estar_r = 1.0 / ((1.0 - pa.poisson * pa.poisson) / E1r + (1.0 - pb.poisson * pb.poisson) / E2r);
                    double Rstar_r = (Rp1 * Rp2) / std::max(Rp1 + Rp2, 1e-12);
                    double k_hr = (4.0 / 3.0) * Estar_r * std::sqrt(Rstar_r);
                    double kt_r = 0.5 * ((E1r * Rp1) / (2.0 * (1.0 + pa.poisson)) + (E2r * Rp2) / (2.0 * (1.0 + pb.poisson)));
                    double m_eff_r = 1.0 / (pa.inv_mass + pb.inv_mass);
                    int n_inc = static_cast<int>(incA.size() + incB.size());
                    double m_eff_v = m_eff_r / std::max(n_inc, 1);  // each vertex carries its mass share
                    double e_r = std::min(pa.restitution, pb.restitution);
                    e_r = std::min(0.9999, std::max(1e-6, e_r));
                    double loge_r = std::log(e_r);
                    double pi2 = 3.141592653589793 * 3.141592653589793;
                    double ct_v = cfg.tangential_damping * std::sqrt(kt_r * m_eff_v);
                    double mu_r = std::min(pa.mu, pb.mu);
                    auto apply_penalty = [&](const Vec3& wv, const Vec3& cp, double d,
                                             const Particle& p_own, const Particle& p_other,
                                             int i_own, int i_other) {
                        Vec3 u = (wv - cp) * (1.0 / std::max(d, 1e-15));  // surface -> interior
                        Vec3 n_push = u * -1.0;                            // separation direction
                        Vec3 v_own = p_own.vel + cross(p_own.omega, wv - p_own.tf.pos);
                        Vec3 v_oth = p_other.vel + cross(p_other.omega, cp - p_other.tf.pos);
                        Vec3 vrel = v_own - v_oth;
                        double vn = dot(vrel, n_push);  // <0 = penetrating deeper
                        double kn_eff_v = 1.5 * k_hr * std::sqrt(std::max(d, 1e-12));
                        double cn_v = std::abs(2.0 * loge_r * std::sqrt(kn_eff_v * m_eff_v) /
                            std::sqrt(loge_r * loge_r + pi2));
                        double fn_v = k_hr * d * std::sqrt(d);
                        if (vn < 0.0) fn_v += -cn_v * vn;
                        if (fn_v < 0.0) fn_v = 0.0;
                        Vec3 fn_vec = n_push * fn_v;
                        forces[i_own] += fn_vec;
                        forces[i_other] -= fn_vec;
                        torques[i_own] += cross(wv - p_own.tf.pos, fn_vec);
                        torques[i_other] += cross(cp - p_other.tf.pos, fn_vec * -1.0);
                        // Tangential friction at the same pair of points
                        Vec3 vt = vrel - n_push * vn;
                        Vec3 ft_v = vt * (-ct_v);
                        double ftv_norm = norm(ft_v);
                        double ftv_max = mu_r * fn_v;
                        if (ftv_norm > ftv_max && ftv_norm > 1e-14) {
                            ft_v = ft_v * (ftv_max / ftv_norm);
                        }
                        forces[i_own] += ft_v;
                        forces[i_other] -= ft_v;
                        torques[i_own] += cross(wv - p_own.tf.pos, ft_v);
                        torques[i_other] += cross(cp - p_other.tf.pos, ft_v * -1.0);
                    };
                    for (std::size_t k = 0; k < incA.size(); ++k) {
                        apply_penalty(incA[k], incA_cp[k], incA_d[k], pa, pb, i, j);
                    }
                    for (std::size_t k = 0; k < incB.size(); ++k) {
                        apply_penalty(incB[k], incB_cp[k], incB_d[k], pb, pa, j, i);
                    }
                    contact_counts[i] += 1;
                    contact_counts[j] += 1;
                    contacts++;
                    // Containment fully covers this pair; skip the loop pipeline.
                    continue;
                }

                if (segments.empty()) {
                    continue;
                }

                t0 = std::chrono::steady_clock::now();
                std::vector<std::vector<int>> comps;
                if (cfg.split_contacts) {
                    label_segments_by_contact_greedy(segments, tol, comps);
                } else {
                    comps.resize(1);
                    comps[0].reserve(segments.size());
                    for (int s = 0; s < static_cast<int>(segments.size()); ++s) {
                        comps[0].push_back(s);
                    }
                }
                total_segments += segments.size();
                total_comps += comps.size();

                // Loop tracing (telemetry only): the per-vertex containment
                // penalty above fully covers this pair's contact force.
                // Edge-edge touches before any vertex actually crosses have a
                // brief force-free grace period (the vertex crosses within
                // ~a mean edge of travel).
                for (const auto& comp : comps) {
                    std::vector<std::pair<Vec3, Vec3>> segs;
                    segs.reserve(comp.size());
                    for (int idx : comp) {
                        segs.push_back(segments[idx]);
                    }
                    std::vector<Vec3> V;
                    std::vector<std::pair<int, int>> edges;
                    snap_endpoints(segs, tol, V, edges);
                    std::vector<Vec3> V2;
                    std::vector<std::pair<int, int>> edges2;
                    split_t_junctions(V, edges, tol, V2, edges2);
                    std::vector<std::vector<int>> loops_vids = extract_loops(V2, edges2);
                    total_loops += loops_vids.size();
                }
                t_rebuild += std::chrono::steady_clock::now() - t0;
            }
        }

        // ============== Particle-Wall Contact Detection ==============
        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            Particle& p = particles[i];
            const Mesh& pmesh = meshes[p.mesh_index];

            for (int w = 0; w < static_cast<int>(walls.size()); ++w) {
                const Wall& wall = walls[w];

                // Broad phase: check if particle bounding sphere overlaps with wall bounding sphere
                Vec3 dpos = p.tf.pos - wall.tf.pos;
                double dist2 = dot(dpos, dpos);
                double rsum = p.radius + wall.mesh.radius;
                if (dist2 > rsum * rsum) {
                    continue;
                }

                std::vector<std::array<Vec3, 3>> wall_tris = transform_tris(wall.mesh, wall.tf);

                double tol = pmesh.mean_edge > 0.0 ? (pmesh.mean_edge * 0.1) : (pmesh.bbox_diag * 1e-2);
                if (tol <= 0.0) tol = 1e-6;

                // ============== Wall contact: per-vertex penalty (v8) ==============
                // Same contact element as particle-particle: every vertex that
                // crossed a wall plane is one penalty spring acting exactly at
                // that vertex. Coplanar wall triangles passing the AABB test
                // are grouped so a vertex in the shared footprint is charged
                // once (a split plane wall is TWO coplanar triangles -- without
                // grouping every contact force would be doubled). The old
                // clip-loop pipeline drove the force by the deepest vertex but
                // applied it at a depth-weighted average point -- different
                // material points, the spin-mismatch injection fixed in the
                // pair contact; after the integrator stopped dissipating
                // numerically it left the e=0.3 wall case oscillating forever.
                struct WallPlaneGroup {
                    Vec3 n; double d;
                    std::vector<std::array<Vec3, 3>> footprints;
                };
                std::vector<WallPlaneGroup> groups;
                for (std::size_t t = 0; t < wall_tris.size(); ++t) {
                    const auto& wtri = wall_tris[t];
                    Vec3 wmn{std::min({wtri[0].x, wtri[1].x, wtri[2].x}),
                             std::min({wtri[0].y, wtri[1].y, wtri[2].y}),
                             std::min({wtri[0].z, wtri[1].z, wtri[2].z})};
                    Vec3 wmx{std::max({wtri[0].x, wtri[1].x, wtri[2].x}),
                             std::max({wtri[0].y, wtri[1].y, wtri[2].y}),
                             std::max({wtri[0].z, wtri[1].z, wtri[2].z})};
                    Vec3 pmn = p.tf.pos - Vec3{p.radius, p.radius, p.radius};
                    Vec3 pmx = p.tf.pos + Vec3{p.radius, p.radius, p.radius};
                    if (pmx.x < wmn.x || pmn.x > wmx.x ||
                        pmx.y < wmn.y || pmn.y > wmx.y ||
                        pmx.z < wmn.z || pmn.z > wmx.z) {
                        continue;
                    }
                    Vec3 pn = tri_normal(wtri);
                    if (dot(pn, p.tf.pos - wtri[0]) < 0.0) {
                        pn = pn * -1.0;  // push toward the particle side
                    }
                    double pd = -dot(pn, wtri[0]);
                    const double pq = 1e6;
                    bool merged = false;
                    for (auto& g : groups) {
                        if (std::llround(g.n.x * pq) == std::llround(pn.x * pq) &&
                            std::llround(g.n.y * pq) == std::llround(pn.y * pq) &&
                            std::llround(g.n.z * pq) == std::llround(pn.z * pq) &&
                            std::llround(g.d * pq) == std::llround(pd * pq)) {
                            g.footprints.push_back(wtri);
                            merged = true;
                            break;
                        }
                    }
                    if (!merged) {
                        WallPlaneGroup g;
                        g.n = pn; g.d = pd;
                        g.footprints.push_back(wtri);
                        groups.push_back(g);
                    }
                }

                double Rp = std::max(p.equiv_radius, 1e-12);
                double Ep = std::max(p.young, 1e-12);
                double nu_p = p.poisson;
                double Estar = Ep / (2.0 * (1.0 - nu_p * nu_p));
                double k_hw = (4.0 / 3.0) * Estar * std::sqrt(Rp);
                double kt_w = (Ep * Rp) / (2.0 * (1.0 + nu_p));
                double e_w = std::min(0.9999, std::max(1e-6, p.restitution));
                double loge_w = std::log(e_w);
                double pi2 = 3.141592653589793 * 3.141592653589793;
                double mu_w = std::min(p.mu, wall.mu);

                for (const auto& g : groups) {
                    // Count crossings first so the damping uses each vertex share
                    // of the particle mass.
                    int n_cross = 0;
                    for (const auto& v : pmesh.vertices) {
                        Vec3 wv = quat_rotate(p.tf.rot, v) + p.tf.pos;
                        double s = dot(g.n, wv) + g.d;
                        if (s > -1e-12) continue;
                        Vec3 p_proj = wv - g.n * s;  // drop onto the wall plane
                        for (const auto& fp : g.footprints) {
                            if (point_in_tri(p_proj, fp, g.n)) { ++n_cross; break; }
                        }
                    }
                    if (n_cross == 0) continue;
                    double m_eff_v = p.mass / static_cast<double>(n_cross);
                    double ct_v = cfg.tangential_damping * std::sqrt(kt_w * m_eff_v);
                    bool any_force = false;
                    for (const auto& v : pmesh.vertices) {
                        Vec3 wv = quat_rotate(p.tf.rot, v) + p.tf.pos;
                        double s = dot(g.n, wv) + g.d;
                        if (s > -1e-12) continue;
                        Vec3 p_proj = wv - g.n * s;
                        bool inside = false;
                        for (const auto& fp : g.footprints) {
                            if (point_in_tri(p_proj, fp, g.n)) { inside = true; break; }
                        }
                        if (!inside) continue;
                        double d_v = -s;
                        Vec3 r_v = wv - p.tf.pos;
                        Vec3 v_v = p.vel + cross(p.omega, r_v);
                        double vn_v = dot(v_v, g.n);  // <0 = penetrating deeper
                        double kn_eff_v = 1.5 * k_hw * std::sqrt(std::max(d_v, 1e-12));
                        double cn_v = std::abs(2.0 * loge_w * std::sqrt(kn_eff_v * m_eff_v) /
                            std::sqrt(loge_w * loge_w + pi2));
                        double fn_v = k_hw * d_v * std::sqrt(d_v);
                        if (vn_v < 0.0) fn_v += -cn_v * vn_v;
                        if (fn_v < 0.0) fn_v = 0.0;
                        Vec3 vt_v = v_v - g.n * vn_v;
                        Vec3 ft_v = vt_v * (-ct_v);
                        double ftv_norm = norm(ft_v);
                        double ftv_max = mu_w * fn_v;
                        if (ftv_norm > ftv_max && ftv_norm > 1e-14) {
                            ft_v = ft_v * (ftv_max / ftv_norm);
                        }
                        Vec3 f_v = g.n * fn_v + ft_v;
                        forces[i] += f_v;
                        torques[i] += cross(r_v, f_v);
                        any_force = true;
                    }
                    if (any_force) {
                        contact_counts[i] += 1;
                        contacts++;
                        if (cfg.contact_debug && contacts < 3) {
                            std::cout << "  wall_contact: particle=" << i << " wall=" << w
                                      << " n_cross=" << n_cross
                                      << " F=" << forces[i].x << "," << forces[i].y << "," << forces[i].z
                                      << "\n";
                        }
                    }
                }
            }
        }

        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            integrate_particle(particles[i], forces[i], torques[i], cfg.gravity, cfg.dt);
        }

        if (step % cfg.output_interval == 0) {
            std::ostringstream oss;
            oss << cfg.vtk_prefix << "_" << std::setw(6) << std::setfill('0') << step << ".vtk";
            std::filesystem::path out_path = std::filesystem::path(cfg.output_dir) / oss.str();
            auto t_out0 = std::chrono::steady_clock::now();
            write_vtk_particles(out_path.string(), meshes, particles, forces, torques, contact_counts);

            // Output particle state to txt
            std::ostringstream oss_txt;
            oss_txt << cfg.vtk_prefix << "_" << std::setw(6) << std::setfill('0') << step << ".txt";
            std::filesystem::path txt_path = std::filesystem::path(cfg.output_dir) / oss_txt.str();
            write_particle_state_txt(txt_path.string(), particles, meshes, inits);

            t_output += std::chrono::steady_clock::now() - t_out0;

            auto step_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - step_t0).count();
            auto tri_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_tri).count();
            auto rebuild_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_rebuild).count();
            auto out_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_output).count();
            std::cout << "step=" << step
                  << " contacts=" << contacts
                  << " step_ms=" << step_ms
                  << " tri_ms=" << tri_ms
                  << " rebuild_ms=" << rebuild_ms
                  << " output_ms=" << out_ms
                  << " tri_checks=" << tri_checks
                  << " segments=" << total_segments
                  << " comps=" << total_comps
                  << " loops=" << total_loops
                  << "\n";
        }
    }

    std::cout << "done\n";
    return 0;
}
