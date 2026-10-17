#include "mex.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

const int MAX_DAXIAO = 32;
const int MAX_NODES = 2 * MAX_DAXIAO;

struct Params {
    int n1, K, daxiao, vmax, n2, num_nodes, threads;
};

struct Row { int u, v, t; };

struct Matrix {
    int data[MAX_DAXIAO][3];
    int rows;
    Matrix() : rows(0) { std::memset(data, 0, sizeof(data)); }
    void init(int n) { rows = n; }
    int get(int r, int c) const { return data[r][c]; }
    void set(int r, int c, int value) { data[r][c] = value; }
};

// A canonical row needs only 13 bits: (u << 7) | (v << 1) | t.
// Unused array entries are zero, so hashing/equality are deterministic.
struct CanonicalKey {
    std::array<std::uint16_t, MAX_DAXIAO> code;
    CanonicalKey() { code.fill(0); }
    bool operator==(const CanonicalKey& b) const { return code == b.code; }
    bool operator<(const CanonicalKey& b) const { return code < b.code; }
};

struct CanonicalHash {
    std::size_t operator()(const CanonicalKey& key) const {
        std::size_t h = static_cast<std::size_t>(1469598103934665603ULL);
        for (std::size_t i = 0; i < key.code.size(); ++i) {
            h ^= static_cast<std::size_t>(key.code[i]);
            h *= static_cast<std::size_t>(1099511628211ULL);
        }
        return h;
    }
};
typedef std::unordered_set<CanonicalKey, CanonicalHash> ResultSet;
typedef std::array<std::vector<Row>, MAX_NODES> CandidateBuckets;

inline int strict_mod(int x, int n) { int r = x % n; return r < 0 ? r + n : r; }
inline int get_next_vertex(int x, int n) {
    return (x / n) * n + ((x % n + 1) % n);
}

struct PathSystem {
    int degree[MAX_NODES], mate[MAX_NODES], size[MAX_NODES];
    bool used_edge[MAX_NODES]; // Polygon edges, distinct from corner degrees.
    int closed_loops_count;
    int excess; // Sum over components of max(component vertex count - 3, 0).
    bool has_intra_t1;
    unsigned cross_t_mask;
    Matrix current_A;

    explicit PathSystem(int n) : closed_loops_count(0), excess(0), has_intra_t1(false), cross_t_mask(0) {
        current_A.init(n);
        std::memset(degree, 0, sizeof(degree));
        std::memset(mate, 0, sizeof(mate));
        std::memset(size, 0, sizeof(size));
        std::memset(used_edge, 0, sizeof(used_edge));
    }

    bool add_segment(int u, int v, int vmax) {
        if (u == v || degree[u] >= 2 || degree[v] >= 2) return false;
        if (degree[u] == 0 && degree[v] == 0) {
            if (vmax < 2) return false;
            degree[u] = degree[v] = 1;
            mate[u] = v; mate[v] = u;
            size[u] = size[v] = 2;
            return true;
        }
        if (degree[u] == 0 && degree[v] == 1) return add_segment(v, u, vmax);
        if (degree[u] == 1 && degree[v] == 0) {
            int z = mate[u], length = size[u] + 1;
            if (length > vmax) return false;
            excess += std::max(length - 3, 0) - std::max(size[u] - 3, 0);
            degree[u] = 2; degree[v] = 1;
            mate[z] = v; mate[v] = z;
            size[z] = size[v] = length;
            return true;
        }
        if (degree[u] == 1 && degree[v] == 1) {
            int z = mate[u], w = mate[v];
            if (z == v) {
                if (size[u] < 3 || size[u] > vmax) return false;
                degree[u] = degree[v] = 2;
                ++closed_loops_count;
            }
            else {
                int length = size[u] + size[v];
                if (length > vmax) return false;
                excess += (length - 3) - std::max(size[u] - 3, 0) - std::max(size[v] - 3, 0);
                degree[u] = degree[v] = 2;
                mate[z] = w; mate[w] = z;
                size[z] = size[w] = length;
            }
            return true;
        }
        return false;
    }

    // Failure may modify this copy; callers always discard a failed copy.
    bool add_row(const Row& row, int row_idx, const Params& p) {
        if (row.u == row.v || used_edge[row.u] || used_edge[row.v]) return false;
        used_edge[row.u] = used_edge[row.v] = true;
        current_A.set(row_idx, 0, row.u);
        current_A.set(row_idx, 1, row.v);
        current_A.set(row_idx, 2, row.t);
        if ((row.u < p.n1) == (row.v < p.n1)) {
            if (row.t == 1) has_intra_t1 = true;
        }
        else {
            cross_t_mask |= (1U << row.t);
        }
        const int un = get_next_vertex(row.u, p.n1);
        const int vn = get_next_vertex(row.v, p.n1);
        if (row.t == 1)
            return add_segment(row.u, row.v, p.vmax) && add_segment(un, vn, p.vmax);
        return add_segment(un, row.v, p.vmax) && add_segment(row.u, vn, p.vmax);
    }
};

// Exactly the synchronous D_n action and optional exchange of the two blocks.
// The same reflection acts on BOTH blocks; t is unchanged. Group order: 4*n.
CanonicalKey get_canonical(const Matrix& a, const Params& p) {
    CanonicalKey best;
    bool initialized = false;
    for (int swap_blocks = 0; swap_blocks < 2; ++swap_blocks) {
        for (int reflect = 0; reflect < 2; ++reflect) {
            for (int d = 0; d < p.n1; ++d) {
                std::uint16_t by_u[MAX_NODES] = {};
                for (int r = 0; r < p.daxiao; ++r) {
                    int x[2] = { a.get(r, 0), a.get(r, 1) };
                    for (int c = 0; c < 2; ++c) {
                        const int block = (x[c] / p.n1) ^ swap_blocks;
                        const int local = x[c] % p.n1;
                        x[c] = block * p.n1 + strict_mod(d + (reflect ? -local : local), p.n1);
                    }
                    if (x[0] > x[1]) std::swap(x[0], x[1]);
                    by_u[x[0]] = static_cast<std::uint16_t>((x[0] << 7) | (x[1] << 1) | a.get(r, 2));
                }
                CanonicalKey candidate;
                int count = 0;
                // A perfect matching uses each edge once, hence distinct first columns.
                for (int u = 0; u < p.num_nodes; ++u)
                    if (by_u[u] != 0) candidate.code[count++] = by_u[u];
                if (!initialized || candidate < best) { best = candidate; initialized = true; }
            }
        }
    }
    return best;
}

void dfs(const PathSystem& state, int row_idx, const CandidateBuckets& buckets,
    ResultSet& results, std::atomic<bool>& abort_search, const Params& p) {
    if (abort_search.load(std::memory_order_relaxed)) return;
    if (state.closed_loops_count > p.n2) return;
    // Component excess cannot decrease under extension, joining or closure.
    // At a valid leaf its value is 2*n - 3*(n-K) = 3*K-n.
    if (state.excess > 3 * p.K - p.n1) return;
    const int remaining = p.daxiao - row_idx;
    // One new segment can close at most one cycle; each row adds two segments.
    if (state.closed_loops_count + 2 * remaining < p.n2) return;
    if (remaining == 0) {
        if (state.closed_loops_count != p.n2) return;
        // Seed guarantees connectedness. Orientation is inconsistent exactly when
        // an intra-block t=1 occurs or cross-block rows contain both t values.
        if (!state.has_intra_t1 && state.cross_t_mask != 3U) return;
        for (int e = 0; e < p.num_nodes; ++e)
            if (!state.used_edge[e] || state.degree[e] != 2) return;
        CanonicalKey key = get_canonical(state.current_A, p);
        // Deduplicate immediately: never retain the full list of raw discoveries.
        std::exception_ptr insertion_error;
#pragma omp critical(two_tile_results)
        {
            // OpenMP requires exceptions to be caught inside this structured block.
            try { results.insert(key); }
            catch (...) { insertion_error = std::current_exception(); }
        }
        if (insertion_error) std::rethrow_exception(insertion_error);
        return;
    }

    int u = 0;
    while (u < p.num_nodes && state.used_edge[u]) ++u;
    if (u == p.num_nodes) return;
    // Every completion must pair this smallest unused edge with a larger unused
    // edge. Branching only on bucket[u] is complete and gives one row order.
    for (std::size_t i = 0; i < buckets[u].size(); ++i) {
        const Row& candidate = buckets[u][i];
        if (state.used_edge[candidate.v]) continue;
        PathSystem next_state = state;
        if (next_state.add_row(candidate, row_idx, p))
            dfs(next_state, row_idx + 1, buckets, results, abort_search, p);
    }
}

int read_integer(const mxArray* value, int low, int high, const char* name) {
    if (!mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value) ||
        mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument(std::string(name) + " must be a full real double scalar.");
    const double x = mxGetScalar(value);
    if (!std::isfinite(x) || x != std::floor(x) || x < low || x > high)
        throw std::invalid_argument(std::string(name) + " must be an integer in [" +
            std::to_string(low) + ", " + std::to_string(high) + "].");
    return static_cast<int>(x);
}

std::vector<Row> read_rows(const mxArray* value, const Params& p, bool seed) {
    const char* name = seed ? "firstE" : "combinations1";
    if (!mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value) ||
        mxGetNumberOfDimensions(value) != 2 || mxGetN(value) != 3)
        throw std::invalid_argument(std::string(name) + " must be a full real double M-by-3 matrix.");
    const mwSize count = mxGetM(value);
    // This also keeps the OpenMP loop index within signed-int range.
    if (count > static_cast<mwSize>(std::numeric_limits<int>::max()))
        throw std::invalid_argument(std::string(name) + " has too many rows.");
    const double* data = mxGetPr(value);
    std::vector<Row> rows;
    rows.reserve(count);
    for (mwSize i = 0; i < count; ++i) {
        const double u = data[i], v = data[i + count], t = data[i + 2 * count];
        if (!std::isfinite(u) || !std::isfinite(v) || !std::isfinite(t) ||
            u != std::floor(u) || v != std::floor(v) || (t != 0 && t != 1) ||
            u < 0 || v >= p.num_nodes || u >= v)
            throw std::invalid_argument(std::string(name) + " rows require integer 0 <= u < v < 2*n and t in {0,1}.");
        if (seed && (u != 0 || v < p.n1))
            throw std::invalid_argument("Each firstE row must pair edge 0 with an edge of the other block.");
        if (!seed && u == 0)
            throw std::invalid_argument("combinations1 must exclude edge 0, which is already used by the seed.");
        rows.push_back(Row{ static_cast<int>(u), static_cast<int>(v), static_cast<int>(t) });
    }
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        return a.u != b.u ? a.u < b.u : (a.v != b.v ? a.v < b.v : a.t < b.t);
        });
    rows.erase(std::unique(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        return a.u == b.u && a.v == b.v && a.t == b.t;
        }), rows.end());
    return rows;
}

mxArray* solve(int nrhs, const mxArray* prhs[]) {
    if (nrhs != 7 && nrhs != 8)
        throw std::invalid_argument("Expected 7 inputs: n, K, daxiao, vmax, firstE, combinations1, n2; optional 8th input: Threads.");
    Params p;
    p.n1 = read_integer(prhs[0], 3, MAX_DAXIAO, "n");
    p.K = read_integer(prhs[1], 1, p.n1 - 1, "K");
    p.daxiao = read_integer(prhs[2], p.n1, p.n1, "daxiao");
    p.vmax = read_integer(prhs[3], 3 * p.K + 3 - p.n1, 3 * p.K + 3 - p.n1, "vmax");
    p.n2 = read_integer(prhs[6], p.n1 - p.K, p.n1 - p.K, "n2");
    p.num_nodes = 2 * p.n1;
    p.threads = 1;
#ifdef _OPENMP
    const int default_threads = std::max(1, static_cast<int>(0.4 * omp_get_num_procs() + 0.5));
    p.threads = std::max(1, std::min(omp_get_max_threads(), default_threads));
#endif
    if (nrhs == 8) p.threads = read_integer(prhs[7], 1, 1024, "Threads");
    std::vector<Row> seeds = read_rows(prhs[4], p, true);
    std::vector<Row> combinations = read_rows(prhs[5], p, false);
    CandidateBuckets buckets;
    for (std::size_t i = 0; i < combinations.size(); ++i)
        buckets[combinations[i].u].push_back(combinations[i]);
    p.threads = std::max(1, std::min(p.threads, static_cast<int>(seeds.size())));
    ResultSet results;
    std::atomic<bool> abort_search(false);
    std::exception_ptr worker_error;
    const bool feasible = p.vmax >= 3 && 3 * p.n2 <= p.num_nodes;
    if (feasible) {
        const int seed_count = static_cast<int>(seeds.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(p.threads)
        for (int i = 0; i < seed_count; ++i) {
            if (abort_search.load(std::memory_order_relaxed)) continue;
            try {
                PathSystem state(p.daxiao);
                if (state.add_row(seeds[i], 0, p))
                    dfs(state, 1, buckets, results, abort_search, p);
            }
            catch (...) {
                abort_search.store(true, std::memory_order_relaxed);
#pragma omp critical(two_tile_error)
                { if (!worker_error) worker_error = std::current_exception(); }
            }
        }
    }
    if (worker_error) std::rethrow_exception(worker_error);

    // Only pointers are sorted; the compact canonical keys are not duplicated.
    std::vector<const CanonicalKey*> ordered;
    ordered.reserve(results.size());
    for (ResultSet::const_iterator it = results.begin(); it != results.end(); ++it)
        ordered.push_back(&*it);
    std::sort(ordered.begin(), ordered.end(), [](const CanonicalKey* a, const CanonicalKey* b) {
        return *a < *b;
        });
    const mwSize count = static_cast<mwSize>(ordered.size());
    const mwSize per_solution = static_cast<mwSize>(3 * p.daxiao);
    if (count > std::numeric_limits<mwSize>::max() / per_solution / sizeof(double))
        throw std::length_error("Output array size overflows the addressable size.");
    const mwSize dims[3] = { static_cast<mwSize>(p.daxiao), 3, count };
    mxArray* output = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    if (!output) throw std::bad_alloc();
    double* data = mxGetPr(output);
    for (mwSize k = 0; k < count; ++k) {
        for (int r = 0; r < p.daxiao; ++r) {
            const unsigned code = ordered[k]->code[r];
            data[k * per_solution + r] = static_cast<double>(code >> 7);
            data[k * per_solution + p.daxiao + r] = static_cast<double>((code >> 1) & 63U);
            data[k * per_solution + 2 * p.daxiao + r] = static_cast<double>(code & 1U);
        }
    }
    return output;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    try {
        if (nlhs > 1) throw std::invalid_argument("At most one output is supported.");
        mxArray* output = solve(nrhs, prhs);
        if (nlhs == 0) mxDestroyArray(output);
        else plhs[0] = output;
    }
    catch (const std::bad_alloc&) {
        mexErrMsgTxt("Insufficient memory for the distinct solutions and MATLAB output array.");
    }
    catch (const std::exception& error) {
        mexErrMsgTxt(error.what());
    }
    catch (...) {
        mexErrMsgTxt("Unexpected C++ exception during enumeration.");
    }
}
