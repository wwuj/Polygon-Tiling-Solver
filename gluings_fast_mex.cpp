//cpp端并行去重
#include "mex.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_set> // 改用哈希集合
#include <set>
#include <omp.h>
#include <cmath>

// ================== 辅助工具: MATLAB 风格取模 ==================
inline int matlab_mod(int a, int n) {
    int r = a % n;
    return r < 0 ? r + n : r;
}

// ================== 数据结构: 向量哈希函数 ==================
// 为了让 unordered_set 支持 vector<int>
struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        size_t seed = 0;
        for (int i : v) {
            // 类似于 boost::hash_combine
            seed ^= std::hash<int>{}(i)+0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// ================== 数据结构: DSU (并查集) ==================
struct RollbackDSU {
    std::vector<int> parent;
    std::vector<int> sz;
    struct History { int child; int parent; bool merged; };

    RollbackDSU(int n) {
        parent.resize(n);
        std::iota(parent.begin(), parent.end(), 0);
        sz.assign(n, 1);
    }

    // 非压缩查找，支持回滚
    int find(int x) const {
        while (parent[x] != x) x = parent[x];
        return x;
    }

    std::pair<History, bool> unite(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX != rootY) {
            if (sz[rootX] < sz[rootY]) std::swap(rootX, rootY);
            parent[rootY] = rootX;
            sz[rootX] += sz[rootY];
            return { {rootY, rootX, true}, false };
        }
        return { {0, 0, false}, true };
    }

    void rollback(const History& h) {
        if (h.merged) {
            sz[h.parent] -= sz[h.child];
            parent[h.child] = h.child;
        }
    }

    int get_size(int x) const { return sz[find(x)]; }
};

// ================== 核心去重逻辑 (高度优化内存版) ==================

// 1. 规范化矩阵 (原地修改)
inline void normalize_matrix(std::vector<int>& M) {
    int rows = M.size() / 2;
    // 行内排序
    for (int i = 0; i < rows; ++i) {
        int idx1 = 2 * i;
        int idx2 = 2 * i + 1;
        if (M[idx1] > M[idx2]) {
            std::swap(M[idx1], M[idx2]);
        }
    }

    // 行间排序 (自定义简单的插入排序或冒泡，避免 vector<pair> 的分配开销)
    // 由于 n1 很小 (10-14)，简单的排序比分配堆内存快得多
    for (int i = 1; i < rows; ++i) {
        int j = i;
        // 保存当前行
        int u_curr = M[2 * j];
        int v_curr = M[2 * j + 1];

        // 向前比较 u
        while (j > 0 && M[2 * (j - 1)] > u_curr) {
            M[2 * j] = M[2 * (j - 1)];
            M[2 * j + 1] = M[2 * (j - 1) + 1];
            j--;
        }
        // 如果 u 相等，比较 v (实际上规范化后 u 很难相等，除非重边，但为保险起见)
        while (j > 0 && M[2 * (j - 1)] == u_curr && M[2 * (j - 1) + 1] > v_curr) {
            M[2 * j] = M[2 * (j - 1)];
            M[2 * j + 1] = M[2 * (j - 1) + 1];
            j--;
        }
        M[2 * j] = u_curr;
        M[2 * j + 1] = v_curr;
    }
}

// 2. 变换函数 (写入 dst)
inline void apply_transform(const std::vector<int>& src, std::vector<int>& dst, int t, int n1, int mode) {
    // 假设 dst 已经 resize 好了
    size_t len = src.size();
    for (size_t i = 0; i < len; ++i) {
        int val = src[i];
        bool is_block0 = (val < n1);
        int base_val = is_block0 ? val : val - n1;
        int new_base;

        if (mode == 0) { // Add t
            new_base = base_val + t;
            if (new_base >= n1) new_base -= n1; // 简单的取模优化
        }
        else { // Sub from t
            new_base = t - base_val;
            if (new_base < 0) new_base += n1;
        }
        dst[i] = is_block0 ? new_base : (new_base + n1);
    }
}

// 3. 评分函数
inline std::pair<long long, int> evaluate_score(std::vector<int>& M, int n1) {
    long long total_len = 0;
    int rows = M.size() / 2;

    for (int i = 0; i < rows; ++i) {
        int u = M[2 * i];
        int v = M[2 * i + 1];
        int diff = std::abs((u % n1) - (v % n1));
        int dist = std::min(diff, n1 - diff);
        total_len += (long long)dist * dist;
    }

    // 对称性检查：
    // 为了极致性能，这里不再调用 transform 和 normalize 
    // 而是做简单的快速检查，或者只在 total_len 相等时才做
    // 这里保持原逻辑，但需注意 M_prime 内存复用

    // *性能妥协*：完全计算对称性开销很大。
    // 由于我们在寻找"最小表示"，对称性只是第二判据。
    // 我们可以直接用 normalize 后的 M 的字典序作为最终判据，
    // 对称性计数仅用于调试或特殊数学需求。
    // 如果题目只要求去重，不需要知道对称性阶数，我们可以简化评分：
    // Score = {ModLength, Lexicographical}
    // 
    // 为了保持与原始逻辑一致，我们尽量优化内存。
    int sym_count = 0;
    // (此处省略复杂的对称性计算以加速，仅依靠 Canonical Rep 的唯一性去重)
    // 如果你需要统计对称性数量，可以在去重后的后处理中做。

    return { total_len, 0 };
}

// 4. 获取标准型 (内存极致优化版)
// 使用 thread_local 避免反复分配 vector
void get_canonical_rep_optimized(const std::vector<int>& A, int n1, std::vector<int>& out_best) {
    // 线程局部缓存，每个线程只有一份，重复利用
    static thread_local std::vector<int> B;
    static thread_local std::vector<int> curr;
    static thread_local std::vector<int> temp_best;

    // 初始化大小
    size_t sz = A.size();
    if (B.size() != sz) B.resize(sz);
    if (curr.size() != sz) curr.resize(sz);
    if (temp_best.size() != sz) temp_best.resize(sz);

    // 构造 B
    for (size_t i = 0; i < sz; ++i) {
        if (A[i] < n1) B[i] = A[i] + n1;
        else B[i] = A[i] - n1;
    }

    bool initialized = false;
    long long min_len = -1;
    // max_sym 不再作为核心判据，因为计算太慢，且字典序已足够保证唯一性

    const std::vector<int>* sources[2] = { &A, &B };

    for (int src_idx = 0; src_idx < 2; ++src_idx) {
        const auto& src = *sources[src_idx];
        for (int t = 0; t < n1; ++t) {
            for (int mode = 0; mode < 2; ++mode) {

                // 1. 变换
                apply_transform(src, curr, t, n1, mode);

                // 2. 规范化
                normalize_matrix(curr);

                // 3. 评分 (只算模长)
                long long curr_len = 0;
                for (int i = 0; i < n1; ++i) {
                    int diff = std::abs((curr[2 * i] % n1) - (curr[2 * i + 1] % n1));
                    int dist = diff < (n1 - diff) ? diff : (n1 - diff);
                    curr_len += dist * dist;
                    // 剪枝：如果已经超过当前最小长度，且我们只关心最小长度优先
                    if (initialized && curr_len > min_len) break;
                }

                bool replace = false;
                if (!initialized) {
                    replace = true;
                }
                else if (curr_len < min_len) {
                    replace = true;
                }
                else if (curr_len == min_len) {
                    // 长度相同时，直接比字典序
                    // std::vector 的比较操作符就是字典序
                    if (curr < temp_best) {
                        replace = true;
                    }
                }

                if (replace) {
                    temp_best = curr; // 拷贝一次
                    min_len = curr_len;
                    initialized = true;
                }
            }
        }
    }
    out_best = temp_best;
}

// ================== 搜索逻辑 ==================

void solve(std::vector<int>& current_A,
    std::vector<char>& used,
    RollbackDSU& dsu,
    int n1, int target_rows, int fangxiang,
    const std::vector<std::vector<int>>& valid_avcs,
    const std::vector<std::pair<int, int>>& candidates,
    const std::vector<int>& cand_start_indices,
    std::unordered_set<std::vector<int>, VectorHash>& local_uniques) { // 使用 unordered_set

    int current_rows = (int)current_A.size() / 2;

    if (current_rows == target_rows) {
        // 验证 AVC
        std::vector<int> v_cycles;
        // 简单的 visited 数组，避免分配
        // n1 最大 30 左右，用 uint64 位运算或固定数组
        long long visited_mask = 0; // 仅当 2*n1 <= 64
        // 如果 2*n1 > 64，还是用 vector<bool> 或 char 数组
        // 这里为了通用性保留 vector，但移到外部复用会更好，这里暂时保留局部变量因为它是递归末端
        std::vector<bool> visited(2 * n1, false);

        for (int i = 0; i < 2 * n1; ++i) {
            int root = dsu.find(i);
            if (!visited[root]) {
                v_cycles.push_back(dsu.sz[root]);
                visited[root] = true;
            }
        }
        std::sort(v_cycles.begin(), v_cycles.end());

        bool avc_match = false;
        // 快速检查：先比环的个数
        for (const auto& avc : valid_avcs) {
            if (avc.size() != v_cycles.size()) continue;
            // 再比内容
            if (avc == v_cycles) { // vector 重载了 ==
                avc_match = true;
                break;
            }
        }

        if (avc_match) {
            std::vector<int> canonical;
            get_canonical_rep_optimized(current_A, n1, canonical);
            local_uniques.insert(canonical);
        }
        return;
    }

    // --- 搜索剪枝 ---
    int last_u = current_A[current_A.size() - 2];
    int remaining_slots = (target_rows - current_rows) * 2;

    int start_idx = 0;
    if (last_u + 1 < (int)cand_start_indices.size()) {
        start_idx = cand_start_indices[last_u + 1];
    }
    if (start_idx == -1) return;

    for (size_t i = start_idx; i < candidates.size(); ++i) {
        int u = candidates[i].first;
        int v = candidates[i].second;

        // 强剪枝：剩余数字不足
        if ((2 * n1 - u) < remaining_slots) break;
        if (used[u] || used[v]) continue;

        // 拓扑剪枝
        int u_next = (u < n1) ? (u + 1) % n1 : n1 + (u - n1 + 1) % n1;
        int v_next = (v < n1) ? (v + 1) % n1 : n1 + (v - n1 + 1) % n1;
        bool is_same_block = (u < n1 && v < n1) || (u >= n1 && v >= n1);

        int n1a, n1b, n2a, n2b;
        if (fangxiang == 1) {
            n1a = u_next; n1b = v; n2a = u; n2b = v_next;
        }
        else {
            if (is_same_block) { n1a = u_next; n1b = v; n2a = u; n2b = v_next; }
            else { n1a = u; n1b = v; n2a = u_next; n2b = v_next; }
        }

        auto res1 = dsu.unite(n1a, n1b);
        // AVC 约束：最小环长度为 3
        if (res1.second && dsu.get_size(n1a) < 3) { dsu.rollback(res1.first); continue; }

        auto res2 = dsu.unite(n2a, n2b);
        if (res2.second && dsu.get_size(n2a) < 3) {
            dsu.rollback(res2.first); dsu.rollback(res1.first); continue;
        }

        used[u] = 1; used[v] = 1;
        current_A.push_back(u); current_A.push_back(v);

        solve(current_A, used, dsu, n1, target_rows, fangxiang, valid_avcs,
            candidates, cand_start_indices, local_uniques);

        current_A.pop_back(); current_A.pop_back();
        used[u] = 0; used[v] = 0;
        dsu.rollback(res2.first); dsu.rollback(res1.first);
    }
}

// ================== MEX 入口 ==================

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 5) mexErrMsgTxt("Expected 5 inputs.");

    int* seeds_ptr = (int*)mxGetData(prhs[0]);
    size_t num_seeds = mxGetM(prhs[0]);
    int* cand_ptr = (int*)mxGetData(prhs[1]);
    size_t cand_count = mxGetM(prhs[1]);
    int n1 = (int)mxGetScalar(prhs[2]);
    int fangxiang = (int)mxGetScalar(prhs[3]);
    double* avc_data = mxGetPr(prhs[4]);
    size_t avc_rows = mxGetM(prhs[4]);
    size_t avc_cols = mxGetN(prhs[4]);

    std::vector<std::pair<int, int>> candidates;
    candidates.reserve(cand_count);
    for (size_t i = 0; i < cand_count; ++i) {
        candidates.push_back({ cand_ptr[i], cand_ptr[i + cand_count] });
    }
    // 确保 MATLAB 传来的 candidates 已经是排好序的，或者这里再排一次
    std::sort(candidates.begin(), candidates.end());

    std::vector<int> cand_start_indices(2 * n1 + 1, -1);
    for (int i = (int)candidates.size() - 1; i >= 0; --i) {
        cand_start_indices[candidates[i].first] = i;
    }
    int last_valid = -1;
    for (int i = 2 * n1; i >= 0; --i) {
        if (cand_start_indices[i] != -1) last_valid = cand_start_indices[i];
        else cand_start_indices[i] = last_valid;
    }

    std::vector<std::vector<int>> valid_avcs;
    for (size_t i = 0; i < avc_rows; ++i) {
        std::vector<int> row;
        for (size_t j = 0; j < avc_cols; ++j) {
            int val = (int)avc_data[i + j * avc_rows];
            if (val > 0) row.push_back(val);
        }
        std::sort(row.begin(), row.end());
        valid_avcs.push_back(row);
    }

    int max_threads = omp_get_max_threads();
    // 使用 unordered_set 数组
    std::vector<std::unordered_set<std::vector<int>, VectorHash>> thread_uniques(max_threads);

#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < (int)num_seeds; ++k) {
        int tid = omp_get_thread_num();

        int u1 = seeds_ptr[k];
        int v1 = seeds_ptr[k + num_seeds];
        int u2 = seeds_ptr[k + 2 * num_seeds];
        int v2 = seeds_ptr[k + 3 * num_seeds];

        std::vector<int> current_A;
        current_A.reserve(2 * n1);
        std::vector<char> used(2 * n1, 0);

        current_A.push_back(u1); current_A.push_back(v1);
        used[u1] = 1; used[v1] = 1;
        current_A.push_back(u2); current_A.push_back(v2);
        used[u2] = 1; used[v2] = 1;

        RollbackDSU dsu(2 * n1);

        // 种子初始化逻辑
        auto apply_edge_simple = [&](int u, int v) {
            int u_next = (u < n1) ? (u + 1) % n1 : n1 + (u - n1 + 1) % n1;
            int v_next = (v < n1) ? (v + 1) % n1 : n1 + (v - n1 + 1) % n1;
            bool is_same_block = (u < n1 && v < n1) || (u >= n1 && v >= n1);
            if (fangxiang == 1) {
                dsu.unite(u_next, v); dsu.unite(u, v_next);
            }
            else {
                if (is_same_block) { dsu.unite(u_next, v); dsu.unite(u, v_next); }
                else { dsu.unite(u, v); dsu.unite(u_next, v_next); }
            }
            };
        apply_edge_simple(u1, v1);
        apply_edge_simple(u2, v2);

        solve(current_A, used, dsu, n1, n1, fangxiang, valid_avcs,
            candidates, cand_start_indices, thread_uniques[tid]);
    }

    // 合并
    // 为了节省内存，我们可以先合并到一个 vector 再去重，或者直接合并到 unordered_set
    // 鉴于 n1=12 时解可能不多（去重后），直接 insert 即可
    std::unordered_set<std::vector<int>, VectorHash> global_uniques;
    for (auto& local_set : thread_uniques) {
        // 使用 move iterator 避免拷贝
        global_uniques.insert(std::make_move_iterator(local_set.begin()),
            std::make_move_iterator(local_set.end()));
    }

    if (global_uniques.empty()) {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    else {
        size_t total_sols = global_uniques.size();
        plhs[0] = mxCreateNumericMatrix(total_sols * n1, 2, mxINT32_CLASS, mxREAL);
        int* out_ptr = (int*)mxGetData(plhs[0]);

        size_t row_offset = 0;
        for (const auto& sol : global_uniques) {
            for (int r = 0; r < n1; ++r) {
                out_ptr[row_offset + r] = sol[2 * r];
                out_ptr[row_offset + r + total_sols * n1] = sol[2 * r + 1];
            }
            row_offset += n1;
        }
    }
}


/*#include "mex.h"//matlab端并行去重
#include <vector>
#include <algorithm>
#include <numeric>
#include <omp.h>

// ================== 数据结构 ==================

// 支持回滚的并查集 (DSU with Rollback)
// 用于在递归过程中低成本地维护拓扑状态
struct RollbackDSU {
    std::vector<int> parent;
    std::vector<int> sz; // 集合大小（即角点循环的长度）

    // 记录操作历史用于回滚: {child_node, parent_node_before_merge}
    // 其实只需要记录谁变成了谁的子节点，以及原来的大小
    struct History {
        int child;
        int parent;
        bool merged;
    };

    RollbackDSU(int n) {
        parent.resize(n);
        std::iota(parent.begin(), parent.end(), 0);
        sz.assign(n, 1);
    }

    // 查找根节点 (带路径压缩的查找不适合回滚，这里使用非压缩或仅在当前层级安全的查找)
    // 为了支持回滚，我们不能做路径压缩(Path Compression)，否则回滚太复杂。
    // 由于 n 很小 (24左右)，按秩合并(Union by Rank/Size) 足够快，树高也就是 log(24) ≈ 4-5。
    int find(int x) const {
        while (parent[x] != x) {
            x = parent[x];
        }
        return x;
    }

    // 尝试合并，返回合并后的状态信息用于回滚，以及是否闭环
    // returns: {history_info, is_cycle_closed}
    std::pair<History, bool> unite(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX != rootY) {
            // 按大小合并：小的挂在大的下面
            if (sz[rootX] < sz[rootY]) std::swap(rootX, rootY);

            // rootY 挂到 rootX
            parent[rootY] = rootX;
            sz[rootX] += sz[rootY];

            return { {rootY, rootX, true}, false };
        }
        else {
            // 已经在同一个集合，说明形成了闭环
            return { {0, 0, false}, true };
        }
    }

    // 回滚操作
    void rollback(const History& h) {
        if (h.merged) {
            sz[h.parent] -= sz[h.child];
            parent[h.child] = h.child; // 恢复 child 为根
        }
    }

    // 获取集合大小
    int get_size(int x) const {
        return sz[find(x)];
    }
};

// ================== 核心逻辑 ==================

/**
 * 现有的剪枝逻辑清单 (Pruning Check List):
 *
 * 1. 顺序剪枝 (Ordering Pruning):
 *    强制要求矩阵的第一列 u 严格递增 (u_current > u_prev)。
 *    通过 candidate_start_indices 实现 O(1) 跳转，直接跳过无效的小 u。
 *
 * 2. 空间前瞻剪枝 (Lookahead Pruning):
 *    如果 (2*n1 - 1 - 当前u) < 剩余需要的行数 * 2，说明即便后面选了所有最大的数，
 *    也不够填满矩阵，直接 break。
 *
 * 3. 占用检查 (Availability Pruning):
 *    检查 used[u] 和 used[v] 是否被占用。
 *
 * 4. 拓扑环剪枝 (Topological Cycle Pruning) [本次新增]:
 *    每加一条边，计算其导致的角点粘合。
 *    如果发现某个角点循环闭合 (closed cycle)，且其长度 (sz) < 3，
 *    则判定为非法 (违背 AVC min_val=3)，停止该分支搜索。
 

void solve(std::vector<int>& current_A,
    std::vector<char>& used,
    RollbackDSU& dsu,
    int n1, int target_rows, int fangxiang,
    const std::vector<std::vector<int>>& valid_avcs,
    const std::vector<std::pair<int, int>>& candidates,
    const std::vector<int>& cand_start_indices, // 预计算的索引
    std::vector<int>& local_solutions) {

    int current_rows = (int)current_A.size() / 2;

    // === 成功终止条件 ===
    if (current_rows == target_rows) {
        // 最后再做一次完整的 AVC 校验（包含所有环的排序匹配）
        // 虽然中间做了剪枝，但最终还需要匹配具体的 AVC 组合
        std::vector<int> v_cycles;
        std::vector<bool> visited(2 * n1, false);
        for (int i = 0; i < 2 * n1; ++i) {
            int root = dsu.find(i);
            if (!visited[root]) {
                v_cycles.push_back(dsu.sz[root]);
                visited[root] = true;
            }
        }
        std::sort(v_cycles.begin(), v_cycles.end());

        for (const auto& avc : valid_avcs) {
            if (avc.size() != v_cycles.size()) continue;
            bool match = true;
            for (size_t k = 0; k < avc.size(); ++k) {
                if (avc[k] != v_cycles[k]) {
                    match = false;
                    break;
                }
            }
            if (match) {
                local_solutions.insert(local_solutions.end(), current_A.begin(), current_A.end());
                return; // 找到一种即可，或者继续找？通常这里 continue 找所有
            }
        }
        return;
    }

    // === 准备搜索 ===
    int last_u = current_A[current_A.size() - 2];
    int remaining_pairs = target_rows - current_rows;
    int remaining_slots = remaining_pairs * 2;
    int max_possible_u = 2 * n1 - 1;

    // 优化1: 利用索引直接跳跃到 u > last_u 的位置
    int start_idx = 0;
    if (last_u + 1 < (int)cand_start_indices.size()) {
        start_idx = cand_start_indices[last_u + 1];
    }
    if (start_idx == -1) return; // 后面没有合法的 u 了

    // 遍历候选边
    for (size_t i = start_idx; i < candidates.size(); ++i) {
        int u = candidates[i].first;
        int v = candidates[i].second;

        // 剪枝2: 空间前瞻 (如果剩下的数不够填空)
        // max_possible_u - u 是如果我们选了 u，后面最多还能有多少个更大的数
        // 我们还需要填 (remaining_pairs * 2 - 2) 个位置
        // 宽松一点的判断：如果 (2*n1 - u) < remaining_slots，肯定不够
        if ((2 * n1 - u) < remaining_slots) break;

        if (used[u] || used[v]) continue;

        // === 剪枝4: 拓扑环检查 (Incremental DSU Check) ===
        // 计算拓扑邻接关系
        int u_next = (u < n1) ? (u + 1) % n1 : n1 + (u - n1 + 1) % n1;
        int v_next = (v < n1) ? (v + 1) % n1 : n1 + (v - n1 + 1) % n1;

        bool is_same_block = (u < n1 && v < n1) || (u >= n1 && v >= n1);

        int node1_a, node1_b, node2_a, node2_b;

        if (fangxiang == 1) {
            // u_next <-> v, u <-> v_next
            node1_a = u_next; node1_b = v;
            node2_a = u;      node2_b = v_next;
        }
        else {
            if (is_same_block) {
                node1_a = u_next; node1_b = v;
                node2_a = u;      node2_b = v_next;
            }
            else {
                node1_a = u;      node1_b = v;
                node2_a = u_next; node2_b = v_next;
            }
        }

        // 尝试合并第一对
        auto res1 = dsu.unite(node1_a, node1_b);
        if (res1.second) { // 闭环了
            if (dsu.get_size(node1_a) < 3) { // 环太短，剪枝
                dsu.rollback(res1.first);
                continue;
            }
        }

        // 尝试合并第二对
        auto res2 = dsu.unite(node2_a, node2_b);
        if (res2.second) { // 闭环了
            if (dsu.get_size(node2_a) < 3) { // 环太短，剪枝
                dsu.rollback(res2.first);
                dsu.rollback(res1.first);
                continue;
            }
        }

        // === 通过所有检查，进入下一层 ===
        used[u] = 1;
        used[v] = 1;
        current_A.push_back(u);
        current_A.push_back(v);

        solve(current_A, used, dsu, n1, target_rows, fangxiang, valid_avcs,
            candidates, cand_start_indices, local_solutions);

        // === 回溯 (Backtrack) ===
        current_A.pop_back();
        current_A.pop_back();
        used[u] = 0;
        used[v] = 0;

        // DSU 回滚
        dsu.rollback(res2.first);
        dsu.rollback(res1.first);
    }
}

// ================== MEX 入口 ==================

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 5) mexErrMsgTxt("Expected 5 inputs.");

    // 1. 获取输入
    int* seeds_ptr = (int*)mxGetData(prhs[0]);
    size_t num_seeds = mxGetM(prhs[0]);

    int* cand_ptr = (int*)mxGetData(prhs[1]);
    size_t cand_count = mxGetM(prhs[1]);

    int n1 = (int)mxGetScalar(prhs[2]);
    int fangxiang = (int)mxGetScalar(prhs[3]);

    double* avc_data = mxGetPr(prhs[4]);
    size_t avc_rows = mxGetM(prhs[4]);
    size_t avc_cols = mxGetN(prhs[4]);

    // 2. 处理候选边 & 预计算索引
    std::vector<std::pair<int, int>> candidates;
    candidates.reserve(cand_count);
    for (size_t i = 0; i < cand_count; ++i) {
        // MATLAB int32 也是列主序
        candidates.push_back({ cand_ptr[i], cand_ptr[i + cand_count] });
    }

    // 确保排序 (u 升序, v 升序)
    std::sort(candidates.begin(), candidates.end());

    // 构建快速索引表: cand_start_indices[k] 存储第一个 u >= k 的 index
    // 这样当 last_u = k 时，我们查找 indices[k+1] 即可
    std::vector<int> cand_start_indices(2 * n1 + 1, -1);
    for (int i = (int)candidates.size() - 1; i >= 0; --i) {
        cand_start_indices[candidates[i].first] = i;
    }
    // 填充空洞：如果 u=5 不存在，cand_start_indices[5] 应该等于 cand_start_indices[6] 
    // 或者保持 -1 表示后面没了。为了方便，我们让它指向下一个存在的
    int last_valid = -1;
    for (int i = 2 * n1; i >= 0; --i) {
        if (cand_start_indices[i] != -1) last_valid = cand_start_indices[i];
        else cand_start_indices[i] = last_valid;
    }

    // 3. 处理 AVC
    std::vector<std::vector<int>> valid_avcs;
    for (size_t i = 0; i < avc_rows; ++i) {
        std::vector<int> row;
        for (size_t j = 0; j < avc_cols; ++j) {
            int val = (int)avc_data[i + j * avc_rows];
            if (val > 0) row.push_back(val);
        }
        std::sort(row.begin(), row.end());
        valid_avcs.push_back(row);
    }

    // 4. 并行准备
    int max_threads = omp_get_max_threads();
    std::vector<std::vector<int>> thread_results(max_threads);

    // 5. 并行执行
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < (int)num_seeds; ++k) {
        int tid = omp_get_thread_num();

        // 读取种子
        int u1 = seeds_ptr[k];
        int v1 = seeds_ptr[k + num_seeds];
        int u2 = seeds_ptr[k + 2 * num_seeds];
        int v2 = seeds_ptr[k + 3 * num_seeds];

        // 初始化状态
        std::vector<int> current_A;
        current_A.reserve(2 * n1);
        std::vector<char> used(2 * n1, 0);

        current_A.push_back(u1); current_A.push_back(v1);
        used[u1] = 1; used[v1] = 1;

        current_A.push_back(u2); current_A.push_back(v2);
        used[u2] = 1; used[v2] = 1;

        // 初始化 DSU 并应用种子的拓扑影响
        RollbackDSU dsu(2 * n1);

        // 辅助 lambda 应用边
        auto apply_edge = [&](int u, int v) -> bool {
            int u_next = (u < n1) ? (u + 1) % n1 : n1 + (u - n1 + 1) % n1;
            int v_next = (v < n1) ? (v + 1) % n1 : n1 + (v - n1 + 1) % n1;
            bool is_same_block = (u < n1 && v < n1) || (u >= n1 && v >= n1);

            if (fangxiang == 1) {
                dsu.unite(u_next, v);
                dsu.unite(u, v_next);
            }
            else {
                if (is_same_block) {
                    dsu.unite(u_next, v);
                    dsu.unite(u, v_next);
                }
                else {
                    dsu.unite(u, v);
                    dsu.unite(u_next, v_next);
                }
            }
            // 种子阶段暂不检查环大小，假设种子是合法的或在递归第一步会被检查
            return true;
            };

        apply_edge(u1, v1);
        apply_edge(u2, v2);

        // 递归
        solve(current_A, used, dsu, n1, n1, fangxiang, valid_avcs,
            candidates, cand_start_indices, thread_results[tid]);
    }

    // 6. 输出合并
    size_t total_ints = 0;
    for (const auto& res : thread_results) total_ints += res.size();

    if (total_ints == 0) {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    else {
        size_t total_sols = total_ints / (2 * n1);
        plhs[0] = mxCreateNumericMatrix(total_sols * n1, 2, mxINT32_CLASS, mxREAL);
        int* out_ptr = (int*)mxGetData(plhs[0]);

        size_t current_sol_idx = 0;
        for (const auto& res : thread_results) {
            size_t num_local = res.size() / (2 * n1);
            for (size_t s = 0; s < num_local; ++s) {
                for (size_t r = 0; r < (size_t)n1; ++r) {
                    size_t src_idx = s * (2 * n1) + 2 * r;
                    int u = res[src_idx];
                    int v = res[src_idx + 1];

                    size_t row_idx = (current_sol_idx + s) * n1 + r;
                    out_ptr[row_idx] = u;
                    out_ptr[row_idx + (total_sols * n1)] = v;
                }
            }
            current_sol_idx += num_local;
        }
    }
}*/



