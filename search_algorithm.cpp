#include "mex.h"
#include <vector>
#include <algorithm>
#include <cstring>
#include <omp.h>

// ================= 参数结构体 =================
struct Params {
    int n1;
    int K;
    int daxiao;
    int vmax;
    int n2;          // 新增 n2
    int num_nodes;
};

// ================= 常量定义 =================
const int MAX_DAXIAO = 32;
const int MAX_NODES = 64;

// ================= 基础结构 =================

struct Row {
    int u, v, t;
};

// 矩阵结构体 (用于存储解和同构去重)
struct Matrix {
    int data[MAX_DAXIAO][3];
    int rows;

    Matrix() : rows(0) {}

    void init(int r) { rows = r; }

    inline void set(int r, int c, int val) {
        data[r][c] = val;
    }
    inline int get(int r, int c) const {
        return data[r][c];
    }

    // 比较函数：用于排序和去重
    bool operator<(const Matrix& other) const {
        if (rows != other.rows) return rows < other.rows;
        return std::memcmp(data, other.data, rows * 3 * sizeof(int)) < 0;
    }
    bool operator==(const Matrix& other) const {
        if (rows != other.rows) return false;
        return std::memcmp(data, other.data, rows * 3 * sizeof(int)) == 0;
    }
};

// ================= 核心状态：顶点库/链系统 =================
// 不使用图论搜索，而是维护顶点串联状态
struct PathSystem {
    // degree[i]: 顶点 i 的当前度数 (0, 1, 2)
    // 0: 未使用
    // 1: 链的端点
    // 2: 链的内部点 (或环上的点)
    int degree[MAX_NODES];

    // mate[i]: 如果 degree[i] == 1, mate[i] 存储该链另一端的顶点编号
    // 如果 degree[i] != 1, mate[i] 的值无意义(或不维护)
    int mate[MAX_NODES];

    // size[i]: 如果 degree[i] == 1, size[i] 存储以 i 为端点的整条链的顶点个数
    // 新增一条边 [u,v] 时，顶点个数为 2
    int size[MAX_NODES];

    // 当前形成的闭环个数
    int closed_loops_count;

    // 当前使用的矩阵
    Matrix current_A;

    PathSystem(int daxiao) {
        current_A.init(daxiao);
        std::memset(degree, 0, sizeof(degree));
        std::memset(mate, 0, sizeof(mate));
        std::memset(size, 0, sizeof(size));
        closed_loops_count = 0;
    }

    // 尝试添加连接 (u, v)
    // 返回 true 表示成功添加且符合所有约束
    // 返回 false 表示违反约束（如度数超标、长度超标等）
    bool add_segment(int u, int v, int vmax) {
        // 1. 基础度数检查：每个点最多连接2条边
        if (degree[u] >= 2 || degree[v] >= 2) return false;
        if (u == v) return false; // 不允许自环（长度1，小于最小环长3）

        // Case A: 两个点都是新点 (0 - 0)
        // 形成新链: [u, v], 顶点数 2
        if (degree[u] == 0 && degree[v] == 0) {
            degree[u] = 1; mate[u] = v; size[u] = 2;
            degree[v] = 1; mate[v] = u; size[v] = 2;
            return true;
        }

        // Case B: 一个新点，一个端点 (1 - 0)
        // 延长链条
        if (degree[u] == 1 && degree[v] == 0) {
            int z = mate[u]; // u 所在的链的另一端
            int new_len = size[u] + 1; // 增加一个顶点 v

            if (new_len > vmax) return false; // 串联起来没有闭环，长度也要 <= vmax

            // 更新状态
            degree[u] = 2; // u 变为内部
            degree[v] = 1; // v 成为新端点

            mate[z] = v;
            mate[v] = z;
            size[z] = new_len;
            size[v] = new_len;
            return true;
        }
        if (degree[u] == 0 && degree[v] == 1) {
            return add_segment(v, u, vmax); // 对称处理
        }

        // Case C: 两个都是端点 (1 - 1)
        // 合并两条链 或 闭合环
        if (degree[u] == 1 && degree[v] == 1) {
            int z = mate[u];
            int w = mate[v];

            if (z == v) {
                // === 形成闭环 ===
                // u 和 v 是同一条链的两端，连起来就闭合了
                int loop_len = size[u]; // 闭环的顶点数等于链的顶点数 (因为加入边(u,v)只是闭合，不加新点)

                // 闭环约束：个数 >= 3 且 <= vmax
                if (loop_len < 3 || loop_len > vmax) return false;

                degree[u] = 2;
                degree[v] = 2;
                closed_loops_count++;
                return true;
            }
            else {
                // === 合并两条不同的链 ===
                // 链1: z...u, 链2: v...w
                // 合并后: z...u-v...w
                int new_len = size[u] + size[v]; // 顶点数相加

                if (new_len > vmax) return false; // 串联长度约束

                degree[u] = 2;
                degree[v] = 2;

                mate[z] = w;
                mate[w] = z;
                size[z] = new_len;
                size[w] = new_len;
                return true;
            }
        }

        return false; // Should not reach here
    }
};

// ================= 辅助函数 =================

// 严格取模：保证结果非负
inline int strict_mod(int a, int n) {
    int res = a % n;
    return res < 0 ? res + n : res;
}

// 顶点变换逻辑：
// 0~n1-1 (内圈): x+1 -> 0 if x=n1-1
// n1~2n1-1 (外圈): x+1 -> n1 if x=2n1-1
inline int get_next_vertex(int x, int n1) {
    if (x < n1) {
        return (x + 1 == n1) ? 0 : x + 1;
    }
    else {
        return (x + 1 == 2 * n1) ? n1 : x + 1;
    }
}

// 矩阵排序与标准化 (用于去重)
void sort_and_normalize_matrix(Matrix& M, int daxiao) {
    // 行内排序
    for (int i = 0; i < daxiao; ++i) {
        if (M.data[i][0] > M.data[i][1]) {
            std::swap(M.data[i][0], M.data[i][1]);
        }
    }
    // 矩阵行排序
    for (int i = 0; i < daxiao - 1; ++i) {
        for (int j = 0; j < daxiao - i - 1; ++j) {
            bool swap = false;
            int* rA = M.data[j];
            int* rB = M.data[j + 1];
            if (rA[0] > rB[0]) swap = true;
            else if (rA[0] == rB[0]) {
                if (rA[1] > rB[1]) swap = true;
                else if (rA[1] == rB[1] && rA[2] > rB[2]) swap = true;
            }
            if (swap) {
                int tmp[3]; std::memcpy(tmp, rA, 12);
                std::memcpy(rA, rB, 12); std::memcpy(rB, tmp, 12);
            }
        }
    }
}

// 计算同构代表元
Matrix get_canonical(const Matrix& A, const Params& p) {
    Matrix sortedA = A;
    sort_and_normalize_matrix(sortedA, p.daxiao);

    Matrix best;
    bool init = false;

    // 构建 B (互补矩阵)
    Matrix B; B.init(p.daxiao);
    for (int i = 0; i < p.daxiao; ++i) {
        int m = sortedA.get(i, 0);
        int n = sortedA.get(i, 1);
        int t = sortedA.get(i, 2);
        int m_new = (m < p.n1) ? m + p.n1 : m - p.n1;
        int n_new = (n < p.n1) ? n + p.n1 : n - p.n1;
        B.set(i, 0, m_new); B.set(i, 1, n_new); B.set(i, 2, t);
    }
    sort_and_normalize_matrix(B, p.daxiao);

    const Matrix* sources[2] = { &sortedA, &B };
    Matrix T; T.init(p.daxiao);

    for (int s = 0; s < 2; ++s) {
        const Matrix& M = *sources[s];
        for (int d = 0; d < p.n1; ++d) {
            // 平移
            for (int i = 0; i < p.daxiao; ++i) {
                int c[2] = { M.get(i,0), M.get(i,1) };
                int nc[2];
                for (int k = 0; k < 2; ++k) {
                    if (c[k] < p.n1) nc[k] = strict_mod(c[k] + d, p.n1);
                    else nc[k] = p.n1 + strict_mod(c[k] - p.n1 + d, p.n1);
                }
                T.set(i, 0, nc[0]); T.set(i, 1, nc[1]); T.set(i, 2, M.get(i, 2));
            }
            sort_and_normalize_matrix(T, p.daxiao);
            if (!init || T < best) { best = T; init = true; }

            // 反射
            for (int i = 0; i < p.daxiao; ++i) {
                int c[2] = { M.get(i,0), M.get(i,1) };
                int nc[2];
                for (int k = 0; k < 2; ++k) {
                    if (c[k] < p.n1) nc[k] = strict_mod(d - c[k], p.n1);
                    else nc[k] = p.n1 + strict_mod(d - (c[k] - p.n1), p.n1);
                }
                T.set(i, 0, nc[0]); T.set(i, 1, nc[1]); T.set(i, 2, M.get(i, 2));
            }
            sort_and_normalize_matrix(T, p.daxiao);
            if (T < best) best = T;
        }
    }
    return best;
}

// ================= 搜索逻辑 =================

// 使用按值传递 State (PathSystem) 来简化回溯
// 因为 State 较小 (~1KB)，且深度有限，直接拷贝比复杂的回溯逻辑更安全且高效
void dfs(PathSystem state, int row_idx, const std::vector<Row>& combinations,
    std::vector<Matrix>& local_results, const Params& p) {

    // 终局判断
    if (row_idx == p.daxiao) {
        // 增加条件：最后达到 daxiao 时，闭环个数必须等于 n2
        if (state.closed_loops_count == p.n2) {
            // 这里还可以加上 t 的求和校验等原始逻辑
            int sum_t = 0;
            bool has_intra = false;
            for (int i = 0; i < p.daxiao; ++i) {
                sum_t += state.current_A.get(i, 2);
                int u = state.current_A.get(i, 0);
                int v = state.current_A.get(i, 1);
                if ((u < p.n1) == (v < p.n1)) has_intra = true;
            }

            bool valid = true;
            if (sum_t == 0) valid = false;
            else if (sum_t == p.daxiao && !has_intra) valid = false;

            if (valid) {
                local_results.push_back(get_canonical(state.current_A, p));
            }
        }
        return;
    }

    // 剪枝：当前闭环数已经超过 n2，直接返回
    if (state.closed_loops_count > p.n2) return;

    // 剪枝：剩余行数不足以凑齐 n2 个闭环 (估算：每个闭环至少3条边，粗略估计)
    // 这里的剪枝比较难精确，保持基础的闭环数不超过即可

    int last_u = state.current_A.get(row_idx - 1, 0);
    int remaining_rows = p.daxiao - row_idx;

    for (const auto& cand : combinations) {
        // 顺序约束：第一列必须非降序
        if (cand.u <= last_u) continue;

        // 简单的剩余行数剪枝
        if ((p.num_nodes - 1 - cand.u) < 2 * remaining_rows - 1) break;

        // 顶点已经被度数占满的快速检查 (degree < 2)
        // 注意：cand.u 和 cand.v 是矩阵中的数，不是所有涉及的顶点
        // 需要根据 t 具体判断涉及哪些顶点

        // 拷贝当前状态准备尝试
        PathSystem next_state = state;

        // 记录矩阵行
        next_state.current_A.set(row_idx, 0, cand.u);
        next_state.current_A.set(row_idx, 1, cand.v);
        next_state.current_A.set(row_idx, 2, cand.t);

        // 解析需要添加的两条线段
        int u1, v1, u2, v2;
        if (cand.t == 1) {
            // 串联 [m, n] 和 [m+1, n+1]
            u1 = cand.u;
            v1 = cand.v;
            u2 = get_next_vertex(cand.u, p.n1);
            v2 = get_next_vertex(cand.v, p.n1);
        }
        else { // t == 0
            // 串联 [m+1, n] 和 [m, n+1]
            u1 = get_next_vertex(cand.u, p.n1);
            v1 = cand.v;
            u2 = cand.u;
            v2 = get_next_vertex(cand.v, p.n1);
        }

        // 尝试添加第一条线段
        if (next_state.add_segment(u1, v1, p.vmax)) {
            // 尝试添加第二条线段
            if (next_state.add_segment(u2, v2, p.vmax)) {
                // 两条都成功添加，进入下一层
                dfs(next_state, row_idx + 1, combinations, local_results, p);
            }
        }
        // 如果失败，next_state 被丢弃，循环继续
    }
}

// ================= MEX 入口 =================

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // 增加 n2，参数变为 7 个
    if (nrhs != 7) mexErrMsgTxt("Expected 7 inputs: n1, K, daxiao, vmax, firstE, combinations1, n2");

    Params p;
    p.n1 = (int)mxGetScalar(prhs[0]);
    p.K = (int)mxGetScalar(prhs[1]);
    p.daxiao = (int)mxGetScalar(prhs[2]);
    p.vmax = (int)mxGetScalar(prhs[3]);
    p.num_nodes = 2 * p.n1;
    // 读取 n2
    p.n2 = (int)mxGetScalar(prhs[6]);

    double* fE_ptr = mxGetPr(prhs[4]);
    size_t fE_rows = mxGetM(prhs[4]);
    std::vector<Row> firstE(fE_rows);
    for (size_t i = 0; i < fE_rows; ++i) {
        firstE[i].u = (int)fE_ptr[i + 0 * fE_rows];
        firstE[i].v = (int)fE_ptr[i + 1 * fE_rows];
        firstE[i].t = (int)fE_ptr[i + 2 * fE_rows];
    }

    double* c1_ptr = mxGetPr(prhs[5]);
    size_t c1_rows = mxGetM(prhs[5]);
    std::vector<Row> combinations1(c1_rows);
    for (size_t i = 0; i < c1_rows; ++i) {
        combinations1[i].u = (int)c1_ptr[i + 0 * c1_rows];
        combinations1[i].v = (int)c1_ptr[i + 1 * c1_rows];
        combinations1[i].t = (int)c1_ptr[i + 2 * c1_rows];
    }

    std::vector<Matrix> global_results;

#pragma omp parallel
    {
        std::vector<Matrix> local_results;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < (int)fE_rows; ++i) {
            PathSystem s(p.daxiao);

            // 初始化第一行
            Row r0 = firstE[i];
            s.current_A.set(0, 0, r0.u);
            s.current_A.set(0, 1, r0.v);
            s.current_A.set(0, 2, r0.t);

            // 处理第一行的两条线段
            int u1, v1, u2, v2;
            if (r0.t == 1) {
                u1 = r0.u; v1 = r0.v;
                u2 = get_next_vertex(r0.u, p.n1); v2 = get_next_vertex(r0.v, p.n1);
            }
            else {
                u1 = get_next_vertex(r0.u, p.n1); v1 = r0.v;
                u2 = r0.u; v2 = get_next_vertex(r0.v, p.n1);
            }

            // 对于第一行，理论上总是合法的（初始），但还是检查一下
            if (s.add_segment(u1, v1, p.vmax) && s.add_segment(u2, v2, p.vmax)) {
                dfs(s, 1, combinations1, local_results, p);
            }
        }

#pragma omp critical
        {
            global_results.insert(global_results.end(), local_results.begin(), local_results.end());
        }
    }

    // 去重
    if (!global_results.empty()) {
        std::sort(global_results.begin(), global_results.end());
        auto last = std::unique(global_results.begin(), global_results.end());
        global_results.erase(last, global_results.end());
    }

    // 输出
    int count = global_results.size();
    const mwSize dims[3] = { (mwSize)p.daxiao, 3, (mwSize)count };
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* out_ptr = mxGetPr(plhs[0]);

    int k = 0;
    for (const auto& mat : global_results) {
        for (int c = 0; c < 3; ++c) {
            for (int r = 0; r < p.daxiao; ++r) {
                out_ptr[r + c * p.daxiao + k * p.daxiao * 3] = (double)mat.get(r, c);
            }
        }
        k++;
    }
}