// =========================================================================
//                  文件: build_mex.cpp (完整且已修正)
// =========================================================================
#include "mex.h"
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <sstream>
#include <string.h>
#include <functional> // 用于 std::greater

// 定义图相关的基本类型
using Edge = std::pair<int, int>;
using EdgeList = std::vector<Edge>;
// 用于在 unordered_set 中存储 Edge
struct EdgeHash {
    std::size_t operator()(const Edge& e) const {
        return std::hash<int>()(e.first) ^ (std::hash<int>()(e.second) << 1);
    }
};
// 定义高级类型
using ResultSet = std::vector<EdgeList>;
using ForbiddenSet = std::unordered_set<Edge, EdgeHash>;

// ===================================================================
// 并查集 (DSU) 数据结构 (带撤销功能)
// ===================================================================
struct DSU {
    std::vector<int> parent;
    std::vector<int> component_size;
    struct UndoInfo {
        int u, v, old_parent_v, old_size_u;
    };
    std::vector<UndoInfo> history;

    DSU(int n) {
        parent.resize(n);
        std::iota(parent.begin(), parent.end(), 0);
        component_size.assign(n, 1);
    }

    int find(int i) {
        if (parent[i] == i) return i;
        return find(parent[i]);
    }

    void unite(int i, int j) {
        int root_i = find(i);
        int root_j = find(j);
        if (root_i != root_j) {
            if (component_size[root_i] < component_size[root_j]) std::swap(root_i, root_j);
            history.push_back({ root_i, root_j, parent[root_j], component_size[root_i] });
            parent[root_j] = root_i;
            component_size[root_i] += component_size[root_j];
        }
        else {
            history.push_back({ -1, -1, -1, -1 });
        }
    }

    void undo() {
        if (!history.empty()) {
            UndoInfo last = history.back();
            history.pop_back();
            if (last.u != -1) {
                parent[last.v] = last.old_parent_v;
                component_size[last.u] = last.old_size_u;
            }
        }
    }
};

// ===================================================================
// 核心辅助函数
// ===================================================================
std::string hash_edges_internal(const EdgeList& edges) {
    std::stringstream ss;
    EdgeList sorted_edges = edges;
    std::sort(sorted_edges.begin(), sorted_edges.end());
    for (const auto& edge : sorted_edges) {
        ss << edge.first << "," << edge.second << ";";
    }
    return ss.str();
}

EdgeList pingfanyiout_internal(const EdgeList& edges, int n1) {
    if (edges.empty()) return {};
    auto normalize_edge = [](Edge e) {
        if (e.first > e.second) std::swap(e.first, e.second);
        return e;
        };
    EdgeList best_edges = edges;
    for (auto& edge : best_edges) {
        edge = normalize_edge(edge);
    }
    std::sort(best_edges.begin(), best_edges.end());

    EdgeList current_edges(edges.size());
    for (int i = 0; i < n1; ++i) {
        for (size_t j = 0; j < edges.size(); ++j) {
            current_edges[j] = normalize_edge({ (edges[j].first + i) % n1, (edges[j].second + i) % n1 });
        }
        std::sort(current_edges.begin(), current_edges.end());
        if (current_edges < best_edges) {
            best_edges = current_edges;
        }
    }
    for (int i = 0; i < n1; ++i) {
        for (size_t j = 0; j < edges.size(); ++j) {
            int u = (n1 - edges[j].first + n1) % n1;
            int v = (n1 - edges[j].second + n1) % n1;
            current_edges[j] = normalize_edge({ (u + i) % n1, (v + i) % n1 });
        }
        std::sort(current_edges.begin(), current_edges.end());
        if (current_edges < best_edges) {
            best_edges = current_edges;
        }
    }
    return best_edges;
}

// ===================================================================
// 辅助函数 (MATLAB <-> C++)
// ===================================================================
static EdgeList mxArrayToEdgeList(const mxArray* arr) {
    double* data = mxGetPr(arr);
    mwSize m = mxGetM(arr);
    EdgeList edges;
    edges.reserve(m);
    for (mwSize i = 0; i < m; ++i) {
        int u = static_cast<int>(data[i]);
        int v = static_cast<int>(data[i + m]);
        edges.push_back({ std::min(u, v), std::max(u, v) });
    }
    return edges;
}

static mxArray* edgeListToMxArray(const EdgeList& edge_list) {
    mxArray* edges_mx = mxCreateDoubleMatrix(edge_list.size(), 2, mxREAL);
    double* outData = mxGetPr(edges_mx);
    for (size_t j = 0; j < edge_list.size(); ++j) {
        outData[j] = static_cast<double>(edge_list[j].first);
        outData[j + edge_list.size()] = static_cast<double>(edge_list[j].second);
    }
    return edges_mx;
}

static std::vector<std::vector<int>> mxArrayToVectorOfVectors(const mxArray* arr) {
    double* data = mxGetPr(arr);
    mwSize num_partitions = mxGetM(arr);
    mwSize n2 = mxGetN(arr);
    std::vector<std::vector<int>> all_partitions;
    all_partitions.reserve(num_partitions);
    for (mwSize i = 0; i < num_partitions; ++i) {
        std::vector<int> partition;
        partition.reserve(n2);
        for (mwSize j = 0; j < n2; ++j) {
            partition.push_back(static_cast<int>(data[i + j * num_partitions]));
        }
        all_partitions.push_back(partition);
    }
    return all_partitions;
}

// ===================================================================
// 核心算法: recursive_build (完整优化版)
// ===================================================================
void recursive_build(
    EdgeList& current_edges, const EdgeList& candidate_pairs, size_t start_index,
    ResultSet& results, int n1, int n2, int daxiao,
    std::vector<char>& used_nodes, ForbiddenSet& forbidden_pairs, DSU& dsu,
    const std::vector<int>& target_partition_asc
) {
    // 剪枝 1: 基础数量剪枝
    if (current_edges.size() + (candidate_pairs.size() - start_index) < daxiao) {
        return;
    }

    // 终止条件
    if (current_edges.size() == daxiao) {
        std::vector<int> final_partition_sizes;
        for (int i = 0; i < n1; ++i) {
            if (dsu.parent[i] == i) {
                final_partition_sizes.push_back(dsu.component_size[i]);
            }
        }
        if (final_partition_sizes.size() == n2) {
            std::sort(final_partition_sizes.begin(), final_partition_sizes.end());
            if (final_partition_sizes == target_partition_asc) {
                results.push_back(current_edges);
            }
        }
        return;
    }

    // 剪枝 2: 基于当前划分状态的剪枝
    bool should_prune = false;
    std::vector<int> current_partition_sizes;
    for (int k = 0; k < n1; ++k) {
        if (dsu.parent[k] == k) {
            current_partition_sizes.push_back(dsu.component_size[k]);
        }
    }
    std::sort(current_partition_sizes.begin(), current_partition_sizes.end(), std::greater<int>());
    std::vector<int> target_partition_desc = target_partition_asc;
    std::sort(target_partition_desc.begin(), target_partition_desc.end(), std::greater<int>());

    if (!current_partition_sizes.empty()) {
        if (current_partition_sizes[0] > target_partition_desc[0]) {
            should_prune = true;
        }
        if (!should_prune) {
            if (current_partition_sizes[0] == target_partition_desc[0]) {
                for (size_t k = 1; k < current_partition_sizes.size(); ++k) {
                    if (k < target_partition_desc.size() && current_partition_sizes[k] > target_partition_desc[k]) {
                        should_prune = true;
                        break;
                    }
                }
            }
            else {
                if (current_partition_sizes.size() > 1 && target_partition_desc.size() > 1) {
                    if (current_partition_sizes[1] > target_partition_desc[1]) {
                        if (current_partition_sizes[0] + current_partition_sizes[1] > target_partition_desc[0]) {
                            should_prune = true;
                        }
                    }
                }
            }
        }
    }
    if (should_prune) {
        return;
    }

    // --- 深度优先搜索主循环 ---
    for (size_t i = start_index; i < candidate_pairs.size(); ++i) {
        const auto& e = candidate_pairs[i];
        int u = e.first;
        int v = e.second;

        // 约束条件: 新边的第一个顶点 u 必须严格大于上一条边的第一个顶点。
        if (!current_edges.empty()) {
            if (u <= current_edges.back().first) {
                continue;
            }
        }

        // 剪枝 3: 精确的前瞻剪枝 (Look-ahead Pruning)
        int edges_still_needed = daxiao - static_cast<int>(current_edges.size()) - 1;
        if (edges_still_needed > 0) {
            int potential_nodes_count = (n1 - 1) - u;
            int used_nodes_count_in_pool = 1; // 计入当前的 v
            for (const auto& prev_edge : current_edges) {
                if (prev_edge.second > u) {
                    used_nodes_count_in_pool++;
                }
            }
            int nodes_available = potential_nodes_count - used_nodes_count_in_pool;
            int max_future_edges = (nodes_available > 0) ? (nodes_available / 2) : 0;
            if (max_future_edges < edges_still_needed) {
                return;
            }
        }

        if (used_nodes[u] || used_nodes[v] || forbidden_pairs.count(e)) {
            continue;
        }

        // --- 做出选择 ---
        current_edges.push_back(e);
        used_nodes[u] = 1;
        used_nodes[v] = 1;
        dsu.unite(u, (v + 1) % n1);
        dsu.unite((u + 1) % n1, v);
        int u_m1 = (u - 1 + n1) % n1, v_p1 = (v + 1) % n1;
        int u_p1 = (u + 1) % n1, v_m1 = (v - 1 + n1) % n1;
        Edge f1 = { std::min(u_m1, v_p1), std::max(u_m1, v_p1) };
        Edge f2 = { std::min(u_p1, v_m1), std::max(u_p1, v_m1) };
        forbidden_pairs.insert(f1);
        forbidden_pairs.insert(f2);

        // --- 进入下一层递归 ---
        recursive_build(current_edges, candidate_pairs, i + 1, results, n1, n2, daxiao,
            used_nodes, forbidden_pairs, dsu, target_partition_asc);

        // --- 回溯 (撤销选择) ---
        forbidden_pairs.erase(f1);
        forbidden_pairs.erase(f2);
        dsu.undo();
        dsu.undo();
        used_nodes[u] = 0;
        used_nodes[v] = 0;
        current_edges.pop_back();
    }
}

// ===================================================================
// 【重要】顶层搜索函数 (必须在 mexFunction 之前定义)
// ===================================================================
void search_all_partitions(
    const EdgeList& initial_edges, const EdgeList& candidate_pairs,
    ResultSet& results, int n1, int n2, int daxiao,
    std::vector<char>& base_used_nodes, ForbiddenSet& base_forbidden_pairs, DSU& base_dsu,
    const std::vector<std::vector<int>>& all_target_partitions
) {
    for (const auto& target_partition_asc : all_target_partitions) {
        EdgeList current_edges = initial_edges;
        recursive_build(current_edges, candidate_pairs, 0, results, n1, n2, daxiao,
            base_used_nodes, base_forbidden_pairs, base_dsu, target_partition_asc);
    }
}

// ===================================================================
// MEX 网关函数
// ===================================================================
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("MyToolbox:main:nrhs", "At least one argument (command) is required.");
    }
    char command[64];
    mxGetString(prhs[0], command, sizeof(command));

    if (strcmp(command, "build") == 0) {
        if (nrhs != 7) {
            mexErrMsgIdAndTxt("MyToolbox:build:nrhs", "build command requires 6 arguments: initial_edges, all_pairs, n1, n2, daxiao, all_valid_partitions_matrix");
        }
        EdgeList initial_edges = mxArrayToEdgeList(prhs[1]);
        EdgeList all_pairs = mxArrayToEdgeList(prhs[2]);
        int n1 = static_cast<int>(mxGetScalar(prhs[3]));
        int n2 = static_cast<int>(mxGetScalar(prhs[4]));
        int daxiao = static_cast<int>(mxGetScalar(prhs[5]));
        std::vector<std::vector<int>> all_target_partitions = mxArrayToVectorOfVectors(prhs[6]);

        ResultSet results_cpp;
        std::vector<char> used_nodes(n1, 0);
        ForbiddenSet forbidden_pairs;
        DSU dsu(n1);

        for (const auto& edge : initial_edges) {
            used_nodes[edge.first] = 1;
            used_nodes[edge.second] = 1;
            dsu.unite(edge.first, (edge.second + 1) % n1);
            dsu.unite((edge.first + 1) % n1, edge.second);
            int u_m1 = (edge.first - 1 + n1) % n1, v_p1 = (edge.second + 1) % n1;
            int u_p1 = (edge.first + 1) % n1, v_m1 = (edge.second - 1 + n1) % n1;
            forbidden_pairs.insert({ std::min(u_m1, v_p1), std::max(u_m1, v_p1) });
            forbidden_pairs.insert({ std::min(u_p1, v_m1), std::max(u_p1, v_m1) });
        }
        dsu.history.clear();

        // 现在调用 search_all_partitions 是安全的，因为它已经在上面被定义了
        search_all_partitions(initial_edges, all_pairs, results_cpp, n1, n2, daxiao,
            used_nodes, forbidden_pairs, dsu, all_target_partitions);

        plhs[0] = mxCreateCellMatrix(results_cpp.size(), 1);
        for (size_t i = 0; i < results_cpp.size(); ++i) {
            mxSetCell(plhs[0], i, edgeListToMxArray(results_cpp[i]));
        }
    }
    else if (strcmp(command, "hash") == 0) {
        if (nrhs != 2) mexErrMsgIdAndTxt("MyToolbox:hash:nrhs", "hash command needs 1 extra arg: edge_list.");
        EdgeList edges = mxArrayToEdgeList(prhs[1]);
        std::string hash_str = hash_edges_internal(edges);
        plhs[0] = mxCreateString(hash_str.c_str());
    }
    else if (strcmp(command, "pingfan") == 0) {
        if (nrhs != 3) mexErrMsgIdAndTxt("MyToolbox:pingfan:nrhs", "pingfan command needs 2 extra args: edge_list, n1.");
        EdgeList edges = mxArrayToEdgeList(prhs[1]);
        int n1 = static_cast<int>(mxGetScalar(prhs[2]));
        EdgeList result_edges = pingfanyiout_internal(edges, n1);
        plhs[0] = edgeListToMxArray(result_edges);
    }
    else {
        mexErrMsgIdAndTxt("MyToolbox:main:unknownCommand", "Unknown command. Valid are 'build', 'hash', 'pingfan'.");
    }
}