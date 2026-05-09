#include "mex.h"
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <sstream>
#include <string.h>
#include <functional>

// ===================================================================
// 数据结构定义：增加 Type 字段
// ===================================================================

struct Edge {
    int u;
    int v;
    int type; // 0 或 1

    // 排序运算符：优先 u，其次 v，最后 type
    bool operator<(const Edge& other) const {
        if (u != other.u) return u < other.u;
        if (v != other.v) return v < other.v;
        return type < other.type;
    }
    
    bool operator==(const Edge& other) const {
        return u == other.u && v == other.v && type == other.type;
    }
};

using EdgeList = std::vector<Edge>;

// 哈希函数需要包含 type
struct EdgeHash {
    std::size_t operator()(const Edge& e) const {
        // 简单的组合哈希
        std::size_t h1 = std::hash<int>()(e.u);
        std::size_t h2 = std::hash<int>()(e.v);
        std::size_t h3 = std::hash<int>()(e.type);
        return h1 ^ (h2 << 1) ^ (h3 << 2); 
    }
};

using ResultSet = std::vector<EdgeList>;
using ForbiddenSet = std::unordered_set<Edge, EdgeHash>;

// ===================================================================
// 并查集 (DSU)
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
// 核心辅助函数：哈希与规范化
// ===================================================================
std::string hash_edges_internal(const EdgeList& edges) {
    std::stringstream ss;
    EdgeList sorted_edges = edges;
    std::sort(sorted_edges.begin(), sorted_edges.end());
    for (const auto& edge : sorted_edges) {
        ss << edge.u << "," << edge.v << "," << edge.type << ";";
    }
    return ss.str();
}

// 【修改】pingfanyiout：只变换 u,v，type 跟随
EdgeList pingfanyiout_internal(const EdgeList& edges, int n1) {
    if (edges.empty()) return {};
    
    // 内部 lambda：确保 u < v，type 不变
    auto normalize_edge = [](Edge e) {
        if (e.u > e.v) std::swap(e.u, e.v);
        return e;
    };

    EdgeList best_edges = edges;
    for (auto& edge : best_edges) {
        edge = normalize_edge(edge);
    }
    std::sort(best_edges.begin(), best_edges.end());

    EdgeList current_edges(edges.size());

    // 1. 旋转
    for (int i = 0; i < n1; ++i) {
        for (size_t j = 0; j < edges.size(); ++j) {
            int new_u = (edges[j].u + i) % n1;
            int new_v = (edges[j].v + i) % n1;
            current_edges[j] = normalize_edge({ new_u, new_v, edges[j].type });
        }
        std::sort(current_edges.begin(), current_edges.end());
        if (current_edges < best_edges) {
            best_edges = current_edges;
        }
    }

    // 2. 翻转 + 旋转
    for (int i = 0; i < n1; ++i) {
        for (size_t j = 0; j < edges.size(); ++j) {
            int u_mir = (n1 - edges[j].u + n1) % n1;
            int v_mir = (n1 - edges[j].v + n1) % n1;
            int new_u = (u_mir + i) % n1;
            int new_v = (v_mir + i) % n1;
            current_edges[j] = normalize_edge({ new_u, new_v, edges[j].type });
        }
        std::sort(current_edges.begin(), current_edges.end());
        if (current_edges < best_edges) {
            best_edges = current_edges;
        }
    }
    return best_edges;
}

// ===================================================================
// MATLAB <-> C++ 转换辅助
// ===================================================================
static EdgeList mxArrayToEdgeList(const mxArray* arr) {
    double* data = mxGetPr(arr);
    mwSize m = mxGetM(arr); // 行数
    mwSize n = mxGetN(arr); // 列数，应该是 3
    EdgeList edges;
    edges.reserve(m);
    for (mwSize i = 0; i < m; ++i) {
        int u = static_cast<int>(data[i]);
        int v = static_cast<int>(data[i + m]);
        // 处理第三列 type，如果没有第三列默认为0（防崩溃），但逻辑上应该有
        int t = (n >= 3) ? static_cast<int>(data[i + 2 * m]) : 0;
        edges.push_back({ std::min(u, v), std::max(u, v), t });
    }
    return edges;
}

static mxArray* edgeListToMxArray(const EdgeList& edge_list) {
    // 输出 3 列
    mxArray* edges_mx = mxCreateDoubleMatrix(edge_list.size(), 3, mxREAL);
    double* outData = mxGetPr(edges_mx);
    for (size_t j = 0; j < edge_list.size(); ++j) {
        outData[j] = static_cast<double>(edge_list[j].u);
        outData[j + edge_list.size()] = static_cast<double>(edge_list[j].v);
        outData[j + 2 * edge_list.size()] = static_cast<double>(edge_list[j].type);
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
// 递归构建与搜索
// ===================================================================
void recursive_build(
    EdgeList& current_edges, const EdgeList& candidate_pairs, size_t start_index,
    ResultSet& results, int n1, int n2, int daxiao,
    std::vector<char>& used_nodes, ForbiddenSet& forbidden_pairs, DSU& dsu,
    const std::vector<int>& target_partition_asc
) {
    if (current_edges.size() + (candidate_pairs.size() - start_index) < daxiao) {
        return;
    }

    // === 终止条件 ===
    if (current_edges.size() == daxiao) {
        // 【新增检查】必须至少含有一个 Type 1 的边
        bool has_type1 = false;
        for (const auto& e : current_edges) {
            if (e.type == 1) {
                has_type1 = true;
                break;
            }
        }
        if (!has_type1) return;

        // 检查分区大小
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

    // === 剪枝：当前分区状态 ===
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

    // === 搜索循环 ===
    for (size_t i = start_index; i < candidate_pairs.size(); ++i) {
        const auto& e = candidate_pairs[i];
        int u = e.u;
        int v = e.v;
        int type = e.type;

        // 约束：顺序 (只比较 u)
        if (!current_edges.empty()) {
            if (u < current_edges.back().u) continue; 
            if (u == current_edges.back().u && v <= current_edges.back().v) continue; 
            // 这里严格控制顺序以避免重复
        }

        // Look-ahead 剪枝 (基于 u 的节点可用性)
        int edges_still_needed = daxiao - static_cast<int>(current_edges.size()) - 1;
        if (edges_still_needed > 0) {
             int potential = 0;
             // 简单估算剩余节点能提供的边的上限
             int nodes_left = n1 - 1 - u;
             if (nodes_left / 2 < edges_still_needed) return; 
        }

        if (used_nodes[u] || used_nodes[v] || forbidden_pairs.count(e)) {
            continue;
        }

        // === 应用选择 ===
        current_edges.push_back(e);
        used_nodes[u] = 1;
        used_nodes[v] = 1;

        // 【修改】DSU 连接逻辑
        if (type == 1) {
            // Type 1: 衍生图 mod([u,v,1]) 和 mod([u+1,v+1,1])
            // 意味着连接节点 u 和 v，以及 u+1 和 v+1
            // 注意：DSU操作的是节点索引
            dsu.unite(u, v);
            dsu.unite((u + 1) % n1, (v + 1) % n1);
        } else {
            // Type 0: 衍生图 mod([u+1,v,0]) 和 mod([u,v+1,0])
            // 意味着连接节点 u+1 和 v，以及 u 和 v+1
            dsu.unite(u, (v + 1) % n1);
            dsu.unite((u + 1) % n1, v);
        }

        // 【修改】计算矛盾边并加入 Forbidden Set
        Edge f1, f2;
        if (type == 0) {
            // mod([u-1, v+1, 0])
            int u_m1 = (u - 1 + n1) % n1;
            int v_p1 = (v + 1) % n1;
            f1 = { std::min(u_m1, v_p1), std::max(u_m1, v_p1), 0 };

            // mod([u+1, v-1, 0])
            int u_p1 = (u + 1) % n1;
            int v_m1 = (v - 1 + n1) % n1;
            f2 = { std::min(u_p1, v_m1), std::max(u_p1, v_m1), 0 };
        } else { // type == 1
            // mod([u-1, v-1, 1])
            int u_m1 = (u - 1 + n1) % n1;
            int v_m1 = (v - 1 + n1) % n1;
            f1 = { std::min(u_m1, v_m1), std::max(u_m1, v_m1), 1 };

            // mod([u+1, v+1, 1])
            int u_p1 = (u + 1) % n1;
            int v_p1 = (v + 1) % n1;
            f2 = { std::min(u_p1, v_p1), std::max(u_p1, v_p1), 1 };
        }

        forbidden_pairs.insert(f1);
        forbidden_pairs.insert(f2);

        // 递归
        recursive_build(current_edges, candidate_pairs, i + 1, results, n1, n2, daxiao,
            used_nodes, forbidden_pairs, dsu, target_partition_asc);

        // === 回溯 (撤销) ===
        forbidden_pairs.erase(f1);
        forbidden_pairs.erase(f2);
        dsu.undo(); // 撤销两个 unite 操作
        dsu.undo();
        used_nodes[u] = 0;
        used_nodes[v] = 0;
        current_edges.pop_back();
    }
}

// ===================================================================
// 顶层搜索入口
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
// MEX 接口
// ===================================================================
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 1) mexErrMsgIdAndTxt("MyToolbox:main", "Command required.");
    char command[64];
    mxGetString(prhs[0], command, sizeof(command));

    if (strcmp(command, "build") == 0) {
        if (nrhs != 7) mexErrMsgIdAndTxt("MyToolbox:build", "Args: cmd, init, all_pairs, n1, n2, daxiao, parts");
        
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

        // 初始化状态
        for (const auto& edge : initial_edges) {
            used_nodes[edge.u] = 1;
            used_nodes[edge.v] = 1;
            int u = edge.u; 
            int v = edge.v;
            int type = edge.type;

            // 初始化 DSU
            if (type == 1) {
                dsu.unite(u, v);
                dsu.unite((u + 1) % n1, (v + 1) % n1);
            } else {
                dsu.unite(u, (v + 1) % n1);
                dsu.unite((u + 1) % n1, v);
            }

            // 初始化 Forbidden
            if (type == 0) {
                int u_m1 = (u - 1 + n1) % n1; int v_p1 = (v + 1) % n1;
                int u_p1 = (u + 1) % n1; int v_m1 = (v - 1 + n1) % n1;
                forbidden_pairs.insert({ std::min(u_m1, v_p1), std::max(u_m1, v_p1), 0 });
                forbidden_pairs.insert({ std::min(u_p1, v_m1), std::max(u_p1, v_m1), 0 });
            } else {
                int u_m1 = (u - 1 + n1) % n1; int v_m1 = (v - 1 + n1) % n1;
                int u_p1 = (u + 1) % n1; int v_p1 = (v + 1) % n1;
                forbidden_pairs.insert({ std::min(u_m1, v_m1), std::max(u_m1, v_m1), 1 });
                forbidden_pairs.insert({ std::min(u_p1, v_p1), std::max(u_p1, v_p1), 1 });
            }
        }
        dsu.history.clear(); // 清除初始历史，避免回溯到初始状态之前

        search_all_partitions(initial_edges, all_pairs, results_cpp, n1, n2, daxiao,
            used_nodes, forbidden_pairs, dsu, all_target_partitions);

        plhs[0] = mxCreateCellMatrix(results_cpp.size(), 1);
        for (size_t i = 0; i < results_cpp.size(); ++i) {
            mxSetCell(plhs[0], i, edgeListToMxArray(results_cpp[i]));
        }
    }
    else if (strcmp(command, "hash") == 0) {
        EdgeList edges = mxArrayToEdgeList(prhs[1]);
        std::string hash_str = hash_edges_internal(edges);
        plhs[0] = mxCreateString(hash_str.c_str());
    }
    else if (strcmp(command, "pingfan") == 0) {
        EdgeList edges = mxArrayToEdgeList(prhs[1]);
        int n1 = static_cast<int>(mxGetScalar(prhs[2]));
        EdgeList result_edges = pingfanyiout_internal(edges, n1);
        plhs[0] = edgeListToMxArray(result_edges);
    }
}