#ifndef BATCHENUM_H
#define BATCHENUM_H

#include "Topo/TopoLogical.h"
#include "DirectedGraph.h"
#include "QuerySharingGraph.h"

class BatchEnum {
public:
    static void Search(DirectedGraph& G, std::vector<Path>& P, std::unordered_map<Query, std::vector<Path>>& R, Path& path, Query& q, QuerySharingGraph& Psi) {
    // 1. 获取当前路径的最后一个顶点 v'
    int v_prime = path.vertices.back();
    // 2. 将当前路径 p 添加到结果集 P 中
    P.push_back(path);

    // 3. 检查路径长度
    if (path.vertices.size() - 1 == q.k) {
        return;
    }

    // 4. 遍历当前顶点 v' 的所有邻居 v''
    for (int v_double_prime : G.getNeighbors(v_prime)) {
        // 5. 检查邻居 v'' 是否满足跳数约束 1->0！！！！！！！！！！！
        if (0 <= q.k - path.vertices.size()) {
            // 6. 如果邻居 v'' 已经在当前路径 p 中，则跳过
            if (std::find(path.vertices.begin(), path.vertices.end(), v_double_prime) != path.vertices.end()) {
                continue;
            }

            // 7. 如果存在一个查询 q'，其目标顶点 q'.v 等于 v''，则将 p 与 R[q'] 的路径合并并添加到 P 中
            bool found = false;
            for (auto& q_prime : Psi.getInNeighbors(q)) {
                if (q_prime.t == v_double_prime) {
                    for (auto& r_path : R[q_prime]) {
                        Path combinedPath = path;
                        combinedPath.vertices.insert(combinedPath.vertices.end(), r_path.vertices.begin(), r_path.vertices.end());
                        if (G.isSimplePath(combinedPath)) {
                            P.push_back(combinedPath);
                        }
                    }
                    found = true;
                    break;
                }
            }

            // 8. 否则，递归调用 Search 函数继续搜索
            if (!found) {
                path.vertices.push_back(v_double_prime);
                Search(G, P, R, path, q, Psi);
                path.vertices.pop_back();
            }
        }
    }
    }

    // Batch查询算法
    static void BatchEnum_notclass(DirectedGraph& G, std::vector<Query>& Q, double gamma) {
        std::cout << "Batch Enumeration" << std::endl;
        std::vector<Path> Output;
        // 多源BFS 1-2
        std::unordered_set<int> S, T;
        for (const Query& q : Q) {
            S.insert(q.s);
            T.insert(q.t);
        }
        // std::vector<Query> Q_r;
        // for (auto& q : Q) {
        //     Query q_r(q.t, q.s, q.k);
        //     Q_r.push_back(q_r);
        // }

        G.initializeBFSIndex(S, T);
        DirectedGraph G_r;
        G_r = G.getReverseGraph();
        G_r.initializeBFSIndex(T, S);
        
        // 初始化缓存和聚类查询
        std::unordered_map<Query, std::vector<Path>> R; // 缓存类型待定,查询和对应的路径集合
        std::vector<QueryGroup> Cs = G.hierarchicalClustering(Q, gamma);
        // std::vector<QueryGroup> Cs_r = G_r.hierarchicalClustering(Q_r, gamma);

        // 构建查询共享图
        for (auto& C : Cs) {
            QuerySharingGraph Psi;
            QuerySharingGraph Psir;
            QueryGroup C_r;
            for (auto& q : C.queries) {
                Query q_r(q.t, q.s, q.k);
                C_r.queries.push_back(q_r);
            }
            Psi.DetectCommonQuery(G, C, Psi);
            Psir.DetectCommonQuery_r(G_r, C_r, Psir);
            Psi.printEdges();
            Psir.printEdges();

            // 查询
            std::unordered_set<Query> processedQueries;
            // 按照拓扑顺序
            // 正向图
            TopoLogical Topo(Psi);
            TopoLogical Topo_r(Psir);
            auto topoOrder = Topo.getorder();
            auto topoOrder_r = Topo_r.getorder();

            // 交替执行正向和反向图的遍历
            size_t i = 0, j = 0;
            while (i < topoOrder.size() || j < topoOrder_r.size()) {
            // while (i < 10 || j < 10) {
                // for (auto& r : R) {
                //     std::cout << "R: " << r.first.s << "->" << r.first.t << " k=" << r.first.k << std::endl;
                //     // for (auto& p : r.second) {
                //     //     std::cout << "Path: ";
                //     //     for (auto& vertex : p.vertices) {
                //     //         std::cout << vertex << " -> ";
                //     //     }
                //     //     std::cout << "\b\b" << std::endl; // 删除最后的 " -> "
                //     // }
                // }
                // 正向图
                if (i < topoOrder.size()) {
                    auto &q = topoOrder[i];
                    // std::cout << "q:" << q.s << "->" << q.t << " k=" << q.k << std::endl;
                    
                    if (q.t == 0) { // HC-s 路径查询
                        std::vector<Path> P;
                        Path initialPath;
                        initialPath.vertices.push_back(q.s);
                        Search(G, P, R, initialPath, q, Psi);
                        R.insert({q, P});
                    }

                    if (q.t != 0 && R.find({q.s, 0, ceil(static_cast<double>(q.k)/2)}) != R.end() && R.find({q.t, 0, floor(static_cast<double>(q.k)/2)}) != R.end()) { // HC-s-t 路径查询
                        processedQueries.insert({q.s, 0, ceil(static_cast<double>(q.k)/2)});
                        processedQueries.insert({q.t, 0, floor(static_cast<double>(q.k)/2)});
                        std::vector<Path> forwardPaths;
                        std::vector<Path> backwardPaths;
                        for (auto& p1 : R[{q.s, 0, ceil(static_cast<double>(q.k)/2)}]) {
                            forwardPaths.push_back(p1);
                        }
                        for (auto& p2 : R[{q.t, 0, floor(static_cast<double>(q.k)/2)}]) {
                            backwardPaths.push_back(p2);
                        }
                        std::vector<Path> combinedPaths;
                        G.connectPaths(forwardPaths, backwardPaths, combinedPaths);
                        for (const auto& combinedPath : combinedPaths) {
                            if (G.isSimplePath(combinedPath) && G.isPathSatisfyConstraints(combinedPath, q.k)) {
                                Output.push_back(combinedPath);
                                processedQueries.insert(q);
                            }
                        }
                    }

                    int tmp = 0;
                    for (auto& q_prime : Psi.getInNeighbors(q)) {
                        for (auto& q_prime_two : Psi.getNeighbors(q_prime)) {
                            if (processedQueries.find(q_prime_two) == processedQueries.end()) {
                                tmp = 1;
                                break;
                            }
                        }
                        if (tmp != 1) {
                            R.erase(q_prime);
                        }
                    }
                }
                // 反向图
                if (j < topoOrder_r.size()) {
                    auto &q = topoOrder_r[j];
                    
                    if (q.t == 0) { // HC-s 路径查询
                        std::vector<Path> P;
                        Path initialPath;
                        initialPath.vertices.push_back(q.s);
                        Search(G_r, P, R, initialPath, q, Psir);
                        R.insert({q, P});
                    }

                    int tmp = 0;
                    for (const auto& q_prime : Psi.getInNeighbors(q)) {
                        for (auto& q_prime_two : Psi.getNeighbors(q_prime)) {
                            if (processedQueries.find(q_prime_two) == processedQueries.end()) {
                                tmp = 1;
                                break;
                            }
                        }
                        if (tmp != 1) {
                            R.erase(q_prime);
                        }
                    }
                }
                i++;
                j++;
            }

            i = 0;
            while (i < topoOrder.size()) {
                auto &q = topoOrder[i];
                if (processedQueries.find(q) == processedQueries.end() && q.t != 0 && R.find({q.s, 0, ceil(static_cast<double>(q.k)/2)}) != R.end() && R.find({q.t, 0, floor(static_cast<double>(q.k)/2)}) != R.end()) { // HC-s-t 路径查询
                        processedQueries.insert({q.s, 0, ceil(static_cast<double>(q.k)/2)});
                        processedQueries.insert({q.t, 0, floor(static_cast<double>(q.k)/2)});
                        std::vector<Path> forwardPaths;
                        std::vector<Path> backwardPaths;
                        for (auto& p1 : R[{q.s, 0, ceil(static_cast<double>(q.k)/2)}]) {
                            forwardPaths.push_back(p1);
                        }
                        for (auto& p2 : R[{q.t, 0, floor(static_cast<double>(q.k)/2)}]) {
                            backwardPaths.push_back(p2);
                        }
                        std::vector<Path> combinedPaths;
                        G.connectPaths(forwardPaths, backwardPaths, combinedPaths);
                        for (const auto& combinedPath : combinedPaths) {
                            // std::cout << "combinedPaths: " << combinedPath.vertices.size() << std::endl;
                            if (G.isSimplePath(combinedPath) && G.isPathSatisfyConstraints(combinedPath, q.k)) {
                                Output.push_back(combinedPath);
                                processedQueries.insert(q);
                            }
                        }
                    }
                    i++;
            }
        }
        std::cout << "\n";
        for (const auto& path : Output) {
            std::cout << "Path: ";
            for (const auto& vertex : path.vertices) {
                std::cout << vertex << " ";
            }
            std::cout << std::endl; // 删除最后的 " -> "
        }
    }
};
#endif // BATCHENUM_H