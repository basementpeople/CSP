#include "QuerySharingGraph.h"

// HC-s查询
// struct Query {
//     int s;
//     int t; t = 0 表示 HC-s 查询
//     int k;
// };

// QuerySharingGraph::QuerySharingGraph() {
//     // 初始化代码
//     m = 0;
//     n = 0;
// }
QuerySharingGraph::QuerySharingGraph() : DirectedGraph(), m(0), n(0) {}

bool QuerySharingGraph::addEdge(const Query& q1, const Query& q2) {
    if (q1 != q2 && adj[q1].find(q2) == adj[q1].end()) {
        edges.insert({q1, q2});
        adj[q1].insert(q2);
        m++;
        return true;
    } else {
        return false;
    }    
}

bool QuerySharingGraph::isSubquery(Query& q1, Query& q2, DirectedGraph& G) {
   // 判断 q1 是否是 q2 的子查询
   // 这里需要根据具体需求实现
   // q1.k <= q2.k - dist (q2.s, q1.s)
   if (q1.k <= q2.k - G.dist(q2.s, q1.s)) {
       return true;
   }
   return false;
}

// void QuerySharingGraph::DetectCommonQuery(DirectedGraph& G, std::vector<Query>& Q, QuerySharingGraph& Psi) {
//    std::unordered_map<int, std::set<Query>> ME;
//    std::unordered_map<int, Query> MQ;
//    std::queue<int> que;

//    // Step 1: Initialize
//    for (const auto& q : Q) {
//        Psi.addEdge({q.s, 0, static_cast<int> (ceil(q.k / 2.0))}, q);
//        ME[q.s].insert({q.s, 0, static_cast<int> (ceil(q.k / 2.0))});
//        que.push(q.s);
//    }

//    int kmax = 0;
//    for (const auto& q : Q) {
//        kmax = std::max(kmax, static_cast<int>(ceil(q.k / 2.0)));
//    }

//    // Step 2: Process each k from 0 to kmax
//    for (int k = 0; k <= kmax; ++k) {
//     std::cout << "k = " << k << std::endl;
//        for (auto& v : G.traverseAllVertices()) { // 假设 G 是一个有向图，这里需要遍历所有节点
//            std::set<Query> SQ;

//            // Step 3: Collect queries with k = kmax - k
//            for (const auto& q : ME[v]) {
//             std::cout << "q = (" << q.s << ", " << q.t << ", " << q.k << ")" << std::endl;
//                if (q.k == kmax - k) {
//                    SQ.insert(q);
//                    ME[v].erase(q);
//                }
//            }

//            // Step 4: Handle SQ
//            if (SQ.empty()) {
//             std::cout << "empty SQ" << std::endl;
//                continue;
//            } else if (SQ.size() == 1) {
//                MQ[v] = *SQ.begin();
//            } else {
//                for (const auto& q : SQ) {
//                    Psi.addEdge({v, 0, kmax - k}, q);
//                }
//                MQ[v] = {v, 0, kmax - k}; // ??
//            }

//            // Step 5: Update ME and MQ
//            std::cout << "打印第五步，点：" << v << std::endl;
//            for (const auto& v_prime : G.getNeighbors(v)) {
//                std::cout << "  邻居：" << v_prime << std::endl;
//                if (MQ.find(v_prime) != MQ.end() && !isSubquery(MQ[v], MQ[v_prime], G)) {
//                    Psi.addEdge(MQ[v_prime], MQ[v]);
//                    std::cout << "add edge (" << MQ[v_prime].s << ", " << MQ[v_prime].t << ", " << MQ[v_prime].k << ") -> (" << MQ[v].s << ", " << MQ[v].t << ", " << MQ[v].k << ")" << std::endl;
//                } else {
//                    ME[v_prime].insert(MQ[v]);
//                    std::cout << "insert MQ[" << v << "] = (" << MQ[v].s << ", " << MQ[v].t << ", " << MQ[v].k << ")" << std::endl;
//                }
//            }
//         }
//     }
// }

// 可以说这只是正向图的共享构建，因为如果k未奇数，就不能两两相连
void QuerySharingGraph::DetectCommonQuery(DirectedGraph& G, QueryGroup& Q, QuerySharingGraph& Psi) {
   std::unordered_map<int, std::set<Query>> ME;
   std::unordered_map<int, Query> MQ;
   std::queue<int> que;

   // Step 1: Initialize
   for (const auto& q : Q.queries) {
       Psi.addEdge({q.s, 0, static_cast<int> (ceil(q.k / 2.0))}, q);
       ME[q.s].insert({q.s, 0, static_cast<int> (ceil(q.k / 2.0))});
       que.push(q.s);
       std::cout << "dw" << q.s << std::endl;
   }

   int kmax = 0;
   for (const auto& q : Q.queries) {
       kmax = std::max(kmax, static_cast<int>(ceil(q.k / 2.0)));
   }

   // Step 2: Process each k from 0 to kmax
   for (int k = 0; k <= kmax; ++k) {
    int flag = que.size();
    std::cout << "k = " << k << std::endl;
    std::cout << "flag = " << flag << std::endl;
       for (int i = 0; i < flag; ++i) { // 假设 G 是一个有向图，这里需要遍历所有节点
           int v = que.front();
           que.pop();
           std::cout << "v = " << v << std::endl;
           std::set<Query> SQ;

           // Step 3: Collect queries with k = kmax - k
           for (const auto& q : ME[v]) {
            std::cout << "q = (" << q.s << ", " << q.t << ", " << q.k << ")" << std::endl;
            SQ.insert(q);
            ME[v].erase(q);
           }

           // Step 4: Handle SQ
           if (SQ.empty()) {
            std::cout << "empty SQ" << std::endl;
               continue;
           } else if (SQ.size() == 1) {
               MQ[v] = *SQ.begin();
           } else {
               for (const auto& q : SQ) {
                    // 扔掉k=1的边
                    if (kmax - k == 1) {
                        continue;
                        }
                    Psi.addEdge({v, 0, kmax - k}, q);
               }
               MQ[v] = {v, 0, kmax - k}; // ??
           }

           // Step 5: Update ME and MQ
           std::cout << "打印第五步，点：" << v << std::endl;
           for (const auto& v_prime : G.getNeighbors(v)) {
               std::cout << "  邻居：" << v_prime << std::endl;
               if (MQ.find(v_prime) != MQ.end() && isSubquery(MQ[v_prime], MQ[v], G) && MQ[v_prime] != MQ[v]) {
                   // 扔掉k=1的边
                   if (MQ[v_prime].k == 1) {
                       continue;
                   }
                   Psi.addEdge(MQ[v_prime], MQ[v]);
                   std::cout << "add edge (" << MQ[v_prime].s << ", " << MQ[v_prime].t << ", " << MQ[v_prime].k << ") -> (" << MQ[v].s << ", " << MQ[v].t << ", " << MQ[v].k << ")" << std::endl;
               } else {
                   ME[v_prime].insert(MQ[v]);
                   que.push(v_prime);
                   std::cout << "insert MQ[" << v << "] = (" << MQ[v].s << ", " << MQ[v].t << ", " << MQ[v].k << ")" << std::endl;
               }
           }
        }
    }
}

// 反向图的共享构建
void QuerySharingGraph::DetectCommonQuery_r(DirectedGraph& G, QueryGroup& Q, QuerySharingGraph& Psi) {
   std::unordered_map<int, std::set<Query>> ME;
   std::unordered_map<int, Query> MQ;
   std::queue<int> que;

   // Step 1: Initialize
   for (const auto& q : Q.queries) {
       Psi.addEdge({q.s, 0, static_cast<int> (floor(q.k / 2.0))}, q);
       ME[q.s].insert({q.s, 0, static_cast<int> (floor(q.k / 2.0))});
       que.push(q.s);
       std::cout << "dw" << q.s << std::endl;
   }

   int kmax = 0;
   for (const auto& q : Q.queries) {
       kmax = std::max(kmax, static_cast<int>(floor(q.k / 2.0)));
   }

   // Step 2: Process each k from 0 to kmax
   for (int k = 0; k <= kmax; ++k) {
    int flag = que.size();
    std::cout << "k = " << k << std::endl;
    std::cout << "flag = " << flag << std::endl;
       for (int i = 0; i < flag; ++i) { // 假设 G 是一个有向图，这里需要遍历所有节点
           int v = que.front();
           que.pop();
           std::cout << "v = " << v << std::endl;
           std::set<Query> SQ;

           // Step 3: Collect queries with k = kmax - k
           for (const auto& q : ME[v]) {
            std::cout << "q = (" << q.s << ", " << q.t << ", " << q.k << ")" << std::endl;
            SQ.insert(q);
            ME[v].erase(q);
           }

           // Step 4: Handle SQ
           if (SQ.empty()) {
            std::cout << "empty SQ" << std::endl;
               continue;
           } else if (SQ.size() == 1) {
               MQ[v] = *SQ.begin();
           } else {
               for (const auto& q : SQ) {
                    // 扔掉k=1的边
                    if (kmax - k == 1) {
                        continue;
                        }
                    Psi.addEdge({v, 0, kmax - k}, q);
               }
               MQ[v] = {v, 0, kmax - k}; // ??
           }

           // Step 5: Update ME and MQ
           std::cout << "打印第五步，点：" << v << std::endl;
           for (const auto& v_prime : G.getNeighbors(v)) {
               std::cout << "  邻居：" << v_prime << std::endl;
               if (MQ.find(v_prime) != MQ.end() && isSubquery(MQ[v_prime], MQ[v], G) && MQ[v_prime] != MQ[v]) {
                   // 扔掉k=1的边
                   if (MQ[v_prime].k == 1) {
                       continue;
                   }
                   Psi.addEdge(MQ[v_prime], MQ[v]);
                   std::cout << "add edge (" << MQ[v_prime].s << ", " << MQ[v_prime].t << ", " << MQ[v_prime].k << ") -> (" << MQ[v].s << ", " << MQ[v].t << ", " << MQ[v].k << ")" << std::endl;
               } else {
                   ME[v_prime].insert(MQ[v]);
                   que.push(v_prime);
                   std::cout << "insert MQ[" << v << "] = (" << MQ[v].s << ", " << MQ[v].t << ", " << MQ[v].k << ")" << std::endl;
               }
           }
        }
    }
}

void QuerySharingGraph::printEdges() {
    std::cout << "Edges in the Query Sharing Graph:" << std::endl;
    for (const auto& [q1, q2] : edges) {
        std::cout << "  Edge: (" << q1.s << ", " << q1.t << ", " << q1.k << ") -> (" << q2.s << ", " << q2.t << ", " << q2.k << ")" << std::endl;
    }
}


// std::vector<int> QuerySharingGraph::getNeighborsWithConstraint(int v, int k, DirectedGraph& G) {
//     std::cout << "Get neighbors of " << v << " with constraint k = " << k << std::endl;
//     std::set<int> neighborsSet;

//     for (const auto& neighbor : G.getNeighbors(v)) {
//         std::cout << "  Neighbor: " << neighbor << std::endl;
//         neighborsSet.insert(neighbor);
//         for (const auto& neigh : G.getNeighbors(neighbor))
//         {
//             if (G.dist(v, neigh) <= k) {
//                 neighborsSet.insert(neigh);
//                 std::cout << "  Neighbor: " << neigh << std::endl;
//                 if (G.dist(v, neigh) <= k - 1) {
//                     // 递归获取邻居并插入
//                     std::vector<int> subNeighbors = getNeighborsWithConstraint(neigh, k, G);
//                     neighborsSet.insert(subNeighbors.begin(), subNeighbors.end());
//                 }
//             }
//         }
//     }
//     // 将 set 转换为 vector
//     std::vector<int> neighbors(neighborsSet.begin(), neighborsSet.end());
//     return neighbors;
// }

std::vector<int> QuerySharingGraph::getNeighborsWithConstraint(int v, int k, DirectedGraph& G) {
    std::cout << "Get neighbors of " << v << " with constraint k = " << k << std::endl;
    std::set<int> neighborsSet;

    for (const auto& neighbor : G.getNeighbors(v)) {
        if (G.dist(v, neighbor) <= k) {
            std::cout << "  Neighbor: " << neighbor << std::endl;
            neighborsSet.insert(neighbor);
            if (k > 1) {
                // 递归获取邻居并插入
                std::vector<int> subNeighbors = getNeighborsWithConstraint(neighbor, k - 1, G);
                neighborsSet.insert(subNeighbors.begin(), subNeighbors.end());
            }
        }
    }
    
    // 将 set 转换为 vector
    std::vector<int> neighbors(neighborsSet.begin(), neighborsSet.end());
    return neighbors;
}

std::unordered_set<Query> QuerySharingGraph::getNodes()
{
    std::unordered_set<Query> nodes;
    for (auto &entry : adj)
    {
        nodes.insert(entry.first);
    }
    return nodes;
}

// 获取邻居
std::vector<Query> QuerySharingGraph::getNeighbors(Query vertex) const {
    auto it = adj.find(vertex);
    if (it != adj.end()) {
        return std::vector<Query>(it->second.begin(), it->second.end());
    }
    return {};
}

// 获取顶点的所有入边点
std::vector<Query> QuerySharingGraph::getInNeighbors(Query& vertex) {
    std::vector<Query> inNeighbors;

    // 遍历所有顶点及其邻接列表
    for (const auto& [source, neighbors] : adj) {
        if (neighbors.find(vertex) != neighbors.end()) {
            inNeighbors.push_back(source);
        }
    }

    return inNeighbors;
}



void QuerySharingGraph::Search(DirectedGraph& G, std::vector<Path>& P, std::unordered_map<Query, std::vector<Path>>& R, Path& path, Query& q, QuerySharingGraph& Psi) {
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
        // 5. 检查邻居 v'' 是否满足跳数约束
        if (1 <= q.k - path.vertices.size()) {
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
void QuerySharingGraph::BatchEnum(DirectedGraph& G, std::vector<Query>& Q, double gamma) {
    std::vector<Path> Output;
    // 多源BFS 1-2
    std::unordered_set<int> S, T;
    for (const Query& q : Q) {
        S.insert(q.s);
        T.insert(q.t);
    }

    initializeBFSIndex(S, T);

    // 初始化缓存和聚类查询
    std::unordered_map<Query, std::vector<Path>> R; // 缓存类型待定,查询和对应的路径集合
    std::vector<QueryGroup> Cs = hierarchicalClustering(Q, gamma);

    // 构建查询共享图
    for (auto& C : Cs) {
        QuerySharingGraph Psi;
        QuerySharingGraph Psir;
        DirectedGraph Gr;
        Gr = G.getReverseGraph();
        Psi.DetectCommonQuery(G, C, Psi);
        Psir.DetectCommonQuery_r(Gr, C, Psir);

        // 查询
        std::unordered_set<Query> processedQueries;
        // 按照拓扑顺序
        // 正向图
        TopoLogical Topo(Psi);
        for (auto& q : Topo.getorder()) {
            processedQueries.insert(q);
            if (q.t == 0) { // HC-s 路径查询
                std::vector<Path> P;
                Path initialPath;
                initialPath.vertices.push_back(q.s);
                Search(G, P, R, initialPath, q, Psi);
                R.insert({q, P});
            }
            if (q.t != 0 && R.find({q.s, 0, ceil((q.k)/2)}) != R.end() && R.find({q.t, 0, floor((q.k)/2)}) != R.end()) { // HC-s-t 路径查询
                std::vector<Path> forwardPaths;
                std::vector<Path> backwardPaths;
                for (auto& p1 : R[{q.s, 0, ceil((q.k)/2)}]) {
                    forwardPaths = {p1};
                }
                for (auto& p2 : R[{q.t, 0, floor((q.k)/2)}]) {
                    backwardPaths = {p2};
                }
                std::vector<Path> combinedPaths;
                G.connectPaths(forwardPaths, backwardPaths, combinedPaths);
                for (const auto& combinedPath : combinedPaths) {
                    if (G.isSimplePath(combinedPath)) {
                        Output.push_back(combinedPath);
                    }
                }
            }

            int tmp = 0;
            for (const auto& q_prime : Psi.getInNeighbors(q)) {
                for (auto& q_prime_two : Psi.getNeighbors(q_prime)) {
                    if (processedQueries.find(q_prime_two) == processedQueries.end()) {
                        tmp = 1;
                        break;
                    }
                    if (tmp != 1) {
                        R.erase(q_prime);
                    }
                }
            }
        }

        // 反向图
        TopoLogical Topo_r(Psir);
        for (auto& q : Topo_r.getorder()) {
            processedQueries.insert(q);
            if (q.t == 0) { // HC-s 路径查询
                std::vector<Path> P;
                Path initialPath;
                initialPath.vertices.push_back(q.s);
                Search(Gr, P, R, initialPath, q, Psir);
                R.insert({q, P});
            }
            if (q.t != 0 && R.find({q.s, 0, ceil((q.k)/2)}) != R.end() && R.find({q.t, 0, floor((q.k)/2)}) != R.end()) { // HC-s-t 路径查询
                std::vector<Path> forwardPaths;
                std::vector<Path> backwardPaths;
                for (auto& p1 : R[{q.s, 0, ceil((q.k)/2)}]) {
                    forwardPaths = {p1};
                }
                for (auto& p2 : R[{q.t, 0, floor((q.k)/2)}]) {
                    backwardPaths = {p2};
                }
                std::vector<Path> combinedPaths;
                G.connectPaths(forwardPaths, backwardPaths, combinedPaths);
                for (const auto& combinedPath : combinedPaths) {
                    if (G.isSimplePath(combinedPath)) {
                        Output.push_back(combinedPath);
                    }
                }
            }

            int tmp = 0;
            for (const auto& q_prime : Psi.getInNeighbors(q)) {
                for (auto& q_prime_two : Psi.getNeighbors(q_prime)) {
                    if (processedQueries.find(q_prime_two) == processedQueries.end()) {
                        tmp = 1;
                        break;
                    }
                    if (tmp != 1) {
                        R.erase(q_prime);
                    }
                }
            }
        }
    }
}