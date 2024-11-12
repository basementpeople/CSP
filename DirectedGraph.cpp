#include <vector>
#include <unordered_set>
#include <algorithm>

#include "DirectedGraph.h"
// 构造函数
// DirectedGraph::DirectedGraph(std::string path) : Graph(path) {}
DirectedGraph::DirectedGraph(std::string path) {
    readFromFile(path);
}

DirectedGraph::DirectedGraph() : Graph() {}

// 从文件读取有向图
void DirectedGraph::readFromFile(std::string fileName) {
    m = 0;
    n = 0;
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Error: File not found." << std::endl;
        exit(1);
    }
    int sum = 0;
    std::string line;
    while (std::getline(file, line)) {
        sum++;
        int from, to;
        if (sscanf(line.c_str(), "%d %d", &from, &to) == 2) {
            if (0 <= from && 0 <= to && from != to) {
                addNode(from);
                addNode(to);
                addEdge(from, to);
            } else {
                continue;
            }
        } else {
        	std::cout << sum << std::endl;
            std::cerr << "Error: Incorrect line format" << std::endl;
        }
    }
    file.close();
}

// 添加边
bool DirectedGraph::addEdge(int from, int to) {
    if (from != to && adj[from].find(to) == adj[from].end()) {
        adj[from].insert(to);
        m++;
        return true;
    } else {
        return false;
    }
}

// 深度优先搜索遍历所有顶点
std::vector<int> DirectedGraph::dfsTraversal(int vertex, std::unordered_set<int>& visited) {
    visited.insert(vertex);
    std::vector<int> result;
    result.push_back(vertex);

    for (int neighbor : getNeighbors(vertex)) {
        if (visited.find(neighbor) == visited.end()) {
            std::vector<int> subResult = dfsTraversal(neighbor, visited);
            result.insert(result.end(), subResult.begin(), subResult.end());
        }
    }

    return result;
}

// 遍历所有顶点的入口函数
std::vector<int> DirectedGraph::traverseAllVertices() {
    std::unordered_set<int> visited;
    std::vector<int> allVertices;

    // 遍历所有顶点，确保所有连通分量都被访问到
    for (const auto& entry : adj) {
        int vertex = entry.first;
        if (visited.find(vertex) == visited.end()) {
            std::vector<int> subResult = dfsTraversal(vertex, visited);
            allVertices.insert(allVertices.end(), subResult.begin(), subResult.end());
        }
    }

    return allVertices;
}

// 获取反向图
DirectedGraph DirectedGraph::getReverseGraph() const {
    DirectedGraph reverseGraph;
    for (const auto& entry : adj) {
        int from = entry.first;
        for (const auto& to : entry.second) {
            reverseGraph.addNode(to);
            reverseGraph.addNode(from);
            reverseGraph.addEdge(to, from);
        }
    }
    reverseGraph.m = m;
    reverseGraph.n = n;
    // 计算反向图的度数和最小度数，有什么用呢
    // reverseGraph.computeDegrees();
    // reverseGraph.computeMinimumDegree();
    // reverseGraph.statistic();
    return reverseGraph;
}

// 多源BFS
void DirectedGraph::initializeBFSIndex(const std::unordered_set<int>& S, const std::unordered_set<int>& T) {

    // 从集合S中的每个顶点执行BFS
    for (const auto& s : S) {
        bfsFromSource(s, distG[s]);
    }
            for (const auto& [vertex, distances] : distG) {
            std::cout << "Vertex: " << vertex << " -> Distances: ";
            for (const auto& distance : distances) {
                std::cout << distance.first << ":" << distance.second << std::endl;
            }
            std::cout << std::endl;
            }

    // 从集合T中的每个顶点执行BFS
    DirectedGraph reverseGraph; 
	reverseGraph = getReverseGraph();
    for (const auto& t : T) {
        reverseGraph.bfsFromSource(t, distGr[t]);
    }
            for (const auto& [vertex, distances] : distGr) {
            std::cout << "Vertex: " << vertex << " -> Distances: ";
            for (const auto& distance : distances) {
                std::cout << distance.first << ":" << distance.second << std::endl;
            }
            std::cout << std::endl;
            }
}

// BFS实现
void DirectedGraph::bfsFromSource(int source, std::unordered_map<int, int>& distances) {
    std::queue<int> q;
    q.push(source);
    distances[source] = 0;

    while (!q.empty()) {
        int v = q.front();
        q.pop();

        for (const auto& neighbor : adj[v]) {
            if (distances.find(neighbor) == distances.end()) {
                distances[neighbor] = distances[v] + 1;
                q.push(neighbor);
            }
        }
    }

    // 将未被访问到的点的 distance 设为无穷大
    for (const auto& vertex : adj) {
        if (distances.find(vertex.first) == distances.end()) {
            distances[vertex.first] = INT_MAX;
        }
    }
}

int DirectedGraph::dist (int s1, int s2) {
    if (s1 == s2) { return 0; }

    std::queue<int> q;
    q.push(s1);
    int distance = 0;
    std::unordered_map<int, int> distances;
    distances[s1] = 0;

    while (!q.empty()) {
        int v = q.front();
        q.pop();

        for (const auto& neighbor : adj[v]) {
            if (neighbor == s2) {
                return distances[v] + 1;
            }
            if (distances.find(neighbor) == distances.end()) {
                distances[neighbor] = distances[v] + 1;
                q.push(neighbor);
            }
        }
    }

    // 将未被访问到的点的 distance 设为无穷大
    return INT_MAX;
}

// 路径是否存在
bool DirectedGraph::isPathPossible(int s, int t, int maxDistance) {
    if (distG.find(s) == distG.end() || distGr.find(t) == distGr.end()) {
        return false;
    }

    for (const auto& [v, distFromS] : distG[s]) {
        if (distGr[t].find(v) != distGr[t].end()) {
            int totalDistance = distFromS + distGr[t][v];
            if (totalDistance <= maxDistance) {
                return true;
            }
        }
    }

    return false;
}

// 调用distG
const std::unordered_map<int, std::unordered_map<int, int>>& DirectedGraph::getDistG() const {
    return distG;
}

// 调用distGr
const std::unordered_map<int, std::unordered_map<int, int>>& DirectedGraph::getDistGr() const {
    return distGr;
}

// 并行多源BFS
// void DirectedGraph::parallelInitializeBFSIndex(const std::unordered_set<int>& S, const std::unordered_set<int>& T) {
//     std::vector<std::thread> threads;

//     // 从集合S中的每个顶点执行BFS
//     for (const auto& s : S) {
//         threads.emplace_back(&DirectedGraph::bfsFromSource, this, s, std::ref(distG[s]));
//     }

//     // 从集合T中的每个顶点执行BFS
//     DirectedGraph reverseGraph;
//     reverseGraph = getReverseGraph();
//     for (const auto& t : T) {
//         threads.emplace_back(&DirectedGraph::bfsFromSource, &reverseGraph, t, std::ref(distGr[t]));
//     }

//     for (auto& thread : threads) {
//         thread.join();
//     }
// }

// ？？？？？？？？？？？分割

// 检查路径是否为简单路径
bool DirectedGraph::isSimplePath(const Path& path) {
    std::unordered_set<int> visited;
    for (int vertex : path.vertices) {
        if (visited.find(vertex) != visited.end()) {
            return false;
        }
        visited.insert(vertex);
    }
    return true;
}

// 正向搜索2
void DirectedGraph::forwardSearch(const DirectedGraph& graph, int s, int t, int k, int K, std::vector<Path>& paths, Path& currentPath) {
    if (currentPath.vertices.size() - 1 <= k) {  // 路径长度小于等于 k
        paths.push_back(currentPath);
//        std::cout << "Forward Path from " << s << " to " << t << ": ";
//        for (int vertex : currentPath.vertices) {
//            std::cout << vertex << " ";
//        }
//        std::cout << std::endl;
    }
    if (currentPath.vertices.size() - 1 >= k) {  // 路径长度达到 k
        // std::cout << "离开1" << std::endl;
        return;
    }
    int lastVertex = currentPath.vertices.back();
    for (int neighbor : graph.getNeighbors(lastVertex)) {
    	// std::cout << "循环" << std::endl;
        if (std::find(currentPath.vertices.begin(), currentPath.vertices.end(), neighbor) == currentPath.vertices.end()) {
            // 剪枝：检查邻居节点到目标节点 t 的距离是否满足条件
            if (distGr[t].find(neighbor) != distGr[t].end() && distGr[t][neighbor] <= K - currentPath.vertices.size()) {
                currentPath.vertices.push_back(neighbor);
                forwardSearch(graph, s, t, k, K, paths, currentPath);
                currentPath.vertices.pop_back();
            }
        }
    }
}

// 反向搜索
void DirectedGraph::backwardSearch(const DirectedGraph& graph, int s, int t, int k, int K, std::vector<Path>& paths, Path& currentPath) {
    if (currentPath.vertices.size() - 1 <= k) {  // 路径长度小于等于 k
        paths.push_back(currentPath);
//        std::cout << "Backward Path from " << t << " to " << s << ": ";
//        for (int vertex : currentPath.vertices) {
//            std::cout << vertex << " ";
//        }
//        std::cout << std::endl;
    }
    if (currentPath.vertices.size() - 1 >= k) {  // 路径长度达到 k
        return;
    }
    int lastVertex = currentPath.vertices.back();
    for (int neighbor : graph.getReverseNeighbors(lastVertex)) {
        if (std::find(currentPath.vertices.begin(), currentPath.vertices.end(), neighbor) == currentPath.vertices.end()) {
            // 剪枝：检查邻居节点到起始节点 s 的距离是否满足条件
            if (distG[s].find(neighbor) != distG[s].end() && distG[s][neighbor] <= K - currentPath.vertices.size()) {
                currentPath.vertices.push_back(neighbor);
                backwardSearch(graph, s, t, k, K, paths, currentPath);
                currentPath.vertices.pop_back();
            }
        }
        // std::cout << "itr" << std::endl;
    }
}

// 连接路径 检查重复路径
// 连接路径
void DirectedGraph::connectPaths(const std::vector<Path>& forwardPaths, const std::vector<Path>& backwardPaths, std::vector<Path>& resultPaths) {
    std::unordered_set<std::string> uniquePaths;

    for (const Path& forwardPath : forwardPaths) {
        for (const Path& backwardPath : backwardPaths) {
            if (forwardPath.vertices.back() == backwardPath.vertices.back()) {
                Path connectedPath;
                connectedPath.vertices = forwardPath.vertices;
                connectedPath.vertices.pop_back();  // 去除重复的中间顶点
                connectedPath.vertices.insert(connectedPath.vertices.end(), backwardPath.vertices.rbegin(), backwardPath.vertices.rend());

                DirectedGraph graph;
                if (graph.isSimplePath(connectedPath)) {
                    // 将路径转换为字符串以便于比较
                    std::string pathStr;
                    for (int vertex : connectedPath.vertices) {
                        pathStr += std::to_string(vertex) + ",";
                    }

                    // 检查路径是否已经存在
                    if (uniquePaths.find(pathStr) == uniquePaths.end()) {
                        resultPaths.push_back(connectedPath);
                        uniquePaths.insert(pathStr);
                    }
                }
            }
        }
    }
}

// 主算法
void DirectedGraph::BasicEnum(const std::vector<Query>& queries) {
    std::unordered_set<int> S, T;
    for (const Query& q : queries) {
        S.insert(q.s);
        T.insert(q.t);
    }

    initializeBFSIndex(S, T);

    for (const Query& q : queries) {
        std::vector<Path> forwardPaths, backwardPaths, resultPaths;

        Path initialForwardPath, initialBackwardPath;
        initialForwardPath.vertices.push_back(q.s);
        initialBackwardPath.vertices.push_back(q.t);
        
        std::cout << "ceil:" << ceil(static_cast<double>(q.k)/2) << std::endl;
        std::cout << "floor:" << floor(static_cast<double>(q.k)/2) << std::endl;

        forwardSearch(*this, q.s, q.t, ceil(static_cast<double>(q.k)/2), q.k, forwardPaths, initialForwardPath);
        backwardSearch(*this, q.s, q.t, floor(static_cast<double>(q.k)/2), q.k, backwardPaths, initialBackwardPath);

        connectPaths(forwardPaths, backwardPaths, resultPaths);

        for (const Path& path : resultPaths) {
            std::cout << "Path from " << q.s << " to " << q.t << ": ";
            for (int vertex : path.vertices) {
                std::cout << vertex << " ";
            }
            std::cout << std::endl;
        }
        std::cout << q.s << std::endl;
    }
    std::cout << "结束" << std::endl;
}

// 获取邻居
std::vector<int> DirectedGraph::getNeighbors(int vertex) const {
    auto it = adj.find(vertex);
    if (it != adj.end()) {
        return std::vector<int>(it->second.begin(), it->second.end());
    }
    return {};
}

// 获取反向邻居
std::vector<int> DirectedGraph::getReverseNeighbors(int vertex) const {
    DirectedGraph reverseGraph;
	reverseGraph = getReverseGraph();
    return reverseGraph.getNeighbors(vertex);
}

// 计算 Hop-Constrained Neighbors H-C s邻居
HopConstrainedNeighbors DirectedGraph::getHopConstrainedNeighbors(int s, int t, int k) {
    std::unordered_map<int, int> forwardDistances, backwardDistances;
    bfsFromSource(s, forwardDistances);

    DirectedGraph reverseGraph;
    reverseGraph = getReverseGraph();
    reverseGraph.bfsFromSource(t, backwardDistances);
    // bfsFromSource(t, backwardDistances);

    HopConstrainedNeighbors hcn;
    for (const auto& entry : forwardDistances) {
        if (entry.second <= k) {
            hcn.forwardReachable.insert(entry.first);
        }
    }

    for (const auto& entry : backwardDistances) {
        if (entry.second <= k) {
            hcn.backwardReachable.insert(entry.first);
        }
    }

    return hcn;
}

// 计算两个查询之间的相似度
double DirectedGraph::querySimilarity(const Query& qA, const Query& qB) {
    HopConstrainedNeighbors hcnA = getHopConstrainedNeighbors(qA.s, qA.t, qA.k);
    HopConstrainedNeighbors hcnB = getHopConstrainedNeighbors(qB.s, qB.t, qB.k);
    std::cout<< qA.s << " " << qB.s << std::endl;
        // 打印 hcnA 和 hcnB 的内容
    std::cout << "hcnA forwardReachable: ";
    for (const auto& node : hcnA.forwardReachable) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    std::cout << "hcnA backwardReachable: ";
    for (const auto& node : hcnA.backwardReachable) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    std::cout << "hcnB forwardReachable: ";
    for (const auto& node : hcnB.forwardReachable) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    std::cout << "hcnB backwardReachable: ";
    for (const auto& node : hcnB.backwardReachable) {
        std::cout << node << " ";
    }
    std::cout << std::endl;


    std::unordered_set<int> intersectionForward, intersectionBackward;
    std::vector<int> sortedA(hcnA.forwardReachable.begin(), hcnA.forwardReachable.end());
    std::vector<int> sortedB(hcnB.forwardReachable.begin(), hcnB.forwardReachable.end());
    std::sort(sortedA.begin(), sortedA.end());
    std::sort(sortedB.begin(), sortedB.end());
    std::set_intersection(sortedA.begin(), sortedA.end(),
                          sortedB.begin(), sortedB.end(),
                          std::inserter(intersectionForward, intersectionForward.begin()));

    std::vector<int> sortedC(hcnA.backwardReachable.begin(), hcnA.backwardReachable.end());
    std::vector<int> sortedD(hcnB.backwardReachable.begin(), hcnB.backwardReachable.end());
    std::sort(sortedC.begin(), sortedC.end());
    std::sort(sortedD.begin(), sortedD.end());
    std::set_intersection(sortedC.begin(), sortedC.end(),
                          sortedD.begin(), sortedD.end(),
                          std::inserter(intersectionBackward, intersectionBackward.begin()));

    int sizeIntersectionForward = intersectionForward.size();
    int sizeIntersectionBackward = intersectionBackward.size();

    int minSizeForward = std::min(hcnA.forwardReachable.size(), hcnB.forwardReachable.size());
    int minSizeBackward = std::min(hcnA.backwardReachable.size(), hcnB.backwardReachable.size());
    
    std::cout<< minSizeForward << " "
    		<< "sizeIntersectionForward" << sizeIntersectionForward << " "
    		 << minSizeBackward << " "
    		 << sizeIntersectionBackward << std::endl;
    if (sizeIntersectionForward == 0 || sizeIntersectionBackward == 0) {
    	return 0;
	}

    double similarity = 2.0 / (minSizeForward / sizeIntersectionForward + minSizeBackward / sizeIntersectionBackward);
    std::cout << "1" << std::endl;
    return similarity;
}

// 计算两个查询组之间的相似度
double DirectedGraph::groupSimilarity(const QueryGroup& groupA, const QueryGroup& groupB) {
    double totalSimilarity = 0.0;
    for (const Query& qA : groupA.queries) {
        for (const Query& qB : groupB.queries) {
            totalSimilarity += querySimilarity(qA, qB);
        }
    }
    return totalSimilarity / (groupA.queries.size() * groupB.queries.size());
}

// 层次聚类算法
std::vector<QueryGroup> DirectedGraph::hierarchicalClustering(const std::vector<Query>& queries, double threshold) {
    std::vector<QueryGroup> groups;

    // 初始化每个查询为一个单独的组
    for (const Query& q : queries) {
        HopConstrainedNeighbors hcn = getHopConstrainedNeighbors(q.s, q.t, q.k);
        QueryGroup group;
        group.addQuery(q, hcn);
        groups.push_back(group);
    }

    while (groups.size() >= 1) {
        double maxSimilarity = -1.0;
        int bestPair[2] = {-1, -1};

        // 找到最相似的两个组
        for (size_t i = 0; i < groups.size(); ++i) {
            for (size_t j = i + 1; j < groups.size(); ++j) {
                double similarity = groupSimilarity(groups[i], groups[j]);
                if (similarity > maxSimilarity) {
                    maxSimilarity = similarity;
                    bestPair[0] = i;
                    bestPair[1] = j;
                    std::cout << maxSimilarity << std::endl;
                }
            }
        }
        if (maxSimilarity < threshold) {
            break;
        }
        // 合并最相似的两个组
        groups[bestPair[0]].queries.insert
            (groups[bestPair[0]].queries.end(), groups[bestPair[1]].queries.begin(), groups[bestPair[1]].queries.end());
        groups[bestPair[0]].hcn.forwardReachable.insert
            (groups[bestPair[1]].hcn.forwardReachable.begin(), groups[bestPair[1]].hcn.forwardReachable.end());
        groups[bestPair[0]].hcn.backwardReachable.insert
            (groups[bestPair[1]].hcn.backwardReachable.begin(), groups[bestPair[1]].hcn.backwardReachable.end());

        groups.erase(groups.begin() + bestPair[1]);
    }

    return groups;
}

// void DirectedGraph::Search(DirectedGraph& G, std::vector<Path>& P, std::unordered_map<Query, std::vector<Path>>& R, Path& path, Query& q, QuerySharingGraph& Psi) {
//     // 1. 获取当前路径的最后一个顶点 v'
//     int v_prime = path.vertices.back();
//     // 2. 将当前路径 p 添加到结果集 P 中
//     P.push_back(path);

//     // 3. 检查路径长度
//     if (path.vertices.size() - 1 == q.k) {
//         return;
//     }

//     // 4. 遍历当前顶点 v' 的所有邻居 v''
//     for (int v_double_prime : G.getNeighbors(v_prime)) {
//         // 5. 检查邻居 v'' 是否满足跳数约束
//         if (1 <= q.k - path.vertices.size()) {
//             // 6. 如果邻居 v'' 已经在当前路径 p 中，则跳过
//             if (std::find(path.vertices.begin(), path.vertices.end(), v_double_prime) != path.vertices.end()) {
//                 continue;
//             }

//             // 7. 如果存在一个查询 q'，其目标顶点 q'.v 等于 v''，则将 p 与 R[q'] 的路径合并并添加到 P 中
//             bool found = false;
//             for (auto& q_prime : Psi.getInNeighbors(q)) {
//                 if (q_prime.t == v_double_prime) {
//                     for (auto& r_path : R[q_prime]) {
//                         Path combinedPath = path;
//                         combinedPath.vertices.insert(combinedPath.vertices.end(), r_path.vertices.begin(), r_path.vertices.end());
//                         if (isSimplePath(combinedPath)) {
//                             P.push_back(combinedPath);
//                         }
//                     }
//                     found = true;
//                     break;
//                 }
//             }

//             // 8. 否则，递归调用 Search 函数继续搜索
//             if (!found) {
//                 path.vertices.push_back(v_double_prime);
//                 Search(G, P, R, path, q, Psi);
//                 path.vertices.pop_back();
//             }
//         }
//     }
// }

// // Batch查询算法
// void DirectedGraph::BatchEnum(DirectedGraph& G, std::vector<Query>& Q, double gamma) {
//     std::vector<Path> Output;
//     // 多源BFS 1-2
//     std::unordered_set<int> S, T;
//     for (const Query& q : Q) {
//         S.insert(q.s);
//         T.insert(q.t);
//     }

//     initializeBFSIndex(S, T);

//     // 初始化缓存和聚类查询
//     std::unordered_map<Query, std::vector<Path>> R; // 缓存类型待定,查询和对应的路径集合
//     std::vector<QueryGroup> Cs = hierarchicalClustering(Q, gamma);

//     // 构建查询共享图
//     for (auto& C : Cs) {
//         QuerySharingGraph Psi;
//         QuerySharingGraph Psir;
//         DirectedGraph Gr;
//         Gr = G.getReverseGraph();
//         Psi.DetectCommonQuery(G, C, Psi);
//         Psir.DetectCommonQuery_r(Gr, C, Psir);

//         // 查询
//         std::unordered_set<Query> processedQueries;
//         // 按照拓扑顺序
//         // 正向图
//         TopoLogical Topo(Psi);
//         for (auto& q : Topo.getorder()) {
//             processedQueries.insert(q);
//             if (q.t == 0) { // HC-s 路径查询
//                 std::vector<Path> P;
//                 Path initialPath;
//                 initialPath.vertices.push_back(q.s);
//                 Search(G, P, R, initialPath, q, Psi);
//                 R.insert({q, P});
//             }
//             if (q.t != 0 && R.find({q.s, 0, ceil((q.k)/2)}) != R.end() && R.find({q.t, 0, floor((q.k)/2)}) != R.end()) { // HC-s-t 路径查询
//                 std::vector<Path> forwardPaths;
//                 std::vector<Path> backwardPaths;
//                 for (auto& p1 : R[{q.s, 0, ceil((q.k)/2)}]) {
//                     forwardPaths = {p1};
//                 }
//                 for (auto& p2 : R[{q.t, 0, floor((q.k)/2)}]) {
//                     backwardPaths = {p2};
//                 }
//                 std::vector<Path> combinedPaths;
//                 connectPaths(forwardPaths, backwardPaths, combinedPaths);
//                 for (const auto& combinedPath : combinedPaths) {
//                     if (isSimplePath(combinedPath)) {
//                         Output.push_back(combinedPath);
//                     }
//                 }
//             }

//             int tmp = 0;
//             for (const auto& q_prime : Psi.getInNeighbors(q)) {
//                 for (auto& q_prime_two : Psi.getNeighbors(q_prime)) {
//                     if (processedQueries.find(q_prime_two) == processedQueries.end()) {
//                         tmp = 1;
//                         break;
//                     }
//                     if (tmp != 1) {
//                         R.erase(q_prime);
//                     }
//                 }
//             }
//         }

//         // 反向图
//         TopoLogical Topo_r(Psir);
//         for (auto& q : Topo_r.getorder()) {
//             processedQueries.insert(q);
//             if (q.t == 0) { // HC-s 路径查询
//                 std::vector<Path> P;
//                 Path initialPath;
//                 initialPath.vertices.push_back(q.s);
//                 Search(Gr, P, R, initialPath, q, Psir);
//                 R.insert({q, P});
//             }
//             if (q.t != 0 && R.find({q.s, 0, ceil((q.k)/2)}) != R.end() && R.find({q.t, 0, floor((q.k)/2)}) != R.end()) { // HC-s-t 路径查询
//                 std::vector<Path> forwardPaths;
//                 std::vector<Path> backwardPaths;
//                 for (auto& p1 : R[{q.s, 0, ceil((q.k)/2)}]) {
//                     forwardPaths = {p1};
//                 }
//                 for (auto& p2 : R[{q.t, 0, floor((q.k)/2)}]) {
//                     backwardPaths = {p2};
//                 }
//                 std::vector<Path> combinedPaths;
//                 connectPaths(forwardPaths, backwardPaths, combinedPaths);
//                 for (const auto& combinedPath : combinedPaths) {
//                     if (isSimplePath(combinedPath)) {
//                         Output.push_back(combinedPath);
//                     }
//                 }
//             }

//             int tmp = 0;
//             for (const auto& q_prime : Psi.getInNeighbors(q)) {
//                 for (auto& q_prime_two : Psi.getNeighbors(q_prime)) {
//                     if (processedQueries.find(q_prime_two) == processedQueries.end()) {
//                         tmp = 1;
//                         break;
//                     }
//                     if (tmp != 1) {
//                         R.erase(q_prime);
//                     }
//                 }
//             }
//         }
//     }
// }