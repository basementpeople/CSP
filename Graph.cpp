#include "Graph.h"

// 1. 构造函数和析构函数
Graph::Graph(const std::string path) {
       readFromFile(path);
       computeDegrees();
       statistic();
}

// 这个构造函数的意义
Graph::Graph() {
    m = 0;
    n = 0;
}

// 深拷贝构造函数，1有没有必要
Graph::Graph(const Graph &graph) {
    // 拷贝邻接表
    adj = graph.adj; // unordered_map 自动进行深拷贝

    // 拷贝节点的度
    degrees = graph.degrees; // unordered_map 的深拷贝

    // 拷贝有相同度的节点的集合的向量
    orderedNodes = graph.orderedNodes; // 深拷贝向量和集合

    // 拷贝最小度
    minimumDegree = graph.minimumDegree;

    // 拷贝最大度
    Dmax = graph.Dmax;

    // 拷贝图的节点数
    n = graph.n;

    // 拷贝图的边数
    m = graph.m;
}

Graph::~Graph() {}

// 社区节点集合转Graph，有没有必要
Graph Graph::getGraph(std::unordered_set<int> &subVertices) {
    // 创建一个新的 Graph 对象
    Graph subGraph;

    // 添加顶点
    for (const auto &vertex : subVertices)
    {
        if (adj.find(vertex) != adj.end())
        {
            subGraph.addNode(vertex);
        }
    }

    // 添加边
    for (const auto &vertex : subVertices)
    {
        if (adj.find(vertex) != adj.end())
        {
            std::unordered_set<int> neighbors = adj[vertex];
            for (const auto &neighbor : neighbors)
            {
                if (subVertices.find(neighbor) != subVertices.end())
                {
                    subGraph.addEdge(vertex, neighbor);
                }
            }
        }
    }
    subGraph.m /= 2;
    // 计算新的 Graph 的度数、最小度数等统计信息
    subGraph.computeDegrees();
    // subGraph.computeMinimumDegree();
    subGraph.statistic();

    return subGraph;
}

// 2 Greedy算法
std::vector<int> Graph::isQuerySetConnected(query_nodes queryNodes, std::unordered_map<int, int> degree) {
    if (queryNodes.empty()) return {}; // 如果查询节点为空，直接返回空

    std::unordered_set<int> visited; // 存储已访问的节点
    auto startIt = queryNodes.begin(); // 从查询节点中的任意一个开始
    int count = 0; // 用于记录访问到的查询节点数
    std::vector<int> tree;

    // 深度优先搜索
    dfs(*startIt, visited, queryNodes, count, degree, tree);

    // 如果访问的查询节点数等于查询节点集的大小，则说明连通
    if (count == queryNodes.size()) {
        return tree; // 返回连通子树
    } else {
        return {}; // 查询节点不连通，返回空
    }
}

void Graph::dfs(int current, std::unordered_set<int>& visited, const query_nodes& queryNodes, 
                int& count, const std::unordered_map<int, int>& degree, std::vector<int>& tree) {
    visited.insert(current); // 标记当前节点已访问
    tree.push_back(current); // 将当前节点加入连通子树

    // 如果当前节点是查询节点，则增加计数
    if (queryNodes.count(current)) {
        count++;
    }

    // 遍历当前节点的所有邻居
    for (const auto& neighbor : getNeighbors(current)) {
        // 如果邻居在查询节点集中且未被访问过，且其度数大于0，则递归访问
        if (!visited.count(neighbor) && degree.at(neighbor) > 0) {
            dfs(neighbor, visited, queryNodes, count, degree, tree);
        }
    }
}



void Graph::removeNode(int min_degree, std::unordered_map<int, std::unordered_set<int>>& resultGraph, std::vector<std::unordered_set<int>>& list) {
    if (list[min_degree].empty()) return; // 确保最低度数列表不为空

    auto it = list[min_degree].begin();
    int nodeToRemove = *it;

    // 从各个结构中移除节点
    list[min_degree].erase(nodeToRemove);
    resultGraph.erase(nodeToRemove);

    // 更新邻居节点的度数
    for (const auto& neighbor : getNeighbors(nodeToRemove)) {
        if (resultGraph.find(neighbor) != resultGraph.end()) {
            resultGraph[neighbor].erase(nodeToRemove);

            // 更新邻居节点的度数
            int neighborDegree = resultGraph[neighbor].size();
            list[resultGraph[neighbor].size() + 1].erase(neighbor); // 移除旧度数
            
            // 如果度数不为0，更新新度数
            if (neighborDegree > 0) {
                list[neighborDegree].insert(neighbor);
            }
        }
    }
}

// Greedy算法，返回图和最大最小度
std::unordered_set<int> Graph::greedy(query_nodes& queryNodes)
{
   std::unordered_map<int, std::unordered_set<int>> resultGraph = adj;
   std::vector<std::unordered_set<int>> list = orderedNodes;
   bool foundQueryNode = false;
   std::unordered_set<int> result_is;
   // 迭代移除度数最低的节点
   int now_degree = INT8_MAX;
    // 找到当前最低度数的节点
    for (int i = 0; i < list.size(); i++) {
        if (i < now_degree && !list[i].empty()) {
            now_degree = i;
            // std::cout << "初始最低度数：" << now_degree << std::endl;
        }
    }
   while (!list.empty()) {
       int min_degree = INT8_MAX;

       // 找到当前最低度数的节点
       if (!list[now_degree].empty()) {
        min_degree = now_degree;
       } else {
        if (!list[now_degree-1].empty()) {
            min_degree = --now_degree;
        } else {
            now_degree += 1;
            continue;
        }
       }
       // 检查是否包含查询节点
        for (auto& pair : list[min_degree]) {
            if (queryNodes.count(pair)) {
                foundQueryNode = true;
            }
        }

       if (foundQueryNode) {
           break;
       }

       // 移除节点
       removeNode(min_degree, resultGraph, list);

        // 检查查询顶点集是否连通
        // if (!isQuerySetConnected(queryNodes)) {
        // 	std::cout << "已断开" << std::endl;
        //         std::unordered_set<int> result;
        //         for (const auto& pair : resultGraph) {
        //         result.insert(pair.first);  // 获取键
        //         }
        //         if (computesubMinimumDegree(result) > computesubMinimumDegree(result_is)){
        //             result_is = result;
        //         }
        //     // 如果不连通，则恢复节点并结束循环
        //     // list[min_degree].insert(nodeToRemove); // 假设 removeNode 中有 nodeToRemove 的记录
        //     // resultGraph[nodeToRemove] = getNeighbors(nodeToRemove); // 恢复节点及其邻居关系
        //     break;
        // }
   }
        std::unordered_set<int> result;

        for (const auto& pair : resultGraph) {
            if (pair.second.size() != 0) {
                result.insert(pair.first);  // 获取键
            }
        }
        result_is = result;

   return result_is;
}

std::unordered_set<int> Graph::greedy_d1(query_nodes& queryNodes) {
    std::cout << "开始执行 Greedy 算法" << std::endl;
        std::vector<std::unordered_set<int>> list = getOrderedNodes();
        std::unordered_map<int, int> degree = getDegrees();

        // 检查查询顶点集是否连通
        std::vector<int> tree;
        tree = isQuerySetConnected(queryNodes, degree);
        if (tree.size() == 0) {
            std::cout << "查询顶点集不连通！" << std::endl;
            return {};
        }

        int node;             // 当前处理的节点
        int lowestDegree = 0; // 当前的最低度数
        int neighborDegree;   // 邻居的度数
        int flag = 0; // 标记是否包含查询节点

        std::unordered_set<int> result_end;
        int k = 0;

        // 主循环
        while (!list.empty())
        {
            if (list[lowestDegree].empty())
            {
                ++lowestDegree; // 如果没有该度数的节点，增加度数
            }
            else
            {
                // 检查是否包含查询节点
                for (auto& pair : list[lowestDegree]) {
                    if (queryNodes.count(pair)) {
                        std::cout << "1:Found query node!" << pair << std::endl;
                        flag = 1;
                        break;
                    }
                }
                if (flag) break; // 包含查询节点，结束搜索
                node = *list[lowestDegree].begin();
                list[lowestDegree].erase(node);
                int debug_1 = degree[node];
                degree[node] = -1;         // 将节点度数设为-1，标记已处理

                // 检查查询顶点集是否连通
                // if (std::find(tree.begin(), tree.end(), node) != tree.end()) {
                //     std::cout << "today" << node << std::endl;
                //     tree = isQuerySetConnected(queryNodes, degree);
                //     if (tree.size() == 0) {
                //         std::cout << "3:Query vertex set not connected！" << std::endl;
                //         degree[node] = 1; // 恢复节点度数,暂恢复为1
                //         // flag = 1;
                //         // degree[node] = debug_1; // 恢复节点度数
                //         // list[lowestDegree].insert(node);
                //         // continue;
                //         // break;
                //         // 不跳出，换一个点删除
                //     }
                // }

                // 检查查询顶点集是否连通
                flag = isConnected(queryNodes, degree);
                if (tree.size() == 0) {
                    std::cout << "3:Query vertex set not connected！" << std::endl;
                    degree[node] = 1; // 恢复节点度数,暂恢复为1
                    // flag = 1;
                    // degree[node] = debug_1; // 恢复节点度数
                    // list[lowestDegree].insert(node);
                    // continue;
                    // break;
                    // 不跳出，换一个点删除
                }

                // 更新所有邻居的度数
                for (int neighbor : getNeighbors(node))
                {
                    neighborDegree = degree[neighbor];
                    if (neighborDegree > lowestDegree)
                    {
                        list[neighborDegree].erase(neighbor);
                        list[neighborDegree - 1].insert(neighbor);
                        degree[neighbor] = neighborDegree - 1;

                        // if (queryNodes.count(neighbor)) {
                        //     std::cout << "2:Found query node!" << neighbor << std::endl;
                        //     degree[node] = 1; // 恢复节点度数,暂恢复为1
                        //     flag = 1;
                        //     break;
                        // }
                    }
                }

                // if (flag) break; // 包含查询节点，结束搜索
            }

            std::unordered_set<int> result;
            for (const auto& pair : degree) {
                if (pair.second != -1 && pair.second != 0) {
                    result.insert(pair.first);  // 获取键
                }
            }

            int minDegree = INT8_MAX;
            for (auto &node : result)
            {
                int degree = 0;
                for (auto &neighbor : getNeighbors(node))
                {
                    if (result.find(neighbor) != result.end())
                        degree++;
                }
                if (degree < minDegree)
                    minDegree = degree;
            }
            if (minDegree > k) {
                k = minDegree;
                result_end = result;
            }

        }

        // std::unordered_set<int> result;
        // for (const auto& pair : degree) {
        //     if (pair.second != -1 && pair.second != 0) {
        //         result.insert(pair.first);  // 获取键
        //     }
        // }
        // return result;
        std::vector<int> tree2;
        tree2 = isQuerySetConnected(queryNodes, degree);
        std::cout << "tree222222222222 " << tree.size() << std::endl;
        if (tree.size() == 0) {
            std::cout << "查询顶点集不连通！" << std::endl;
            return {};
        }

        result_end = getLastResult(queryNodes, result_end);

        std::cout << "得到结果为：" << k << std::endl;
        return result_end;
}

std::unordered_set<int> Graph::greedy_d2(query_nodes& queryNodes) {
    std::vector<std::unordered_set<int>> list = getOrderedNodes();
    std::unordered_map<int, int> degree = getDegrees();

    int node;             // 当前处理的节点
    int lowestDegree = 0; // 当前的最低度数
    int neighborDegree;   // 邻居的度数

    // 主循环
    while (!list.empty())
    {
        if (lowestDegree > 0)
        {
            if (!list[lowestDegree-1].empty())
            {
                --lowestDegree; // 如果没有该度数的节点，减少度数
            }
        }
        if (list[lowestDegree].empty())
        {
            ++lowestDegree; // 如果没有该度数的节点，增加度数
        }
        else
        {
            std::cout << "000000000000000!" << std::endl;
            int flag = 0;
            // 检查是否包含查询节点
            for (auto& pair : list[lowestDegree]) {
                if (queryNodes.count(pair)) {
                    flag = 1;
                    break;
                    std::cout << "Found query node!" << pair << std::endl;
                }
            }
            if (flag) break; // 包含查询节点，结束搜索
            node = *list[lowestDegree].begin();
            list[lowestDegree].erase(node);
            degree[node] = -1;         // 将节点度数设为-1，标记已处理

            // 更新所有邻居的度数
            for (int neighbor : getNeighbors(node))
            {
                neighborDegree = degree[neighbor];
                
                list[neighborDegree].erase(neighbor);
                if (neighborDegree > 1) {
                    list[neighborDegree - 1].insert(neighbor);
                }
                degree[neighbor] = neighborDegree - 1;
                
            }
        }
    }
    std::unordered_set<int> result;
    for (const auto& pair : degree) {
        if (pair.second != -1) {
            result.insert(pair.first);  // 获取键
        }
    }
    return result;
}

std::unordered_set<int> Graph::greedy_d3(query_nodes& queryNodes) {
    std::vector<std::unordered_set<int>> list = getOrderedNodes();
    std::unordered_map<int, int> degree = getDegrees();

    // 检查查询顶点集是否连通
    bool flag = false;
    flag = isConnected(queryNodes, degree);
    if (!flag) { return {}; }

    int node;             // 当前处理的节点
    int lowestDegree = 0; // 当前的最低度数
    int neighborDegree;   // 邻居的度数
    std::unordered_set<int> result_end; // 结果集
    int k = 0;

    // 主循环
    while (!list.empty()) {
        if (list[lowestDegree].empty()) {
            ++lowestDegree; // 如果没有该度数的节点，增加度数
        }
        else {
            // 检查是否包含查询节点
            for (auto& pair : list[lowestDegree]) {
                if (queryNodes.count(pair)) {
                    std::cout << "1:Found query node!" << pair << std::endl;
                    flag = false;
                    break;
                }
            }
            if (!flag) { break; } // 包含查询节点，结束搜索
            node = *list[lowestDegree].begin();
            list[lowestDegree].erase(node);
            int debug_1 = degree[node];
            degree[node] = -1;         // 将节点度数设为-1，标记已处理

            // 检查查询顶点集是否连通
            flag = isConnected(queryNodes, degree);
            if (!flag) {
                std::cout << "3:Query vertex set not connected！" << std::endl;
                break;
            }

            // 更新所有邻居的度数
            for (int neighbor : getNeighbors(node)) {
                neighborDegree = degree[neighbor];
                if (neighborDegree > lowestDegree) {
                    list[neighborDegree].erase(neighbor);
                    list[neighborDegree - 1].insert(neighbor);
                    degree[neighbor] = neighborDegree - 1;
                }
            }
        }

        std::unordered_set<int> result;
        for (const auto& pair : degree) {
            if (pair.second != -1 && pair.second != 0) {
                result.insert(pair.first);  // 获取键
            }
        }

        int minDegree = INT8_MAX;
        for (auto &node : result) {
            int degree = 0;
            for (auto &neighbor : getNeighbors(node)) {
                if (result.find(neighbor) != result.end()) { degree++; }
            }
            if (degree < minDegree) { minDegree = degree; }
        }
        if (minDegree > k) {
            k = minDegree;
            result_end = result;
        }
    }

    flag = isConnected(queryNodes, degree);
    if (!flag) {
        std::cout << "查询顶点集不连通！" << std::endl;
        return {};
    }

    result_end = getLastResult(queryNodes, result_end);

    std::cout << "得到结果为：" << k << std::endl;
    return result_end;
}

std::unordered_set<int> Graph::greedy_end(query_nodes& queryNodes) {
    std::vector<std::unordered_set<int>> list = getOrderedNodes();
    std::unordered_map<int, int> degree = getDegrees();

    // 检查查询顶点集是否连通
    bool flag = false;
    flag = isConnected(queryNodes, degree);
    if (!flag) { return {}; }

    int node;             // 当前处理的节点
    int lowestDegree = 0; // 当前的最低度数
    int neighborDegree;   // 邻居的度数
    std::unordered_set<int> result_end; // 结果集
    int k = 0;

    // 主循环
    while (!list.empty()) {
        if (list[lowestDegree].empty()) {
            ++lowestDegree; // 如果没有该度数的节点，增加度数
        }
        else {
            // 检查是否包含查询节点
            for (auto& pair : list[lowestDegree]) {
                if (queryNodes.count(pair)) {
                    std::cout << "1:Found query node!" << pair << std::endl;
                    flag = false;
                    break;
                }
            }
            if (!flag) { break; } // 包含查询节点，结束搜索
            node = *list[lowestDegree].begin();
            list[lowestDegree].erase(node);
            int debug_1 = degree[node];
            degree[node] = -1;         // 将节点度数设为-1，标记已处理

            // 检查查询顶点集是否连通
            flag = isConnected(queryNodes, degree);
            if (!flag) {
                std::cout << "3:Query vertex set not connected！" << std::endl;
                break;
            }

            // 更新所有邻居的度数
            for (int neighbor : getNeighbors(node)) {
                neighborDegree = degree[neighbor];
                if (neighborDegree > lowestDegree) {
                    list[neighborDegree].erase(neighbor);
                    list[neighborDegree - 1].insert(neighbor);
                    degree[neighbor] = neighborDegree - 1;
                }
                if (degree[neighbor] < lowestDegree) { lowestDegree = degree[neighbor]; }
            }
        }

        std::unordered_set<int> result;
        for (const auto& pair : degree) {
            if (pair.second != -1 && pair.second != 0) {
                result.insert(pair.first);  // 获取键
            }
        }

        int minDegree = INT8_MAX;
        for (auto &node : result) {
            int degree = 0;
            for (auto &neighbor : getNeighbors(node)) {
                if (result.find(neighbor) != result.end()) { degree++; }
            }
            if (degree < minDegree) { minDegree = degree; }
        }
        if (minDegree > k) {
            k = minDegree;
            result_end = result;
        }
    }

    flag = isConnected(queryNodes, degree);
    if (!flag) {
        std::cout << "查询顶点集不连通！" << std::endl;
        return {};
    }

    result_end = getLastResult(queryNodes, result_end);

    std::cout << "得到结果为：" << k << std::endl;
    return result_end;
}


// 针对单个节点的 Greedy 算法
std::unordered_set<int> Graph::greedy(int v0)
{
   std::unordered_map<int, std::unordered_set<int>> resultGraph = adj;
   std::vector<std::unordered_set<int>> list = orderedNodes;
   bool foundQueryNode = false;
   std::unordered_set<int> result_is;
   // 迭代移除度数最低的节点
   while (!list.empty()) {
       int min_degree = INT8_MAX;

       // 找到当前最低度数的节点
       for (int i = 0; i < list.size(); i++) {
           if (i < min_degree && !list[i].empty()) {
               min_degree = i;
           }
       }

       // 检查是否包含查询节点
        for (auto& pair : list[min_degree]) {
            if (v0 == pair) {
                foundQueryNode = true;
                // cout << "Found query node!" << pair << endl;
            }
        }

       if (foundQueryNode) break;

       // 移除节点
       removeNode(min_degree, resultGraph, list);
   }
    std::unordered_set<int> result;
    for (const auto& pair : resultGraph) {
        result.insert(pair.first);  // 获取键
        std::cout << "节点：" << pair.first << "度数：" << pair.second.size() << std::endl;
    }
   return result;
}


// 3 Local Search 算法

// 为什么查询顶点需要唯一
void Graph::search(std::unordered_set<int> H0, int k, std::unordered_set<int>& H)
{
    if(computesubMinimumDegree(H0) == k)
    {
        H = H0;
        std::cout << "得到结果为：" << H.size() << std::endl;
        return ;
    }
    
    for (auto &node : H0)
    {
        for (auto &neighbor : getNeighbors(node))
        {
            if (H0.find(neighbor) == H0.end())
            {
                std::unordered_set<int> H1 = H0;
                H1.insert(neighbor);
                if (computesubMinimumDegree(H1) >= computesubMinimumDegree(H0))
                {
                    search(H1, k, H);
                    std::cout << "结束一次递归：" << neighbor << std::endl;
                    if (H.size() != 0) return ;
                    std::cout << "没有得到结果：" << std::endl;
                }
            }
        }
    }
}

std::unordered_set<int> Graph::baseline_search(int v0, int k)
{
	std::cout << "结果：" << std::endl;
    std::vector<int> visited = std::vector<int>(n+1, 0);
    std::queue<int> queue;
    queue.push(v0);
    std::unordered_set<int> C;

    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();
        C.insert(v);
        std::cout << "C加入: " << v << std::endl;
        visited[v] = 1;

        if (computesubMinimumDegree(C) >= k) {
        	std::cout << "结果：" << computesubMinimumDegree(C) << std::endl;
            return C;  // 找到有效解
        }

        for (auto &w : getNeighbors(v)) {
            std::unordered_set<int> H = C;
            H.insert(w);
            if (C.find(w) == C.end() && degrees[w] >= k && !visited[w] && computesubMinimumDegree(H) >= computesubMinimumDegree(C)) {
                queue.push(w); visited[w] = 1;
            }
        }
    }
    std::cout << "结果：" << computesubMinimumDegree(C) << std::endl;
    std::cout << "k：" << k << std::endl;
    return C;
}

std::unordered_set<int> Graph::baseline_search2(int v0, int k)
{
	std::cout << "结果：" << std::endl;
    std::vector<int> visited = std::vector<int>(n+1, 0);
    std::queue<int> queue;
    queue.push(v0);
    std::unordered_set<int> C;

    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();
        C.insert(v);
        std::cout << "C加入: " << v << std::endl;
        visited[v] = 1;

        if (computesubMinimumDegree(C) >= k) {
        	std::cout << "结果：" << computesubMinimumDegree(C) << std::endl;
            return C;  // 找到有效解
        }

        for (auto &v : C){
            int f = 0;
            for (auto &w : getNeighbors(v)) {
            std::unordered_set<int> H = C;
            H.insert(w);
            if (C.find(w) == C.end() && degrees[w] >= k && !visited[w] && computesubMinimumDegree(H) >= computesubMinimumDegree(C)) {
                queue.push(w); visited[w] = 1;
                f = 1;
            }
            if (f == 1) { break; }
        }
        if (f == 1) { break; }
        }
    }
    std::cout << "结果：" << computesubMinimumDegree(C) << std::endl;
    std::cout << "k：" << k << std::endl;
    return C;
}

std::unordered_set<int> Graph::naiveCandidateGeneration(int v0, int k)
{
    std::vector<int> visited = std::vector<int>(n+1, 0);
    std::queue<int> queue;
    queue.push(v0);
    std::unordered_set<int> C;

    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();
        C.insert(v);
        visited[v] = 1;

        if (computesubMinimumDegree(C) >= k) {
        	std::cout << "结果：" << computesubMinimumDegree(C) << std::endl;
            return C;  // 找到有效解
        }

        for (auto &w : getNeighbors(v)) {
            if (C.find(w) == C.end() && degrees[w] >= k && !visited[w]) {
                queue.push(w); visited[w] = 1;
            }
        }
    }
    std::cout << "结果：" << computesubMinimumDegree(C) << std::endl;
    std::cout << "k：" << k << std::endl;
    return C;
}

// CST框架，确定图G中存在解 , V点E边 m为边数，n为节点数，需初始化
bool Graph::upperBound(int k)
{
   int bond = std::ceil((1 + sqrt(9 + 8 * (m - n))) / 2);
   if (bond < k) return false;
   else return true;
}

// CST框架
std::unordered_set<int> Graph::CSTframework(int v0, int k)
{
    if(upperBound(k)){
        std::unordered_set<int> C = std::unordered_set<int>();
        return naiveCandidateGeneration(v0, k);
    }
    else {
    	return std::unordered_set<int>();
	}
}

std::unordered_set<int> Graph::CSMframework(int v0, double gamma) {
    std::unordered_set<int> H, A, B;
    A.insert(v0);
    for (auto &neighbor : getNeighbors(v0)) {
            B.insert(neighbor);
    }
    int s = 0;

    // s的限制条件未检查 e−γ (b |E|−|V |  (δ(G[H])+1)/2−1 c − |H| 
    // exp(-gamma) * (n - ((n - m) / (double)(H.size() + 1))
    // exp(-gamma) * ( (n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1) - H.size())
    // !B.empty() && s <= exp(-gamma) * ( (n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1) - H.size())
    while (!B.empty()) { 
    	// std::cout << " succeed in!!" << exp(-gamma) * ( std::ceil((n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1)) - H.size()) << std::endl;
        int v = 0;
        int mintmp = 0;
        for (auto &node : B) {
            int tmp = 0;
            for (auto &neighbor : getNeighbors(node)) {
                if (A.find(neighbor) != A.end())
                tmp++;
            }
            if (tmp > mintmp) {
                mintmp = tmp;
                v = node;
            }
        }

        A.insert(v);
        B.erase(v);
        s++;
        if (computesubMinimumDegree(A) > computesubMinimumDegree(H)) {
            H = A;
            s = 0;
            if (computesubMinimumDegree(H) == std::min(degrees[v0], static_cast<int>(floor((1 + std::sqrt(9 + 8 * (m - n))/2))))) {
                return H;
            }
        }

        for (auto &neighbor : getNeighbors(v)){
            if (B.find(neighbor) == B.end() && A.find(neighbor) == A.end() && degrees[neighbor] > computesubMinimumDegree(H)) {
                B.insert(neighbor);
                std::cout << "neighbor " << neighbor << std::endl;
            }
            if (degrees[neighbor] <= computesubMinimumDegree(H)) {break;}
        }
        
        std::cout << "s   " << s << std::endl;
        std::cout << "H   " << computesubMinimumDegree(H) << std::endl;
        if (computesubMinimumDegree(H) == 1 || computesubMinimumDegree(H) == 2) {
            continue;
         }
        if (s > exp(-gamma) * ( floor ((m - n) / ((computesubMinimumDegree(H) + 1)/2 - 1)) - H.size())) {
        	std::cout << "公式: " << exp(-gamma) * ( floor ((m - n) / ((computesubMinimumDegree(H) + 1)/2 - 1)) - H.size()) << std::endl;
            break;
        }
        std::cout << "公式: " << exp(-gamma) * ( floor ((m - n) / ((computesubMinimumDegree(H) + 1)/2 - 1)) - H.size()) << std::endl;
        }
        int k = computesubMinimumDegree(H);
        std::cout << "k   " << k << std::endl;

        //第一步，得到了查询顶点的搜索空间，相当于v0，图G
        //有效的输出：H,A
        return maxcore(A, v0);
    // return H;
}


// generateCandidates C = A

// maxcore
std::unordered_set<int> Graph::maxcore(std::unordered_set<int>& C, int v0) {
    Graph new_G;
	new_G = getGraph(C);
    // std::unordered_set<int> queryNodes={v0};
    int queryNodes = v0;
    std::unordered_set<int> solution = {};
    solution = new_G.greedy(queryNodes);
    std::cout <<new_G.computesubMinimumDegree(solution) << " ";
    return solution;
}

// 计算两个查询之间的相似度
double Graph::querySimilarity(const query_nodes& qA, const query_nodes& qB) {
    // qA的邻居
    std::unordered_set<int> neighborsA;
    for (auto &v :qA) {
        neighborsA.insert(getNeighbors(v).begin(), getNeighbors(v).end());
    }
    // qB的邻居
    std::unordered_set<int> neighborsB;
    for (auto &v :qB) {
        neighborsB.insert(getNeighbors(v).begin(), getNeighbors(v).end());
    }
    // 计算两个邻居的交集
    std::unordered_set<int> intersection;
    std::set_intersection(neighborsA.begin(), neighborsA.end(), neighborsB.begin(), neighborsB.end(), std::inserter(intersection, intersection.begin()));
    
    // 计算相似度
    return intersection.size() / (double)(std::min(qA.size(), qB.size()));
}

// 计算两个查询顶点集之间的相似度, 1 - 邻居
double Graph::groupSimilarity(const query_group& groupA, const query_group& groupB) {
    double totalSimilarity = 0.0;
    for (const query_nodes& qA : groupA) {
        for (const query_nodes& qB : groupB) {
            totalSimilarity += querySimilarity(qA, qB);
        }
    }
    return totalSimilarity / (groupA.size() * groupB.size());
}

// 聚类算法
std::vector<query_group> Graph::Clustering(std::vector<query_nodes>& query_groups, int k, double threshold) {
    std::vector<query_group> groups;

    // 初始化每个查询为一个单独的组
    for (const query_nodes& q : query_groups) {
        // HopConstrainedNeighbors hcn = getHopConstrainedNeighbors(q.s, q.t, q.k);
        query_group group;
        group.push_back(q);
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
                    // std::cout << maxSimilarity << std::endl;
                }
            }
        }
        if (maxSimilarity < threshold) {
            break;
        }
        // 合并最相似的两个组
        groups[bestPair[0]].insert(groups[bestPair[0]].end(), 
                           std::make_move_iterator(groups[bestPair[1]].begin()), 
                           std::make_move_iterator(groups[bestPair[1]].end()));
        groups.erase(groups.begin() + bestPair[1]);
    }

    return groups;

}


// 辅助函数，分界线----------------------------------------------------------------------------------
void Graph::readFromFile(const std::string fileName) {
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        std::cerr << "Error: File not found." << std::endl;
        exit(1);
    }
    
    m = 0;
	n = 0;
    std::string line;
    // 测试数据集的边是否适用无向图
    while (std::getline(file, line)) {
        int from, to;
        if (sscanf(line.c_str(), "%d %d", &from, &to) == 2) {
        	// std::cout << from << "  " << to << std::endl;
            if (from == to)
            {
                continue;
            }
            
            addNode(from);
            addNode(to);

            addEdge(from, to);
        }
        else {
            std::cerr << "Error: Incorrect line format" << std::endl;
        }
    }
    file.close();
}

// 如果adj没有这个点，加入，更新n
void Graph::addNode(int node) {
    if (adj.find(node) == adj.end()) {
        adj[node] = std::unordered_set<int>();
        n++;
    }
}

// 如果adj没有这条边，加入，更新m
void Graph::addEdge(int from, int to) {
    if (adj[from].find(to) == adj[from].end()) {
        adj[from].insert(to);
    }
    if (adj[to].find(from) == adj[to].end()) {
        adj[to].insert(from);
    }
    m++;
}

// 计算degrees，orderedNodes，minimumDegree，Dmax
std::unordered_map<int, int> Graph::computeDegrees() {
    degrees.clear();
    orderedNodes.clear();
    orderedNodes.resize(adj.size());

    Dmax = 0;
    minimumDegree = INT32_MAX;

    for (auto &entry : adj)
    {
        degrees[entry.first] = entry.second.size();
        orderedNodes[entry.second.size()].insert(entry.first);
        if (entry.second.size() > Dmax) {
            Dmax = entry.second.size();
        }
        if (entry.second.size() < minimumDegree) {
            minimumDegree = entry.second.size();
        }
    }

    return degrees;
}

// 更新最小度，有必要吗
// int Graph::computeMinimumDegree()
// {
//     for (size_t i = 0; i < orderedNodes.size(); ++i)
//     {
//         if (!orderedNodes[i].empty())
//         {
//             minimumDegree = i;
//             return i;
//         }
//     }
//     return minimumDegree;
// }

// 打印图的相关信息
void Graph::statistic() {
    std::cout << "初始图的信息：" << std::endl;
    // 节点数
    std::cout << "节点数为 " << n << std::endl;
    // 边数
    std::cout << "边数为 " << m << std::endl;
    // 最大度数
    std::cout << "图的最大度数为 " << Dmax << std::endl;
    // 最小度数
    std::cout << "图的最小度数为 " << minimumDegree << std::endl;
}

// 检查查询集是否连通
bool Graph::isConnected(query_nodes queryNodes, std::unordered_map<int, int> degree) {
    if (queryNodes.empty() || degree.empty()) return false; // 如果查询节点为空，直接返回空

    std::unordered_set<int> visited; // 存储已访问的节点
    std::queue<int> q;
    for (int startNode : queryNodes)
    {
        if (degree.find(startNode) != degree.end() && degree.at(startNode) > 0)
        {
            q.push(startNode);
            visited.insert(startNode);
            break;
        }
    }

    if (q.empty())
    {
        return false; // 如果没有有效的起始节点，则直接返回不连通
    }

    while (!q.empty())
    {
        // std::cout<<"BFS"<<std::endl;
        int node = q.front();
        q.pop();
        for (int neighbor : getNeighbors(node))
        {
            if (degree.find(neighbor) != degree.end() && visited.insert(neighbor).second && degree.at(neighbor) > 0)
            {
                q.push(neighbor);
            }
        }
    }

    // 检查所有查询节点是否都在访问集合中
    for (int node : queryNodes)
    {
        if (visited.find(node) == visited.end())
        {
            return false;
        }
    }

    return true;
}

// 获取并返回排序后的邻居节点
std::unordered_set<int> Graph::getNeighbors(int node) {
    std::vector<int> sortedNeighbors = sortNeighbors(node);
    return std::unordered_set<int>(sortedNeighbors.begin(), sortedNeighbors.end());
}

// 对特定节点的邻居节点进行排序
std::vector<int> Graph::sortNeighbors(int node) {
    // 获取节点的邻居节点
    std::unordered_set<int> &neighbors = adj[node];

    // 将邻居节点存储到向量中
    std::vector<int> sortedNeighbors(neighbors.begin(), neighbors.end());
    auto compareByDegree = [this](int a, int b) 
    {
        return degrees[a] > degrees[b]; // 排序
    };

    // 对向量中的元素进行排序
    std::sort(sortedNeighbors.begin(), sortedNeighbors.end(), compareByDegree);

    return sortedNeighbors;
}

// 计算子社区中的最小度，比较常用，但打算更新子社区的数据结构
int Graph::computesubMinimumDegree(std::unordered_set<int> &nodes) {
	if (nodes.size() == 0) return 0;
    int minDegree = INT8_MAX;
    for (auto &node : nodes)
    {
        int degree = 0;
        for (auto &neighbor : getNeighbors(node))
        {
            if (nodes.find(neighbor) != nodes.end())
                degree++;
        }
        if (degree < minDegree)
            minDegree = degree;
    }
    return minDegree;
}

// 确保结果联通
std::unordered_set<int> Graph::getLastResult(query_nodes& queryNodes, std::unordered_set<int>& result_end) {
    if (queryNodes.empty() || result_end.empty()) {
        return {}; // 如果输入集合为空，直接返回空结果
    }

    std::unordered_set<int> result;
    std::queue<int> q;

    // 初始化队列
    for (auto& node : queryNodes) {
        q.push(node);
        result.insert(node); // 避免起始节点重复访问
    }

    while (!q.empty()) {
        int node = q.front();
        q.pop();

        for (auto& neighbor : getNeighbors(node)) {
            if (result.find(neighbor) != result.end()) {
                continue; // 如果已经访问过，跳过
            }
            if (result_end.find(neighbor) != result_end.end()) {
                q.push(neighbor);
                result.insert(neighbor);
            }
        }
    }

    return result;
}


// 低价值的辅助函数，可删除
// 打印adj
void Graph::printAdj() {
    for (const auto& entry : adj)
    {
        int vertex = entry.first;
        const std::unordered_set<int>& neighbors = entry.second;

        std::cout << "Vertex " << vertex << ": ";
        for (int neighbor : neighbors)
        {
            std::cout << neighbor << " ";
        }
        std::cout << std::endl;
    }
}

// 获取图中节点的数量
unsigned int Graph::getNumberOfNodes() {
    return n;
}

// 在 Graph 类中添加一个新方法来找出度数最大的点
int Graph::findMaxDegreeNode() {
    if (degrees.empty()) {
        std::cerr << "Error: No nodes in the graph." << std::endl;
        return -1;
    }

    int maxDegreeNode = -1;
    int maxDegree = -1;

    for (const auto &entry : degrees) {
        if (entry.second > maxDegree) {
            maxDegree = entry.second;
            maxDegreeNode = entry.first;
        }
    }

    return maxDegreeNode;
}

// 获取图中所有节点
std::unordered_set<int> Graph::getNodes() {
    std::unordered_set<int> nodes;
    for (auto &entry : adj) {
        nodes.insert(entry.first);
    }
    return nodes;
}