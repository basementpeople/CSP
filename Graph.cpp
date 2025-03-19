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

Graph::Graph(std::unordered_set<int> &subVertices, Graph &graph) {
    m = 0;
    n = 0;
    // 添加顶点
    for (const auto &vertex : subVertices) {
        if (graph.adj.find(vertex) != graph.adj.end()) {
            addNode(vertex);
        }
    }

    // 添加边
    for (const auto &vertex : subVertices) {
        if (graph.adj.find(vertex) != graph.adj.end()) {
            std::unordered_set<int> neighbors = graph.adj[vertex];
            for (const auto &neighbor : neighbors) {
                if (subVertices.find(neighbor) != subVertices.end()) {
                    addEdge(vertex, neighbor);
                }
            }
        }
    }
    
    // 计算每个节点的度
    computeDegrees();

    // 打印图的相关信息
    statistic();
}

Graph::~Graph() {}

// 社区节点集合转Graph，有必要， 改成构造函数了，需要将所有此函数替换为Graph（）
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

// 2. Greedy算法
Graph Graph::globalsearch(query_nodes& queryNodes) {
    clock_t start = clock(); // 记录开始时间
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
            degree[node] = -1;         // 将节点度数设为-1，标记已处理

            // 检查查询顶点集是否连通
            flag = isConnected(queryNodes, degree);
            if (!flag) {
                std::cout << "2:Query vertex set not connected！" << std::endl;
                break;
            }

            // 更新所有邻居的度数
            for (int neighbor : getNeighbors(node)) {
                neighborDegree = degree[neighbor];
                
                list[neighborDegree].erase(neighbor);
                list[neighborDegree - 1].insert(neighbor);
                degree[neighbor] = neighborDegree - 1;

                if (degree[neighbor] < lowestDegree) { lowestDegree = degree[neighbor]; }
            }
        }

        std::unordered_set<int> result;
        for (const auto& pair : degree) {
            if (pair.second != -1 && pair.second != 0) {
                result.insert(pair.first);  // 获取键
            }
        }
        result = getComponent(queryNodes, result);

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

    std::cout <<"执行到此！" << std::endl;
    flag = isConnected(queryNodes, degree);
    if (!flag) {
        std::cout << "查询顶点集不连通！" << std::endl;
        return {};
    }

    // result_end = getComponent(queryNodes, result_end);
    Graph ans(result_end, *this);
    std::cout << "得到结果为：" << k << std::endl;
    clock_t end = clock(); // 记录结束时间
    std::cout << "花费了" << (double)(end - start) / CLOCKS_PER_SEC << "秒" << std::endl;
    return ans;
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
    // 避免自环
    if (from == to) return;

    // 如果边已存在，则不增加计数
    if (adj[from].find(to) == adj[from].end()) {
        adj[from].insert(to);
        adj[to].insert(from);
        m++;  // 只有第一次插入时计数
    }
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
std::unordered_set<int> Graph::getComponent(query_nodes& queryNodes, std::unordered_set<int>& result_end) {
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