#include "Graph.h"

Graph::Graph(std::string path)
{
    readFromFile(path);
    computeDegrees();
    computeMinimumDegree();
    statistic();
}

Graph::Graph()
{
    m = 0;
    n = 0;
}

Graph Graph::getGraph(std::unordered_set<int>& subVertices)
{
    // 创建一个新的 Graph 对象
    Graph subGraph;

    // 添加顶点
    for (const auto& vertex : subVertices)
    {
        if (adj.find(vertex) != adj.end())
        {
            subGraph.addNode(vertex);
        }
    }

    // 添加边
    for (const auto& vertex : subVertices)
    {
        if (adj.find(vertex) != adj.end())
        {
            std::unordered_set<int> neighbors = adj[vertex];
            for (const auto& neighbor : neighbors)
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
    subGraph.computeMinimumDegree();
    subGraph.statistic();

    return subGraph;
}

Graph Graph::getReverseGraph() {
    Graph reverseGraph;
    for (const auto& pair : adj) {
        reverseGraph.addNode(pair.first);
        for (const int node : pair.second) {
            reverseGraph.addEdge(node, pair.first);
        }
    }
    return reverseGraph;
}

void Graph::addNode(int node)
{
    if (adj.find(node) == adj.end())
    {
        adj[node] = std::unordered_set<int>();
        n++;
    }
}

bool Graph::addEdge(int from, int to)
{
    if (from != to && adj[from].find(to) == adj[from].end())
    {
        adj[from].insert(to);
        m++;
        return true;
    }
    else return false;
}

void Graph::printGraph() {
    for (const auto& pair : adj) {
        std::cout << pair.first << ":";
        for (const int node : pair.second) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
}

/* 从邻接表格式文件中构建无向图 */
void Graph::buildGraphFromAdjFile(std::string fileName) {
    m = 0;
    n = 0;
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << " Error in File Openning." << std::endl;
        exit(1);
    }
    std::unordered_set<int> nodes;
    int sum = 0;
    int tmp = 0;
    std::string line;
    while (getline(file, line)) {
        std::istringstream iss(line);
        int node;
        iss >> node;
        addNode(node);
        int neighbor;
        while (iss >> neighbor) {
            addEdge(node, neighbor);
            addEdge(neighbor, node);
        }
    }
    file.close();
    this->printGraph();
}

/* 从文件中构建有向图 */
void Graph::buildDigraphFromFile(std::string fileName) {
    m = 0;
    n = 0;
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Error in File Openning." << std::endl;
        exit(1);
    }
    std::unordered_set<int> nodes;
    int sum = 0;
    int tmp = 0;
    std::string line;
    while (std::getline(file, line)) {
        sum++;
        int from;
        int to;
        if (sscanf(line.c_str(), "%d %d", &from, &to) == 2) {
            // std::cout << from << " " << to << std::endl;
            if (from >= 0 && to >= 0 && from != to) {}
            else {
                continue;
            }
            addNode(from);
            addNode(to);
            addEdge(from, to);
            nodes.insert(from);
            nodes.insert(to);
            if (from > tmp) {
                tmp = from;
            }
            if (to > from) {
                tmp = to;
            }
        }
        else {
            std::cerr << "Error: Incorrect line format" << std::endl;
        }
    }
    file.close();
}





void Graph::readFromFile(std::string fileName)
{
    m = 0;
    n = 0;
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        std::cerr << "Error: File not found." << std::endl;
        exit(1);
    }
    std::unordered_set<int> nodes;
    int sum = 0;
    int tmp = 0;
    std::string line;
    // 测试数据集的边是否适用无向图
    while (std::getline(file, line))
    {
        sum++;
        int from, to;
        if (sscanf(line.c_str(), "%d %d", &from, &to) == 2)
        {
            // std::cout << from << "  " << to << std::endl;
            if (0 <= from && 0 <= to && from != to)
            {
                ;
            }
            else
            {
                continue;
            }

            addNode(from);
            addNode(to);

            addEdge(from, to);
            addEdge(to, from);


            nodes.insert(from);
            nodes.insert(to);

            if (from > tmp)
                tmp = from;
            if (to > tmp)
                tmp = to;
        }
        else
        {
            std::cerr << "Error: Incorrect line format" << std::endl;
        }
    }
    file.close();
    m /= 2;
}

std::unordered_map<int, int> Graph::computeDegrees()
{
    degrees.clear();
    orderedNodes.clear();
    orderedNodes.resize(adj.size());

    Dmax = 0;

    for (auto& entry : adj)
    {
        degrees[entry.first] = entry.second.size();
        orderedNodes[entry.second.size()].insert(entry.first);
        if (entry.second.size() > Dmax)
            Dmax = entry.second.size();
    }

    return degrees;
}

int Graph::computeMinimumDegree()
{
    for (size_t i = 0; i < orderedNodes.size(); ++i)
    {
        if (!orderedNodes[i].empty())
        {
            minimumDegree = i;
            return i;
        }
    }
    return minimumDegree;
} // 有瑕疵，现在又找不到了

void Graph::statistic()
{
    // 节点数
    std::cout << "节点数为 " << n << std::endl;
    // 边数
    std::cout << "边数为 " << m << std::endl;
    // 最大度数
    std::cout << "最大度数为 " << Dmax << std::endl;
}

Graph::Graph(Graph& graph)
{
    // 拷贝邻接表
    adj = graph.adj; // unordered_map 自动进行深拷贝

    // 拷贝节点的度
    degrees = graph.degrees; // unordered_map 的深拷贝

    // 拷贝节点标签
    //nodeLables = graph.nodeLables; // 深拷贝unordered_map

    // 拷贝边标签
    //edgeLables = graph.edgeLables; // 深拷贝内部unordered_map和集合

    // 拷贝有相同度的节点的集合的向量
    orderedNodes = graph.orderedNodes; // 深拷贝向量和集合

    // 拷贝最小度
    minimumDegree = graph.minimumDegree;
}

Graph::~Graph() {}

int Graph::getNumberOfNodes()
{
    return adj.size();
}

std::unordered_map<int, int> Graph::getDegrees()
{
    return degrees;
}

std::unordered_set<int> Graph::getNeighbors(int node)
{
    // 获取并返回排序后的邻居节点
    std::vector<int> sortedNeighbors = sortNeighbors(node);

    return std::unordered_set<int>(sortedNeighbors.begin(), sortedNeighbors.end());
}

std::vector<int> Graph::sortNeighbors(int node)
{
    // 获取节点的邻居节点
    std::unordered_set<int>& neighbors = adj[node];

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

std::unordered_set<int> Graph::getNodes()
{
    std::unordered_set<int> nodes;
    for (auto& entry : adj)
    {
        nodes.insert(entry.first);
    }
    return nodes;
}

int Graph::computesubMinimumDegree(std::unordered_set<int>& nodes)
{
    if (nodes.size() == 0) return 0;
    int minDegree = INT_MAX;
    for (auto& node : nodes)
    {
        int degree = 0;
        for (auto& neighbor : getNeighbors(node))
        {
            if (nodes.find(neighbor) != nodes.end())
                degree++;
        }
        if (degree < minDegree)
            minDegree = degree;
    }
    return minDegree;
}

// 为什么查询顶点需要唯一
void Graph::search(std::unordered_set<int> H0, int k, std::unordered_set<int>& H)
{
    if (computesubMinimumDegree(H0) == k)
    {
        H = H0;
        std::cout << "得到结果为：" << H.size() << std::endl;
        return;
    }

    for (auto& node : H0)
    {
        for (auto& neighbor : getNeighbors(node))
        {
            if (H0.find(neighbor) == H0.end())
            {
                std::unordered_set<int> H1 = H0;
                H1.insert(neighbor);
                if (computesubMinimumDegree(H1) >= computesubMinimumDegree(H0))
                {
                    search(H1, k, H);
                    std::cout << "结束一次递归：" << neighbor << std::endl;
                    if (H.size() != 0) return;
                    std::cout << "没有得到结果：" << std::endl;
                }
            }
        }
    }
}

std::unordered_set<int> Graph::baseline_search(int v0, int k)
{
    std::cout << "结果：" << std::endl;
    std::vector<int> visited = std::vector<int>(n + 1, 0);
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

        for (auto& w : getNeighbors(v)) {
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
    std::vector<int> visited = std::vector<int>(n + 1, 0);
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

        for (auto& v : C) {
            int f = 0;
            for (auto& w : getNeighbors(v)) {
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
    std::vector<int> visited = std::vector<int>(n + 1, 0);
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

        for (auto& w : getNeighbors(v)) {
            if (C.find(w) == C.end() && degrees[w] >= k && !visited[w]) {
                queue.push(w); visited[w] = 1;
            }
        }
    }
    std::cout << "结果：" << computesubMinimumDegree(C) << std::endl;
    std::cout << "k：" << k << std::endl;
    return C;
}

// 全局搜索，不遍历全局，而使用深搜确定一定会返回解（如果有）
// std::unordered_set<int> Graph::globalsearch(std::unordered_set<int> C, int k)
// {
//     std::unordered_set<int> H = std::unordered_set<int>();
//     search(C, k, H);
//     return H;
// }

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
    if (upperBound(k)) {
        std::unordered_set<int> C = std::unordered_set<int>();
        return naiveCandidateGeneration(v0, k);
    }
    else {
        return std::unordered_set<int>();
    }
}

//CSM
// std::unordered_set<int> Graph::CSMframework(int v0, double gamma) {
//     std::unordered_set<int> H, A, B;
//     A.insert(v0);
//     for (auto &neighbor : getNeighbors(v0)) {
//         B.insert(neighbor);
//     }
//     int s = 0;

//     // s的限制条件未检查 e−γ (b |E|−|V |  (δ(G[H])+1)/2−1 c − |H| 
//     // exp(-gamma) * (n - ((n - m) / (double)(H.size() + 1))
//     // exp(-gamma) * ( (n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1) - H.size())
//     // !B.empty() && s <= exp(-gamma) * ( (n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1) - H.size())
//     int x =0 ;
//     while (!B.empty() && x++<=10)
// 	//&& s <= exp(-gamma) * ( std::ceil((n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1)) - H.size()))
//     { 
//     	// std::cout << " succeed in!!" << exp(-gamma) * ( std::ceil((n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1)) - H.size()) << std::endl;
//         std::unordered_set<int> V;
//         int mintmp = 0;
//         for (auto &node : B)
//         {
//             int tmp = 0;
//             for (auto &neighbor : getNeighbors(node)) 
//             {
//                 if (A.find(neighbor) != A.end())
//                 tmp++;
//             }
//             if (tmp > mintmp) 
//             {
//                 mintmp = tmp;
//                 // V.insert(node);
//                 // std::cout << " node" << node << std::endl;
//             }
//         }
//         for (auto &node : B)
//         {
//             int tmp = 0;
//             for (auto &neighbor : getNeighbors(node)) 
//             {
//                 if (A.find(neighbor) != A.end())
//                 tmp++;
//             }
//             if (tmp == mintmp) 
//             {
//                 V.insert(node);
//                 std::cout<< "tmp " << tmp << std::endl;
//                 std::cout << " node " << node << std::endl;
//             }
//         }
//         for (auto &v : V) {
//         A.insert(v);
//         B.erase(v);
//         s++;
//         }
//         if (computesubMinimumDegree(A) > computesubMinimumDegree(H)) {
//             H = A;
//             s = 0;
//             if (computesubMinimumDegree(H) == std::min(degrees[v0], static_cast<int>((1 + std::sqrt(9 + 8 * (n - m))/2)))) {
//                 return H;
//             }
//         }
//         for (auto &v : V){
//         for (auto &neighbor : getNeighbors(v))
//         {
//             if (B.find(neighbor) == B.end() && degrees[neighbor] >= computesubMinimumDegree(H)) {
//                 B.insert(neighbor);
//                 std::cout << "neighbor " << neighbor << std::endl;
//             }
//         }
//         }
//         std::cout << "s   " << s << std::endl;
//         std::cout << "H   " << computesubMinimumDegree(H) << std::endl;
//         }
//         int k = computesubMinimumDegree(H);
//         std::cout << "k   " << k << std::endl;
//         std::unordered_set<int> C = generateCandidates(H, k);
//         return maxcore(C, v0, k);
//     //return H;
// }


std::unordered_set<int> Graph::CSMframework(int v0, double gamma) {
    std::unordered_set<int> H, A, B;
    A.insert(v0);
    for (auto& neighbor : getNeighbors(v0)) {
        B.insert(neighbor);
    }
    int s = 0;

    // s的限制条件未检查 e−γ (b |E|−|V |  (δ(G[H])+1)/2−1 c − |H| 
    // exp(-gamma) * (n - ((n - m) / (double)(H.size() + 1))
    // exp(-gamma) * ( (n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1) - H.size())
    // !B.empty() && s <= exp(-gamma) * ( (n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1) - H.size())
    int x = 0;
    while (!B.empty() && x++ <= 10000)
        //&& ((s <= exp(-gamma) * ( std::ceil((n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1)) - H.size())) || (computesubMinimumDegree(H) == 1)))
    {
        // std::cout << " succeed in!!" << exp(-gamma) * ( std::ceil((n - m) / ((computesubMinimumDegree(H) + 1)/2 - 1)) - H.size()) << std::endl;
        int v = 0;
        int mintmp = 0;
        for (auto& node : B) {
            int tmp = 0;
            for (auto& neighbor : getNeighbors(node)) {
                if (A.find(neighbor) != A.end())
                    tmp++;
            }
            if (tmp > mintmp) {
                mintmp = tmp;
                v = node;
                // std::cout << " node" << node << std::endl;
            }
        }

        A.insert(v);
        B.erase(v);
        s++;
        if (computesubMinimumDegree(A) > computesubMinimumDegree(H)) {
            H = A;
            s = 0;
            if (computesubMinimumDegree(H) == std::min(degrees[v0], static_cast<int>(floor((1 + std::sqrt(9 + 8 * (n - m)) / 2))))) {
                return H;
            }
        }

        for (auto& neighbor : getNeighbors(v)) {
            if (B.find(neighbor) == B.end() && A.find(neighbor) == A.end() && degrees[neighbor] > computesubMinimumDegree(H)) {
                B.insert(neighbor);
                std::cout << "neighbor " << neighbor << std::endl;
            }
            if (degrees[neighbor] <= computesubMinimumDegree(H)) { break; }
        }

        std::cout << "s   " << s << std::endl;
        std::cout << "H   " << computesubMinimumDegree(H) << std::endl;
    }
    int k = computesubMinimumDegree(H);
    std::cout << "k   " << k << std::endl;

    //第一步，得到了查询顶点的搜索空间，相当于v0，图G
    //有效的输出：H,A
    //std::unordered_set<int> C = generateCandidates(H, k);
    return maxcore(A, v0);
    // return H;
}


// generateCandidates
std::unordered_set<int> Graph::generateCandidates(std::unordered_set<int>& H, int k) {
    std::unordered_set<int> C = H;
    for (auto& v : H) {
        for (auto& neighbor : getNeighbors(v)) {
            if (degrees[neighbor] < k) {
                C.erase(neighbor);
            }
        }
    }
    return C;
}

// maxcore
std::unordered_set<int> Graph::maxcore(std::unordered_set<int>& C, int v0) {
    Graph new_G;
    new_G = getGraph(C);
    std::unordered_set<int> queryNodes = { v0 };
    std::unordered_set<int> solution = {};
    solution = new_G.greedy(queryNodes);
    std::cout << new_G.computesubMinimumDegree(solution) << " ";
    return solution;
}

bool Graph::isQuerySetConnected(const std::unordered_set<int>& queryNodes, const std::unordered_map<int, std::unordered_set<int>>& graph) {
    if (queryNodes.empty()) return true;

    std::unordered_set<int> visited;
    std::queue<int> q;

    // 从查询节点中的任意一个开始遍历
    auto startIt = queryNodes.begin();
    q.push(*startIt);
    visited.insert(*startIt);

    while (!q.empty()) {
        int current = q.front();
        q.pop();

        // 遍历当前节点的所有邻居
        for (const auto& neighbor : graph.at(current)) {
            // 如果邻居也在查询节点集中且未被访问过，则加入队列
            if (queryNodes.count(neighbor) && !visited.count(neighbor)) {
                q.push(neighbor);
                visited.insert(neighbor);
            }
        }
    }

    // 如果所有查询节点都被访问过，则说明它们是连通的
    return visited.size() == queryNodes.size();
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
std::unordered_set<int> Graph::greedy(std::unordered_set<int>& queryNodes)
{
    std::unordered_map<int, std::unordered_set<int>> resultGraph = adj;
    std::vector<std::unordered_set<int>> list = orderedNodes;
    bool foundQueryNode = false;
    std::unordered_set<int> result_is;
    // 迭代移除度数最低的节点
    while (!list.empty()) {
        int min_degree = INT_MAX;

        // 找到当前最低度数的节点
        for (int i = 0; i < list.size(); i++) {
            if (i < min_degree && !list[i].empty()) {
                min_degree = i;
            }
        }

        // 检查是否包含查询节点
        for (auto& pair : list[min_degree]) {
            if (queryNodes.count(pair)) {
                foundQueryNode = true;
                // cout << "Found query node!" << pair << endl;
            }
        }

        if (foundQueryNode) break;

        // 移除节点
        removeNode(min_degree, resultGraph, list);

        // 检查查询顶点集是否连通
        if (!isQuerySetConnected(queryNodes, resultGraph)) {
            std::cout << "已断开" << std::endl;
            std::unordered_set<int> result;
            for (const auto& pair : resultGraph) {
                result.insert(pair.first);  // 获取键
            }
            if (computesubMinimumDegree(result) > computesubMinimumDegree(result_is)) {
                result_is = result;
            }
            // 如果不连通，则恢复节点并结束循环
            // list[min_degree].insert(nodeToRemove); // 假设 removeNode 中有 nodeToRemove 的记录
            // resultGraph[nodeToRemove] = getNeighbors(nodeToRemove); // 恢复节点及其邻居关系
            break;
        }
        std::unordered_set<int> result;
        for (const auto& pair : resultGraph) {
            result.insert(pair.first);  // 获取键
        }
        if (computesubMinimumDegree(result) > computesubMinimumDegree(result_is)) {
            result_is = result;
        }
    }
    return result_is;
}


//std::unordered_set<int> Graph::greedy(std::unordered_set<int>& queryNodes)
//{
//    std::unordered_map<int, std::unordered_set<int>> resultGraph = adj;
//    std::vector<std::unordered_set<int>> list = orderedNodes;
//    bool foundQueryNode = false;
//    std::unordered_set<int> result_is;
//    int result_is_minDegree = INT_MAX;
//    // 迭代移除度数最低的节点
//    while (!list.empty()) {
//        int min_degree = INT_MAX;
//
//        // 找到当前最低度数的节点
//        for (int i = 0; i < list.size(); i++) {
//            if (i < min_degree && !list[i].empty()) {
//                min_degree = i;
//            }
//        }
//
//        // 检查是否包含查询节点
//        for (auto& pair : list[min_degree]) {
//            if (queryNodes.count(pair)) {
//                foundQueryNode = true;
//                // cout << "Found query node!" << pair << endl;
//            }
//        }
//
//        if (foundQueryNode) break;
//
//        // 移除节点
//        removeNode(min_degree, resultGraph, list);
//
//        // 检查查询顶点集是否连通
//        if (!isQuerySetConnected(queryNodes, resultGraph)) {
//            std::cout << "已断开" << std::endl;
//            std::unordered_set<int> result;
//            for (const auto& pair : resultGraph) {
//                result.insert(pair.first);  // 获取键
//            }
//            if (computesubMinimumDegree(result) > computesubMinimumDegree(result_is)) {
//                result_is = result;
//            }
//            // 如果不连通，则恢复节点并结束循环
//            // list[min_degree].insert(nodeToRemove); // 假设 removeNode 中有 nodeToRemove 的记录
//            // resultGraph[nodeToRemove] = getNeighbors(nodeToRemove); // 恢复节点及其邻居关系
//            break;
//        }
//        std::unordered_set<int> result;
//        for (const auto& pair : resultGraph) {
//            result.insert(pair.first);  // 获取键
//        }
//        int result_minDegree = computesubMinimumDegree(result)
//            if (result_minDegree > result_is_minDegree) {
//                result_is = result;
//                result_is_minDegree = result_minDegree;
//            }
//    }
//    return result_is;
//}

/* 多源BFS函数 */
void Graph::multiSourceBFS(Graph& graph, const std::vector<int> source, std::unordered_map<int, std::unordered_map<int, int>>& distIndex) {

    std::queue<int> q;
    std::unordered_set<int> visited;

    for (const auto& s : source) {
        q.push(s);
        for (int i = 0; i <= graph.getNumberOfNodes(); i ++) {
            distIndex[s][i] = INT_MAX;
        }
        distIndex[s][s] = 0;
    }

    std::cout << graph.getNumberOfNodes() << std::endl;

    while (!q.empty()) {
        int curr = q.front();
        q.pop();

        std::cout << curr << ":" << std::endl;
        
        for (const auto& neighbor : graph.getNeighbors(curr)) {
            if (visited.count(neighbor)) {
                for (const auto& s : source) {
                    if (distIndex[s][curr] != INT_MAX && distIndex[s][neighbor] > distIndex[s][curr] + 1) {
                        distIndex[s][neighbor] = distIndex[s][curr] + 1;
                    }
                }
            }
            else if (!visited.count(neighbor) && std::find(source.begin(), source.end(), neighbor) == source.end()) {
                q.push(neighbor);
                visited.insert(neighbor);
                for (const auto& s : source) {
                    if (distIndex[s][curr] != INT_MAX && distIndex[s][neighbor] > distIndex[s][curr] + 1) {
                        distIndex[s][neighbor] = distIndex[s][curr] + 1;
                    }
                }
            }
            else if (std::find(source.begin(), source.end(), neighbor) != source.end() && !visited.count(neighbor)) {
                for (const auto& s : source) {
                    if (distIndex[s][curr] != INT_MAX && distIndex[s][neighbor] > distIndex[s][curr] + 1) {
                        distIndex[s][neighbor] = distIndex[s][curr] + 1;
                    }
                    if (s != neighbor) {
                        for (int i = 0; i <= n; i++) {
                            if (distIndex[s][neighbor] != INT_MAX && distIndex[neighbor][i] != INT_MAX && distIndex[s][i] > distIndex[s][neighbor] + distIndex[neighbor][i]) {
                                distIndex[s][i] = distIndex[s][neighbor] + distIndex[neighbor][i];
                            }
                        }
                    }
                }
                visited.insert(neighbor);
            }
        }
        
        graph.printDistIndex(source, distIndex);

    }

}

/* 测试用打印distIndex函数 */
void Graph::printDistIndex(std::vector<int> source, std::unordered_map<int, std::unordered_map<int, int>>& distIndex) {
    for (const auto& s : source) {
        std::cout << "From source " << s << ":\n";
        for (int i = 0; i < this->getNumberOfNodes(); i++) {
            std::cout << "  To node " << i << " : Distance = " << (distIndex[s][i] == INT_MAX ? -1 : distIndex[s][i]) << std::endl;
        }
    }
}


/* BasicEnum中的search过程 */
void Graph::basicEnumSearch(std::vector<std::vector<int>> &P, std::vector<int>& p, int v,const int k, std::unordered_map<int, std::unordered_map<int, int>>& distIndex) {

    //将末结点赋值给v'
    int v_prime = p.back();
    P.push_back(p);

    // 递归结束条件
    if (p.size() == k || k == 0) {
        return;
    }
    //递归搜索
    for (const auto& neighbor : this->getNeighbors(v_prime)) {
        if (std::find(p.begin(), p.end(), neighbor) == p.end() && p.size() + (distIndex[v][neighbor]) < k) {
            p.push_back(neighbor);
            this->basicEnumSearch(P, p, v, k, distIndex);
            p.pop_back();
        }
    }
}

/* 测试用打印P函数 */
void Graph::printP(std::vector<std::vector<int>>& P) {
    for (const auto& path : P) {
        for (int node : path) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
}

/* BasicEnum算法 */
void Graph::basicEnum(const Graph &graph, const std::vector<std::pair<int,int>> &query, const int k){
    
    std::vector<int> source;  //Q.s
    std::vector<int> target;  //Q.t

    for (auto &q : query)    {
        source.push_back(q.first);
        target.push_back(q.second);
    } 
}

