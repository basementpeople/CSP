#include "TreeIndex.h"
// struct Item {
//     std::pair<int, int> priority; // 优先级 {22, 11}
//     int value;                    // 值 1

//     // 构造函数
//     Item(std::pair<int, int> p, int v) : priority(p), value(v) {}
// };
struct Item {
    std::tuple<int, int, int> priority; // 优先级 {22, 11, 5}
    int value;                          // 值 1

    // 构造函数
    Item(std::tuple<int, int, int> p, int v) : priority(p), value(v) {}
};
// struct CompareItem {
//     bool operator()(const Item& a, const Item& b) const {
//         // 先比较第一个优先级
//         if (a.priority.first != b.priority.first) {
//             return a.priority.first < b.priority.first;
//         }
//         // 如果第一个优先级相同，比较第二个优先级
//         return a.priority.second < b.priority.second;
//     }
// };

// struct CompareItem {
//     bool operator()(const Item& a, const Item& b) const {
//         // 自定义比较逻辑，逐个比较 tuple 中的元素
//         if (std::get<0>(a.priority) != std::get<0>(b.priority)) {
//             return std::get<0>(a.priority) < std::get<0>(b.priority); // 比较第一个优先级
//         } else if (std::max(0,std::get<1>(a.priority) - std::get<2>(a.priority)) != std::max(0,std::get<1>(b.priority) - std::get<2>(b.priority))) {
//             return std::max(0,std::get<1>(a.priority) - std::get<2>(a.priority)) < std::max(0,std::get<1>(b.priority) - std::get<2>(b.priority)); // 比较第二个优先级
//         } else {
//             return std::get<1>(a.priority) < std::get<1>(b.priority); // 比较第三个优先级
//         }
//     }
// };
struct CompareItem {
    bool operator()(const Item& a, const Item& b) const {
        // 自定义比较逻辑，逐个比较 tuple 中的元素
        if (std::get<0>(a.priority) != std::get<0>(b.priority)) {
            return std::get<0>(a.priority) < std::get<0>(b.priority); // 比较第一个优先级
        } else if (std::get<1>(a.priority) - std::get<2>(a.priority) != std::get<1>(b.priority) - std::get<2>(b.priority)) {
            return std::get<1>(a.priority) - std::get<2>(a.priority) < std::get<1>(b.priority) - std::get<2>(b.priority); // 比较第二个优先级
        } else {
            return std::get<1>(a.priority) < std::get<1>(b.priority); // 比较第三个优先级
        }
    }
};

class UnionFind {
    public:
        UnionFind(int n) : parent(n), rank(n, 0), setNum(n) {
            for (int i = 0; i < n; ++i) {
                parent[i] = i;
            }
        }
    
        // 递归地实现，简洁易懂
        // int find(int x) {
        //     if (parent[x] != x) {
        //         parent[x] = find(parent[x]);
        //     }
        //     return parent[x];
        // }

        // 非递归的实现，避免递归调用栈溢出
        int find(int x) {
            if (x < 0 || x >= parent.size()) {
                throw std::out_of_range("Index out of range in UnionFind::find");
            }
        
            // Step 1: Find the root node
            int root = x;
            while (parent[root] != root) {
                root = parent[root];
            }
        
            // Step 2: Path compression
            while (parent[x] != root) {
                int next = parent[x];       // Save the next node in the path
                parent[x] = root;           // Directly connect the current node to the root
                x = next;                   // Move to the next node
            }
        
            return root;
        }
    
        void unite(int x, int y) {
            if (x < 0 || x >= parent.size() || y < 0 || y >= parent.size()) {
                throw std::out_of_range("Index out of range in UnionFind::unite");
            }
            int rootX = find(x);
            int rootY = find(y);
            if (rootX != rootY) {
                if (rank[rootX] < rank[rootY]) {
                    parent[rootX] = rootY;
                } else {
                    parent[rootY] = rootX;
                    if (rank[rootX] == rank[rootY]) {
                        rank[rootX]++;
                    }
                }
                setNum--;  // 更新集合数量
            }
        }
    
        int getSetNum() const {
            return setNum;
        }

    private:
        std::vector<int> parent;
        std::vector<int> rank;
        int setNum;  // 记录集合数量
};

TreeIndex::TreeIndex(Graph &graph)
{
    count_port = graph.getN();
    setNeighbors(graph);
    setDegrees(graph);

    // 不用管，没太主要的作用
    computeCoreIndex(graph);

    // 划分了每个shell中的不同连通分量
    identifyAndStoreComponents(graph);

    // 构建父子连通分量关系:一个子可能有多个不同的父,不同的父可能有相同的子
    buildParentChildRelationships(graph);

}

void TreeIndex::setNeighbors(Graph &graph) {
    adj = graph.getAdj();
}

void TreeIndex::setDegrees(Graph &graph) {
    degrees = graph.getDegrees();
}

std::unordered_set<int> TreeIndex::getNeighbors(int node)
{
    // 获取并返回排序后的邻居节点
    std::unordered_set<int> &neighbors = adj[node];

    return neighbors;
}

// 相当于做了一个抽象，点对应核心度，中间插了个核心ID保证连续
// void TreeIndex::computeCoreIndex(Graph &graph)
// {
//     shell_count = 0;
//     coreIndex = CoreGroup::coreGroupsAlgorithm(graph);
    
//     // 计算每个核心的最小度数并存储在 coreMinimumDegree 中
//     std::set<int> coreIndexSet;  // 用于存储核心度的集合
//     for (auto &pair : coreIndex)
//     {
//         coreIndexSet.insert(pair.second);
//     }
//     int node = 0;
//     while (!coreIndexSet.empty())
//     {
//         coreMinimumDegree[node] = *coreIndexSet.begin();
//         for (auto &pair : coreIndex)
//         {
//             if (pair.second == *coreIndexSet.begin())
//             {
//                 pair.second = node;
//             }
//         }
//         coreIndexSet.erase(coreIndexSet.begin());
//         node++;
//     }

//     shell_count = node;
//     std::cout << "核心层数：" << shell_count << std::endl;
// }
void TreeIndex::computeCoreIndex(Graph &graph)
{
    shell_count = 0;
    coreIndex = CoreGroup::coreGroupsAlgorithm(graph);
    
    // 计算每个核心的最小度数并存储在 coreMinimumDegree 中
    std::set<int> coreIndexSet;  // 用于存储核心度的集合
    for (auto &pair : coreIndex)
    {
        coreIndexSet.insert(pair.second);
    }
    int node = 0;
    while (!coreIndexSet.empty())
    {
        coreMinimumDegree[node] = *coreIndexSet.begin();

        // 1.16 优化，如果核心度大于当前记录的层数
        if (coreMinimumDegree[node] > shell_count) {
            shell_count = coreMinimumDegree[node];
        }

        for (auto &pair : coreIndex)
        {
            if (pair.second == *coreIndexSet.begin())
            {
                pair.second = node;
            }
        }
        coreIndexSet.erase(coreIndexSet.begin());
        node++;
    }

    std::cout << "核心（Shell）层数： " << shell_count << std::endl;

    // 打印最大k
    // for (const auto &pair : coreMinimumDegree)
    // {
    //     std::cout << "核心 " << pair.first << " 的最小度数: " << pair.second << std::endl;
    // }
}

// 划分了每个shell中的不同连通分量
void TreeIndex::identifyAndStoreComponents(Graph &graph)
{
    std::unordered_map<int, std::unordered_set<int>> nodesInShell; // 核心ID  对应  点集

    for (const auto &pair : coreIndex)
    {

        nodesInShell[pair.second].insert(pair.first);
    }

    // 根据每一层k-shell的节点，找到每一层的连通分量
    for (const auto &shell : nodesInShell)
    {
        int shellIndex = shell.first;
        const auto &nodes = shell.second;

        // BFS
        std::unordered_set<int> visited;
        for (int node : nodes)
        {

            if (visited.count(node) == 0)
            {
                std::queue<int> queue;
                std::unordered_set<int> component;
                queue.push(node);
                visited.insert(node);

                while (!queue.empty())
                {
                    int currentNode = queue.front();
                    queue.pop();
                    component.insert(currentNode);

                    for (int neighbor : graph.getNeighbors(currentNode))
                    {
                        // 只找未访问过的同一层邻居
                        if (nodes.count(neighbor) > 0 && visited.count(neighbor) == 0)
                        {
                            queue.push(neighbor);
                            visited.insert(neighbor);
                        }
                    }
                }

                int componentId = nextComponentId++;
                layerToComponentToNodes[shellIndex][componentId] = component;
                ComponentToNodes[componentId] = component;
                for (int nodeId : component) { 
                    nodeToComponentId[nodeId] = componentId;
                }
            }
        }
    }
}

void TreeIndex::printComponents() const
{
    std::cout << "Layer to Component to Nodes Mapping:" << std::endl;
    for (const auto &layer : layerToComponentToNodes)
    {
        std::cout << "Shell Index: " << layer.first << " the k is : " << coreMinimumDegree.at(layer.first) << std::endl;
        for (const auto &comp : layer.second)
        {
            std::cout << "  Component ID: " << comp.first << ", Nodes: ";
            for (int node : comp.second)
            {
                std::cout << node << " ";
            }
            std::cout << std::endl;
        }
    }
}

// 建立了不同shell连通分量之间的父子关系
void TreeIndex::buildParentChildRelationships(Graph &graph)
{
    // 父子索引的处理
    int size = layerToComponentToNodes.size();
    for (int currentLevel = 0; currentLevel < size; ++currentLevel)
    {
        auto &currentLayerComponents = layerToComponentToNodes[currentLevel];

        for (auto &component : currentLayerComponents)
        {
            int currentComponentId = component.first;
            auto &nodesInCurrentComponent = component.second;

            for (int node : nodesInCurrentComponent)
            {
                for (int neighbor : graph.getNeighbors(node))
                {
                    if (coreMinimumDegree[coreIndex[node]] >= coreMinimumDegree[coreIndex[neighbor]])
                    {
                        continue;
                    }
                    if (nodeToComponentId.find(neighbor) != nodeToComponentId.end() && nodeToComponentId[neighbor] != currentComponentId)
                    {
                        int neighborComponentId = nodeToComponentId[neighbor];
                        connectedComponentParent[currentLevel][currentComponentId] = neighborComponentId;
                        ComponentParent[currentComponentId].insert(neighborComponentId);
                    }
                }
            }
            // 顶层分量的父分量是-1
            if (connectedComponentParent[currentLevel].find(currentComponentId) == connectedComponentParent[currentLevel].end() && ComponentParent[currentComponentId].empty())
            {
                connectedComponentParent[currentLevel][currentComponentId] = -1;
                ComponentParent[currentComponentId].insert(-1);
            }
        }
    }
    for (int currentLevel = size - 1; currentLevel >= 0; currentLevel--)
    {

        auto &currentLayerComponents = layerToComponentToNodes[currentLevel];

        for (auto &component : currentLayerComponents)
        {
            int currentComponentId = component.first;
            auto &nodesInCurrentComponent = component.second;

            for (int node : nodesInCurrentComponent)
            {
                // 再打印父索引加入了几个分量 通过这样来判断是否有父子关系遗漏
                for (int neighbor : graph.getNeighbors(node))
                {
                    if (coreMinimumDegree[coreIndex[node]] <= coreMinimumDegree[coreIndex[neighbor]])
                    {
                        continue;
                    }
                    if (nodeToComponentId.find(neighbor) != nodeToComponentId.end() && nodeToComponentId[neighbor] != currentComponentId)
                    {
                        int neighborComponentId = nodeToComponentId[neighbor];
                        connectedComponentChildren[currentLevel][currentComponentId].insert(neighborComponentId);
                        ComponentChildren[currentComponentId].insert(neighborComponentId);
                    }
                }
            }
            // 底层分量的子分量是-1
            if (connectedComponentChildren[currentLevel].find(currentComponentId) == connectedComponentChildren[currentLevel].end())
            {
                connectedComponentChildren[currentLevel][currentComponentId].insert(-1);
                ComponentChildren[currentComponentId].insert(-1);
            }
        }
    }

}


std::unordered_set<int> TreeIndex::findTopComponents(std::vector<int> &queryNodes, int k)
{
    std::unordered_set<int> topComponents;
    std::unordered_set<int> queryComponents;
    std::queue<int> bfsQueue;
    std::unordered_set<int> visited;

    // 初始化：将查询节点的组件ID添加到队列和已访问集合中
    for (int node : queryNodes)
    {
        int componentId = nodeToComponentId[node];

        if (visited.find(componentId) == visited.end()) {

            queryComponents.insert(componentId);
            bfsQueue.push(componentId);
            visited.insert(componentId);

        }
    }

    while (!bfsQueue.empty())
    {
        int currentComponentId = bfsQueue.front();
        bfsQueue.pop();

        // 探测所有父分量
        for (int parentComponentId : ComponentParent[currentComponentId])
        {
            //std::cout<<"the parent compId is "<<parentComponentId<<std::endl;
            if (parentComponentId == -1) {
                // 如果父分量为-1，保存当前分量
                topComponents.insert(currentComponentId);
            }
            else {
                if (visited.find(parentComponentId) == visited.end()) {
                    bfsQueue.push(parentComponentId);
                    visited.insert(parentComponentId);

                }
            }
        }

        // 探测所有子分量，意义何在
        // 1.9 意义在于主动找到Q的点大于等于k的子分量，并将其的顶层分量加入结果集，这一步到底对不对呢
        for (int childComponentId : ComponentChildren[currentComponentId])
        {
            if (visited.find(childComponentId) == visited.end())
            {
                // 确保所有子分量的节点核心等级必须大于等于k
                bool eligible = true;
                for (int childNode : ComponentToNodes[childComponentId])
                {
                    if (coreMinimumDegree[coreIndex[childNode]] < k)
                    {
                        eligible = false;
                        break;
                    }
                }
                if (eligible)
                {
                    bfsQueue.push(childComponentId);
                    visited.insert(childComponentId);
                }
            }
        }
    }

    return topComponents;
}

// 
std::unordered_set<int> TreeIndex::findKCoreSubgraph(std::vector<int> &queryNodes, int k)
{
    // 初始化最小核心索引为最大整数，以便后续比较
    int minCoreIndex = INT_MAX;

    // 遍历查询节点，确定最小核心索引
    for (int node : queryNodes)
    {
        // 如果节点不在核心索引中，抛出异常
        // if (coreIndex.find(node) == coreMinimumDegree.end())
        if (coreIndex.find(node) == coreIndex.end())
        {
            throw std::runtime_error("Error: Node does not exist in core index.");
        }
        std::cout << "the coreIndex of node " << node << " is " << coreMinimumDegree[coreIndex[node]]<< std::endl;
        // 更新最小核心索引
        minCoreIndex = std::min(minCoreIndex, coreMinimumDegree[coreIndex[node]]);
    }

    // 找到最大的公共shell，如果这个shell的核心度>=k,不改变k；如果核心度<k，则k=核心度
    // int k_tmp = findCommenShell();
    // if (k_tmp < k) {
    //     k = k_tmp;
    // }

    // 确定k值
    if (k == -1) {
    k = minCoreIndex;
    }
    std::cout << "the k is : " << k << std::endl;
    // -------------------------- 以上的解决方案暗含了这样一个定理，即：核心分解后，核心度就是社区k的大小
    // -------------------------- 但事实上，核心度不一定就是最小度，因为可能会有两个集中的子图被连在一起，但核心度都很高，实际的k值只有2

    // 这里我们增加一个函数，去寻找最大的公共连通分量 1.11 不在这里加
    // k = getKfromIndex(queryNodes);
    // std::cout << "the k is : " << k << std::endl;
     
    // if (k == -1) {
    //     k = minCoreIndex;
    // }

    // 开始计时
    // auto start = std::chrono::high_resolution_clock::now();
    // 用于存储所有顶层分量的ID
    std::unordered_set<int> topComponentIds;

    // 找到所有父分量为-1的顶层分量
    topComponentIds = findTopComponents(queryNodes, k);

    // 中间计时点
    // auto mid = std::chrono::high_resolution_clock::now();
    // 计算中间耗时
    // auto midduration = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start);
    // 输出找到父索引所需时间
    // std::cout << "找到父索引所需时间为 :" << midduration.count() << " 毫秒" << std::endl;
    // 输出顶层分量的大小
    // std::cout << "the topComponentIds size is :" << topComponentIds.size() << std::endl;

    // 用于存储最终结果的节点集合
    std::unordered_set<int> resultNodes;
    // 用于存储最终结果
    std::unordered_set<int> res;
    // 用于广度优先搜索的队列
    std::queue<int> componentQueue;

    // 将顶层分量的ID加入队列
    for (int compId : topComponentIds)
    {
        componentQueue.push(compId);
    }

    // 广度优先搜索所有分量
    while (!componentQueue.empty())
    {
        int currentId = componentQueue.front();
        componentQueue.pop();

        // 如果当前分量已经在结果集中，则跳过
        if (resultNodes.find(currentId) != resultNodes.end())
            continue;

        // 将当前分量添加到结果集中
        resultNodes.insert(currentId);

        // 遍历当前分量的所有子分量
        for (int childCompId : ComponentChildren[currentId])
        {
            // 1.10 目的在于获得最小社区
            // if (visitedComponents.find(childCompId) != visitedComponents.end()) {

            // 获取子分量中的一个节点
            int childNode;
            for (int node : ComponentToNodes[childCompId])
            {
                childNode = node;
                break;
            }

            // 如果子分量的节点的核心度大于等于k，将子分量加入队列
            if (coreMinimumDegree[coreIndex[childNode]] >= k)
            {
                componentQueue.push(childCompId);
            }

            // }

        }
    }

    // 将结果集中的所有连通分量的点添加到最终结果中
    for (int compId : resultNodes)
    {
        for (int node : ComponentToNodes[compId])
        {
            res.insert(node);
        }
    }

    // 结束计时
    // auto end = std::chrono::high_resolution_clock::now();
    // 计算找到子索引所需时间
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid);
    // 输出找到子索引所需时间
    //std::cout << "找到子索引所需时间为 : " << duration.count() << " 毫秒" << std::endl;

    // 返回最终结果
    return res;
}

std::unordered_set<int> TreeIndex::getChildShell(int shell) {
    return ComponentChildren[shell];
}

std::unordered_set<int> TreeIndex::getParentShell(int shell) {
    return ComponentParent[shell];
}

int TreeIndex::findCommenShell(std::vector<int>& queryNodes) {
    int componentID = -1;
    std::unordered_map<int, query_nodes> connectShell;

    // 遍历每个顶点对应连通分量的所有子分量
    for (auto& p : queryNodes) {
        std::queue<int> que;
        std::unordered_set<int> visited;
        que.push(nodeToComponentId[p]);
        visited.insert(nodeToComponentId[p]);

        while (!que.empty()) {
            componentID = que.front();
            que.pop();

            connectShell[p].insert(componentID);
            for (auto &childrenSet: getChildShell(componentID)) {
                if (visited.find(childrenSet) == visited.end() && childrenSet!= -1) { // 防止耗时过长
                    visited.insert(childrenSet);
                    que.push(childrenSet);
                }
            }
        }

        que.push(nodeToComponentId[p]);
        while (!que.empty()) {
            componentID = que.front();
            que.pop();

            connectShell[p].insert(componentID);
            for (auto &parentSet: getParentShell(componentID)) {
                if (visited.find(parentSet) == visited.end() && parentSet!= -1) { // 防止耗时过长
                    visited.insert(parentSet);
                    que.push(parentSet);
                }
            }
        }

    }

    // 找到chilren中所有元素的并集
    query_nodes result;  // 并集
    std::unordered_multiset<int> countSum;
    for (auto& node : connectShell) {
        countSum.insert(node.second.begin(), node.second.end());
        result.insert(node.second.begin(), node.second.end());
    }

    // 找到result中元素个数等于query.size()的元素
    int minSize = queryNodes.size();
    for (auto& c : result) {
        if (countSum.count(c) != minSize) {
            result.erase(c);
        }
        // if (countSum.count(c) == minSize) {
        //     std::cout << "存在" << c << std::endl;
        // }
    }

    // 找到result中核心度最高的元素
    int maxK = -1;
    for (auto& c : result) {
        int tmp = coreMinimumDegree[coreIndex[*ComponentToNodes[c].begin()]];
        // std::cout << "coreIndex: " << tmp << std::endl;
        if (tmp > maxK) {
            maxK = tmp;
            // std::cout << "maxK: " << maxK << std::endl;
        }
    }

    // std::cout << "result: " << maxK << std::endl;
    return maxK;
}

std::unordered_set<int> TreeIndex::findKCoreSubgraph_d1(std::vector<int> &queryNodes)
{
    // 初始化最小核心索引为最大整数，以便后续比较
    int minCoreIndex = INT_MAX;

    // 遍历查询节点，确定最小核心索引
    for (int node : queryNodes)
    {
        // 如果节点不在核心索引中，抛出异常
        // if (coreIndex.find(node) == coreMinimumDegree.end())
        if (coreIndex.find(node) == coreIndex.end())
        {
            throw std::runtime_error("Error: Node does not exist in core index.");
        }
        std::cout << "the coreIndex of node " << node << " is " << coreMinimumDegree[coreIndex[node]]<< std::endl;
        // 更新最小核心索引
        minCoreIndex = std::min(minCoreIndex, coreMinimumDegree[coreIndex[node]]);
    }

    // 找到最大的公共shell，如果这个shell的核心度>=k,不改变k；如果核心度<k，则k=核心度
    int k = findCommenShell(queryNodes);
    if (k >= minCoreIndex || k == -1) {
        k = minCoreIndex;
    }

    // 确定k值
    std::cout << "the k is : " << k << std::endl;
    // -------------------------- 以上的解决方案暗含了这样一个定理，即：核心分解后，核心度就是社区k的大小
    // -------------------------- 但事实上，核心度不一定就是最小度，因为可能会有两个集中的子图被连在一起，但核心度都很高，实际的k值只有2

    // 开始计时
    // auto start = std::chrono::high_resolution_clock::now();
    // 用于存储所有顶层分量的ID
    std::unordered_set<int> topComponentIds;

    // 找到所有父分量为-1的顶层分量
    topComponentIds = findTopComponents(queryNodes, k);

    // 中间计时点
    // auto mid = std::chrono::high_resolution_clock::now();
    // 计算中间耗时
    // auto midduration = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start);
    // 输出找到父索引所需时间
    // std::cout << "找到父索引所需时间为 :" << midduration.count() << " 毫秒" << std::endl;
    // 输出顶层分量的大小
    // std::cout << "the topComponentIds size is :" << topComponentIds.size() << std::endl;

    // 用于存储最终结果的节点集合
    std::unordered_set<int> resultNodes;
    // 用于存储最终结果
    std::unordered_set<int> res;
    // 用于广度优先搜索的队列
    std::queue<int> componentQueue;

    // 将顶层分量的ID加入队列
    for (int compId : topComponentIds)
    {
        componentQueue.push(compId);
    }

    // 广度优先搜索所有分量
    while (!componentQueue.empty())
    {
        int currentId = componentQueue.front();
        componentQueue.pop();

        // 如果当前分量已经在结果集中，则跳过
        if (resultNodes.find(currentId) != resultNodes.end())
            continue;

        // 将当前分量添加到结果集中
        resultNodes.insert(currentId);

        // 遍历当前分量的所有子分量
        for (int childCompId : ComponentChildren[currentId])
        {
            // 1.10 目的在于获得最小社区
            // if (visitedComponents.find(childCompId) != visitedComponents.end()) {

            // 获取子分量中的一个节点
            int childNode;
            for (int node : ComponentToNodes[childCompId])
            {
                childNode = node;
                break;
            }

            // 如果子分量的节点的核心度大于等于k，将子分量加入队列
            if (coreMinimumDegree[coreIndex[childNode]] >= k)
            {
                componentQueue.push(childCompId);
            }

            // }

        }
    }

    // 将结果集中的所有连通分量的点添加到最终结果中
    for (int compId : resultNodes)
    {
        for (int node : ComponentToNodes[compId])
        {
            res.insert(node);
        }
    }

    // 结束计时
    // auto end = std::chrono::high_resolution_clock::now();
    // 计算找到子索引所需时间
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid);
    // 输出找到子索引所需时间
    //std::cout << "找到子索引所需时间为 : " << duration.count() << " 毫秒" << std::endl;

    // 返回最终结果
    return res;
}

std::unordered_set<int> TreeIndex::checkResult(std::vector<int> &queryNodes, std::unordered_set<int> &result, Graph &graph) {
    std::vector<int> tree2;
    std::unordered_set<int> querySet(queryNodes.begin(), queryNodes.end());

    std::unordered_map<int, int> d2 = {};
    for (auto& node : graph.getDegrees()) {
        d2[node.first] = 0;
    }
    for (auto& node : result) {
        d2[node] = 1;
    }

    tree2 = graph.isQuerySetConnected(querySet, d2);
    std::cout << "辅助函数已启用 " << tree2.size() << std::endl;
    if (tree2.size() == 0) {
        std::cout << "查询顶点集不连通！" << std::endl;
        // 这里我们增加一个函数，去寻找最大的公共连通分量
        int k = getKfromIndex(queryNodes);
        std::cout << "the k is : " << k << std::endl;
        if (k == -1) {
            k = 1;
        }
        result = findKCoreSubgraph(queryNodes, k);
    }

    return result;
}

void TreeIndex::findSubShells(int node, int wander) {
    int shellIndex = nodeToComponentId[node];
    int wanderIndex = nodeToComponentId[wander];
    int flag = 0;

    std::vector<int> queryShell;
    std::vector<int> visited;
    visited.push_back(shellIndex);
    visited.push_back(wanderIndex);
    for (auto& comp : ComponentChildren[shellIndex]) {
        // std::cout << "the comp is " << comp << std::endl;
        if (comp == wanderIndex) {
            flag = 1;
            break;
        }
        queryShell.push_back(comp);
        visited.push_back(comp);
    }
    if (flag == 1) {
        std::cout << "the ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    } 

    while (!queryShell.empty()) {
        int currentComp = queryShell.back();
        queryShell.pop_back();

        // std::cout << "the currentComp is " << currentComp << std::endl;
        // std::cout << "the childComp is " << std::endl;
        for (int childComp : ComponentChildren[currentComp]) {
            if (childComp == wanderIndex) {
                flag = 1;
                break;
            }
            // std::cout << childComp << " " << std::endl;
            if (std::find(visited.begin(), visited.end(), childComp) == visited.end()) {
                queryShell.push_back(childComp);
                visited.push_back(childComp);
            }
            }
            if (flag == 1) {
                std::cout << "the ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            } 
    }

    if (flag == 1) {
        std::cout << "the ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    }

}

int TreeIndex::getKfromIndex(const std::vector<int>& query) {
    int componentID = -1;
    std::unordered_map<int, query_nodes> children;

    // 不记录已访问的连通分量ID，遍历每个顶点对应连通分量的所有子分量
    for (auto& p : query) {
        std::queue<int> que;
        std::unordered_set<int> visited;
        que.push(nodeToComponentId[p]);
        visited.insert(nodeToComponentId[p]);

        while (!que.empty()) {
            componentID = que.front();
            que.pop();

            children[p].insert(componentID);
            for (auto &childrenSet: getConnectedComponentChildren(componentID)) {
                if (visited.find(childrenSet) == visited.end() && childrenSet!= -1) { // 防止耗时过长
                    visited.insert(childrenSet);
                    que.push(childrenSet);
                }
            }
        }
    }

    // 找到chilren中所有元素的并集
    query_nodes result;  // 并集
    std::unordered_multiset<int> countSum;
    for (auto& child : children) {
        // std::cout << "child: " << child.first << std::endl;
        // std::cout << "child.size: " << child.second.size() << std::endl;
        countSum.insert(child.second.begin(), child.second.end());
        result.insert(child.second.begin(), child.second.end());
    }

    // 找到result中元素个数少于query.size()的元素
    int minSize = query.size();
    // std::cout << "minSize: " << minSize << std::endl;
    for (auto& c : result) {
        if (countSum.count(c) == minSize) {
            std::cout << "countSum: " << c << std::endl;
        }
        if (countSum.count(c) != minSize) {
            result.erase(c);
        }
        // else 
        // std::cout << "countSum: " << coreMinimumDegree[coreIndex[*ComponentToNodes[c].begin()]] << std::endl;
    }
    
    // 找到result中核心度最高的元素
    int maxK = -1;
    for (auto& c : result) {
        // int tmp = coreMinimumDegree[c]; // 这是核心id对应的核心度！
        int tmp = coreMinimumDegree[coreIndex[*ComponentToNodes[c].begin()]];
        std::cout << "coreIndex: " << tmp << std::endl;
        if (tmp > maxK) {
            maxK = tmp;
            std::cout << "c " << c << std::endl;
        }
    }
    // std::cout << "coreIndex: " << maxK << std::endl;
    return maxK;
}

TreeIndex::TreeIndex(Graph &graph, std::string datasetName)
{
    // computeCoreIndex(graph);

    printCoreIndex();

    computeCoreComposition(graph);

    // printConnectedComponentNodes();
    // printConnectedComponentChildren();
    // printConnectedComponentParent();

    // computeNodeNeighbors(graph);
}

void TreeIndex::computeCoreComposition(Graph &graph)
{

    std::unordered_map<int, int> highestGroup;

    std::vector<std::stack<int>> coreGroup(graph.getNumberOfNodes());

    for (int i = 0; i < graph.getNumberOfNodes(); ++i)
    {
        coreGroup[i] = std::stack<int>();
    }

    for (const auto &ci : coreIndex)
    {
        coreGroup[coreGroup.size() - ci.second - 1].push(ci.first);
    }

    int connectedComponentIndex = 0;

    int coreIndex = coreGroup.size() - 1;
    for (std::stack<int> &group : coreGroup)
    {
        if (!group.empty())
        {
            connectedComponentChildren[coreIndex] = std::unordered_map<int, std::unordered_set<int>>();
            connectedComponentNodes[coreIndex] = std::unordered_map<int, std::unordered_set<int>>();
            connectedComponentParent[coreIndex] = std::unordered_map<int, int>();
        }

        while (!group.empty())
        {
            int node = group.top();
            // std::cout<<node<<" endl"<<std::endl;
            group.pop();

            for (int neighbor : graph.getNeighbors(node))
            {
                if (this->coreIndex[node] == this->coreIndex[neighbor])
                {
                    if (!nodeGroup.count(neighbor))
                    {
                        continue;
                    }
                    else if (nodeGroup.find(node) == nodeGroup.end())
                    {
                        // Branch 1
                        connectedComponentNodes[coreIndex][nodeGroup[neighbor]].insert(node);
                        nodeGroup[node] = nodeGroup[neighbor];
                        highestGroup[node] = nodeGroup[neighbor];
                    }
                    else if (nodeGroup[node] != nodeGroup[neighbor])
                    {
                        // Branch 2
                        int oldGroup = nodeGroup[neighbor];
                        connectedComponentNodes[coreIndex][nodeGroup[node]].insert(
                            connectedComponentNodes[coreIndex][oldGroup].begin(),
                            connectedComponentNodes[coreIndex][oldGroup].end());
                        connectedComponentChildren[coreIndex][nodeGroup[node]].insert(
                            connectedComponentChildren[coreIndex][oldGroup].begin(),
                            connectedComponentChildren[coreIndex][oldGroup].end());

                        for (int oldGroupNode : connectedComponentNodes[coreIndex][oldGroup])
                        {
                            nodeGroup[oldGroupNode] = nodeGroup[node];
                        }

                        for (auto &oldHighestNode : highestGroup)
                        {
                            if (oldHighestNode.second == oldGroup)
                            {
                                highestGroup[oldHighestNode.first] = nodeGroup[node];
                            }
                        }

                        for (int newParentGroup : this->connectedComponentChildren[coreIndex][oldGroup])
                        {
                            this->connectedComponentParent[coreIndex + 1][newParentGroup] = nodeGroup[node];
                        }

                        connectedComponentNodes[coreIndex].erase(oldGroup);
                        connectedComponentChildren[coreIndex].erase(oldGroup);
                        connectedComponentParent[coreIndex].erase(oldGroup);
                    }
                }
                else
                {
                    if (nodeGroup.find(neighbor) == nodeGroup.end())
                    {
                        continue;
                    }
                    else if (nodeGroup.find(node) == nodeGroup.end())
                    {
                        if (connectedComponentNodes[coreIndex].find(highestGroup[neighbor]) != connectedComponentNodes[coreIndex].end())
                        {
                            // Branch 3
                            connectedComponentNodes[coreIndex][highestGroup[neighbor]].insert(node);
                            nodeGroup[node] = highestGroup[neighbor];
                            highestGroup[node] = highestGroup[neighbor];
                        }
                        else if (connectedComponentNodes[coreIndex + 1].find(highestGroup[neighbor]) != connectedComponentNodes[coreIndex + 1].end())
                        {
                            // Branch 4
                            connectedComponentNodes[coreIndex][connectedComponentIndex].insert(node);
                            connectedComponentChildren[coreIndex][connectedComponentIndex].insert(highestGroup[neighbor]);
                            nodeGroup[node] = connectedComponentIndex;
                            highestGroup[node] = connectedComponentIndex;
                            connectedComponentParent[coreIndex + 1][highestGroup[neighbor]] = connectedComponentIndex;

                            int oldHighest = highestGroup[neighbor];
                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == oldHighest)
                                {
                                    highestGroup[oldHighestNode.first] = connectedComponentIndex;
                                }
                            }
                            connectedComponentIndex++;
                        }
                        else
                        {
                            // Branch 5

                            int highestLevel = coreIndex + 2;

                            while (connectedComponentNodes[highestLevel].find(highestGroup[neighbor]) == connectedComponentNodes[highestLevel].end())
                            {
                                highestLevel++; // 不断向上查找，直到找到一个包含邻居最高组的级别
                            }

                            while (connectedComponentNodes[coreIndex + 1].find(highestGroup[neighbor]) == connectedComponentNodes[highestLevel].end())
                            {
                                if (connectedComponentNodes.find(highestLevel - 1) == connectedComponentNodes.end())
                                {
                                    connectedComponentNodes[highestLevel - 1] = std::unordered_map<int, std::unordered_set<int>>();
                                }
                                connectedComponentNodes[highestLevel - 1][connectedComponentIndex] = std::unordered_set<int>();
                                connectedComponentChildren[coreIndex][connectedComponentIndex] = std::unordered_set<int>();
                                connectedComponentChildren[coreIndex][connectedComponentIndex].insert(highestGroup[neighbor]);
                                connectedComponentParent[highestLevel][highestGroup[neighbor]] = connectedComponentIndex;

                                // 更新highestGroup映射
                                for (auto &oldHighestNode : highestGroup)
                                {
                                    if (oldHighestNode.second == highestGroup[neighbor])
                                    {
                                        highestGroup[oldHighestNode.first] = connectedComponentIndex;
                                    }
                                }

                                highestLevel--;
                                connectedComponentIndex++;
                            }

                            // 创建最初的连接组
                            connectedComponentNodes[coreIndex][connectedComponentIndex] = std::unordered_set<int>({node});
                            connectedComponentChildren[coreIndex][connectedComponentIndex] = std::unordered_set<int>({highestGroup[neighbor]});
                            nodeGroup[node] = connectedComponentIndex;
                            highestGroup[node] = connectedComponentIndex;
                            connectedComponentParent[coreIndex + 1][highestGroup[neighbor]] = connectedComponentIndex;

                            // 更新highestGroup映射
                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == highestGroup[neighbor])
                                {
                                    highestGroup[oldHighestNode.first] = connectedComponentIndex;
                                }
                            }
                            connectedComponentIndex++;
                        }
                    }
                    else if (nodeGroup[node] != highestGroup[neighbor])
                    {
                        if (connectedComponentNodes[coreIndex].find(highestGroup[neighbor]) != connectedComponentNodes[coreIndex].end() && nodeGroup[node] != highestGroup[neighbor])
                        {
                            // Branch 6
                            int oldGroup = highestGroup[neighbor];
                            connectedComponentNodes[coreIndex][nodeGroup[node]].insert(
                                connectedComponentNodes[coreIndex][oldGroup].begin(),
                                connectedComponentNodes[coreIndex][oldGroup].end());
                            connectedComponentChildren[coreIndex][nodeGroup[node]].insert(
                                connectedComponentChildren[coreIndex][oldGroup].begin(),
                                connectedComponentChildren[coreIndex][oldGroup].end());

                            for (int oldGroupNode : connectedComponentNodes[coreIndex][oldGroup])
                            {
                                nodeGroup[oldGroupNode] = nodeGroup[node];
                            }

                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == oldGroup)
                                {
                                    highestGroup[oldHighestNode.first] = nodeGroup[node];
                                }
                            }

                            for (int newParentGroup : this->connectedComponentChildren[coreIndex][oldGroup])
                            {
                                connectedComponentParent[coreIndex + 1][newParentGroup] = nodeGroup[node];
                            }

                            connectedComponentNodes[coreIndex].erase(oldGroup);
                            connectedComponentChildren[coreIndex].erase(oldGroup);
                            connectedComponentParent[coreIndex].erase(oldGroup);
                        }
                        else if (connectedComponentNodes[coreIndex + 1].find(highestGroup[neighbor]) != connectedComponentNodes[coreIndex + 1].end())
                        {
                            connectedComponentChildren[coreIndex][nodeGroup[node]].insert(highestGroup[neighbor]);
                            connectedComponentParent[coreIndex + 1][highestGroup[neighbor]] = nodeGroup[node];

                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == highestGroup[neighbor])
                                {
                                    highestGroup[oldHighestNode.first] = nodeGroup[node];
                                }
                            }
                        }
                        else
                        { // branch 8
                            int highestLevel = coreIndex + 2;
                            while (connectedComponentNodes[highestLevel].find(highestGroup[neighbor]) == connectedComponentNodes[highestLevel].end())
                            {
                                highestLevel++;
                            }

                            while (connectedComponentNodes[coreIndex + 1].find(highestGroup[neighbor]) == connectedComponentNodes[coreIndex + 1].end())
                            {
                                connectedComponentNodes[highestLevel - 1][connectedComponentIndex] = std::unordered_set<int>();
                                connectedComponentChildren[coreIndex][connectedComponentIndex] = std::unordered_set<int>();
                                connectedComponentChildren[coreIndex][connectedComponentIndex].insert(highestGroup[neighbor]);
                                connectedComponentParent[highestLevel][highestGroup[neighbor]] = connectedComponentIndex;

                                for (auto &oldHighestNode : highestGroup)
                                {
                                    if (oldHighestNode.second == highestGroup[neighbor])
                                    {
                                        highestGroup[oldHighestNode.first] = connectedComponentIndex;
                                    }
                                }

                                highestLevel--;
                                connectedComponentIndex++;
                            }

                            connectedComponentChildren[coreIndex][nodeGroup[node]].insert(highestGroup[neighbor]);
                            connectedComponentParent[coreIndex + 1][highestGroup[node]] = nodeGroup[node];

                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == highestGroup[neighbor])
                                {
                                    highestGroup[oldHighestNode.first] = nodeGroup[node];
                                }
                            }
                        }
                    }
                }
            }

            // If it has no neighbors
            if (nodeGroup.find(node) == nodeGroup.end())
            {
                connectedComponentNodes[coreIndex][connectedComponentIndex] = std::unordered_set<int>{node};
                connectedComponentChildren[coreIndex][connectedComponentIndex] = std::unordered_set<int>();
                nodeGroup[node] = connectedComponentIndex;
                highestGroup[node] = connectedComponentIndex;
                connectedComponentIndex++;
            }
        }

        coreIndex--;
    }
}

void TreeIndex::computeCoreCompositionByLcy(Graph &graph)
{
    //
}

void printHighestGroup(const std::unordered_map<int, int> &highestGroup)
{
    for (const auto &pair : highestGroup)
    {
        std::cout << "Node: " << pair.first << ", Group: " << pair.second << std::endl;
    }
}

void TreeIndex::computeNodeNeighbors(Graph &graph)
{
    for (auto &node : graph.getNodes())
    {
        for (int neighbor : graph.getNeighbors(node))
        {
            int shellIndex = coreIndex[neighbor];
            nodeNeighbors[node][shellIndex].insert(neighbor);
        }
    }
}

int TreeIndex::getMinimumCoreIndex(std::vector<int> queryNodes)
{
    std::unordered_set<int> queryIndexes;
    for (int queryNode : queryNodes)
    {
        queryIndexes.insert(coreIndex[queryNode]);
    }
    int k = *std::max_element(queryIndexes.begin(), queryIndexes.end());

    std::unordered_set<int> Q(queryNodes.begin(), queryNodes.end());
    std::unordered_set<int> C;
    for (int queryNode : queryNodes)
    {
        if (coreIndex[queryNode] == k)
        {
            C.insert(nodeGroup[queryNode]);
            Q.erase(queryNode);
        }
    }

    while (!Q.empty() || C.size() != 1)
    {
        std::unordered_set<int> cFirst;
        for (int connectedComponentID : C)
        {
            cFirst.insert(connectedComponentParent[k][connectedComponentID]);
        }

        for (int queryNode : queryNodes)
        {
            if (coreIndex[queryNode] == k - 1)
            {
                cFirst.insert(nodeGroup[queryNode]);
                Q.erase(queryNode);
            }
        }

        C = cFirst;
        --k;
    }

    return k;
}

std::unordered_set<int> TreeIndex::getNeighbors(int node, int coreIndex)
{
    std::unordered_set<int> neighbors;

    for (auto &[shellIndex, nodes] : nodeNeighbors[node])
    {
        if (shellIndex >= coreIndex)
        {
            neighbors.insert(nodes.begin(), nodes.end());
        }
    }

    return neighbors;
}

std::unordered_set<int> TreeIndex::getNeighbors(int node, int coreIndex, std::unordered_set<int> &subcore)
{
    std::unordered_set<int> neighbors;
    for (auto &[shellIndex, nodes] : nodeNeighbors[node])
    {
        if (shellIndex >= coreIndex)
        {
            for (int x : nodes)
            {
                if (subcore.count(x))
                {
                    neighbors.insert(x);
                }
            }
        }
    }

    return neighbors;
}

int TreeIndex::getCoreMinimumDegree(int coreIndex)
{
    return coreMinimumDegree[coreIndex];
}

// 获取 map 中所有的键
std::unordered_set<int> getKeySet(const std::unordered_map<int, std::unordered_set<int>> &map)
{
    std::unordered_set<int> keys;
    for (const auto &pair : map)
    {
        keys.insert(pair.first);
    }
    return keys;
}

// 获取指定核心中的节点数量
int TreeIndex::getNumberOfNodes(int coreIndex)
{
    int currentCoreIndex = coreIndex;
    std::unordered_set<int> currentConnectedComponents = getKeySet(connectedComponentNodes[currentCoreIndex]);

    int numberOfNodes = 0;
    while (!currentConnectedComponents.empty())
    {
        std::unordered_set<int> childrenConnectedComponents;
        for (int connectedComponent : currentConnectedComponents)
        {
            numberOfNodes += connectedComponentNodes[currentCoreIndex][connectedComponent].size();
            if (connectedComponentChildren[currentCoreIndex].count(connectedComponent))
            {
                for (int child : connectedComponentChildren[currentCoreIndex][connectedComponent])
                {
                    childrenConnectedComponents.insert(child);
                }
            }
        }

        currentConnectedComponents = childrenConnectedComponents;
        ++currentCoreIndex;
    }

    return numberOfNodes;
}

std::unordered_set<int> TreeIndex::getCore(std::vector<int> queryNodes)
{
    std::unordered_set<int> core;
    std::unordered_set<int> queryIndexes;
    for (int queryNode : queryNodes)
    {
        queryIndexes.insert(coreIndex[queryNode]);
    }
    int k = *std::max_element(queryIndexes.begin(), queryIndexes.end());

    std::unordered_set<int> Q(queryNodes.begin(), queryNodes.end());
    std::unordered_set<int> C;
    for (int queryNode : queryNodes)
    {
        if (coreIndex[queryNode] == k)
        {
            int cc = nodeGroup[queryNode];
            C.insert(cc);
            Q.erase(queryNode);
        }
    }

    while (!Q.empty() || C.size() != 1)
    {
        std::unordered_set<int> cFirst;
        for (int connectedComponentID : C)
        {
            int parentID = connectedComponentParent[k][connectedComponentID];
            cFirst.insert(parentID);
        }

        for (int queryNode : queryNodes)
        {
            if (coreIndex[queryNode] == k - 1)
            {
                int cc = nodeGroup[queryNode];
                cFirst.insert(cc);
                Q.erase(queryNode);
            }
        }

        C = cFirst;
        --k;
    }

    while (connectedComponentChildren.count(k))
    {
        std::unordered_set<int> newC;
        for (int cc : C)
        {
            newC.insert(cc);
            if (connectedComponentChildren[k].count(cc))
            {
                for (int ccc : connectedComponentChildren[k][cc])
                {
                    newC.insert(ccc);
                }
            }
            if (connectedComponentNodes[k].count(cc))
            {
                core.insert(connectedComponentNodes[k][cc].begin(), connectedComponentNodes[k][cc].end());
            }
        }
        C = newC;
        ++k;
    }

    for (int queryNode : queryNodes)
    {
        core.insert(queryNode);
    }

    return core;
}

int TreeIndex::getCoreIndex(int node)
{
    int index = coreIndex[node];
    return coreMinimumDegree[index];
}

// 得到  点  对应的  连通分量ID
int TreeIndex::getComponent(int node)
{
    auto it = nodeGroup.find(node);
    if (it != nodeGroup.end())
    {
        return it->second;
    }
    else
    {
        return -1; // Or throw an exception or handle error as appropriate
    }
}

// Get the parent component ID for a given component
int TreeIndex::getParentComponent(int componentId)
{
    for (const auto &level : connectedComponentParent)
    {
        auto it = level.second.find(componentId);
        if (it != level.second.end())
        {
            return it->second;
        }
    }
    return 0; // Or handle the case where the component has no parent
}

std::unordered_set<int> TreeIndex::getConnectedComponentChildren(int componentId)
{
    std::unordered_set<int> childrenComponents;

    // 遍历所有核心层
    for (const auto &layer : connectedComponentChildren)
    {
        auto it = layer.second.find(componentId);
        if (it != layer.second.end())
        {
            // 如果找到了给定连通分量ID，加入其所有子连通分量到结果集
            childrenComponents.insert(it->second.begin(), it->second.end());
        }
    }

    return childrenComponents;
}

// Get all nodes in a specified component
std::unordered_set<int> TreeIndex::getNodesInComponent(int componentId)
{
    for (const auto &level : connectedComponentNodes)
    {
        auto it = level.second.find(componentId);
        if (it != level.second.end())
        {
            return it->second;
        }
    }
    return std::unordered_set<int>(); // Return empty if not found
}

// 打印核心指标
void TreeIndex::printCoreIndex()
{
    std::cout << "Core Index:" << std::endl;
    for (const auto &pair : coreIndex)
    {
        std::cout << "Node: " << pair.first << ", Core: " << pair.second << std::endl;
    }
}

// 打印每个核心的最小度数
void TreeIndex::printCoreMinimumDegree()
{
    std::cout << "Core Minimum Degree:" << std::endl;
    for (const auto &pair : coreMinimumDegree)
    {
        std::cout << "Core: " << pair.first << ", Minimum Degree: " << pair.second << std::endl;
    }
}

// 打印连通分量的父子关系
void TreeIndex::printConnectedComponentParent()
{
    std::cout << "Connected Component Parent:" << std::endl;
    for (const auto &outer_pair : connectedComponentParent)
    {
        int parent = outer_pair.first;
        std::cout << "Core Index: " << parent << std::endl;
        const auto &children = outer_pair.second;
        for (const auto &inner_pair : children)
        {
            int child = inner_pair.first;
            std::cout << "    Component ID : " << child << ", Parent: " << inner_pair.second << std::endl;
        }
    }
}

// 打印连通分量的孩子节点集合
void TreeIndex::printConnectedComponentChildren()
{
    std::cout << "Connected Component Children:" << std::endl;
    for (const auto &outer_pair : connectedComponentChildren)
    {
        int parent = outer_pair.first;
        std::cout << "Core Index: " << parent << std::endl;
        const auto &children = outer_pair.second;

        for (const auto &inner_pair : children)
        {
            int child = inner_pair.first;
            const auto &child_set = inner_pair.second;
            std::cout << "    Component ID: " << child << ", Children: ";
            for (int c : child_set)
            {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }
    }
}

// 打印连通分量的节点集合
void TreeIndex::printConnectedComponentNodes()
{
    std::cout << "Connected Component Nodes:" << std::endl;
    for (const auto &outer_pair : connectedComponentNodes)
    {
        int connectedComponent = outer_pair.first;
        std::cout << "Core Index: " << connectedComponent << std::endl;
        const auto &nodes = outer_pair.second;
        for (const auto &inner_pair : nodes)
        {
            if (inner_pair.second.size() == 0)
            {
                continue;
            }
            int coreIndex = inner_pair.first;
            const auto &node_set = inner_pair.second;
            std::cout << " Component ID:: " << coreIndex << ", Nodes: ";
            for (int node : node_set)
            {
                std::cout << node << " ";
            }
            std::cout << std::endl;
        }
    }
}

// 打印节点分组
void TreeIndex::printNodeGroup()
{
    std::cout << "Node Group:" << std::endl;
    for (const auto &pair : nodeGroup)
    {
        std::cout << "Node: " << pair.first << ", Group: " << pair.second << std::endl;
    }
}

// 打印节点邻居信息
void TreeIndex::printNodeNeighbors()
{
    std::cout << "Node Neighbors:" << std::endl;
    for (const auto &outer_pair : nodeNeighbors)
    {
        int node = outer_pair.first;
        std::cout << "Node: " << node << std::endl;
        const auto &core_neighbors = outer_pair.second;
        for (const auto &inner_pair : core_neighbors)
        {
            int coreIndex = inner_pair.first;
            const auto &neighbor_set = inner_pair.second;
            std::cout << "    Core Index: " << coreIndex << ", Neighbors: ";
            for (int neighbor : neighbor_set)
            {
                std::cout << neighbor << " ";
            }
            std::cout << std::endl;
        }
    }
}

void TreeIndex::printQueryNodesCoreIndex(const std::vector<int> &queryNodes)
{
    std::cout << "Core index for query nodes:" << std::endl;
    for (int node : queryNodes)
    {
        std::cout << "Node " << node << ": Core Index = " << coreMinimumDegree.at(coreIndex.at(node)) << " Component is " << nodeToComponentId.at(node) << std::endl;
    }
}

// 获得连通分量对应的核心度
int TreeIndex::getCoreFromComponent(int componentId) {
    for (const auto &level : connectedComponentNodes) {
        auto it = level.second.find(componentId);
        if (it != level.second.end()) {
            // 如果找到了对应的连通分量ID，返回其对应的核心度
            return coreMinimumDegree.at(level.first);
        }
    }
    // 如果没有找到对应的连通分量，返回-1或抛出异常
    return -1; // 表示未找到
}

std::unordered_set<int> TreeIndex::greedyStep(std::vector<int>& queryNodes, int k) {
    std::unordered_set<int> H_min_star(queryNodes.begin(), queryNodes.end());
    std::unordered_set<int> visited;
    // std::priority_queue<std::pair<std::pair<int, int>, int>, 
    //                     std::vector<std::pair<std::pair<int, int>, int>>, 
    //                     std::greater<>> pq;
    // 定义优先队列
    std::priority_queue<Item, std::vector<Item>, CompareItem> pq;
    std::unordered_map<int, std::tuple<int, int, int>> pq_sustain; // 维护优先队列
    std::unordered_map<int, int> min_degree;

    // 初始化
    for (int node : queryNodes) {
        // pq.push({{0, 0}, node}); // 初始优先级为 (0, 0)
        visited.insert(node);
        min_degree[node] = 0;
    }
    for (int i = 0; i < count_port; i++) {
        pq_sustain[i] = {0, 0, 0};
    }

    // 最小度数 mu_star
    int mu_star = INT8_MAX;

    // 初始化并查集
    UnionFind uf(count_port);

    // 将查询节点加入并查集 计算当前最小度数 mu_star 维护度数列表
    for (int node : queryNodes) {
        uf.find(node); // 确保查询节点在并查集中
        // 更新并查集
        for (int neighbor : getNeighbors(node)) {
            if (H_min_star.find(neighbor) != H_min_star.end()) {
                uf.unite(node, neighbor);
                min_degree[node]++;
            }
        }
        mu_star = std::min(mu_star, min_degree[node]);
        std::cout << node << "------------uf.find(node)----------: " << uf.find(node) << std::endl;
        std::cout << "chushi----------: " << mu_star << std::endl;
    }

    for (int node : queryNodes) {
        // 遍历当前节点的邻居
        for (int neighbor : getNeighbors(node)) {

                // 计算连接分数 p'(u)
                std::unordered_set<int> components;
                for (int neighbor_neighbor : getNeighbors(neighbor)) {
                    if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                        components.insert(uf.find(neighbor_neighbor));
                    }
                }
                int connection_score = components.size() - 1;  //max(0,x)
                connection_score = std::max(0, connection_score);

                // 计算最小度数分数 p''(u)
                int degree = 0;  // 当前节点的度数
                for (int neighbor_neighbor : getNeighbors(neighbor)) {
                    if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                        degree++;
                    }
                }

                // 增益效应：统计当前节点的邻居中，哪些邻居的加入可以使当前节点的度数更接近 k
                int degree_gain = 0;
                for (int neighbor_neighbor : getNeighbors(neighbor)) {
                    if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                        int neighbor_degree = 0;
                        // for (int nn : getNeighbors(neighbor_neighbor)) {
                        //     if (H_min_star.find(nn) != H_min_star.end()) {
                        //         neighbor_degree++;
                        //     }
                        // }
                        neighbor_degree = min_degree[neighbor_neighbor];
                        if (neighbor_degree < k) {
                            degree_gain++;
                        }
                    }
                }

                // 惩罚效应：计算当前节点需要多少额外的邻居才能达到目标最小度数 k
                int degree_penalty = std::max(0, k - degree);
                // degree_penalty = degree_penalty / 2;

                // 最小度数分数 p''(u) = 增益效应 - 惩罚效应
                int degree_score = std::max(0, degree_gain - degree_penalty);

                // 综合分数 p(u) = (p'(u), p''(u))
                pq.push({{connection_score, degree_gain, degree_penalty}, neighbor});
                pq_sustain[neighbor] = {connection_score, degree_gain, degree_penalty};
            
        }
    }

    while (!pq.empty()) {

        auto top = pq.top();
        pq.pop();
        int node = top.value;

        // 跳过已经访问过的节点
        if (visited.find(node) != visited.end()) {
            continue;
        }
        if (getCoreIndex(node) < k) {
            visited.insert(node);
            continue;
        }

        int connection_score = std::get<0>(top.priority);
        int degree_gain = std::get<1>(top.priority);
        int degree_penalty = std::get<2>(top.priority);

        // 避免重复
        if (std::get<0>(pq_sustain[node]) != connection_score || 
            std::get<1>(pq_sustain[node]) != degree_gain ||
            std::get<2>(pq_sustain[node]) != degree_penalty ) {
            continue;
        }

        H_min_star.insert(node);
        if (connection_score == 1) {
            std::cout << "charu: " << node << "  " << connection_score << "  " << degree_gain << "  " << degree_penalty << std::endl;
        }

        visited.insert(node);
        min_degree[node] = 0;

        // 更新并查集 维护度数列表
        std::vector<int> tmp;
        tmp.push_back(node);
        for (int neighbor : getNeighbors(node)) {
            if (H_min_star.find(neighbor) != H_min_star.end()) {
                uf.unite(node, neighbor);
                min_degree[node]++;
                min_degree[neighbor]++;
                tmp.push_back(neighbor);
            }
            // std::cout << neighbor << "  min_degree[neighbor]: " << min_degree[neighbor] << std::endl;
        }
        // std::cout << node << "  min_degree[node]: " << min_degree[node] << std::endl;
        mu_star = INT8_MAX;
        for (auto& x : H_min_star) {
            mu_star = std::min(mu_star, min_degree[x]);
        }

        // 判断是否跳出循环
        int count = count_port - uf.getSetNum() + 1;
        // std::cout << "count: " << count << "         ";
        // std::cout << "H_min_star.size(): " << H_min_star.size() << std::endl;
        // std::cout << "mu_star: " << mu_star << "         ";
        // std::cout << "k: " << k << std::endl;
        if (H_min_star.size() == count && mu_star >= k) {
        // if (H_min_star.size() == count) {    
            break;
        }

        // 遍历当前节点的邻居 以及必要节点的邻居
        for (auto& node : tmp) {
            for (int neighbor : getNeighbors(node)) {
                if (visited.find(neighbor) == visited.end()) {
    
                    // 计算连接分数 p'(u)
                    std::unordered_set<int> components;
                    for (int neighbor_neighbor : getNeighbors(neighbor)) {
                        if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                            components.insert(uf.find(neighbor_neighbor));
                            // std::cout << "uf: " << uf.find(neighbor_neighbor) << std::endl;
                        }
                    }
                    int connection_score = components.size() - 1;
    
                    // 计算最小度数分数 p''(u)
                    int degree = 0;  // 当前节点的度数
                    for (int neighbor_neighbor : getNeighbors(neighbor)) {
                        if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                            degree++;
                        }
                    }
    
                    // 增益效应：统计当前节点的邻居中，哪些邻居的加入可以使当前节点的度数更接近 k
                    int degree_gain = 0;
                    for (int neighbor_neighbor : getNeighbors(neighbor)) {
                        if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                            int neighbor_degree = 0;
                            // for (int nn : getNeighbors(neighbor_neighbor)) {
                            //     if (H_min_star.find(nn) != H_min_star.end()) {
                            //         neighbor_degree++;
                            //     }
                            // }
                            neighbor_degree = min_degree[neighbor_neighbor];
                            if (neighbor_degree < k) {
                                degree_gain++;
                            }
                        }
                    }
    
                    // 惩罚效应：计算当前节点需要多少额外的邻居才能达到目标最小度数 k
                    int degree_penalty = std::max(0, k - degree);
                    // degree_penalty = degree_penalty / 2;
                    // std::cout << "degree_gain: " << degree_gain << std::endl;
                    // std::cout << "degree_penalty: " << degree_penalty << std::endl;
    
                    // 最小度数分数 p''(u) = 增益效应 - 惩罚效应
                    int degree_score = std::max(0, degree_gain - degree_penalty);
                    // std::cout << "degree_score: " << degree_score << std::endl;
    
                    // 综合分数 p(u) = (p'(u), p''(u))
                    pq.push({{connection_score, degree_gain, degree_penalty}, neighbor});
                    pq_sustain[neighbor] = {connection_score, degree_gain, degree_penalty};
                }
            }
        }
    }

    return H_min_star;
}

std::unordered_set<int> TreeIndex::greedyStep_simply(std::vector<int>& queryNodes, int k, std::unordered_set<int>& realnodes) {
    std::unordered_set<int> H_min_star(queryNodes.begin(), queryNodes.end());
    std::unordered_set<int> visited;
    // std::priority_queue<std::pair<std::pair<int, int>, int>, 
    //                     std::vector<std::pair<std::pair<int, int>, int>>, 
    //                     std::greater<>> pq;
    // 定义优先队列
    std::priority_queue<Item, std::vector<Item>, CompareItem> pq;
    std::unordered_map<int, std::tuple<int, int, int>> pq_sustain; // 维护优先队列
    std::unordered_map<int, int> min_degree;

    // 初始化
    for (int node : queryNodes) {
        // pq.push({{0, 0}, node}); // 初始优先级为 (0, 0)
        visited.insert(node);
        min_degree[node] = 0;
    }
    for (int i = 0; i < count_port; i++) {
        pq_sustain[i] = {0, 0, 0};
    }

    // 最小度数 mu_star
    int mu_star = INT8_MAX;

    // 初始化并查集
    UnionFind uf(count_port);

    // 将查询节点加入并查集 计算当前最小度数 mu_star 维护度数列表
    for (int node : queryNodes) {
        uf.find(node); // 确保查询节点在并查集中
        // 更新并查集
        for (int neighbor : getNeighbors(node)) {
            if (H_min_star.find(neighbor) != H_min_star.end()) {
                uf.unite(node, neighbor);
                min_degree[node]++;
            }
        }
        mu_star = std::min(mu_star, min_degree[node]);
        std::cout << node << "------------uf.find(node)----------: " << uf.find(node) << std::endl;
        std::cout << "最小度----------: " << mu_star << std::endl;
    }

    for (int node : queryNodes) {
        // 遍历当前节点的邻居
        for (int neighbor : getNeighbors(node)) {
            if (realnodes.find(neighbor) != realnodes.end()) {
                // 计算连接分数 p'(u)
                std::unordered_set<int> components;
                for (int neighbor_neighbor : getNeighbors(neighbor)) {
                    if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                        components.insert(uf.find(neighbor_neighbor));
                    }
                }
                int connection_score = components.size() - 1;  //max(0,x)
                connection_score = std::max(0, connection_score);

                // 计算最小度数分数 p''(u)
                int degree = 0;  // 当前节点的度数
                for (int neighbor_neighbor : getNeighbors(neighbor)) {
                    if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                        degree++;
                    }
                }

                // 增益效应：统计当前节点的邻居中，哪些邻居的加入可以使当前节点的度数更接近 k
                int degree_gain = 0;
                for (int neighbor_neighbor : getNeighbors(neighbor)) {
                    if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                        int neighbor_degree = 0;
                        // for (int nn : getNeighbors(neighbor_neighbor)) {
                        //     if (H_min_star.find(nn) != H_min_star.end()) {
                        //         neighbor_degree++;
                        //     }
                        // }
                        neighbor_degree = min_degree[neighbor_neighbor];
                        if (neighbor_degree < k) {
                            degree_gain++;
                        }
                    }
                }

                // 惩罚效应：计算当前节点需要多少额外的邻居才能达到目标最小度数 k
                int degree_penalty = std::max(0, k - degree);
                // degree_penalty = degree_penalty / 2;

                // 最小度数分数 p''(u) = 增益效应 - 惩罚效应
                int degree_score = std::max(0, degree_gain - degree_penalty);

                // 综合分数 p(u) = (p'(u), p''(u))
                pq.push({{connection_score, degree_gain, degree_penalty}, neighbor});
                pq_sustain[neighbor] = {connection_score, degree_gain, degree_penalty};
            } 
        }
    }

    while (!pq.empty()) {

        auto top = pq.top();
        pq.pop();
        int node = top.value;

        // 跳过已经访问过的节点
        if (visited.find(node) != visited.end()) {
            continue;
        }
        if (getCoreIndex(node) < k) {
            visited.insert(node);
            continue;
        }

        int connection_score = std::get<0>(top.priority);
        int degree_gain = std::get<1>(top.priority);
        int degree_penalty = std::get<2>(top.priority);

        // 避免重复
        if (std::get<0>(pq_sustain[node]) != connection_score || 
            std::get<1>(pq_sustain[node]) != degree_gain ||
            std::get<2>(pq_sustain[node]) != degree_penalty ) {
            continue;
        }

        H_min_star.insert(node);
        if (connection_score == 1) {
            std::cout << "charu: " << node << "  " << connection_score << "  " << degree_gain << "  " << degree_penalty << std::endl;
        }

        visited.insert(node);
        min_degree[node] = 0;

        // 更新并查集 维护度数列表
        std::vector<int> tmp;
        tmp.push_back(node);
        for (int neighbor : getNeighbors(node)) {
            if (realnodes.find(neighbor) != realnodes.end()) {
                if (H_min_star.find(neighbor) != H_min_star.end()) {
                    uf.unite(node, neighbor);
                    min_degree[node]++;
                    min_degree[neighbor]++;
                    tmp.push_back(neighbor);
                }
                // std::cout << neighbor << "  min_degree[neighbor]: " << min_degree[neighbor] << std::endl;
            }
        }
        // std::cout << node << "  min_degree[node]: " << min_degree[node] << std::endl;
        mu_star = INT8_MAX;
        for (auto& x : H_min_star) {
            mu_star = std::min(mu_star, min_degree[x]);
        }

        // 判断是否跳出循环
        int count = count_port - uf.getSetNum() + 1;
        // std::cout << "count: " << count << "         ";
        // std::cout << "H_min_star.size(): " << H_min_star.size() << std::endl;
        // std::cout << "mu_star: " << mu_star << "         ";
        // std::cout << "k: " << k << std::endl;
        if (H_min_star.size() == count && mu_star >= k) {
        // if (H_min_star.size() == count) {    
            break;
        }

        // 遍历当前节点的邻居 以及必要节点的邻居
        for (auto& node : tmp) {
            for (int neighbor : getNeighbors(node)) {
                if (realnodes.find(neighbor) != realnodes.end()) {
                    if (visited.find(neighbor) == visited.end()) {
    
                        // 计算连接分数 p'(u)
                        std::unordered_set<int> components;
                        for (int neighbor_neighbor : getNeighbors(neighbor)) {
                            if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                                components.insert(uf.find(neighbor_neighbor));
                                // std::cout << "uf: " << uf.find(neighbor_neighbor) << std::endl;
                            }
                        }
                        int connection_score = components.size() - 1;
        
                        // 计算最小度数分数 p''(u)
                        int degree = 0;  // 当前节点的度数
                        for (int neighbor_neighbor : getNeighbors(neighbor)) {
                            if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                                degree++;
                            }
                        }
        
                        // 增益效应：统计当前节点的邻居中，哪些邻居的加入可以使当前节点的度数更接近 k
                        int degree_gain = 0;
                        for (int neighbor_neighbor : getNeighbors(neighbor)) {
                            if (H_min_star.find(neighbor_neighbor) != H_min_star.end()) {
                                int neighbor_degree = 0;
                                // for (int nn : getNeighbors(neighbor_neighbor)) {
                                //     if (H_min_star.find(nn) != H_min_star.end()) {
                                //         neighbor_degree++;
                                //     }
                                // }
                                neighbor_degree = min_degree[neighbor_neighbor];
                                if (neighbor_degree < k) {
                                    degree_gain++;
                                }
                            }
                        }
        
                        // 惩罚效应：计算当前节点需要多少额外的邻居才能达到目标最小度数 k
                        int degree_penalty = std::max(0, k - degree);
                        // degree_penalty = degree_penalty / 2;
                        // std::cout << "degree_gain: " << degree_gain << std::endl;
                        // std::cout << "degree_penalty: " << degree_penalty << std::endl;
        
                        // 最小度数分数 p''(u) = 增益效应 - 惩罚效应
                        int degree_score = std::max(0, degree_gain - degree_penalty);
                        // std::cout << "degree_score: " << degree_score << std::endl;
        
                        // 综合分数 p(u) = (p'(u), p''(u))
                        pq.push({{connection_score, degree_gain, degree_penalty}, neighbor});
                        pq_sustain[neighbor] = {connection_score, degree_gain, degree_penalty};
                    }
                }
            }
            
        }
    }

    return H_min_star;
}

// 斯坦纳树的近似算法，prim
std::vector<int> TreeIndex::steinerTree(std::unordered_set<int>& H_min_star, const std::vector<int>& terminals) {
    int n = count_port; // 图的节点数
    std::vector<bool> isTerminal(n, false); // 标记是否为终端节点
    for (int v : terminals) isTerminal[v] = true;

    // 最小生成树（MST）的 Prim 算法
    std::vector<int> mstParent(n, -1); // MST 中的父节点
    std::vector<bool> inMST(n, false); // 标记是否在 MST 中
    std::vector<int> key(n, INT_MAX); // 每个节点的键值（最小边权重）
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq; // 优先队列

    // 从任意一个终端节点开始
    pq.push({0, terminals[0]});
    key[terminals[0]] = 0;

    while (!pq.empty()) {
        int u = pq.top().second; // 当前节点
        pq.pop();

        if (inMST[u]) continue; // 如果已经在 MST 中，跳过
        inMST[u] = true;

        // 遍历当前节点的所有邻居
        for (int v : getNeighbors(u)) {
            if (H_min_star.find(v) != H_min_star.end()) {
                if (!inMST[v] && key[v] > 1) { // 权重为1（无权图）
                    key[v] = 1; // 更新键值
                    pq.push({key[v], v}); // 将邻居加入优先队列
                    mstParent[v] = u; // 设置父节点
                }
            }
        }
    }

    // 提取包含终端节点的子树
    std::vector<int> steinerTree;
    for (int v : terminals) {
        while (v != -1) { // 从终端节点向上追溯到根
            steinerTree.push_back(v);
            v = mstParent[v];
        }
    }
    sort(steinerTree.begin(), steinerTree.end()); // 排序
    steinerTree.erase(unique(steinerTree.begin(), steinerTree.end()), steinerTree.end()); // 去重
    return steinerTree;
}

std::unordered_set<int> TreeIndex::connectionStep(std::unordered_set<int>& H_min_star, std::vector<int>& queryNodes, int k) {
    std::unordered_set<int> result;
    std::vector<int> r1 = steinerTree(H_min_star, queryNodes);
    std::unordered_set<int> r2(r1.begin(), r1.end());
    // result = r2;
    result = greedyStep_simply(r1, k, H_min_star);
    if (result.size() )
    return result;
}

std::unordered_set<int> TreeIndex::greedyConnection(std::vector<int>& queryNodes, int k) {
    // Step 1: Greedy Step
    std::unordered_set<int> H_min_star = greedyStep(queryNodes, k);

    // Step 2: Connection Step
    std::unordered_set<int> result = connectionStep(H_min_star, queryNodes, k);
    // std::unordered_set<int> result = H_min_star;
    return result;
}