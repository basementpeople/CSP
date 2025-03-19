#include "TreeIndex.h"
struct Item {
    std::tuple<int, int, int> priority; // 优先级 {22, 11, 5}
    int value;                          // 值 1

    // 构造函数
    Item(std::tuple<int, int, int> p, int v) : priority(p), value(v) {}
};

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

// 1 构造函数和析构函数
TreeIndex::TreeIndex(Graph &graph) {
    
    // 与Graph相同的数据成员，有必要继承吗
    // 图的节点数，图的邻接表，图中节点的度
    n = graph.getN(); 
    adj = graph.getAdj();
    degrees = graph.getDegrees();

    // 计算shell_count，coreIndex，coreMinimumDegree
    computeCoreIndex(graph);

    // 划分了每个shell中的不同连通分量
    identifyAndStoreComponents(graph);

    // 构建父子连通分量关系:一个子可能有多个不同的父,不同的父可能有相同的子
    buildParentChildRelationships(graph);

}

// 2 shell索引解决CSP
std::unordered_set<int> TreeIndex::shellsearch(std::vector<int> &queryNodes) {

    // 遍历查询节点，确定最小核心索引
    // 这步最重要
    int k = coreMinimumDegree[coreIndex[queryNodes[0]]];
    if (queryNodes.size() != 1) {
        k = findCommenK_2(queryNodes);
    }
    std::cout << "the k is @@@@@@@@@: " << k << std::endl;

    // 用于存储所有顶层分量的ID
    std::unordered_set<int> topComponentIds;
    // 找到所有父分量为-1的顶层分量
    topComponentIds = findTopComponents(queryNodes, k);

    // 用于存储最终结果的节点集合
    std::unordered_set<int> resultNodes;
    // 用于存储最终结果
    std::unordered_set<int> res;
    // 用于广度优先搜索的队列
    std::queue<int> componentQueue;

    // 将顶层分量的ID加入队列
    for (int compId : topComponentIds) {
        componentQueue.push(compId);
    }

    // 广度优先搜索所有分量
    while (!componentQueue.empty()) {
        int currentId = componentQueue.front();
        componentQueue.pop();

        // 如果当前分量已经在结果集中，则跳过
        if (resultNodes.find(currentId) != resultNodes.end())
            continue;

        // 将当前分量添加到结果集中
        resultNodes.insert(currentId);

        // 遍历当前分量的所有子分量
        for (int childCompId : ComponentChildren[currentId]) {
            // 获取子分量中的一个节点
            int childNode;
            for (int node : ComponentToNodes[childCompId]) {
                childNode = node;
                break;
            }

            // 如果子分量的节点的核心度大于等于k，将子分量加入队列
            if (coreMinimumDegree[coreIndex[childNode]] >= k) {
                componentQueue.push(childCompId);
            }
        }
    }

    // 将结果集中的所有连通分量的点添加到最终结果中
    for (int compId : resultNodes) {
        for (int node : ComponentToNodes[compId]) {
            res.insert(node);
        }
    }

    // 返回最终结果
    return res;
}


std::unordered_set<int> TreeIndex::getNeighbors(int node) {
    // 获取并返回邻居节点
    std::unordered_set<int> &neighbors = adj[node];
    return neighbors;
}





void TreeIndex::printComponents() const {
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

// 
std::unordered_set<int> TreeIndex::findKCoreSubgraph(std::vector<int> &queryNodes, int k) {
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
    // int k_tmp = findCommenK();
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

// 
int TreeIndex::findCommenK(std::vector<int>& queryNodes) {
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
    }

    // 找到result中核心度最高的元素
    int maxK = -1;
    for (auto& c : result) {
        int tmp = coreMinimumDegree[coreIndex[*ComponentToNodes[c].begin()]];
        std::cout << "coreIndex: " << tmp << std::endl;
        if (tmp > maxK) {
            maxK = tmp;
            // std::cout << "maxK: " << maxK << std::endl;
        }
    }

    std::cout << "result: " << maxK << std::endl;
    return maxK;
}

int TreeIndex::findCommenK_2(std::vector<int>& queryNodes) {
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
    }

    query_nodes result_;
    for (auto & node : result) {
        for (auto & x : ComponentToNodes[node]) {
            result_.insert(x);
        }
    }
    // for (auto & node : adj) {
    //     result_.insert(node.first);
    // }

    int maxK = INT32_MAX;
    std::vector<int> r1 = steinerTree(result_, queryNodes);
    for (auto & node : r1) {
        std::cout << "r1: " << node << std::endl;
        int tmp = coreMinimumDegree[coreIndex[node]];
        // std::cout << "coreIndex: " << tmp << std::endl;
        if (tmp < maxK) {
            maxK = tmp;
            // std::cout << "maxK: " << maxK << std::endl;
        }
    }

    std::cout << "result: " << maxK << std::endl;
    return maxK;
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


void printHighestGroup(const std::unordered_map<int, int> &highestGroup) {
    for (const auto &pair : highestGroup)
    {
        std::cout << "Node: " << pair.first << ", Group: " << pair.second << std::endl;
    }
}

void TreeIndex::computeNodeNeighbors(Graph &graph) {
    for (auto &node : graph.getNodes())
    {
        for (int neighbor : graph.getNeighbors(node))
        {
            int shellIndex = coreIndex[neighbor];
            nodeNeighbors[node][shellIndex].insert(neighbor);
        }
    }
}

int TreeIndex::getMinimumCoreIndex(std::vector<int> queryNodes) {
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

std::unordered_set<int> TreeIndex::getNeighbors(int node, int coreIndex) {
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

std::unordered_set<int> TreeIndex::getNeighbors(int node, int coreIndex, std::unordered_set<int> &subcore) {
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

int TreeIndex::getCoreMinimumDegree(int coreIndex) {
    return coreMinimumDegree[coreIndex];
}

// 获取 map 中所有的键
std::unordered_set<int> getKeySet(const std::unordered_map<int, std::unordered_set<int>> &map) {
    std::unordered_set<int> keys;
    for (const auto &pair : map)
    {
        keys.insert(pair.first);
    }
    return keys;
}

// 获取指定核心中的节点数量
int TreeIndex::getNumberOfNodes(int coreIndex) {
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

std::unordered_set<int> TreeIndex::getCore(std::vector<int> queryNodes) {
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

int TreeIndex::getCoreIndex(int node) {
    int index = coreIndex[node];
    return coreMinimumDegree[index];
}

// 得到  点  对应的  连通分量ID
int TreeIndex::getComponent(int node) {
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
int TreeIndex::getParentComponent(int componentId) {
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

std::unordered_set<int> TreeIndex::getConnectedComponentChildren(int componentId) {
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
std::unordered_set<int> TreeIndex::getNodesInComponent(int componentId) {
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
void TreeIndex::printCoreIndex() {
    std::cout << "Core Index:" << std::endl;
    for (const auto &pair : coreIndex)
    {
        std::cout << "Node: " << pair.first << ", Core: " << pair.second << std::endl;
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
    for (int i = 0; i < n; i++) {
        pq_sustain[i] = {0, 0, 0};
    }

    // 最小度数 mu_star
    int mu_star = INT8_MAX;

    // 初始化并查集
    UnionFind uf(n);

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
        int count = n - uf.getSetNum() + 1;
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
    for (int i = 0; i < n; i++) {
        pq_sustain[i] = {0, 0, 0};
    }

    // 最小度数 mu_star
    int mu_star = INT8_MAX;

    // 初始化并查集
    UnionFind uf(n);

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
        int count = n - uf.getSetNum() + 1;
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
    std::cout << "开始" << std::endl;
    // 图的节点数 如果 int n = n  n就会变成0
    std::cout << "开始" << n << std::endl;
    std::vector<bool> isTerminal(n, false); // 标记是否为终端节点
    for (int v : terminals) isTerminal[v] = true;

    // 最小生成树（MST）的 Prim 算法
    std::vector<int> mstParent(n, -1); // MST 中的父节点
    std::vector<bool> inMST(n, false); // 标记是否在 MST 中
    std::vector<int> key(n, INT_MAX); // 每个节点的键值（最小边权重）
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq; // 优先队列

    std::cout << "开始" << std::endl;
    // 从任意一个终端节点开始
    pq.push({0, terminals[0]});
    key[terminals[0]] = 0;

    while (!pq.empty()) {
        int u = pq.top().second; // 当前节点
        std::cout << "node" << u << std::endl;
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

    std::cout << "这里" << std::endl;
    // 提取包含终端节点的子树
    std::vector<int> steinerTree;
    for (int v : terminals) {
        while (v != -1) { // 从终端节点向上追溯到根
            steinerTree.push_back(v);
            v = mstParent[v];
        }
    }
    std::cout << "哪里" << std::endl;
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

// 辅助函数，固定最后-----------------------------------------------------------------------------
// 计算shell_count，coreIndex，coreMinimumDegree
void TreeIndex::computeCoreIndex(Graph &graph) {
    shell_count = 0;
    coreIndex = CoreGroup::coreGroupsAlgorithm(graph);
    
    std::set<int> coreIndexSet;  // 用于存储核心度的集合
    for (auto &pair : coreIndex) {
        coreIndexSet.insert(pair.second);
    }
    int node = 0;
    while (!coreIndexSet.empty()) {
        coreMinimumDegree[node] = *coreIndexSet.begin();

        // 1.16 优化，如果核心度大于当前记录的层数
        if (coreMinimumDegree[node] > shell_count) {
            shell_count = coreMinimumDegree[node];
        }

        for (auto &pair : coreIndex) {
            if (pair.second == *coreIndexSet.begin()) {
                pair.second = node;
            }
        }
        coreIndexSet.erase(coreIndexSet.begin());
        node++;
    }
}

// 划分了每个shell中的不同连通分量
void TreeIndex::identifyAndStoreComponents(Graph &graph) {
    std::unordered_map<int, std::unordered_set<int>> nodesInShell; // 核心ID  对应  点集

    for (const auto &pair : coreIndex) {

        nodesInShell[pair.second].insert(pair.first);
    }

    // 根据每一层k-shell的节点，找到每一层的连通分量
    for (const auto &shell : nodesInShell) {
        int shellIndex = shell.first;
        const auto &nodes = shell.second;

        // BFS
        std::unordered_set<int> visited;
        for (int node : nodes) {

            if (visited.count(node) == 0) {
                std::queue<int> queue;
                std::unordered_set<int> component;
                queue.push(node);
                visited.insert(node);

                while (!queue.empty()) {
                    int currentNode = queue.front();
                    queue.pop();
                    component.insert(currentNode);

                    for (int neighbor : graph.getNeighbors(currentNode)) {
                        // 只找未访问过的同一层邻居
                        if (nodes.count(neighbor) > 0 && visited.count(neighbor) == 0) {
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

// 建立了不同shell连通分量之间的父子关系
void TreeIndex::buildParentChildRelationships(Graph &graph) {
    // 父子索引的处理
    int size = layerToComponentToNodes.size();
    for (int currentLevel = 0; currentLevel < size; ++currentLevel) {
        auto &currentLayerComponents = layerToComponentToNodes[currentLevel];

        for (auto &component : currentLayerComponents) {
            int currentComponentId = component.first;
            auto &nodesInCurrentComponent = component.second;

            for (int node : nodesInCurrentComponent) {
                for (int neighbor : graph.getNeighbors(node)) {
                    if (coreMinimumDegree[coreIndex[node]] >= coreMinimumDegree[coreIndex[neighbor]]) {
                        continue;
                    }
                    if (nodeToComponentId.find(neighbor) != nodeToComponentId.end() && nodeToComponentId[neighbor] != currentComponentId) {
                        int neighborComponentId = nodeToComponentId[neighbor];
                        connectedComponentParent[currentLevel][currentComponentId] = neighborComponentId;
                        ComponentParent[currentComponentId].insert(neighborComponentId);
                    }
                }
            }
            // 顶层分量的父分量是-1
            if (connectedComponentParent[currentLevel].find(currentComponentId) == connectedComponentParent[currentLevel].end() && ComponentParent[currentComponentId].empty()) {
                connectedComponentParent[currentLevel][currentComponentId] = -1;
                ComponentParent[currentComponentId].insert(-1);
            }
        }
    }

    for (int currentLevel = size - 1; currentLevel >= 0; currentLevel--) {
        auto &currentLayerComponents = layerToComponentToNodes[currentLevel];

        for (auto &component : currentLayerComponents) {
            int currentComponentId = component.first;
            auto &nodesInCurrentComponent = component.second;

            for (int node : nodesInCurrentComponent) {
                for (int neighbor : graph.getNeighbors(node)) {
                    if (coreMinimumDegree[coreIndex[node]] <= coreMinimumDegree[coreIndex[neighbor]]) {
                        continue;
                    }
                    if (nodeToComponentId.find(neighbor) != nodeToComponentId.end() && nodeToComponentId[neighbor] != currentComponentId) {
                        int neighborComponentId = nodeToComponentId[neighbor];
                        connectedComponentChildren[currentLevel][currentComponentId].insert(neighborComponentId);
                        ComponentChildren[currentComponentId].insert(neighborComponentId);
                    }
                }
            }
            // 底层分量的子分量是-1
            if (connectedComponentChildren[currentLevel].find(currentComponentId) == connectedComponentChildren[currentLevel].end()) {
                connectedComponentChildren[currentLevel][currentComponentId].insert(-1);
                ComponentChildren[currentComponentId].insert(-1);
            }
        }
    }
}

// 找到顶层分量，包括多余的顶层分量
std::unordered_set<int> TreeIndex::findTopComponents(std::vector<int> &queryNodes, int k) {
    std::unordered_set<int> topComponents;
    std::unordered_set<int> queryComponents;
    std::queue<int> bfsQueue;
    std::unordered_set<int> visited;

    // 初始化：将查询节点的组件ID添加到队列和已访问集合中
    for (int node : queryNodes) {
        int componentId = nodeToComponentId[node];

        if (visited.find(componentId) == visited.end()) {

            queryComponents.insert(componentId);
            bfsQueue.push(componentId);
            visited.insert(componentId);

        }
    }

    while (!bfsQueue.empty()) {
        int currentComponentId = bfsQueue.front();
        bfsQueue.pop();

        // 探测所有父分量
        for (int parentComponentId : ComponentParent[currentComponentId]) {
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
        // 3.19 这一步会加入无关的点，但会使结果大小和globalsearch相同
        for (int childComponentId : ComponentChildren[currentComponentId]) {
            if (visited.find(childComponentId) == visited.end()) {
                // 确保所有子分量的节点核心等级必须大于等于k
                bool eligible = true;
                for (int childNode : ComponentToNodes[childComponentId]) {
                    if (coreMinimumDegree[coreIndex[childNode]] < k) {
                        eligible = false;
                        break;
                    }
                }
                if (eligible) {
                    bfsQueue.push(childComponentId);
                    visited.insert(childComponentId);
                }
            }
        }
    }

    return topComponents;
}

// 找到查询顶点集的k值，即最小的shell层
// int TreeIndex::findCommenShell(std::vector<int> queryNodes) {

//     // dfs
//     std::stack<int> line;

//     int next_end = -1;
//     int k_ans = INT32_MAX;
//     int k_end = 0;
//     for (int i = 0; i < queryNodes.size()-1; i++) {
//         std::unordered_set<int> visited;
//         std::vector<int> tree;

//         std::cout << "cishi：" << i << std::endl;
//         if (next_end == -1) {
//             line.push(nodeToComponentId[queryNodes[0]]);
//             visited.insert(nodeToComponentId[queryNodes[0]]);
//         } else {
//             line.push(next_end);
//             visited.insert(next_end);
//         }
//         int need = nodeToComponentId[queryNodes[i+1]];

//         int next = need;
//         int k = INT32_MAX;

//         while (!line.empty()) {
//             int tmp = line.top();
//             line.pop();
//             tree.push_back(tmp);

//             // for (auto node : ComponentToNodes[tmp]) {
//             //     if (k > coreMinimumDegree[coreIndex[node]]) {
//             //         k = coreMinimumDegree[coreIndex[node]];
//             //         next = tmp;
//             //     }
//             //     break;
//             // }
    
//             if (tmp == need) {
//                 std::cout << "============：" << k << std::endl;
//                 for (auto node : tree) {
//                     // 如果最小核心度大于k
//                     if (coreMinimumDegree[coreIndex[*ComponentToNodes[node].begin()]] > k) {
//                         k = coreMinimumDegree[coreIndex[*ComponentToNodes[node].begin()]];
//                         next = node;
//                     }
//                 }
//                 if (k_end < k) {
//                     k_end = k;
//                     next_end = next;
//                 }
//             }
            
//             bool flag_1 = false, flag_2 = false;
//             for (auto node : getChildShell(tmp)) {
//                 if (visited.find(node) == visited.end()) {
//                     line.push(node);
//                     visited.insert(node);
//                     flag_1 = true;
//                 }
//             }
//             for (auto node : getParentShell(tmp)) {
//                 if (visited.find(node) == visited.end()) {
//                     line.push(node);
//                     visited.insert(node);
//                     flag_2 = true;
//                 }
//             }
//             // 刷新k
//             if (!flag_1 && !flag_2) {
//                 std::cout << "刷新！！！！！！！！！" <<std::endl;
//                 int k = INT32_MAX;
//                 tree.pop_back();
//             }

//         }
//         k_ans = std::min(k_ans, k_end);
//     }

//     return k_ans;
// }

// 我们先默认第一次找到的结果就是最好的
// void TreeIndex::dfs(int com_1, int com_2, std::unordered_set<int>& visited, std::unordered_set<int>& result, int depth = 0) {
//     if (depth > 1000) { // 限制递归深度
//         std::cerr << "Depth limit exceeded" << std::endl;
//         return;
//     }

//     std::cerr << "Visiting node " << com_1 << " at depth " << depth << std::endl;
//     // std::cout << "调用dfs" << std::endl;
//     if (com_1 == com_2) {
//         std::cout << "jieshu" << std::endl;
//         result = visited;
//         return ;
//     }

//     // 先看父节点
//     for (auto node : getParentShell(com_1)) {
//         if (visited.find(node) == visited.end() && node != -1) {
//             visited.insert(node);
//             dfs(node, com_2, visited, result, depth+1);
//             visited.erase(node);
//         }
//     }

//     for (auto node : getChildShell(com_1)) {
//         if (visited.find(node) == visited.end() && node != -1) {
//             visited.insert(node);
//             dfs(node, com_2, visited, result, depth+1);
//             visited.erase(node);
//         }
//     }
// }

void TreeIndex::dfs(int com_1, int com_2, std::unordered_set<int>& visited, std::vector<int>& currentPath, std::vector<std::vector<int>>& allPaths) {
    if (visited.find(com_1) != visited.end()) {
        return; // 检测到环路，终止搜索
    }

    visited.insert(com_1);
    currentPath.push_back(com_1);

    if (com_1 == com_2) {
        allPaths.push_back(currentPath); // 找到一条路径，存储起来
        currentPath.pop_back(); // 回溯，移除当前节点
        visited.erase(com_1); // 回溯，移除当前节点
        return; // 从当前路径回退，继续搜索其他路径
    }

    for (auto node : getParentShell(com_1)) {
        if (visited.find(node) == visited.end()) {
            dfs(node, com_2, visited, currentPath, allPaths);
        }
    }

    for (auto node : getChildShell(com_1)) {
        if (visited.find(node) == visited.end()) {
            dfs(node, com_2, visited, currentPath, allPaths);
        }
    }

    visited.erase(com_1); // 回溯
    currentPath.pop_back(); // 回溯
}

int TreeIndex::findCommenShell(std::vector<int> queryNodes) {
    int k = INT32_MAX;
    int next = nodeToComponentId[queryNodes[0]];
    for (int i = 0; i < queryNodes.size()-1; i++) {
        std::unordered_set<int> visited;
        // std::unordered_set<int> result;
        // visited.insert(next);

        std::vector<int> result;
        std::vector<std::vector<int>> allPaths;
        int tmp_k = 0;
        int tmp_next = 0;
        dfs(next, nodeToComponentId[queryNodes[i+1]], visited, result, allPaths);
        for (auto &x : allPaths) {
            for (auto &node : result) {
                for (auto subNode : ComponentToNodes[node]) {
                    if (tmp_k < coreMinimumDegree[coreIndex[subNode]]) {
                        tmp_k = coreMinimumDegree[coreIndex[subNode]];
                        tmp_next = subNode;
                    }
                    break;
                }
            }  
        }
        if (k > tmp_k) {
            k = tmp_k;
            next = tmp_next;
        }

        // for (auto &node : result) {
        //     for (auto subNode : ComponentToNodes[node]) {
        //         if (k > coreMinimumDegree[coreIndex[subNode]]) {
        //             k = coreMinimumDegree[coreIndex[subNode]];
        //             next = subNode;
        //         }
        //         break;
        //     }
        // }
    }

    return k;
}