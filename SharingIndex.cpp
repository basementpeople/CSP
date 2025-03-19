#include "SharingIndex.h"

void SharingIndex::getKCoreToQuery(const std::vector<std::vector<int>>  &group, Graph &graph, std::string path) {
    q_count = 0;
    std::vector<int> q_;
    int k_min = INT_MAX;

    KCoreToQuery.clear();
    for (std::vector<int> q : group) {
        QueryCode[q_count] = q;
        q_.insert(q_.end(), q.begin(), q.end());

        // 初始化最小核心索引为最大整数，以便后续比较
        int minCoreIndex = INT_MAX;

        // 遍历查询节点，确定最小核心索引
        for (int node : q) {
            // 如果节点不在核心索引中，抛出异常
            // if (coreIndex.find(node) == coreMinimumDegree.end())
            if (coreIndex.find(node) == coreIndex.end()) {
                throw std::runtime_error("Error: Node does not exist in core index.");
            }
            std::cout << "the coreIndex of node " << node << " is " << coreMinimumDegree[coreIndex[node]]<< std::endl;
            // 更新最小核心索引
            minCoreIndex = std::min(minCoreIndex, coreMinimumDegree[coreIndex[node]]);
        }

        // 找到最大的公共shell，如果这个shell的核心度>=k,不改变k；如果核心度<k，则k=核心度
        int k = findCommenK(q);
        if (k >= minCoreIndex || k == -1) {
            k = minCoreIndex;
        }

        // 确定k值
        std::cout << "the k is : " << k << std::endl;
        if (k < k_min) {
            k_min = k;
        }
        KCoreToQuery[k].insert(q_count);
        q_count++;
    }
    // 用于存储所有顶层分量的ID
    std::unordered_set<int> topComponentIds;

    for (auto node : q_) {
        std::cout << node << " ";
    }

    // 找到所有父分量为-1的顶层分量
    topComponentIds = findTopComponents(q_, k_min);

    // 用于存储最终结果的节点集合
    std::unordered_set<int> resultNodes;
    // 用于存储已经访问过的节点集合
    std::unordered_set<int> visited;
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
        if (visited.find(currentId) != visited.end())
            continue;

        visited.insert(currentId);

        // 获取shell中的一个节点
        int currentNode;
        for (int node : ComponentToNodes[currentId]) {
            currentNode = node;
            break;
        }

        // 如果子分量的节点的核心度大于等于k
        int k_tmp = coreMinimumDegree[coreIndex[currentNode]];
        if (k_tmp >= k_min) {
            // 将当前分量添加到结果集中
            kToResult[k_tmp].insert(currentId);
        }

        // 遍历当前分量的所有子分量
        for (int childCompId : ComponentChildren[currentId])
        {
            // 获取子分量中的一个节点
            int childNode;
            for (int node : ComponentToNodes[childCompId])
            {
                childNode = node;
                break;
            }

            // 如果子分量的节点的核心度大于等于k，将子分量加入队列
            int k_tmp = coreMinimumDegree[coreIndex[childNode]];
            if (k_tmp >= k_min)
            {
                componentQueue.push(childCompId);
            }

        }
    }

    // 得到结果
    for (int i = shell_count; i > 0; i--) {
        if (!kToResult[i].empty()) {
            res.insert(kToResult[i].begin(), kToResult[i].end());
            for (int compId : res) {
                for (int tmp : ComponentToNodes[compId]) {
                    resultNodes.insert(tmp);
                }
            }
        }
        if (KCoreToQuery[i].size() > 0) {
            for (int node : KCoreToQuery[i]) {
                queryToKcore[node] = i;
                queryToResult[node] = getLastResult(QueryCode[node], resultNodes, graph);
            }
        }
    }

    // 打印结果
    for (int i = 0; i < q_count; ++i) {
        std::cout << "查询顶点集" << i << ": " << std::endl;
        for (int node : QueryCode[i]) {
            std::cout << node << " ";
        }
        std::cout << std::endl;

        std::cout << "查询顶点集" << i << "的K-core: " << queryToKcore[i] << std::endl;
        std::cout << "查询顶点集" << i << "的K-core对应的子社区大小: " << queryToResult[i].size() << std::endl;
    }

    // 把结果存到csv文件
    std::ofstream outFile(path);
    if (!outFile.is_open()) {
        std::cerr << "无法打开 results.csv 文件" << std::endl;
        return;
    }
    // 写入 CSV 头部
    outFile << "TestNumber,QueryNodes,SharingIndexSize,SharingIndexMinDegree\n";

    for (int i = 0; i < q_count; ++i) {
        std::cout << "查询顶点集" << i << ": " << std::endl;
        for (int node : QueryCode[i]) {
            std::cout << node << " ";
        }
        std::cout << std::endl;

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        queryNodesStr += "{";
        for (int node : QueryCode[i]) {
            queryNodesStr += std::to_string(node) + " 、";
        }
        queryNodesStr += "}";
        if (!queryNodesStr.empty()) {
            queryNodesStr.pop_back(); // 移除最后一个逗号
        }

        // 将结果写入 CSV 文件
        outFile << q_count << "," 
                << queryNodesStr << "," 
                << queryToResult[i].size() << "," 
                << queryToKcore[i] << "\n";
    }

    // std::cout << "k_min: " << k_min << std::endl;
}

std::unordered_set<int> SharingIndex::getLastResult(std::vector<int>& queryNodes, std::unordered_set<int>& result_end, Graph &graph) {
    std::unordered_set<int> result;
    std::queue<int> q;
    for (auto& node : queryNodes) {
        q.push(node);
    }

    while (!q.empty()) {
        int node = q.front();
        q.pop();
        // std::cout << "加入节点：" << node << std::endl;
        result.insert(node);

        for (auto& neighbor : graph.getNeighbors(node)) {
            if (result_end.find(neighbor) != result_end.end() && result.find(neighbor) == result.end()) {
                q.push(neighbor);
                result.insert(neighbor);
                }
        }
    }
    return result;
}
// void SharingIndex::getMinPortToQuery(const query_group &group) {
//     for (query_nodes q : group) {
//         query_nodes tmp;
//         for (auto& p : q) {
//             std::cout << "p: " << p << std::endl;
//             // 找到q对应的连通分量
//             tmp.insert(nodeToComponentId[p]);
//             std::cout << nodeToComponentId[p] << std::endl;
//         }
//         if (tmp.size() == 1) {
//             // 找到q的最小分量
//             std::cout << "1-------" << std::endl;
//             int componentID = *tmp.begin();
//             MinPortToQuery.push_back({componentID, q});
//         }
//         else {
//             std::cout << "2-------" << std::endl;
//             int componentID = findCommonChildren(q);
//             // 如果找不到q的最小分量，说明q的子集都不在同一个连通分量中，需要将q拆分为多个子集
//             if (componentID == -1) {
//                 std::cout << "拆分子集" << std::endl;
//                 // continue;
//                 componentID = 1;
//             }
//             // 找到q的最小分量
//             // std::cout << "componentID: " << componentID << " " << coreMinimumDegree[coreIndex[*ComponentToNodes[MinPortToQuery.back().first].begin()]] <<std::endl;
//             MinPortToQuery.push_back({componentID, q});
//         }
//         // 得到核心度
//         queryToKcore[q] = coreMinimumDegree[coreIndex[*ComponentToNodes[MinPortToQuery.back().first].begin()]]; 
//         // std::cout << "MinPortToQuery: " <<  MinPortToQuery.back().first << std::endl;
//         // std::cout << "queryToKcore: " <<  queryToKcore[q] << std::endl;
//     }
// }

int SharingIndex::findCommonChildren(const query_nodes &query) {
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
            for (auto &childrenSet: getChildShell(componentID)) {
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
        countSum.insert(child.second.begin(), child.second.end());
        result.insert(child.second.begin(), child.second.end());
    }

    // 找到result中元素个数等于query.size()的元素
    int minSize = query.size();
    for (auto& c : result) {
        if (countSum.count(c) != minSize) {
            result.erase(c);
        }
        // if (countSum.count(c) == minSize) {
            // std::cout << "存在" << c << std::endl;
        // }
    }

    // 找到result中核心度最高的元素
    int maxCore = -1;
    for (auto& c : result) {
        // int tmp = coreMinimumDegree[c]; // 这是核心id对应的核心度！
        int tmp = coreMinimumDegree[coreIndex[*ComponentToNodes[c].begin()]];
        // std::cout << "coreIndex: " << tmp << std::endl;
        if (tmp > maxCore) {
            maxCore = c;
            // std::cout << "maxCore: " << maxCore << std::endl;
        }
    }

    // std::cout << "result: " << maxCore << std::endl;
    return maxCore;
}

void SharingIndex::getSharingComponents()
{
    std::unordered_set<int> queryComponents;
    std::queue<int> bfsQueue;
    std::unordered_set<int> visited;

    // 初始化：将查询节点的组件ID添加到队列和已访问集合中
    for (auto& node : MinPortToQuery)
    {
        int componentId = node.first;
        
        queryComponents.insert(componentId);
        bfsQueue.push(componentId);
        visited.insert(componentId);
    }

    while (!bfsQueue.empty())
    {
        int currentComponentId = bfsQueue.front();
        bfsQueue.pop();

        // 将当前连通分量添加到子图中
        subGraph.insert(currentComponentId);

        // 探测所有父分量
        for (int parentComponentId : ComponentParent[currentComponentId])
        {
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
    }
}

// SharingIndex SharingIndex::buildSubgraphIndex() {
//     // 创建子图的TreeIndex
//     SharingIndex subgraphIndex = SharingIndex(*this);
//     // 检查这种复制方法

//     // 构建子图的核心索引
//     std::unordered_map<int, int> subgraphCoreIndex;
//     for (int node : subGraph) {
//         // 1.9 应该是 连通分量的ID 对应 核心度
//         int coreIndexValue = coreIndex[*ComponentToNodes[node].begin()];
//         subgraphCoreIndex[node] = coreIndexValue;
//     }

//     // 更新子图的coreIndex 连通分量ID 对应 核心度，有什么用处吗？ 1.9
//     subgraphIndex.coreIndex = subgraphCoreIndex;

//     // 构建子图的连通分量关系，之后需要判断多余的节点
//     for (int node : subGraph) {
//         subgraphIndex.ComponentParent[node] = ComponentParent[node];
//     }
//     for (int node : subGraph) {
//         subgraphIndex.ComponentChildren[node] = ComponentChildren[node];
//     }

//     return subgraphIndex;
// }

// void SharingIndex::findCommunities() {
//     std::queue<int> componentQueue; // 用于广度优先搜索的队列
//     std::unordered_set<int> visited; // 记录已经访问过的分量

//     for (auto& node : topComponents) {
//         componentQueue.push(node);
//     }

//     // 广度优先搜索所有分量
//     while (!componentQueue.empty()) {
//         int currentId = componentQueue.front();
//         componentQueue.pop();
//         componentToResult[currentId].insert(currentId);
//         visited.erase(currentId);

//         // 遍历当前分量的所有子分量
//         for (int child : ComponentChildren[currentId]) {
//             if (subGraph.find(child) != subGraph.end()) {

//                 componentToResult[child].insert(componentToResult[currentId].begin(), componentToResult[currentId].end());
                
//                 if (child != -1 && visited.find(child) == visited.end()) {
//                     componentQueue.push(child);
//                     visited.insert(child);
//                 }
//             }
//         }
//     }


//     for (auto& pair: MinPortToQuery) {
//         int componentId = pair.first;
//         query_nodes query = pair.second;

//         for (auto &node : componentToResult[componentId]) {
//             queryToResult[query].insert(ComponentToNodes[node].begin(), ComponentToNodes[node].end());
//         }
//     }

//     // 打印结果
//     for (auto &node : queryToResult) {
//         std::cout << "Query: ";
//         for (auto &queryNode : node.first) {
//             std::cout << queryNode << " ";
//         }
//         std::cout << "k-core: " << queryToKcore[node.first] << std::endl;
//         std::cout << "size: " << node.second.size() << std::endl;
//     }

// }