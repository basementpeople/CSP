#include "SharingIndex.h"

// batch查找解决CSP
void SharingIndex::batchsearch(query_group  &group, std::string path = "") {
    
    // 重置成员变量
    // 一次调用清空一次
    q_count = 0;  // 查询顶点集的数量
    queryToKcore.clear();  // 清空查询编号到核心度的映射
    queryToResult.clear();  // 清空查询编号到结果节点的映射
    
    // 分成不同的group
    std::vector<query_group> groups;
    for (query_nodes& q : group) {
        query_group group;
        group.push_back(q);
        groups.push_back(group);
    }

    while (groups.size() > 1) {
        double CP = 0;
        int bestPair[2] = {-1, -1};

        // 找到最相似的两个组
        for (size_t i = 0; i < groups.size(); ++i) {
            for (size_t j = i + 1; j < groups.size(); ++j) {
                for (query_nodes& qA : groups[i]) {
                    for (query_nodes& qB : groups[j]) {
                        std::unordered_set<int> x = findTopComponents(qA);
                        std::unordered_set<int> y = findTopComponents(qB);
                        std::unordered_set<int> uni;
                        uni.insert(x.begin(), x.end());
                        uni.insert(y.begin(), y.end());
                        if (uni.size() == x.size()) {
                            CP = 1;
                        } else {
                            CP = 0;
                        }
                        break;
                    }
                    break;
                }

                // std::cout << "此次聚类的相似度为：" << similarity << std::endl;
                if (CP == 1) {
                    bestPair[0] = i;
                    bestPair[1] = j;
                }
            }
        }
        // 合并最相似的两个组
        groups[bestPair[0]].insert(groups[bestPair[0]].end(), 
        std::make_move_iterator(groups[bestPair[1]].begin()), 
        std::make_move_iterator(groups[bestPair[1]].end()));
        groups.erase(groups.begin() + bestPair[1]); 

    }

    
    for (auto group : groups) {
        std::unordered_set<int> q_; // 存储所有的查询顶点
        int k_min = INT_MAX; // 对于这个group，最小的k
        KCoreToQuery.clear();  // 清空核心度到查询编号的映射
        // queryToResult
        // queryToKcore

        for (std::unordered_set<int> q : group) {
            QueryCode[q_count] = q;
            q_.insert(q.begin(), q.end());
    
            int k1 = INT32_MAX;
            for (auto node : q) {
                k1 = std::min(k1, coreMinimumDegree[coreIndex[node]]);
            }
    
            // 确定最小核心索引
            int k = k1;
    
            // 确定k_min值
            if (k < k_min) {
                k_min = k;
            }
            KCoreToQuery[k].insert(q_count);
            q_count++;
        }
        
        searchstep(q_, k_min);
        std::unordered_set<int> res; // 用于存储最终结果的连通分量集合
        std::unordered_set<int> resultNodes; // 用于存储最终结果的节点集合
        // 得到结果
        for (int i = shell_count; i > 0; i--) {

            // k_min不准确，很少出现
            if (i < k_min && KCoreToQuery[i].size() > 0) {
                if (k_min > 1) { k_min--; }
                searchstep(q_, k_min);
            }

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
                    // queryToResult[node] = resultNodes;
                    int k = i;
                    // checkConnected();

                    // 检查Q是否连通
                    std::unordered_map<int, int> degree;
                    degree.clear();
                    for (auto node : Graph::getDegrees()) {
                        if (resultNodes.find(node.first) != resultNodes.end()) {
                            degree.insert(node);
                        }
                    }
                    bool flag = isConnected(QueryCode[node], degree, resultNodes);
                    queryToResult[node] = resultNodes;
            
                    // 如果不连通，k-1，加入新的点再检查
                    if (!flag) {
                        if (k >= k_min) { KCoreToQuery[k-1].insert(node); }
                    }
                    
                }
            }
        }    

    }

    printAndwrite(path);
}

// 检查查询集是否连通
bool SharingIndex::isConnected(query_nodes queryNodes, std::unordered_map<int, int> degree, std::unordered_set<int> &resultNodes) {
    if (queryNodes.empty() || degree.empty()) return false; // 如果查询节点为空，直接返回空

    std::unordered_set<int> visited; // 存储已访问的节点
    std::queue<int> q;
    for (int startNode : queryNodes) {
        if (degree.find(startNode) != degree.end() && degree.at(startNode) > 0) {
            q.push(startNode);
            visited.insert(startNode);
            break;
        }
    }

    if (q.empty()) {
        return false; // 如果没有有效的起始节点，则直接返回不连通
    }

    while (!q.empty()) {
        int node = q.front();
        q.pop();
        for (int neighbor : Graph::getNeighbors(node)) {
            if (degree.find(neighbor) != degree.end() && visited.insert(neighbor).second && degree.at(neighbor) > 0) {
                q.push(neighbor);
            }
        }
    }

    // 检查所有查询节点是否都在访问集合中
    for (int node : queryNodes) {
        if (visited.find(node) == visited.end()) {
            return false;
        }
    }

    resultNodes = visited;
    return true;
}

void SharingIndex::searchstep(std::unordered_set<int> &q_, int k_min) {
    kToResult.clear();  // 清空核心度到分量的映射
    // 找到所有父分量为-1的顶层分量，包括多余的顶层分量，同Treeindex
    std::unordered_set<int> topComponentIds = findTopComponents(q_, k_min);

    std::unordered_set<int> visited; // 用于存储已经访问过的节点集合
    std::queue<int> componentQueue; // 用于广度优先搜索的队列

    // 将顶层分量的ID加入队列
    for (int compId : topComponentIds) {
        componentQueue.push(compId);
    }

    // 广度优先搜索所有分量，更新kToResult
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
        for (int childCompId : ComponentChildren[currentId]) {
            // 获取子分量中的一个节点
            int childNode;
            for (int node : ComponentToNodes[childCompId]) {
                childNode = node;
                break;
            }

            // 如果子分量的节点的核心度大于等于k，将子分量加入队列
            int k_tmp = coreMinimumDegree[coreIndex[childNode]];
            if (k_tmp >= k_min) {
                componentQueue.push(childCompId);
            }

        }
    }
}

// 辅助函数
std::unordered_set<int> SharingIndex::getLastResult(std::vector<int>& queryNodes, std::unordered_set<int>& result_end, Graph &graph) {
    std::unordered_set<int> result;
    std::queue<int> q;
    for (auto& node : queryNodes) {
        q.push(node);
    }

    while (!q.empty()) {
        int node = q.front();
        q.pop();
        result.insert(node);

        for (auto& neighbor : graph.getNeighbors(node)) {
            if (result_end.find(neighbor) != result_end.end() && result.find(neighbor) == result.end()) {
                q.push(neighbor);
                }
        }
    }
    return result;
}

void SharingIndex::printAndwrite(std::string path){

    // 打印结果
    // for (int i = 0; i < q_count; ++i) {
    //     std::cout << "查询顶点集" << i << ": " << std::endl;
    //     for (int node : QueryCode[i]) {
    //         std::cout << node << " ";
    //     }
    //     std::cout << std::endl;

    //     std::cout << "查询顶点集" << i << "的K-core: " << queryToKcore[i] << std::endl;
    //     std::cout << "查询顶点集" << i << "的K-core对应的子社区大小: " << queryToResult[i].size() << std::endl;
    // }

    // 不写入
    if (path.empty()) {
        return ;
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
        // for (int node : QueryCode[i]) {
        //     std::cout << node << " ";
        // }
        // std::cout << std::endl;

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        for (int node : QueryCode[i]) {
            queryNodesStr += std::to_string(node) + " ";
        }

        // 将结果写入 CSV 文件
        outFile << q_count << "," 
                << queryNodesStr << "," 
                << queryToResult[i].size() << "," 
                << queryToKcore[i] << "\n";
    }

}

void SharingIndex::printAndwrite_(std::string path){

    // 不写入
    if (path.empty()) {
        return ;
    }

    // 把结果存到csv文件
    std::ofstream outFile(path);
    if (!outFile.is_open()) {
        std::cerr << "无法打开 results.csv 文件" << std::endl;
        return;
    }
    // 写入 CSV 头部
    outFile << "QueryCode,QueryNodes,ResultSize,K,time,species,time_sum,isConnected\n";

    for (int i = 0; i < q_count_; ++i) {

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        for (int node : QueryCode_[i]) {
            queryNodesStr += std::to_string(node) + " ";
        }

        // 将结果写入 CSV 文件
        outFile << i << "," 
                << queryNodesStr << "," 
                << queryToResult_[i].size() << ","
                << queryToKcore_[i] << "," 
                << time[i] << ","
                << species[i] << ","
                << time_sum[species[i]] << ","
                << flag_[i] << "\n";

        std::cout << "re: ";    
        for (auto node : queryToResult_[i]) {
            std::cout << node << " ";
        }
        std::cout << std::endl; 
    }

}

// 4 batch查找解决MIN_CSP
void SharingIndex::batchMinsearch(query_group& group) {
    for (int i = 0; i < group.size(); ++i) {
        queryToKcore_[i] = i;
        queryToResult_[i] = {i};
        flag_.push_back(false);
        time.push_back(0);
        species.push_back(0);
        time_sum.push_back(0);
    }

    std::cout << "开始！！！" << std::endl;
    // 初始化数据成员
    q_count_ = 0;
    double record = 0;
    int x = 0;
    int count = 0;

    clock_t start = clock(); // 记录开始时间
    // 先聚类，得到多个group
    std::vector<query_group> new_group;
    new_group = Clustering(group, 0.1);
    clock_t end = clock(); // 记录结束时间
    record += (double)(end - start) / CLOCKS_PER_SEC;
    std::cout << "此次聚类花费了" << (double)(end - start) / CLOCKS_PER_SEC << "秒" << std::endl;

    // 对每个group分别做一次greedy算法
    for (auto g : new_group) {
        
        // 前置准备，不计入时间
        batchsearch(g);
        // 动态生成字符串 s，例如 "a1", "a2", "a3" 等
        // std::string s = "a" + std::to_string(counter);
        // printAndwrite(s);
        // 更新计数器
        // counter++;

        // 得到Q和对应的k_min
        // std::unordered_map<std::unordered_set<int>, int, UnorderedSetHash, UnorderedSetEqual> QtoK;
        // for (auto pair : KCoreToQuery) {
        //     for (auto node : pair.second) {
        //         QtoK[QueryCode[node]] = pair.first;
        //     }
        // }
        clock_t start = clock(); // 记录开始时间

        int k_min = INT32_MAX;
        query_nodes Q;
        for (auto pair : KCoreToQuery) {
            if (!pair.second.empty()) {
                if (pair.first < k_min) {
                    k_min = pair.first;
                }
            }
        }

        for (auto nodes : g) {
            Q.insert(nodes.begin(), nodes.end());
        }

        std::unordered_set<int> result;
        for (auto nodes : Q) {
            std::cout << "q  " << nodes << "  "; 
        }
        std::cout << "k_min  " << k_min << "  "; 
        result = greedyStep(Q, k_min);
        std::cout << "未扩展的结果大小：" << result.size() << std::endl;

        // 扩展结果
        // 方法一：如果邻居的核心度大于等于k_min，加入结果

        // 严格
        // for (auto node : Q) {
        //     for (auto neighbor : Graph::getNeighbors(node)) {
        //         if (coreMinimumDegree[coreIndex[neighbor]] >= k_min) {
        //             result.insert(neighbor);
        //         }
        //     }
        // }

        // 宽泛
        // std::queue<int> q;
        // std::unordered_set<int> visited;
        // for (auto node : Q) {
        //     q.push(node);
        //     visited.insert(node);
        // }
        // while (!q.empty()) {
        //     int tmp = q.front();
        //     q.pop();
        //     visited.insert(tmp);

        //     for (auto neighbor : Graph::getNeighbors(tmp)) {
        //         if (visited.find(neighbor) == visited.end()) {
        //             visited.insert(neighbor);
        //             if (result.find(neighbor) == result.end()) {
        //                 if (coreMinimumDegree[coreIndex[neighbor]] >= k_min) {
        //                     result.insert(neighbor);
        //                     q.push(neighbor);
        //                 }
        //             }
        //         }
        //     }
        // }
        // std::cout << "扩展后的结果大小：" << result.size() << std::endl;

        clock_t end = clock();
        double AA = ((double)(end - start) / CLOCKS_PER_SEC);
        record += ((double)(end - start) / CLOCKS_PER_SEC);
        clock_t start_ = clock();
        // 得到不同group的结果 result
        // 对不同的group中的每个Q单独做Connection
        for (auto q : g) {
            clock_t start = clock();
            int k = 0;
            for (auto pair : QueryCode) {
                if (q == pair.second) {
                    k = pair.first;
                }
            }
            k = queryToKcore[k];
// ========
            std::unordered_set<int> res = result;

            if (g.size() != 1) {
                std::queue<int> que;
                std::unordered_set<int> visited;
                for (auto node : Q) {
                    que.push(node);
                    visited.insert(node);
                }
                while (!que.empty()) {
                    int tmp = que.front();
                    que.pop();
                    visited.insert(tmp);
        
                    for (auto neighbor : Graph::getNeighbors(tmp)) {
                        if (visited.find(neighbor) == visited.end()) {
                            visited.insert(neighbor);
                            if (res.find(neighbor) == res.end()) {
                                if (coreMinimumDegree[coreIndex[neighbor]] >= k) {
                                    res.insert(neighbor);
                                    que.push(neighbor);
                                }
                            }
                        }
                    }
                }
                std::cout << "扩展后的结果大小：" << res.size() << std::endl;
            }
            
// ========
            std::unordered_set<int> ans;
            ans = connectionStep(res, q, k);

            // 保存每个Q的结果图和k
            int tmp = 0;
            for (auto node : QueryCode_) {
                if (node.second == q) {
                    tmp = node.first;
                }
            }
            queryToKcore_[tmp] = k;
            queryToResult_[tmp] = ans;

            q_count_++;

            std::unordered_map<int, int> degree_;
            for (auto node : degrees) {
                if (ans.find(node.first) != ans.end()) {
                    degree_.insert(node);
                }
            }
            flag_[tmp] = Graph::isConnected(ans, degree_);
            std::cout << "flag: " << flag_[tmp] << std::endl;
            species[tmp] = count;

            clock_t end = clock(); // 记录开始时间
            time[tmp] = AA + ((double)(end - start) / CLOCKS_PER_SEC);
        }
        clock_t end_ = clock();
        time_sum[count] = (double)(end_ - start) / CLOCKS_PER_SEC;
        count++;
        record += (double)(end_ - start_) / CLOCKS_PER_SEC;
        std::cout << "此次单个类的处理花费了" << (double)(end - start) / CLOCKS_PER_SEC << "秒" << std::endl;

    }
    std::string s = "tt2.csv";
    printAndwrite_(s);
    std::cout << "结束！！！" << std::endl;
    std::cout << "此次batchMinsearch花费了" << record << "秒" << std::endl;
}

// 与batch相关的 聚类算法
// 计算两个查询之间的相似度
double SharingIndex::querySimilarity(query_nodes& qA, query_nodes& qB) {

    // 判断两个集合是否属于同一连通分量
    // 不能用公共父节点判断
    // findCommenParent(qA, qB);

    // qA的邻居
    std::unordered_set<int> neighborsA;
    for (auto &v :qA) {
        neighborsA.insert(Graph::getNeighbors(v).begin(), Graph::getNeighbors(v).end());
    }
    // std::cout << "neighborsA " << neighborsA.size() << std::endl;

    // qB的邻居
    std::unordered_set<int> neighborsB;
    for (auto &v :qB) {
        neighborsB.insert(Graph::getNeighbors(v).begin(), Graph::getNeighbors(v).end());
    }
    // std::cout << "neighborsB " << neighborsB.size() << std::endl;

    // 计算两个邻居的交集
    // 得排序
    // 将 unordered_set 转换为 vector 并排序
    std::vector<int> sortedA(neighborsA.begin(), neighborsA.end());
    std::vector<int> sortedB(neighborsB.begin(), neighborsB.end());
    std::sort(sortedA.begin(), sortedA.end());
    std::sort(sortedB.begin(), sortedB.end());

    std::unordered_set<int> intersection;
    std::set_intersection(sortedA.begin(), sortedA.end(),
                          sortedB.begin(), sortedB.end(),
                          std::inserter(intersection, intersection.begin()));

    // 计算两个邻居的并集
    std::unordered_set<int> unionSet;
    unionSet.insert(neighborsA.begin(), neighborsA.end());
    unionSet.insert(neighborsB.begin(), neighborsB.end());

    // std::cout << "交集：" << intersection.size() << std::endl;
    // std::cout << "并集：" << unionSet.size() << std::endl;

    // 得到cp值
    double CP = 2;
    std::unordered_set<int> x = findTopComponents(qA);
    std::unordered_set<int> y = findTopComponents(qB);
    std::unordered_set<int> uni;
    uni.insert(x.begin(), x.end());
    uni.insert(y.begin(), y.end());
    if (uni.size() == x.size()) {
        CP = 1;
    } else {
        CP = 0;
    }

    if (unionSet.empty()) {
        return 0;
    } else {
        // 计算相似度
        // return intersection.size() / (double)(std::min(qA.size(), qB.size()));
        // std::cout << static_cast<double>(intersection.size()) / unionSet.size() << std::endl;
        return CP*static_cast<double>(intersection.size()) / unionSet.size();
    }
}

// 计算两个查询顶点集之间的相似度, 1 - 邻居
double SharingIndex::groupSimilarity(query_group& groupA, query_group& groupB) {
    double totalSimilarity = 0.0;
    for (query_nodes& qA : groupA) {
        for (query_nodes& qB : groupB) {
            totalSimilarity += querySimilarity(qA, qB);
        }
    }
    return totalSimilarity / (groupA.size() * groupB.size());
}

// 聚类算法
std::vector<query_group> SharingIndex::Clustering(query_group& query_groups, double threshold) {
    std::vector<query_group> groups;
    // 初始化每个查询为一个单独的组
    int i = 0;
    for (const query_nodes& q : query_groups) {
        QueryCode_[i] = q;
        query_group group;
        group.push_back(q);
        groups.push_back(group);
        i++;
    }

    while (groups.size() > 1) {
        double maxSimilarity = -1.0;
        int bestPair[2] = {-1, -1};

        // 找到最相似的两个组
        for (size_t i = 0; i < groups.size(); ++i) {
            for (size_t j = i + 1; j < groups.size(); ++j) {
                double similarity = groupSimilarity(groups[i], groups[j]);
                // std::cout << "此次聚类的相似度为：" << similarity << std::endl;
                if (similarity > maxSimilarity) {
                    maxSimilarity = similarity;
                    bestPair[0] = i;
                    bestPair[1] = j;
                }
            }
        }
        std::cout << "此次聚类的最大相似度为：" << maxSimilarity << std::endl;
        if (maxSimilarity < threshold) {
            break;
        }
        // 合并最相似的两个组
        groups[bestPair[0]].insert(groups[bestPair[0]].end(), 
                           std::make_move_iterator(groups[bestPair[1]].begin()), 
                           std::make_move_iterator(groups[bestPair[1]].end()));
        groups.erase(groups.begin() + bestPair[1]); 
    }

    
    std::cout << "groups数量：" << groups.size() << std::endl;
    return groups;
}