#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <stdio.h>
#include <string.h>
#include <unordered_set>
#include <random>
#include <ctime>
#include <fstream>
#include <cassert>
#include <chrono> // 包含chrono库
#include <algorithm> // 用于std::shuffle
#include <random> // 用于随机数生成
// #include <filesystem> // C++17 for file existence check
#include "Graph.h"
#include "CoreGroup.h"
#include "TreeIndex.h"
#include "SharingIndex.h"

// 洗牌算法
void fisherYatesShuffle(std::vector<int>& vec) {
    // 使用当前时间作为随机数种子
    std::srand(static_cast<unsigned>(std::time(0)));

    for (size_t i = vec.size() - 1; i > 0; --i) {
        // 生成一个0到i之间的随机索引
        size_t j = std::rand() % (i + 1);

        // 交换vec[i]和vec[j]
        std::swap(vec[i], vec[j]);
    }
}

void pro_1_1(Graph& graph, std::string path) {
    // 获取图中所有节点
    std::vector<int> allNodes;
    for (const auto& pair : graph.getAdj()) {
        allNodes.push_back(pair.first);
    }

    int k = 50;
    int num = k;
    std::ofstream outFile(path);
    if (!outFile.is_open()) {
        std::cerr << "无法打开 results.csv 文件" << std::endl;
        return;
    }

    std::vector<std::vector<int>> group;

    // 写入 CSV 头部
    outFile << "TestNumber,QueryNodes,GreedySize,GreedyMinDegree,TreeIndexSize,TreeIndexMinDegree,SameResult\n";

    int sameResultCount1 = 0;
    int sameResultCount2 = 0;

    while (k--) {
        std::cout << "第 " << (num - k) << " 次测试" << std::endl;
        // 随机打乱节点顺序
        fisherYatesShuffle(allNodes);

        // 初始化随机数生成器
        srand(static_cast<unsigned int>(time(0)));

        // 随机选择1到5个节点
        int numQueryNodes = 1 + rand() % 5;  // 随机生成1到5的数量
        std::vector<int> queryNodes(allNodes.begin(), allNodes.begin() + numQueryNodes);

        // 将 vector 转换为 unordered_set
        std::unordered_set<int> querySet(queryNodes.begin(), queryNodes.end());

        // 检查查询顶点集是否连通
        std::unordered_map<int, int> degree = graph.getDegrees();

        bool flag = true;
        flag = graph.isConnected(querySet, degree);
        if (!flag) {
            std::cout << "查询顶点集不连通！" << std::endl;
            k++;
            continue;
        }
        
        // 保存查询顶点集
        group.push_back(queryNodes);

        // 测试Greedy算法
        Graph greedySolution = graph.globalsearch(querySet);
        int greedySize = greedySolution.getN();
        int greedyMinDegree = greedySolution.getminimumDegree();

        // 测试TreeIndex算法
        TreeIndex index = TreeIndex(graph);
        std::unordered_set<int> indexSolution = index.shellsearch(querySet);
        int indexSize = indexSolution.size();
        int indexMinDegree = graph.computesubMinimumDegree(indexSolution);

        // 检查两个结果是否相同
        bool sameResult1 = (greedySize == indexSize);
        if (sameResult1) {
            sameResultCount1++;
        }
        bool sameResult2 = (greedyMinDegree == indexMinDegree);
        if (sameResult2) {
            sameResultCount2++;
        }

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        for (int node : queryNodes) {
            queryNodesStr += std::to_string(node) + " 、";
        }
        if (!queryNodesStr.empty()) {
            queryNodesStr.pop_back(); // 移除最后一个逗号
        }

        // 将结果写入 CSV 文件
        outFile << num-k << "," 
                << queryNodesStr << "," 
                << greedySize << "," 
                << greedyMinDegree << ","
                << indexSize << ","
                << indexMinDegree << ","
                << ((sameResult1 && sameResult2) ? "True" : "False") << "\n";
    }
    // 计算相同率
    double sameRate1 = static_cast<double>(sameResultCount1) / num * 100;
    double sameRate2 = static_cast<double>(sameResultCount2) / num * 100;

    // 将相同率写入 CSV 文件
    outFile << "社区大小相同率," << sameRate1 << "%,,,,,\n";
    outFile << "k相同率," << sameRate2 << "%,,,,,\n";

    outFile.close();
    
    std::cout << "两种算法社区大小相同率为: " << sameRate1 << "%" << std::endl;
    std::cout << "k相同率为: " << sameRate2 << "%" << std::endl;

    for (auto& q : group) {
        std::cout << "查询集: ";
        for (auto& i : q) { 
            std::cout << i << " ";
            }
        std::cout << std::endl;
    }

    SharingIndex dex= SharingIndex(graph);
    // dex.getKCoreToQuery(group, graph, "batch_1.csv");
}

void pro_1_2(Graph& graph, std::string path) {
    // 获取图中所有节点
    std::vector<int> allNodes;
    for (const auto& pair : graph.getAdj()) {
        allNodes.push_back(pair.first);
    }

    int k = 4; // 随机测试100次
    int num = k;
    std::ofstream outFile(path);
    if (!outFile.is_open()) {
        std::cerr << "无法打开 results.csv 文件" << std::endl;
        return;
    }

    // 写入 CSV 头部
    outFile << "TestNumber,QueryNodes,TreeIndexSize,TreeIndexMinDegree,isConnected\n";

    query_group group;

    while (k--) {

        // 打乱allNodes中元素的顺序
        fisherYatesShuffle(allNodes);

        // 初始化随机数生成器
        srand(static_cast<unsigned int>(time(0)));

        // 随机选择节点
        int numQueryNodes = 1 + rand() % 5;  // 随机生成1到5的数量
        std::vector<int> queryNodes(allNodes.begin(), allNodes.begin() + numQueryNodes);

        // 将 vector 转换为 unordered_set
        std::unordered_set<int> querySet(queryNodes.begin(), queryNodes.end());

        // 检查查询顶点集是否连通
        std::unordered_map<int, int> degree = graph.getDegrees();

        bool flag = true;
        flag = graph.isConnected(querySet, degree);
        if (!flag) {
            std::cout << "查询顶点集不连通！" << std::endl;
            k++;
            continue;
        }

        // 保存查询顶点集
        group.push_back(querySet);

        // 测试TreeIndex算法
        TreeIndex index = TreeIndex(graph);
        std::unordered_set<int> indexSolution = index.shellsearch(querySet);
        int indexSize = indexSolution.size();
        int indexMinDegree = graph.computesubMinimumDegree(indexSolution);

        // 检查single的结果是否为连通图
        degree.clear();
        for (auto node : graph.getDegrees()) {
            if (indexSolution.find(node.first) != indexSolution.end()) {
                degree.insert(node);
            }
        }
        flag = index.isConnected(indexSolution, degree);

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        for (int node : queryNodes) {
            queryNodesStr += std::to_string(node) + " ";
        }

        // 将结果写入 CSV 文件
        outFile << num << "," 
                << queryNodesStr << ","
                << indexSize << ","
                << indexMinDegree << "," 
                << flag << "\n";
    }
    
    outFile.close();
    // for (auto& q : group) {
    //     std::cout << "查询集: ";
    //     for (auto& i : q) { 
    //         std::cout << i << " ";
    //         }
    //     std::cout << std::endl;
    // }

    SharingIndex dex= SharingIndex(graph);
    dex.batchsearch(group, "batch_2.csv");

}

void pro_2_2(Graph& graph, std::string path) {
    // 获取图中所有节点
    std::vector<int> allNodes;

    for (const auto& pair : graph.getAdj()) {
        allNodes.push_back(pair.first);
    }

    int k = 1; // 随机测试10次
    std::ofstream outFile(path);
    if (!outFile.is_open()) {
        std::cerr << "无法打开 results.csv 文件" << std::endl;
        return;
    }

    // 写入 CSV 头部
    outFile << "TestNumber,QueryNodes,TreeIndexSize,TreeIndexMinDegree,min_size,k,isConnected,time\n";

    query_group group;
    std::vector<std::vector<int>> queryNodesii = {{3079, 3202, 3598}};
    double less = 0;

    for (int i = 0; i < k; i++) {

        // 打乱allNodes中元素的顺序
        fisherYatesShuffle(allNodes);

        // 初始化随机数生成器
        srand(static_cast<unsigned int>(time(0)));

        // 随机选择节点
        // int numQueryNodes = 1 + rand() % 2;  // 随机生成1到5的数量
        // // int numQueryNodes = 3;
        // std::vector<int> queryNodes(allNodes.begin(), allNodes.begin() + numQueryNodes);
        std::vector<int> queryNodes = queryNodesii[i];

        // 将 vector 转换为 unordered_set
        std::unordered_set<int> querySet(queryNodes.begin(), queryNodes.end());

        // 检查查询顶点集是否连通
        std::unordered_map<int, int> degree = graph.getDegrees();

        bool flag = true;
        flag = graph.isConnected(querySet, degree);
        if (!flag) {
            std::cout << "查询顶点集不连通！" << std::endl;
            continue;
        }

        // 保存查询顶点集
        group.push_back(querySet);

        // 测试TreeIndex算法
        TreeIndex index = TreeIndex(graph);
        std::unordered_set<int> indexSolution = index.shellsearch(querySet);
        int indexSize = indexSolution.size();
        int indexMinDegree = graph.computesubMinimumDegree(indexSolution);
        clock_t start = clock(); // 记录开始时间
        std::unordered_set<int> community = index.greedyConnection(querySet, indexMinDegree);
        clock_t end = clock(); // 记录结束时间

        // 检查single的结果是否为连通图
        degree.clear();
        for (auto node : graph.getDegrees()) {
            if (community.find(node.first) != community.end()) {
                degree.insert(node);
            }
        }
        flag = index.isConnected(community, degree);

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        for (int node : queryNodes) {
            queryNodesStr += std::to_string(node) + " ";
        }

        // 将结果写入 CSV 文件
        outFile << i << "," 
                << queryNodesStr << ","
                << indexSize << ","
                << indexMinDegree << ","
                << community.size() << ","
                << graph.computesubMinimumDegree(community) << ","
                << flag << ","
                << (double)(end - start) / CLOCKS_PER_SEC << "\n";
                less += (double)(end - start) / CLOCKS_PER_SEC;
        
        // std::cout << "re: ";    
        // for (auto node : community) {
        //     std::cout << node << " ";
        // }
        // std::cout << std::endl; 
    }
    
    outFile.close();

    SharingIndex dex= SharingIndex(graph);
    dex.batchMinsearch(group);
    std::cout<< "less: " << less << std::endl;
}

void pro_2_3(Graph& graph, std::string path) {
    // 获取图中所有节点
    std::vector<int> allNodes;

    for (const auto& pair : graph.getAdj()) {
        allNodes.push_back(pair.first);
    }

    int k = 3;
    int num = k;
    std::ofstream outFile(path);
    if (!outFile.is_open()) {
        std::cerr << "无法打开 results.csv 文件" << std::endl;
        return;
    }

    // 写入 CSV 头部
    outFile << "TestNumber,QueryNodes,TreeIndexSize,TreeIndexMinDegree,min_size,k,isConnected\n";

    query_group group = {{785, 1084},
                         {785, 1084, 1722},
                         {785, 1084, 1722, 1117}};
    query_group group_;
    for (int i = 0; i < k; i++) {

        // 将 vector 转换为 unordered_set
        std::unordered_set<int> querySet = group[i];

        // 检查查询顶点集是否连通
        std::unordered_map<int, int> degree = graph.getDegrees();

        bool flag = true;
        flag = graph.isConnected(querySet, degree);
        if (!flag) {
            std::cout << "查询顶点集不连通！" << std::endl;
            continue;
        }
        group_.push_back(querySet);

        // 测试TreeIndex算法
        TreeIndex index = TreeIndex(graph);
        std::unordered_set<int> indexSolution = index.shellsearch(querySet);
        int indexSize = indexSolution.size();
        int indexMinDegree = graph.computesubMinimumDegree(indexSolution);
        std::unordered_set<int> community = index.greedyConnection(querySet, indexMinDegree);
        flag = graph.isConnected(community, degree);

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        for (int node : querySet) {
            queryNodesStr += std::to_string(node) + " ";
        }

        // 将结果写入 CSV 文件
        outFile << i << "," 
                << queryNodesStr << ","
                << indexSize << ","
                << indexMinDegree << ","
                << community.size() << ","
                << graph.computesubMinimumDegree(community) << ","
                << flag << "\n";
    }
    
    outFile.close();

    SharingIndex dex= SharingIndex(graph);
    // dex.batchMinsearch(group_);
}

int main(int argc, char *argv[])
{
    // 测试使用的无向图
    Graph graph("D:\\mySecre\\workspace\\csp_old\\CSP\\dataset\\facebook_combined.txt");


    pro_2_2(graph, "tt1.csv");
    // pro_1_2(graph, "single_2.csv");
    // pro_1_1(graph, "single_1.csv");

    // Graph greedySolution = graph.globalsearch(querySet);
    // int greedySize = greedySolution.getN();
    // int greedyMinDegree = greedySolution.getminimumDegree();
    // std::cout << "11: " << greedySize << "\n" << "22: " << greedyMinDegree << std::endl;
    
    // TreeIndex index = TreeIndex(graph);
    // std::unordered_set<int> indexSolution = index.shellsearch(queryNodes);
    // int indexSize = indexSolution.size();
    // int indexMinDegree = graph.computesubMinimumDegree(indexSolution);
    // std::cout << "11: " << indexSize << "\n" << "22: " << indexMinDegree << std::endl;


    

 /*
    // Graph graph; // 假设你已经加载了图数据
    TreeIndex treeIndex(graph);

    query_nodes queryNodes = {785, 1084, 1722, 1117}; // 查询节点集合 3581 | 1471, 2340 | 1, 6
    std::unordered_set<int> indexSolution = treeIndex.findKCoreSubgraph(queryNodes);
    int indexSize = indexSolution.size();
    int indexMinDegree = graph.computesubMinimumDegree(indexSolution);
    int k = indexMinDegree; // 最小度数约束
    std::cout << "k: " << k <<std::endl;

    std::unordered_set<int> community = treeIndex.greedyConnection(queryNodes, k);

    std::cout << "Community nodes: ";
    for (int node : community) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
    std::cout << "Community size: " << community.size() << std::endl;

    indexMinDegree = graph.computesubMinimumDegree(community);
    k = indexMinDegree;
    std::cout << "k: " << k <<std::endl;
 */

    return 0;
}

