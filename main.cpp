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
        std::unordered_set<int> indexSolution = index.findKCoreSubgraph(queryNodes);
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
    dex.getKCoreToQuery(group, graph, "batch_1.csv");
}

void pro_1_2(Graph& graph, std::string path) {
    // 获取图中所有节点
    std::vector<int> allNodes;
    // for (int i = 0; i < graph.getAdj().size(); ++i) {
    //     // 默认数据集的每个点都有意义，不会出现1，2，70，71~
    //     allNodes.push_back(i);
    // }
    for (const auto& pair : graph.getAdj()) {
        allNodes.push_back(pair.first);
    }

    int k = 10; // 随机测试100次
    int num = k;
    std::ofstream outFile(path);
    if (!outFile.is_open()) {
        std::cerr << "无法打开 results.csv 文件" << std::endl;
        return;
    }

    // 写入 CSV 头部
    outFile << "TestNumber,QueryNodes,TreeIndexSize,TreeIndexMinDegree\n";

    std::vector<std::vector<int>> group;

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
        group.push_back(queryNodes);

        // 测试TreeIndex算法
        TreeIndex index = TreeIndex(graph);
        std::unordered_set<int> indexSolution = index.findKCoreSubgraph(queryNodes);
        int indexSize = indexSolution.size();
        int indexMinDegree = graph.computesubMinimumDegree(indexSolution);

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        queryNodesStr += "{";
        for (int node : queryNodes) {
            queryNodesStr += std::to_string(node) + " 、";
        }
        queryNodesStr += "}";
        if (!queryNodesStr.empty()) {
            queryNodesStr.pop_back(); // 移除最后一个逗号
        }

        // 将结果写入 CSV 文件
        outFile << num << "," 
                << queryNodesStr << ","
                << indexSize << ","
                << indexMinDegree << "\n";
    }
    
    outFile.close();
    for (auto& q : group) {
        std::cout << "查询集: ";
        for (auto& i : q) { 
            std::cout << i << " ";
            }
        std::cout << std::endl;
    }

    SharingIndex dex= SharingIndex(graph);
    dex.getKCoreToQuery(group, graph, "batch_2.csv");

}

void pro_2_2(Graph& graph, std::string path) {
    // 获取图中所有节点
    std::vector<int> allNodes;

    for (const auto& pair : graph.getAdj()) {
        allNodes.push_back(pair.first);
    }

    int k = 10; // 随机测试10次
    int num = k;
    std::ofstream outFile(path);
    if (!outFile.is_open()) {
        std::cerr << "无法打开 results.csv 文件" << std::endl;
        return;
    }

    // 写入 CSV 头部
    outFile << "TestNumber,QueryNodes,TreeIndexSize,TreeIndexMinDegree,test1_size,test2_k\n";

    std::vector<std::vector<int>> group;

    while (k--) {

        // 打乱allNodes中元素的顺序
        fisherYatesShuffle(allNodes);

        // 初始化随机数生成器
        srand(static_cast<unsigned int>(time(0)));

        // 随机选择节点
        int numQueryNodes = 1 + rand() % 5;  // 随机生成1到10的数量
        // int numQueryNodes = 3;
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

        // 测试TreeIndex算法
        TreeIndex index = TreeIndex(graph);
        std::unordered_set<int> indexSolution = index.findKCoreSubgraph(queryNodes);
        int indexSize = indexSolution.size();
        int indexMinDegree = graph.computesubMinimumDegree(indexSolution);
        std::unordered_set<int> community = index.greedyConnection(queryNodes, indexMinDegree);

        // 将查询顶点集转换为字符串
        std::string queryNodesStr;
        queryNodesStr += "{";
        for (int node : queryNodes) {
            queryNodesStr += std::to_string(node) + " 、";
        }
        queryNodesStr += "}";
        if (!queryNodesStr.empty()) {
            queryNodesStr.pop_back(); // 移除最后一个逗号
        }

        // 将结果写入 CSV 文件
        outFile << num << "," 
                << queryNodesStr << ","
                << indexSize << ","
                << indexMinDegree << ","
                << community.size() << ","
                << graph.computesubMinimumDegree(community) << "\n";
    }
    
    outFile.close();
    // for (auto& q : group) {
    //     std::cout << "查询集: ";
    //     for (auto& i : q) { 
    //         std::cout << i << " ";
    //         }
    //     std::cout << std::endl;
    // }

    // SharingIndex dex= SharingIndex(graph);
    // dex.getKCoreToQuery(group, graph, "batch_4.csv");

}

int main(int argc, char *argv[])
{
    // 测试使用的无向图
    Graph graph("D:\\mySecre\\workspace\\csp_old\\CSP\\dataset\\facebook_combined.txt");

    // pro_1_1(graph, "single_1.csv");
    query_nodes querySet = {113, 14};
    // Graph greedySolution = graph.globalsearch(querySet);
    // int greedySize = greedySolution.getN();
    // int greedyMinDegree = greedySolution.getminimumDegree();
    // std::cout << "11: " << greedySize << "\n" << "22: " << greedyMinDegree << std::endl;
    
    std::vector<int> queryNodes = {113, 14};
    TreeIndex index = TreeIndex(graph);
    std::unordered_set<int> indexSolution = index.shellsearch(queryNodes);
    // int indexSize = indexSolution.size();
    // int indexMinDegree = graph.computesubMinimumDegree(indexSolution);
    // std::cout << "11: " << indexSize << "\n" << "22: " << indexMinDegree << std::endl;



    

 /*
    // Graph graph; // 假设你已经加载了图数据
    TreeIndex treeIndex(graph);

    std::vector<int> queryNodes = {1477 ,1101 ,3686 ,1758 ,2467}; // 查询节点集合 3581 | 1471, 2340 | 1, 6
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

