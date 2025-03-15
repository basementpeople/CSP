#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stack>
#include <map>
#include <queue>
#include <cmath>

#define query_nodes std::unordered_set<int>  // 查询顶点集
#define query_group std::vector<query_nodes>  // 多个查询顶点集

class Graph
{
public:

    // 1 构造函数和析构函数
    Graph(std::string path);
    Graph(Graph &graph);
    Graph();
    ~Graph(); 
    Graph getGraph(std::unordered_set<int> &subVertices); // 根据节点集构造图，不规范，后续改成构造函数

    // 2 Greedy算法
    void dfs(int current, std::unordered_set<int>& visited, const query_nodes& queryNodes, 
                int& count, const std::unordered_map<int, int>& degree, std::vector<int>& tree);
    std::vector<int> isQuerySetConnected(query_nodes queryNodes, std::unordered_map<int, int> degree);
    // std::vector<int> dfs(int current, std::unordered_set<int> visited, query_nodes queryNodes, int count, std::unordered_map<int, int> degree);
    // std::vector<int> isQuerySetConnected(query_nodes queryNodes, std::unordered_map<int, int> degree); // 检查查询集是否连通
    bool isConnected(query_nodes queryNodes, std::unordered_map<int, int> degree);
    // 移除节点
    void removeNode(int min_degree, std::unordered_map<int, std::unordered_set<int>>& resultGraph, std::vector<std::unordered_set<int>>& list);
    std::unordered_set<int> greedy(query_nodes& queryNodes); // 贪心算法
    std::unordered_set<int> greedy(int v0); // 单点贪心算法， 返回节点集或图
    std::unordered_set<int> greedy_d1(query_nodes& queryNodes); // 贪心算法,复用核心分解算法
    std::unordered_set<int> greedy_d2(query_nodes& queryNodes); // 1.7号重写
    std::unordered_set<int> greedy_d3(query_nodes& queryNodes); // 1.15号重写
    // 有问题，且目前只保留了basic算法，后续需要补充其他算法
    // 其它算法也只是优化输出的社区在某个大小内，我们的问题不需要这个约束参数

    
    // 3 Local Search 算法
    // Local Search ，但查询顶点H0需唯一
    void search(std::unordered_set<int> H0, int k, std::unordered_set<int>& H);
    std::unordered_set<int> baseline_search(int v0, int k); // 单点Local Search
    std::unordered_set<int> baseline_search2(int v0, int k);
    bool upperBound(int k);
    std::unordered_set<int> naiveCandidateGeneration(int v0, int k);
    std::unordered_set<int> globalsearch(std::unordered_set<int> C, int k);
    std::unordered_set<int> CSTframework(int v0, int k);
    
    // CSM local search
    std::unordered_set<int> CSMframework(int v0, double gamma);
    std::unordered_set<int> generateCandidates(std::unordered_set<int>& H, int k);
    std::unordered_set<int> maxcore(std::unordered_set<int>& C, int v0);
    
    // 4 batch index 算法
    // 计算两个查询之间的相似度
    double querySimilarity(const query_nodes& qA, const query_nodes& qB);
    // 计算两个查询顶点集之间的相似度
    double groupSimilarity(const query_group& groupA, const query_group& groupB);
    // 聚类算法
    std::vector<query_group> Clustering(std::vector<query_nodes>& query_groups, int k, double threshold);

    // 辅助函数，固定最后
    void addNode(int node); // 添加节点到图中
    bool addEdge(int from, int to); // 添加两个节点之间的边
    std::unordered_map<int, int> computeDegrees(); // 计算每个节点的度
    int computeMinimumDegree(); // 计算图中的最小度
    void statistic(); // 统计图的相关信息
    unsigned int getNumberOfNodes(); // 获取图中节点的数量
    std::vector<std::unordered_set<int>> getOrderedNodes(); // 获取图中度数相同的节点的集合
    std::unordered_map<int, int> getDegrees(); // 获取所有节点的度
    std::unordered_set<int> getNeighbors(int node); // 获取特定节点的邻居节点
    std::vector<int> sortNeighbors(int node); // 对特定节点的邻居节点进行排序
    std::unordered_set<int> getNodes(); // 获取图中所有节点
    void printAdj(); // 打印adj
    int findMaxDegreeNode(); // 找出度最大的节点，便于测试
    int computesubMinimumDegree(std::unordered_set<int> &nodes); // 计算子图中的最小度，只包含子图本身的节点
    std::unordered_map<int, std::unordered_set<int>>& getAdj() {
        return adj;
    }
    int getAdjlast() {
        return last;
    }
    int getCount() {
        return last;
    }
    std::unordered_set<int> getLastResult(query_nodes& queryNodes, std::unordered_set<int>& result_end);
    
protected:
    std::unordered_map<int, std::unordered_set<int>> adj; // 图的邻接表表示：节点 ID 和其相邻节点
    int last;

    std::unordered_map<int, int> degrees; // 图中节点的度

    std::vector<std::unordered_set<int>> orderedNodes; // 表示具有相同度的节点的集合的向量,key为degree
    // orderedNodes的数据结构有问题，没有标识度数的键值，需要重定义一个结构，vector<new<int>>

    int minimumDegree; // 图中所有节点的最小度
    
    int Dmax;

    int m;

    int n;

    int count_port;

    // 辅助函数，逻辑上不允许使用时调用
    void readFromFile(std::string fileName); // 从文件中读取图
};
#endif