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

class Graph
{
public:

    // 1 构造函数和析构函数
    Graph(std::string path);
    Graph();
    Graph getGraph(std::unordered_set<int> &subVertices);
    Graph(Graph &graph);
    ~Graph(); 
    
    void addNode(int node); // 添加节点到图中
    bool addEdge(int from, int to); // 添加两个节点之间的边
    std::unordered_map<int, int> computeDegrees(); // 计算每个节点的度
    int computeMinimumDegree(); // 计算图中的最小度
    void statistic();
    
    // 2 CoreGroup中直接使用的函数
    int getNumberOfNodes(); // 获取图中节点的数量
    std::unordered_map<int, int> getDegrees(); // 获取所有节点的度
    std::unordered_set<int> getNeighbors(int node); // 获取特定节点的邻居节点

    // 3 CoreGroup中间接使用的函数
    // std::unordered_map<int, std::unordered_set<int>> adj; 图的邻接表表示：节点 ID 和其相邻节点
    // std::unordered_map<int, int> degrees; 图中节点的度
    std::vector<int> sortNeighbors(int node);

    // 4 TreeIndex中直接使用的函数
    std::unordered_set<int> getNodes(); // 获取图中所有节点
    
    // 5 
    int computesubMinimumDegree(std::unordered_set<int> &nodes);
    void search(std::unordered_set<int> H0, int k, std::unordered_set<int>& H);
    std::unordered_set<int> baseline_search(int v0, int k);
    std::unordered_set<int> baseline_search2(int v0, int k);
    bool upperBound(int k);
    std::unordered_set<int> naiveCandidateGeneration(int v0, int k);
    std::unordered_set<int> globalsearch(std::unordered_set<int> C, int k);
    std::unordered_set<int> CSTframework(int v0, int k);
    
    // CSM
    std::unordered_set<int> CSMframework(int v0, double gamma);
    std::unordered_set<int> generateCandidates(std::unordered_set<int>& H, int k);
    std::unordered_set<int> maxcore(std::unordered_set<int>& C, int v0);
    
    // Greedy算法
    bool isQuerySetConnected(const std::unordered_set<int>& queryNodes, const std::unordered_map<int, std::unordered_set<int>>& graph);
    void removeNode(int min_degree, std::unordered_map<int, std::unordered_set<int>>& resultGraph, std::vector<std::unordered_set<int>>& list);
    std::unordered_set<int> greedy(std::unordered_set<int>& queryNodes);
    std::unordered_set<int> greedy(int v0);
    
private:
    std::unordered_map<int, std::unordered_set<int>> adj; // 图的邻接表表示：节点 ID 和其相邻节点

    std::unordered_map<int, int> degrees; // 图中节点的度

    std::vector<std::unordered_set<int>> orderedNodes; // 表示具有相同度的节点的集合的向量,key为degree
    
    int minimumDegree; // 图中所有节点的最小度
    
    int Dmax;

    int m;

    int n;
    
    void readFromFile(std::string fileName); // 从文件中读取图
    
};
#endif