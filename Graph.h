#ifndef GRAPH_H
#define GRAPH_H

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

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
    Graph getGraph(std::unordered_set<int>& subVertices);
    Graph getReverseGraph();  //构造反向图
    Graph(Graph& graph);
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
    int computesubMinimumDegree(std::unordered_set<int>& nodes);
    void search(std::unordered_set<int> H0, int k, std::unordered_set<int>& H);
    std::unordered_set<int> baseline_search(int v0, int k);
    std::unordered_set<int> baseline_search2(int v0, int k);
    bool upperBound(int k);
    std::unordered_set<int> naiveCandidateGeneration(int v0, int k);
    // std::unordered_set<int> globalsearch(std::unordered_set<int> C, int k);
    std::unordered_set<int> CSTframework(int v0, int k);

    // CSM
    std::unordered_set<int> CSMframework(int v0, double gamma);
    std::unordered_set<int> generateCandidates(std::unordered_set<int>& H, int k);
    std::unordered_set<int> maxcore(std::unordered_set<int>& C, int v0);

    // Greedy算法
    bool isQuerySetConnected(const std::unordered_set<int>& queryNodes, const std::unordered_map<int, std::unordered_set<int>>& graph);
    void removeNode(int min_degree, std::unordered_map<int, std::unordered_set<int>>& resultGraph, std::vector<std::unordered_set<int>>& list);
    std::unordered_set<int> greedy(std::unordered_set<int>& queryNodes);

    /* HC-s-t 用到的函数 */

    void multiSourceBFS(Graph& graph, const std::vector<int> source, std::unordered_map<int, std::unordered_map<int, int>>& distIndex);
    void basicEnumSearch(std::vector<std::vector<int>>& P, std::vector<int>& p, int v, const int k, std::unordered_map<int, std::unordered_map<int, int>>& distIndex);
    void basicEnum(const Graph& graph, const std::vector<std::pair<int, int>>& query, const int k);

    void printP(std::vector<std::vector<int>>& P);
    void printGraph(); //打印图
    void printDistIndex(std::vector<int> source, std::unordered_map<int, std::unordered_map<int, int>>& distIndex); //打印distIndex

    /* 从文件中读取构建图 */
    void buildDigraphFromFile(std::string fileName);  //从文件中读取有向图
    void buildGraphFromAdjFile(std::string fileName); //从邻接表格式文件中读取构建无向图

private:

    std::unordered_map<int, std::unordered_set<int>> adj; // 图的邻接表表示：节点 ID 和其相邻节点

    std::unordered_map<int, int> degrees; // 图中节点的度

    std::vector<std::unordered_set<int>> orderedNodes; // 表示具有相同度的节点的集合的向量,key为degree

    int minimumDegree; // 图中所有节点的最小度

    int Dmax;

    int m;

    int n;

    void readFromFile(std::string fileName); // 从文件中读取无向图

};
#endif
