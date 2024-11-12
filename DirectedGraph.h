#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <thread>
#include <mutex>


#include "Graph.h"
#include "Structures.h"
#include "QueryGroup.h"
// #include "QuerySharingGraph.h"


// 拓扑排序
// #include "Topo/TopoLogical.h"



class DirectedGraph : public Graph {
public:
    // 构造函数
    DirectedGraph(std::string path);

    DirectedGraph();

    // 从文件读取有向图
    void readFromFile(std::string fileName);

    // 添加边
    bool addEdge(int from, int to);

    // 获取所有顶点
    std::vector<int> dfsTraversal(int vertex, std::unordered_set<int>& visited);
    std::vector<int> traverseAllVertices();


    // 获取反向图
    DirectedGraph getReverseGraph() const;
    
    // 辅助聚类查询，构建共享图
    int dist (int s1, int s2);

    // 多源BFS
    void initializeBFSIndex(const std::unordered_set<int>& S, const std::unordered_set<int>& T);
    
    // 优化
    void parallelInitializeBFSIndex(const std::unordered_set<int>& S, const std::unordered_set<int>& T);

    void BasicEnum(const std::vector<Query>& queries); // 1

    std::vector<int> getNeighbors(int vertex) const; // 1

    std::vector<int> getReverseNeighbors(int vertex) const; // 1

    // 路径是否存在，未用到
    bool isPathPossible(int s, int t, int maxDistance);

    bool isSimplePath(const Path& path);
    void connectPaths(const std::vector<Path>& forwardPaths, const std::vector<Path>& backwardPaths, std::vector<Path>& resultPaths);

    // 调用distG
    const std::unordered_map<int, std::unordered_map<int, int>>& getDistG() const;

    // 调用distGr
    const std::unordered_map<int, std::unordered_map<int, int>>& getDistGr() const;


    // 正向搜索
    void forwardSearch(const DirectedGraph& graph, int s, int t, int k, int K, std::vector<Path>& paths, Path& currentPath);
    // 反向搜索
    void backwardSearch(const DirectedGraph& graph, int s, int t, int k, int K, std::vector<Path>& paths, Path& currentPath);
    
    // 聚类算法
    HopConstrainedNeighbors getHopConstrainedNeighbors(int s, int t, int k); // 计算 Hop-Constrained Neighbors
    // 计算两个查询之间的相似度
    double querySimilarity(const Query& qA, const Query& qB);
    // 计算两个查询组之间的相似度
    double groupSimilarity(const QueryGroup& groupA, const QueryGroup& groupB);

    // 层次聚类算法
    std::vector<QueryGroup> hierarchicalClustering(const std::vector<Query>& queries, double threshold);

    // Batch
    // void Search(DirectedGraph& G, std::vector<Path>& P, std::unordered_map<Query, std::vector<Path>>& R, Path& path, Query& q, QuerySharingGraph& Psi);
    // void BatchEnum(DirectedGraph& G, std::vector<Query>& Q, double gamma);
    
private:
    std::unordered_map<int, std::unordered_map<int, int>> distG; // 从源点到各顶点的距离
    std::unordered_map<int, std::unordered_map<int, int>> distGr; // 从目标点到各顶点的距离

    // BFS实现
    void bfsFromSource(int source, std::unordered_map<int, int>& distances);
};

#endif // DIRECTEDGRAPH_H