#ifndef TREEINDEX_H
#define TREEINDEX_H

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <stack>
#include <queue>
#include <set>
#include <map>
#include <climits>
#include <chrono>
#include "Graph.h"
#include "CoreGroup.h"

class TreeIndex 
{
public:
    // 第一个键是壳的层数，内层的键是连接组件ID，值是一个HashSet，包含连接组件中的节点ID
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> connectedComponentNodes;

    // 键是节点的ID，内层的键是壳的核心指数，值是一个HashSet，包含节点在该壳中的邻居节点ID。
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> nodeNeighbors;

    
    int nextComponentId = 0; // 用于生成唯一的连通分量 ID

public:
    TreeIndex(Graph &graph);

    // 建立shellstruct的辅助函数
    void computeCoreIndex(Graph &graph);
    void identifyAndStoreComponents(Graph &graph);
    void buildParentChildRelationships(Graph &graph);
    void setNeighbors(Graph &graph);
    void setDegrees(Graph &graph);

    std::unordered_set<int> getNeighbors(int node);
    
    void printComponents() const;
    
    
    std::unordered_set<int> findTopComponents(std::vector<int> &queryNodes,int k);
    
    std::unordered_set<int> findKCoreSubgraph(std::vector<int>& queryNodes, int k=-1);

    std::unordered_set<int> findKCoreSubgraph_d1(std::vector<int>& queryNodes); // 1.15
    
    int findCommenShell(std::vector<int>& queryNodes);

    // 构造函数，根据图构建树索引
    TreeIndex(Graph &graph, std::string datasetName);

    // 构造函数，从文件中加载树索引
    TreeIndex(const std::string &path);

    // 计算核心组成
    void computeCoreComposition(Graph &graph);

    void computeCoreCompositionByLcy(Graph &graph);

    // 计算节点邻居信息
    void computeNodeNeighbors(Graph &graph);

    // 序列化树索引到文件
    void serializeToFile(const std::string &path);

    // 从文件中反序列化树索引
    void deserializeFromFile(const std::string &path);

    // 获取一组节点中的最小核心编号
    int getMinimumCoreIndex(std::vector<int> queryNodes);

    // 获取节点的邻居节点集合（在指定核心下）
    std::unordered_set<int> getNeighbors(int node, int coreIndex);

    // 获取节点的邻居节点集合（在指定核心下，并且在指定子核心中）
    std::unordered_set<int> getNeighbors(int node, int coreIndex, std::unordered_set<int> &subcore);

    // 获取指定核心的最小度数
    int getCoreMinimumDegree(int coreIndex);

    // 获取指定核心中的节点数量
    int getNumberOfNodes(int coreIndex);

    // 获取一组节点的核心集合
    std::unordered_set<int> getCore(std::vector<int> queryNodes);

    // 获取节点的核心编号
    int getCoreIndex(int node);

    // Get the parent component ID for a given component
    int getParentComponent(int componentId);
    std::unordered_set<int> getConnectedComponentChildren(int componentId);

    // Get all nodes in a specified component
    std::unordered_set<int> getNodesInComponent(int componentId);
    // 打印核心指标
    void printCoreIndex();

    // 打印每个核心的最小度数
    void printCoreMinimumDegree();

    // 打印连通分量的父子关系
    void printConnectedComponentParent();

    // 打印连通分量的孩子节点集合
    void printConnectedComponentChildren();

    // 打印连通分量的节点集合
    void printConnectedComponentNodes();
    // 打印节点分组
    void printNodeGroup();

    // 打印节点邻居信息
    void printNodeNeighbors();
    void printQueryNodesCoreIndex(const std::vector<int>& queryNodes);

// 修改index算法--------------------------------------------------------------------------------    
    // 找到公共子连通分量
    int getKfromIndex(const std::vector<int>& group); // 找到最大公共连通分量

    // 测试函数，查找某个shell的所有子shell
    void findSubShells(int node, int wander);

    // 辅助函数

    std::unordered_set<int> getChildShell(int shell);
    std::unordered_set<int> getParentShell(int shell);

// 整理TreeIndex文件----------------------------------------------------------------------------
    // 得到 点 对应的 连通分量ID  weijiancha
    int getComponent(int node);

    // 获得 连通分量 对应的 核心度  weijiancha
    int getCoreFromComponent(int componentId);

    // 一个新的构造函数，根据连通分量集构建图
    TreeIndex(TreeIndex &treeIndex, std::unordered_set<int> &subGraph);

    std::unordered_set<int> greedyStep(std::vector<int> &queryNodes, int k);
    std::unordered_set<int> greedyStep_simply(std::vector<int> &queryNodes, int k, std::unordered_set<int>& realnodes);
    std::unordered_set<int> connectionStep(std::unordered_set<int> &H_min, std::vector<int> &queryNodes, int k);
    std::unordered_set<int> steinerTreeApproximation(std::unordered_set<int> &H_min, std::vector<int> &queryNodes);
    std::unordered_set<int> greedyConnection(std::vector<int> &queryNodes, int k);
    std::vector<int> steinerTree(std::unordered_set<int>& H_min_star, const std::vector<int>& terminals);


public:
    // 点 对应 连通分量ID
    // findTopComponents、
    std::unordered_map<int, int> nodeToComponentId;

    // 点 对应 核心ID
    // findKCoreSubgraph
    std::unordered_map<int, int> coreIndex;

    // 核心ID 对应 核心度
    // findKCoreSubgraph
    std::unordered_map<int, int> coreMinimumDegree;

    // 连通分量ID 对应 点（集合）
    // findKCoreSubgraph
    std::unordered_map<int, std::unordered_set<int>> ComponentToNodes;

    // 核心ID（对应k-shell） 对应 连通分量ID 对应 连通分量
    // buildParentChildRelationships
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> layerToComponentToNodes;

    // 核心ID（对应k-shell） 对应 连通分量ID 对应 父连通分量ID
    // buildParentChildRelationships
    std::unordered_map<int, std::unordered_map<int, int>> connectedComponentParent;

    // 连通分量ID 对应 父连通分量ID
    // buildParentChildRelationships
    std::unordered_map<int, std::unordered_set<int>> ComponentParent;

    // 核心ID（对应k-shell） 对应 连通分量ID 对应 子连通分量ID（集合）
    // buildParentChildRelationships
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> connectedComponentChildren;

    // 连通分量ID 对应 子连通分量ID
    // buildParentChildRelationships
    std::unordered_map<int, std::unordered_set<int>> ComponentChildren;

    // 节点分组，记录每个节点所属的connectedComponentNodes分组编号
    std::unordered_map<int, int> nodeGroup;

    // 记录分量ID
    std::unordered_set<int> visitedComponents;

    // 记录shell层数， 更准确一些核心度的最大值
    int shell_count;

    // 全局邻居
    std::unordered_map<int, std::unordered_set<int>> adj;

    // 图中节点的度
    std::unordered_map<int, int> degrees;

    // 图中节点数
    int count_port;

};


#endif // TREE_H