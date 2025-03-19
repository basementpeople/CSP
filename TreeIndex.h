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

class TreeIndex : public Graph {
public:
    // 第一个键是壳的层数，内层的键是连接组件ID，值是一个HashSet，包含连接组件中的节点ID
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> connectedComponentNodes;

    // 键是节点的ID，内层的键是壳的核心指数，值是一个HashSet，包含节点在该壳中的邻居节点ID。
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> nodeNeighbors;


public:
    // 1 构造函数和析构函数
    TreeIndex(Graph &graph); // 得到的索引不是树结构

    // 2 shell索引解决CSP
    std::unordered_set<int> shellsearch(std::vector<int> &queryNodes);

    std::unordered_set<int> findKCoreSubgraph(std::vector<int>& queryNodes, int k=-1);


    // 辅助函数，固定最后-----------------------------------
    void computeCoreIndex(Graph &graph); // 计算shell_count，coreIndex，coreMinimumDegree
    void identifyAndStoreComponents(Graph &graph); // 计算layerToComponentToNodes，ComponentToNodes，nodeToComponentId
    void buildParentChildRelationships(Graph &graph); // 计算connectedComponentParent，ComponentParent
    std::unordered_set<int> findTopComponents(std::vector<int> &queryNodes,int k); // 找到顶层分量，包括多余的顶层分量
    int findCommenShell(std::vector<int> queryNodes); // 找到查询顶点集的k值，即最小的shell层
    // void dfs(int com_1, int com_2, std::unordered_set<int>& visited, std::unordered_set<int>& result, int depth);
    void dfs(int com_1, int com_2, std::unordered_set<int>& visited, std::vector<int>& currentPath, std::vector<std::vector<int>>& allPaths);
    
    std::unordered_set<int> getNeighbors(int node);
    
    void printComponents() const;
    
    
    int findCommenK(std::vector<int>& queryNodes);
    int findCommenK_2(std::vector<int>& queryNodes);


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
    int n;

    // 用于生成唯一的连通分量 ID，只在建立索引时调用
    int nextComponentId = 0;

};


#endif // TREE_H