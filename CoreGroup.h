#ifndef COREGROUP_H
#define COREGROUP_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include "Graph.h"

class CoreGroup
{
public:
// 返回节点的核索引（所处的最高k-core）
    static std::unordered_map<int, int> coreGroupsAlgorithm(Graph &graph)
    {
    	// std::cout << "正常进入函数 " << std::endl;
        int maxDegree = 0;                  // 最大度数
        std::unordered_map<int, int> cores; // 存储每个节点的核心度
        std::vector<std::unordered_set<int>> orderedNodes(graph.getNumberOfNodes());

        // 初始化orderedNodes
        for (int i = 0; i < graph.getNumberOfNodes(); ++i)
        {
            orderedNodes[i] = std::unordered_set<int>();
        }

        // 将节点按度数分类
        std::unordered_map<int, int> degrees = graph.getDegrees();
        for (const auto &entry : degrees)
        {
            orderedNodes[entry.second].insert(entry.first);
        }

        int node;             // 当前处理的节点
        int lowestDegree = 0; // 当前的最低度数
        int neighborDegree;   // 邻居的度数

        // 主循环
        while (lowestDegree < graph.getNumberOfNodes())
        {
            if (orderedNodes[lowestDegree].empty())
            {
                ++lowestDegree; // 如果没有该度数的节点，增加度数
            }
            else
            {
                node = *orderedNodes[lowestDegree].begin();
                orderedNodes[lowestDegree].erase(node);
                cores[node] = lowestDegree; // 设置核心度
                // std::cout << "核心度" << ": " << cores[node] << std::endl;
                degrees[node] = -1;         // 将节点度数设为-1，标记已处理

                // 更新所有邻居的度数
                for (int neighbor : graph.getNeighbors(node))
                {
                    neighborDegree = degrees[neighbor];
                    if (neighborDegree > lowestDegree)
                    {
                        orderedNodes[neighborDegree].erase(neighbor);
                        orderedNodes[neighborDegree - 1].insert(neighbor);
                        degrees[neighbor] = neighborDegree - 1;
                    }
                }
            }
        }
        //std::cout << "the cores's size " << cores.size() << std::endl;
        //std::cout << "the node's size " << graph.getNumberOfNodes() << std::endl;
        return cores;
    }
};
#endif