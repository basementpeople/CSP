// Copyright 2023 Your Name

#ifndef TOPO_DIRECTEDCYCLE_H_
#define TOPO_DIRECTEDCYCLE_H_

#include <vector>
#include <stack>
#include "TopoLogical.h"
#include "../QuerySharingGraph.h"

class DirectedCycle {
public:
    // 创建一个检测环对象，检测图G中是否有环
    DirectedCycle(QuerySharingGraph G);

    // 判断当前有向图G中是否有环
    bool gethasCycle();

private:
    std::unordered_map<Query, bool> marked; // 索引代表顶点，值表示当前顶点是否已经被搜索
    bool hasCycle; // 记录图中是否有环
    std::unordered_map<Query, bool> onStack; // 索引代表顶点，使用栈的思想，记录当前顶点有没有已经处于正在搜索的栈上，如果有，则证明有环。

    // 基于深度优先搜索，检测图G中是否有环
    void dfs(QuerySharingGraph G, Query v);
};

#endif  // TOPO_DIRECTEDCYCLE_H_