#ifndef TOPO_DEPTHFIRSTORDER_H_
#define TOPO_DEPTHFIRSTORDER_H_

#include <stack>
#include <vector>
#include "TopoLogical.h"
#include "../QuerySharingGraph.h"

class DepthFirstOrder {
private:
    std::unordered_map<Query, bool> marked; // 索引代表顶点，值表示当前顶点是否已经被搜索
    std::stack<Query> reversePost; // 使用栈，存储顶点序列，打印出栈中的顶点，即是排序后的顶点

public:
    DepthFirstOrder(QuerySharingGraph& G);

    // 基于深度优先搜索，生成顶点线性序列
    void dfs(QuerySharingGraph& G, Query v);

    // 获取顶点线性序列
    std::stack<Query> getreversePost() const;
};

#endif // TOPO_DEPTHFIRSTORDER_H_