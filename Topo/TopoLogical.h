#ifndef TOPO_TOPOLOGICAL_H_
#define TOPO_TOPOLOGICAL_H_


#include <stack>
#include "DirectedCycle.h"
#include "DepthFirstOrder.h"
#include "../QuerySharingGraph.h"

class Query;

class TopoLogical {
private:
    std::stack<Query> order; // 顶点的拓扑排序

public:
    TopoLogical(QuerySharingGraph& G);

    // 判断图G是否有环
    bool isCycle() const;

    // 获取拓扑排序的所有顶点
    std::vector<Query> getorder();
};

#endif // TOPO_TOPOLOGICAL_H_