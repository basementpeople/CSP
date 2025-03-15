#include "TopoLogical.h"
TopoLogical::TopoLogical(QuerySharingGraph& G) {
    // 1. 检测是否有环
    DirectedCycle directedCycle(G);
    if (!directedCycle.gethasCycle()) {
        // 2. 调用顶点排序算法
        DepthFirstOrder depthFirstOrder(G);
        this->order = depthFirstOrder.getreversePost();
    }
}

// 判断图G是否有环
bool TopoLogical::isCycle() const {
    return order.empty();
}

// 获取拓扑排序的所有顶点
std::vector<Query> TopoLogical::getorder() {
    std::vector<Query> result;
    std::stack<Query> temp = order;
    while (!temp.empty()) {
        result.push_back(temp.top());
        temp.pop();
    }
    return result;
}