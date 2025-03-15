#include "DepthFirstOrder.h"
DepthFirstOrder::DepthFirstOrder(QuerySharingGraph& G) {
    // 初始化辅助变量
    // marked.resize(G.getNumberOfNodes(), false); // 默认全部赋值为false
    // 对每一个未搜索过的顶点进行深度优先遍历
    for (auto v : G.getNodes()) {
        if (!marked[v])
            dfs(G, v);
    }
}

// 基于深度优先搜索，生成顶点线性序列
void DepthFirstOrder::dfs(QuerySharingGraph& G, Query v) {
    // 1. 将当前顶点标记为已搜索
    marked[v] = true;
    // 2. 遍历当前顶点的邻接表，对邻接表中未搜索的顶点递归调用深度优先搜索
    for (auto w : G.getNeighbors(v)) {
        if (!marked[w])
            dfs(G, w);
    }
    // 3. 当前顶点v深度优先搜索完毕后，入栈
    reversePost.push(v);
}

// 获取顶点线性序列
std::stack<Query> DepthFirstOrder::getreversePost() const {
        return reversePost;
    }