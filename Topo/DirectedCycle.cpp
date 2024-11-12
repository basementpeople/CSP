#include "DirectedCycle.h"

DirectedCycle::DirectedCycle(QuerySharingGraph G)
{
    // marked.resize(G.getNumberOfNodes(),false);  // 用于标识顶点是否搜索过
    hasCycle = false;
    // onStack.resize(G.getNumberOfNodes()); // 用于标识顶点是否在搜索中

    // 遍历所有顶点，将未搜索过的顶点作为入口，进行深度优先遍历，检测是否有环，一旦检测到有环，则结束；
    // 因为对于不连通图，有很多个子图，也许某个子图存在环，因此，要对每个子图进行深度优先遍历检测，而不能只检测某一个子图。
    for (auto v : G.getNodes()) {
        if (!marked[v])
            dfs(G, v); // 每次搜索一个子图，判断子图内是否有环，如果没环，继续搜索下一个子图（一次搜索后，未搜索的顶点一定在另一个子图中）
    }
}

// 判断当前有向图G中是否有环
bool DirectedCycle::gethasCycle()
{
    return hasCycle;
}

// 基于深度优先搜索，检测图G中是否有环
void DirectedCycle::dfs(QuerySharingGraph G, Query v)
{
    // 1. 当前顶点标记为已搜索
    marked[v] = true;
    // 2. 当前顶点入栈
    onStack[v] = true;
    // 3. 递归深度优先遍历，检查遍历结点是否已经在栈中，如果在栈中，则表明该顶点被两次搜索到，证明有环，则结束
    for (auto w : G.getNeighbors(v)) {
        if (!marked[w])
            dfs(G, w);
        // 如果该顶点已经被搜索过，且如果该顶点在搜索的路径上，则代表又一次搜索到该顶点，证明有环，结束搜索。
        if (onStack[w]) {
            hasCycle = true;
            return;
        }
    }
    // 4. 当前顶点出栈，为下一个节点作为入口，检测是否有环做准备(为什么需要这样，图2.png可以解释)
    onStack[v] = false;
}