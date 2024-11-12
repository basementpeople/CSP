#ifndef QUERY_SHARING_GRAPH_H
#define QUERY_SHARING_GRAPH_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <cmath>
#include "Graph.h"
#include "Structures.h"
#include "QueryGroup.h"
#include "DirectedGraph.h"

// Topo排序
#include "Topo/TopoLogical.h"



class QuerySharingGraph : public DirectedGraph {
public:
    QuerySharingGraph();
    bool addEdge(const Query& q1, const Query& q2);

    bool isSubquery(Query& q1, Query& q2, DirectedGraph& G);

    void DetectCommonQuery(DirectedGraph& G, QueryGroup& Q, QuerySharingGraph& Psi);
    void DetectCommonQuery_r(DirectedGraph& G, QueryGroup& Q, QuerySharingGraph& Psi);

    std::vector<int> getNeighborsWithConstraint(int v, int k, DirectedGraph& G);

    // 打印
    void printEdges();

    // 获得每个点
    std::unordered_set<Query> getNodes();

    // 获得邻居
    std::vector<Query> getNeighbors(Query vertex) const;
    std::vector<Query> getInNeighbors(Query& vertex);

    // Batch算法
    void Search(DirectedGraph& G, std::vector<Path>& P, std::unordered_map<Query, std::vector<Path>>& R, Path& path, Query& q, QuerySharingGraph& Psi);
    void BatchEnum(DirectedGraph& G, std::vector<Query>& Q, double gamma);

private:
    std::set<std::pair<Query, Query>> edges;
    std::unordered_map<Query, std::unordered_set<Query>> adj; // 图的邻接表表示：查询和其相邻查询
    int n; // 图的节点数
    int m; // 图的边数
};

#endif // QUERY_SHARING_GRAPH_H