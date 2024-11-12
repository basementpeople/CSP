#ifndef QUERYGROUP_H
#define QUERYGROUP_H

#include <string>
#include <vector>
#include "Structures.h"

class QueryGroup {
public:
    std::vector<Query> queries;
    HopConstrainedNeighbors hcn;

    // 构造函数
    QueryGroup() = default;

    // 添加查询和HopConstrainedNeighbors
    void addQuery(const Query& q, const HopConstrainedNeighbors& hcn) {
        queries.push_back(q);
        this->hcn.forwardReachable.insert(hcn.forwardReachable.begin(), hcn.forwardReachable.end());
        this->hcn.backwardReachable.insert(hcn.backwardReachable.begin(), hcn.backwardReachable.end());
    }

    // 获取所有查询
    const std::vector<Query>& getQueries() const {
        return queries;
    }

    // 获取HopConstrainedNeighbors
    const HopConstrainedNeighbors& getHCN() const {
        return hcn;
    }

    // 清空查询列表
    void clearQueries() {
        queries.clear();
    }

    // 清空HopConstrainedNeighbors
    void clearHCN() {
        hcn.forwardReachable.clear();
        hcn.backwardReachable.clear();
    }

    // 检查是否包含特定查询
    bool containsQuery(const Query& q) const {
        for (const auto& query : queries) {
            if (query == q) {
                return true;
            }
        }
        return false;
    }

    // 获取查询数量
    size_t getQueryCount() const {
        return queries.size();
    }
};

#endif // QUERYGROUP_H