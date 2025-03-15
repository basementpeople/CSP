#ifndef SHARINGINDEX_H
#define SHARINGINDEX_H

#include <vector>
#include <map>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "TreeIndex.h"

// 自定义哈希函数，用于 std::vector<int>
struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        size_t hash = 0;
        for (const auto& item : v) {
            hash ^= std::hash<int>{}(item) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// 自定义相等性比较函数，用于 std::vector<int>
struct VectorEqual {
    bool operator()(const std::vector<int>& a, const std::vector<int>& b) const {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }
};

// 自定义哈希函数，用于 std::unordered_set<int>
struct SetHash {
    size_t operator()(const std::unordered_set<int>& s) const {
        size_t hash = 0;
        for (const auto& item : s) {
            hash ^= std::hash<int>{}(item) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// 自定义相等性比较函数，用于 std::unordered_set<int>
struct SetEqual {
    bool operator()(const std::unordered_set<int>& a, const std::unordered_set<int>& b) const {
        if (a.size() != b.size()) return false;
        for (const auto& item : a) {
            if (b.find(item) == b.end()) return false;
        }
        return true;
    }
};

// 自定义哈希函数，用于 int
struct IntHash {
    size_t operator()(const int& value) const {
        return std::hash<int>()(value);
    }
};

// 自定义相等性比较函数，用于 int
struct IntEqual {
    bool operator()(const int& a, const int& b) const {
        return a == b;
    }
};

// 自定义哈希函数，用于 std::pair<const int, std::unordered_set<int>>
struct PairEqual {
    bool operator()(const std::pair<const int, std::unordered_set<int>>& a,
                    const std::pair<const int, std::unordered_set<int>>& b) const {
        // 首先比较 pair 的第一个元素
        if (a.first != b.first) return false;
        // 如果第一个元素相等，比较第二个元素（std::unordered_set<int>）
        return a.second == b.second;
    }
};

// 自定义哈希函数，用于 std::pair<const std::unordered_set<int>, std::unordered_set<int>>
struct PairHash {
    size_t operator()(const std::pair<const std::unordered_set<int>, std::unordered_set<int>>& p) const {
        // 为第一个元素（std::unordered_set<int>）计算哈希值
        size_t h = 0;
        for (const int& elem : p.first) {
            h ^= std::hash<int>{}(elem) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        
        // 为第二个元素（std::unordered_set<int>）计算哈希值，并与第一个元素的哈希值组合
        for (const int& elem : p.second) {
            h ^= std::hash<int>{}(elem) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        
        return h;
    }
};


class SharingIndex : public TreeIndex {
    public:
        // SharingIndex(Graph &graph): TreeIndex(graph), communities(VectorHash{}, VectorEqual{}) {};
        SharingIndex(Graph &graph): TreeIndex(graph) {};

        // 1 找到公共子连通分量
        void getMinPortToQuery(const query_group &group); // 找到最大公共连通分量
        int findCommonChildren(const query_nodes &components); // 辅助函数，代价较大的方法

        void getKCoreToQuery(const std::vector<std::vector<int>> &group, Graph &graph, std::string path); // 找到最大K
        std::unordered_set<int> getLastResult(std::vector<int>& queryNodes, std::unordered_set<int>& result_end, Graph &graph);

        // 2 构建共享图
        void getSharingComponents(); // 获得共享图节点
        // checked 包含虚拟分量-1，可以修改返回值
        SharingIndex buildSubgraphIndex();  // 构建共享图

        // 3 查找k社区
        // 从顶层分量开始，查找子分量，并添加顶点到社区中
        void findCommunities(); // 找到多个社区

        // 添加一个新的公共方法来返回 queryToKcore 的引用
        // std::unordered_map<query_nodes, std::unordered_set<int>, SetHash, SetEqual>& getQueryToResult() {
        //     return queryToResult;
        // }
    
    private:
        std::vector<std::pair<int, query_nodes>> MinPortToQuery; // 最大公共连通分量ID 对应 查询顶点集
        std::unordered_set<int> subGraph; // 存储共享图节点的容器 共享图的节点是连通分量
        std::unordered_set<int> topComponents; // 顶层分量
        std::unordered_set<std::vector<int>, VectorHash, VectorEqual> communities; // 存储社区的容器
        std::unordered_map<int, std::unordered_set<int>, IntHash, IntEqual> componentToResult; // 连通分量ID 对应 结果集(连通分量ID集)
        // std::unordered_map<query_nodes, std::unordered_set<int>, SetHash, SetEqual> queryToResult; // 查询顶点集 对应 结果集(顶点集)
        // std::unordered_map<query_nodes, int, SetHash, SetEqual> queryToKcore; // 查询顶点集 对应 k

        // 1.16 修改，不用shell对应顶点集
        int q_count; // 查询顶点集个数 
        std::unordered_map<int, std::vector<int>> QueryCode;
        std::unordered_map<int, std::unordered_set<int>> KCoreToQuery;
        // std::vector<std::unordered_set<int>> KCoreToQuery; // 结构暂定
        std::unordered_map<int, std::unordered_set<int>> kToResult; // k 对应 结果集(连通分量ID集)
        std::unordered_map<int, std::unordered_set<int>> queryToResult; // 查询顶点集code 对应 结果集(顶点集)
        std::unordered_map<int, int> queryToKcore; // 查询顶点集code 对应 k


};

#endif // SHARINGINDEX_H