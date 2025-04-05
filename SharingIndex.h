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

// 3.21 unordered_set
struct UnorderedSetHash {
    std::size_t operator()(const std::unordered_set<int>& set) const {
        std::size_t hash = 0;
        for (const int& elem : set) {
            hash ^= std::hash<int>()(elem) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

struct UnorderedSetEqual {
    bool operator()(const std::unordered_set<int>& a, const std::unordered_set<int>& b) const {
        return a == b;
    }
};

class SharingIndex : public TreeIndex {
public:
    // 1 构造函数同父类
    // SharingIndex(Graph &graph): TreeIndex(graph) {};
    SharingIndex(Graph &graph) : TreeIndex(graph) {};
    SharingIndex() {}; // 暂时不写，没用
    ~SharingIndex() {};

    // 2 batch查找解决CSP
    void batchsearch(query_group &group, std::string path);
    bool isConnected(query_nodes queryNodes, std::unordered_map<int, int> degree, std::unordered_set<int> &resultNodes);
    void searchstep(std::unordered_set<int> &q_, int k_min);

    // 3 与batch相关的 聚类算法，暂时没用到
    // 计算两个查询之间的相似度
    double querySimilarity(query_nodes& qA, query_nodes& qB);
    // 计算两个查询顶点集之间的相似度
    double groupSimilarity(query_group& groupA, query_group& groupB);
    // 聚类算法
    std::vector<query_group> Clustering(query_group& query_groups, double threshold);
    // CP聚类，目的在于划分不同连通分量（总），但事实上如果一个无向图总是连通的。
    // 判断的方法，通过找到顶层分量集，比较得到CP(qA,qB)，1为连通，0为不连通



    // 3 辅助函数
    std::unordered_set<int> getLastResult(std::vector<int>& queryNodes, std::unordered_set<int>& result_end, Graph &graph); // 确保最后的结果连通，没有必要，暂时不用
    void printAndwrite(std::string path); // 写入文件
    void printAndwrite_(std::string path); // 写入文件

    // 4 batch查找解决MIN_CSP
    void batchMinsearch(query_group& group); // 返回Graph集合


protected:
//CSP
    int q_count; // 查询顶点集个数 
    std::unordered_map<int, query_nodes> QueryCode; // 查询顶点集code 对应 查询顶点集
    std::unordered_map<int, std::unordered_set<int>> KCoreToQuery; // k 对应 查询顶点集codes
    std::unordered_map<int, std::unordered_set<int>> kToResult; // k 对应 结果集(连通分量ID集)
    std::unordered_map<int, std::unordered_set<int>> queryToResult; // 查询顶点集code 对应 结果集(顶点集)
    std::unordered_map<int, int> queryToKcore; // 查询顶点集code 对应 k

//MIN_CSP
    int q_count_;
    std::unordered_map<int, query_nodes> QueryCode_;
    std::unordered_map<int, int> queryToKcore_;
    std::unordered_map<int, std::unordered_set<int>> queryToResult_;

    std::vector<double> time;
    std::vector<double> time_sum;
    std::vector<bool> flag_;
    std::vector<int> species;

};

#endif // SHARINGINDEX_H