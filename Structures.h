#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>

// 查询结构体
struct Query {
    int s;  // 源顶点
    int t;  // 目标顶点
    int k;  // 跳数限制

    // 默认构造函数
    Query() : s(0), t(0), k(0) {}

    // 带参数的构造函数
    Query(int s, int t, int k) : s(s), t(t), k(k) {}

    // 重载 < 运算符
    bool operator<(const Query& other) const {
        if (s != other.s) return s < other.s;
        if (t != other.t) return t < other.t;
        return k < other.k;
    }

    // 重载 == 运算符
    bool operator==(const Query& other) const {
        return s == other.s && t == other.t && k == other.k;
    }

    // 重载 != 运算符
    bool operator!=(const Query& other) const {
        return s != other.s || t != other.t || k != other.k;
    }
    
    // 析构函数（如果有动态分配的资源）
    ~Query() = default;  // 使用默认析构函数
};

// 为 Query 类型提供 std::hash 特化
namespace std {
    template<>
    struct hash<Query> {
        size_t operator()(const Query& query) const {
            // 使用位运算和哈希组合来生成哈希值
            return (hash<int>()(query.s) ^ (hash<int>()(query.t) << 1)) ^ (hash<int>()(query.k) << 2);
        }
    };
}

// 路径结构体
struct Path {
    std::vector<int> vertices;  // 路径上的顶点序列
};

// 表示查询的前向和后向可达顶点集
struct HopConstrainedNeighbors {
    std::unordered_set<int> forwardReachable;
    std::unordered_set<int> backwardReachable;
};

#endif // STRUCTURES_H