#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <stdio.h>
#include <string.h>
#include <unordered_set>
#include <random>
#include <ctime>
#include <fstream>
#include <cassert>
// #include <filesystem> // C++17 for file existence check
// #include "Graph.h"
// #include "CoreGroup.h"
// #include "TreeIndex.h"
#include "DirectedGraph.h"
// #include "Structures.h"
#include "QuerySharingGraph.h"
// #include "QueryGroup.h"

int main(int argc, char *argv[])
{
    // 检查是否有输入文件路径
    // if (argc < 2)
    // {
    //     std::cerr << "Usage: " << argv[0] << " <path_to_dataset>" << std::endl;
    //     return 1;
    // }
    // // 获取命令行输入的文件路径
    // std::string dataset_path = argv[1];
    // std::string indexFile = dataset_path + ".index"; // 索引文件路径

//    Graph graph1("/ka/dataset/com-dblp.ungraph.txt"); // 图  CA-AstroPh
    //Graph graph1("C:/Users/86139/Desktop/index_community/output/com-dblp.ungraph.txt"); // 图 email-Eu-core.txt
    //std::cout << "图的最大度数的点：" << graph1.findMaxDegreeNode() << std::endl;
//    Graph graph1("C:/Users/86139/Desktop/index_community/output/example1.txt");
    //CoreGroup::coreGroupsAlgorithm(graph1);
    //std::vector<int> queryNodes={1,2,3}; // 目标集
//    std::unordered_set<int> queryNodes={1,2,3}; // 目标集
//    
//    std::unordered_set<int> H = {};
//    std::unordered_set<int> H0 = {6};
//    H = graph1.greedy(queryNodes); // 调用greedy算法
//    H = graph1.CSTframework(6, 5);
//    H = graph1.CSMframework(8, -2);
//   for (auto node : H) {
//       std::cout << node << " ";
//   }
//    std::cout <<graph1.computesubMinimumDegree(H) << " ";
    // graph1.computesubMinimumDegree(H);
    //  TreeIndex index = TreeIndex(graph1);  // 创建索引
    //  std::unordered_set<int> solution3 = index.findKCoreSubgraph(queryNodes);
    //  //std::unordered_set<int> solution3 = graph1.greedy(queryNodes);
    //  std::cout << "9" << std::endl;
    //  std::cout << solution3.size() << std::endl;
    //  int k = graph1.computesubMinimumDegree(solution3);
    //  std::cout << k << std::endl;
     // 打印找到的社区成员
    //  for (const auto& node : solution3) {
    //  	for (auto& i : node.second){
    //  	std::cout << "社区成员: " << i << std::endl;	
	// 	}
    //  }
    
    
    
    // Greedy算法
    // std::unordered_set<int> queryNodes={3};
//      int queryNodes= 38868;
//      std::unordered_set<int> solution1 = {};
//      solution1 = graph1.greedy(queryNodes);
//      std::cout <<graph1.computesubMinimumDegree(solution1) << " ";
    
    // GreedyDist算法
    
    
    // TreeIndex算法
//    std::vector<int> queryNodes={53213};
//    TreeIndex index = TreeIndex(graph1);
//    std::unordered_set<int> solution3 = index.findKCoreSubgraph(queryNodes);
//    std::cout << solution3.size() << std::endl;
//    std::cout <<graph1.computesubMinimumDegree(solution3) << " ";
//    for (auto& node : solution3) {
//      	std::cout << "社区成员: " << node << " ";	
//    }
    
    // search算法 localsearch的baseline算法
// 	int k = 20;
// 	std::unordered_set<int> queryNode = {3};
// 	std::unordered_set<int> solution4 = {};
// //    graph1.search(queryNode, k, solution4);
//     int q = 3;
//     solution4 = graph1.baseline_search2(q, k);
//     std::cout << solution4.size() << std::endl;
//     for (auto& node : solution4) {
//       	std::cout << "社区成员: " << node << " ";	
//     }
    
    // CST框架算法
//    int k = 28;
//    int queryNode = 3;
//    std::unordered_set<int> solution5 = {};
//    solution5 = graph1.CSTframework(queryNode, k);

	// CSM框架算法
//	int queryNode = 53213;
//	// 84424;
//	
//	int gamma = 2;
//	std::unordered_set<int> solution6 = {};
//	solution6 = graph1.CSMframework(queryNode, gamma);
//	std::cout << "最小度" <<graph1.computesubMinimumDegree(solution6) << std::endl;
//	for (auto& node : solution6) {
//     	std::cout << "社区成员: " << node << std::endl;	
//   }


    // 创建图对象并读取数据
    DirectedGraph graph("C:/Users/86139/Desktop/15_GreedyConnect/demo1/dataset/Batch_.txt");

    // 定义源点集合和目标点集合
   //std::unordered_set<int> S = {0, 2, 5, 4, 9}; // 源点集合
   //std::unordered_set<int> T = {11, 3, 12, 14}; // 目标点集合

//    // 初始化BFS索引
   //graph.initializeBFSIndex(S, T);

    // // 打印从源点1到各顶点的距离
    // for (const auto& [v, dist] : graph.getDistG().at(3)) {
    //     std::cout << "Distance from 3 to " << v << ": " << dist << std::endl;
    // }

    // // 打印从目标点5到各顶点的距离
    // for (const auto& [v, dist] : graph.getDistGr().at(5)) {
    //     std::cout << "Distance from 5 to " << v << ": " << dist << std::endl;
    // }

    // 定义查询列表
   std::vector<Query> queries = {
//        {1, 5, 19},  // 查询从1到5，最多3跳
//        {30, 1412, 2},  // 查询从2到6，最多4跳
//        {3, 5, 2}   // 查询从3到5，最多2跳
    {0, 11, 5},
    {2, 13, 5},
    {5, 12, 5},
    {4, 14, 4},
    {9, 14, 3}
   };

    // 调用BasicEnum方法
    // graph.BasicEnum(queries);
    

//    std::vector<Query> queries = {
//        {1, 5, 3},
//        {2, 6, 4},
//        {3, 7, 2},
//        {4, 8, 3}
//    };
//    graph.BasicEnum(queries);

   double threshold = 0.8;
   std::vector<QueryGroup> clusters = graph.hierarchicalClustering(queries, threshold);
   // graph.BatchEnum(graph, queries, 0.5);

//    for (size_t i = 0; i < clusters.size(); ++i) {
//        std::cout << "Cluster " << i + 1 << ": ";
//        for (const Query& q : clusters[i].queries) {
//            std::cout << "(" << q.s << ", " << q.t << ", " << q.k << ") ";
//        }
//        std::cout << std::endl;
//    }
    
    // // 创建查询
    // std::vector<Query> Q = {
    //     {1, 0, 3},
    //     {2, 0, 4},
    //     {3, 0, 2},
    //     {4, 0, 5}
    // };

    // 创建 QuerySharingGraph
    // QuerySharingGraph Psi;


    // 调用 DetectCommonQuery 函数，共享图构建
    // Psi.DetectCommonQuery(graph, queries, Psi);

    for (auto& cluster : clusters) {
        QuerySharingGraph Psi;
        // Psi.DetectCommonQuery(graph, cluster.queries, Psi);
        Psi.printEdges();
    }

    // 检查结果
    // 打印共享图结果
    // Psi.printEdges();

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
