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
#include <filesystem> // C++17 for file existence check
#include "Graph.h"
#include "CoreGroup.h"
#include "TreeIndex.h"

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

    // Graph graph1("/ka/dataset/CA-AstroPh.txt"); // 图  CA-AstroPh
    Graph graph1("C:/Users/86139/Desktop/index_community/output/CA-AstroPh.txt"); // 图 email-Eu-core.txt
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
//    std::unordered_set<int> queryNodes={3};
//    std::unordered_set<int> solution1 = {};
//    solution1 = graph1.greedy(queryNodes);
//    std::cout <<graph1.computesubMinimumDegree(solution1) << " ";
    
    // GreedyDist算法
    
    
    // TreeIndex算法
//    std::vector<int> queryNodes={84424};
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
	int queryNode = 84424;
	int gamma = 4;
	std::unordered_set<int> solution6 = {};
	solution6 = graph1.CSMframework(queryNode, gamma);
	std::cout << "最小度" <<graph1.computesubMinimumDegree(solution6) << std::endl;
//	for (auto& node : solution6) {
//     	std::cout << "社区成员: " << node << std::endl;	
//   }
    
    return 0;
}
