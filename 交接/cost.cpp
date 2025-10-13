#include "cost.h"
#include <queue>
#include <fstream>
#include <iostream>
#include "blocked_area_analyzer.h"  

CostUpdater::CostUpdater(Grid& grid) : grid(grid) {}

void CostUpdater::updateCost() {
    int directions[4][2] = { {-1, 0}, {1, 0}, {0, -1}, {0, 1} };//方向

    int currentCost = 5;
    //for (int layer = 0; layer < grid.Layers; ++layer) {
    for (int j = 111; j < 274; ++j) {//边界条件
        for (int k = 83; k < 246; ++k) {
            if (grid.grid[1][j][k].isblocked == 1) {
                grid.grid[1][j][k].cost = 5;
            }
        }
    }
    //}
    for (int iteration = 0; iteration < 5; ++iteration) {
        //for (int layer = 0; layer < grid.Layers; ++layer) {
        for (int j = 111; j < 274; ++j) {//边界条件
            for (int k = 83; k < 246; ++k) {

                if (grid.grid[1][j][k].cost == currentCost) {
                    for (const auto& dir : directions) {
                        int newJ = j + dir[0];
                        int newK = k + dir[1];

                        if (newJ >= 0 && newJ < 275 && newK >= 0 && newK < 247) {// 记得改边界条件
                            if (grid.grid[1][newJ][newK].cost == 0) {
                                grid.grid[1][newJ][newK].cost = currentCost - 1;
                            }
                        }
                    }
                }
            }
        }
        //}
        currentCost--;
    }
}

#include <queue>

//void CostUpdater::updateCost() {
//    const int directions[4][2] = { {-1, 0}, {1, 0}, {0, -1}, {0, 1} };
//
//    const int maxCost = 5;
//    const int unitStep = 30; // 每10格降低一个代价
//    const int maxDistance = maxCost * unitStep;
//
//    int maxJ = grid.grid[1].size();
//    int maxK = grid.grid[1][0].size();
//
//    std::queue<std::pair<int, int>> q;
//
//    // Step 1：初始化障碍物为5，加入BFS队列
//    for (int j = 111; j < 274; ++j) {
//        for (int k = 83; k < 246; ++k) {
//            if (grid.grid[1][j][k].isblocked == 1) {
//                grid.grid[1][j][k].cost = maxCost;
//                q.emplace(j, k); // 将障碍点加入队列
//            }
//            else {
//                grid.grid[1][j][k].cost = 0; // 清空非障碍点
//            }
//        }
//    }
//
//    // Step 2：BFS 扩散，记录每个点到障碍点的最小曼哈顿距离
//    std::vector<std::vector<int>> visited(maxJ, std::vector<int>(maxK, -1));
//
//    while (!q.empty()) {
//        auto [j, k] = q.front(); q.pop();
//
//        for (const auto& dir : directions) {
//            int newJ = j + dir[0];
//            int newK = k + dir[1];
//
//            if (newJ >= 0 && newJ < maxJ && newK >= 0 && newK < maxK) {
//                // 如果还没访问过且不是障碍点
//                if (grid.grid[1][newJ][newK].isblocked == 0 && visited[newJ][newK] == -1) {
//                    int dist = std::abs(newJ - j) + std::abs(newK - k);
//                    visited[newJ][newK] = visited[j][k] + 1;
//
//                    int manhattanDist = visited[j][k] + 1;
//                    if (manhattanDist <= maxDistance) {
//                        int cost = maxCost - (manhattanDist / unitStep);
//                        if (grid.grid[1][newJ][newK].cost < cost) {
//                            grid.grid[1][newJ][newK].cost = cost;
//                            q.emplace(newJ, newK);
//                        }
//                    }
//                }
//            }
//        }
//    }
//}



//void CostUpdater::updateCost() {
//    const int cost_blocked = 5000;
//    const int cost_horizontal = 400;
//    const int horizontal_distance = 60;
//
//    int maxJ = grid.grid[1].size();            // 行数
//    int maxK = grid.grid[1][0].size();         // 列数
//
//    // Step 1: 初始化所有 cost 为 0
//    for (int j = 0; j < maxJ; ++j) {
//        for (int k = 0; k < maxK; ++k) {
//            grid.grid[1][j][k].cost = 0;
//        }
//    }
//
//    // Step 2: 遍历所有格子，找阻塞点，处理其左右60格
//    for (int j = 0; j < maxJ; ++j) {
//        for (int k = 0; k < maxK; ++k) {
//            if (grid.grid[1][j][k].isblocked == 1) {
//                // 设置自身为 cost = 5
//                grid.grid[1][j][k].cost = cost_blocked;
//
//                // 向左扩展最多 60 格
//                for (int offset = 1; offset <= horizontal_distance; ++offset) {
//                    int left = k - offset;
//                    if (left < 0) break;
//
//                    if (grid.grid[1][j][left].isblocked != 1 && grid.grid[1][j][left].cost < cost_horizontal) {
//                        grid.grid[1][j][left].cost = cost_horizontal;
//                    }
//                }
//
//                // 向右扩展最多 60 格
//                for (int offset = 1; offset <= horizontal_distance; ++offset) {
//                    int right = k + offset;
//                    if (right >= maxK) break;
//
//                    if (grid.grid[1][j][right].isblocked != 1 && grid.grid[1][j][right].cost < cost_horizontal) {
//                        grid.grid[1][j][right].cost = cost_horizontal;
//                    }
//                }
//            }
//        }
//    }
//}

void CostUpdater::outputWeightsToFile(const std::string& outputFile) {

    std::ofstream outFile(outputFile);

    if (!outFile.is_open()) {
        std::cerr << "无法打开文件 " << outputFile << " 进行写入\n";
        return;
    }

    //for (int i = 0; i <grid.Layers; i++){
    for (int j = 111; j < 274; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
        for (int k = 83; k < 246; k++) {
            outFile << grid.grid[1][j][k].cost << " ";
            outFile << "\n";
        }
    }

    for (int i = 111; i < 274; i++) {
        for (int k = 83; k < 246; k++) {
            outFile << 0 << " ";
            outFile << "\n";
        }
    }
    //}

    // 关闭文件
    outFile.close();
    std::cout << "权重已成功输出到 " << outputFile << "\n";
}



//int main() {
//    Grid grid;
//    grid.setClearance(0.7);
//    grid.setInputScale(10);
//    grid.SetUp("D:\\梁彦诗的项目\\Liang的科研\\混杂的版本\\testcase.kicad_pcb");
//    CostUpdater costUpdater(grid);
//    costUpdater.updateCost();
//    std::string outputFile = "D:\\梁彦诗的项目\\Liang的科研\\混杂的版本\\weights_output.txt";
//    costUpdater.outputWeightsToFile(outputFile);
//    
//    
//    ////输出格子占据情况
//    //for (int i = 0; i <grid.Layers; i++){
//    //for (int j = 36; j < 72; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
//    //    for (int k = 27; k < 66; k++) {
//    //        //if (grid.grid[i][j][k].isoccupied ==1){
//    //        std::cout << grid.grid[1][j][k].isoccupied << " ";
//    //        }
//    //    }
//
//    //    std::cout << std::endl;
//    //}
//    //std::cout << "\n" << std::endl;
//    ////}
//
//    ////输出格子堵塞情况
//    ////for (int i = 0; i <grid.Layers; i++){
//    //for (int j = 9; j < 25; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
//    //    for (int k = 12; k < 28; k++) {
//    //        //if (grid.grid[i][j][k].isoccupied ==1){
//    //        std::cout << grid.grid[1][j][k].isblocked << " ";
//    //        //}
//
//
//    //    }
//
//    //    std::cout << std::endl;
//    //}
//    //std::cout << "\n" << std::endl;
//    ////}
//
//    //
//
//    ////输出格子成本
//    ////for (int i = 0; i <grid.Layers; i++){
//    //for (int j = 9; j < 25; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
//    //    for (int k = 12; k < 28; k++) {
//    //        //if (grid.grid[i][j][k].isoccupied ==1){
//    //        std::cout << grid.grid[1][j][k].cost << " ";
//    //        //}
//    //    }
//
//    //    std::cout << std::endl;
//    //}
//    //std::cout << "\n" << std::endl;
//    //}
//
//
//    
//
//    //for (int i = 0; i <grid.Layers; i++){
//    for (int j = 111; j < 274; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
//        for (int k = 83; k < 246; k++) {
//            std::cout << grid.grid[1][j][k].cost << " ";
//        }
//        std::cout << std::endl;
//    }
//
//    
//
//    return 0;
//}
//
