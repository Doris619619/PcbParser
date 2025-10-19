#include "cost.h"
#include <queue>
#include <fstream>
#include <iostream>
#include "blocked_area_analyzer.h"  
#include "Separate3.h"  

CostUpdater::CostUpdater(Grid& grid) : grid(grid) {}



void CostUpdater::updateCost(int layerId,int start_x, int start_y, int end_x, int end_y) {
    int directions[4][2] = { {-1, 0}, {1, 0}, {0, -1}, {0, 1} };

    int currentCost = 5;
    //for (int layer = 0; layer < grid.Layers; ++layer) {
    for (int j = start_y + 1; j < end_y + 1; ++j) {
        for (int k = start_x; k < end_x; ++k) {
            if (grid.grid[layerId][j][k].isblocked == 1) {
                grid.grid[layerId][j][k].cost = 5;
            }
        }
    }
    //}
    for (int iteration = 0; iteration < 5; ++iteration) {
        //for (int layer = 0; layer < grid.Layers; ++layer) {
        for (int j = start_y + 1; j < end_y + 1; ++j) {
            for (int k = start_x; k < end_x; ++k) {

                if (grid.grid[layerId][j][k].cost == currentCost) {
                    for (const auto& dir : directions) {
                        int newJ = j + dir[0];
                        int newK = k + dir[1];

                        if (newJ >= 0 && newJ < end_y + 2 && newK >= 0 && newK < end_x + 1) {
                            if (grid.grid[layerId][newJ][newK].cost == 0) {
                                grid.grid[layerId][newJ][newK].cost = currentCost - 1;
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




/*
#include <queue>
void CostUpdater::updateCost(int start_x, int start_y, int end_x, int end_y) {
    // 只标记阻塞区域为9000，不进行扩散
    for (int j = start_y; j <= end_y; ++j) {
        for (int k = start_x; k <= end_x; ++k) {
            if (grid.grid[1][j][k].isblocked == 1) {
                grid.grid[1][j][k].cost = 9000; // 阻塞区域被赋值上9000
            }
            else {
                grid.grid[1][j][k].cost = 0; // 非阻塞区域保持为0
            }
        }
    }
}
*/
//
 


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
void CostUpdater::outputWeightsToFile(int layerId,int start_x, int start_y, int end_x, int end_y, const std::string& outputFile) {

    std::ofstream outFile(outputFile);

    if (!outFile.is_open()) {
        std::cerr << "无法打开文件 " << outputFile << " 进行写入\n";
        return;
    }



    for (int i = start_y + 1; i < end_y + 1; i++) {
        for (int k = start_x; k < end_x; k++) {
            outFile << 0 << " ";
            outFile << "\n";
        }
    }


    //for (int i = 0; i <grid.Layers; i++){
    for (int j = start_y + 1; j < end_y + 1; j++) {
        for (int k = start_x; k < end_x; k++) {
            outFile << grid.grid[layerId][j][k].cost << " ";
            outFile << "\n";
        }
    }


    //}

    // �ر��ļ�
    outFile.close();
    std::cout << "权重已成功输出到 " << outputFile << "\n";
}





/*
void CostUpdater::outputWeightsToFile(const std::string& outputFile) {

    std::ofstream outFile(outputFile);

    if (!outFile.is_open()) {
        std::cerr << "无法打开文件 " << outputFile << " 进行写入\n";
        return;
    }


    for (int i = 111; i < 274; i++) {
        for (int k = 83; k < 246; k++) {
            outFile << 0 << " ";
            outFile << "\n";
        }
    }

    //for (int i = 0; i <grid.Layers; i++){
    for (int j = 111; j < 274; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
        for (int k = 83; k < 246; k++) {
            outFile << grid.grid[1][j][k].cost << " ";
            outFile << "\n";
        }
    }


    //}

    // 关闭文件
    outFile.close();
    std::cout << "权重已成功输出到 " << outputFile << "\n";
}
*/









/*
int main() {
    Grid grid;
    grid.setClearance(0.7);
    grid.setInputScale(10);
    grid.SetUp("D:\\梁彦诗的项目\\Liang的科研\\混杂的版本\\testcase.kicad_pcb");
    CostUpdater costUpdater(grid);
    costUpdater.updateCost();
    std::string outputFile = "D:\\梁彦诗的项目\\Liang的科研\\第五步最后9999999999999999999999999999999999999\\文件存储\\weights_output1.txt";
    costUpdater.outputWeightsToFile(outputFile);


    ////输出格子占据情况
    //for (int i = 0; i <grid.Layers; i++){
    //for (int j = 36; j < 72; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
    //    for (int k = 27; k < 66; k++) {
    //        //if (grid.grid[i][j][k].isoccupied ==1){
    //        std::cout << grid.grid[1][j][k].isoccupied << " ";
    //        }
    //    }

    //    std::cout << std::endl;
    //}
    //std::cout << "\n" << std::endl;
    ////}

    ////输出格子堵塞情况
    ////for (int i = 0; i <grid.Layers; i++){
    //for (int j = 9; j < 25; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
    //    for (int k = 12; k < 28; k++) {
    //        //if (grid.grid[i][j][k].isoccupied ==1){
    //        std::cout << grid.grid[1][j][k].isblocked << " ";
    //        //}


    //    }

    //    std::cout << std::endl;
    //}
    //std::cout << "\n" << std::endl;
    ////}

    //

    ////输出格子成本
    ////for (int i = 0; i <grid.Layers; i++){
    //for (int j = 9; j < 25; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
    //    for (int k = 12; k < 28; k++) {
    //        //if (grid.grid[i][j][k].isoccupied ==1){
    //        std::cout << grid.grid[1][j][k].cost << " ";
    //        //}
    //    }

    //    std::cout << std::endl;
    //}
    //std::cout << "\n" << std::endl;
    //}




    //for (int i = 0; i <grid.Layers; i++){
    for (int j = 111; j < 274; j++) { //这里的数字后面读四个坐标的vector,j是y，k是x
        for (int k = 83; k < 246; k++) {
            std::cout << grid.grid[1][j][k].cost << " ";
        }
        std::cout << std::endl;
    }



    return 0;
}

*/




/*
//new main（version 2）
int main() {
    try {
        Grid grid;
        grid.setClearance(0.7);
        grid.setInputScale(10);
        std::string filename = "D:\\梁彦诗的项目\\Liang的科研\\混杂的版本\\testcase.kicad_pcb";
        grid.SetUp(filename);
        CostUpdater costUpdater(grid);
        //costUpdater.updateCost();
        std::string outputFile = "D:\\梁彦诗的项目\\Liang的科研\\第五步最后9999999999999999999999999999999999999\\文件存储\\weights_output1.txt";
        //costUpdater.outputWeightsToFile(outputFile);
        // === 1) 初始化 Grid ===
        //std::string filename = "testcase.kicad_pcb";   // 按需修改
        //Grid grid;
        //grid.SetUp(filename);
        std::cout << "Grid初始化成功! 网格尺寸: " << grid.width << " x " << grid.height << std::endl;

        //std::cout<< std::fixed << std::setprecision(6)<<grid.min_x<<' '<<grid.min_y<<' '<<grid.max_x<<' '<<grid.max_y;

        // === 2) 参数 ===
        int layerId = 1;   // 层序号
        int blockedAreaId = 1;   // 阻塞区域编号
        int expandPixels = 10;   // 像素域扩展像素数（有补偿）
        int kWholePixels = 20;   // 现实域再扩 k 个整像素（与 0.5px 一起）
        double pxSizeMM = 1.0 / static_cast<double>(grid.inputScale);
        // === 3) 读取真实板框（Edge.Cuts 外接矩形） ===
        KiCadParser parser2;
        if (!parser2.parseFile(filename)) {
            throw std::runtime_error("解析 KiCad 文件失败，无法获取板框");
        }
        KiCadParser::Point2D bl, tl, tr, br; // 左下、左上、右上、右下
        if (!parser2.getBoardCorners(bl, tl, tr, br)) {
            throw std::runtime_error("未找到 Edge.Cuts 板框");
        }
        const double boardMinX = bl.x;
        const double boardMinY = bl.y;
        const double boardMaxX = tr.x;
        const double boardMaxY = tr.y;

        //std::cout << std::fixed << std::setprecision(6);
        //std::cout << "[Board] min=(" << boardMinX << "," << boardMinY
        //    << "), max=(" << boardMaxX << "," << boardMaxY << ")\n";

        // === 4) 打印“像素扩展阶段”的像素坐标（有补偿） ===
        // 4.1 阻塞区域像素（1-based）
        auto blockedPixels = GetBlockedAreaPixels_1Based(grid, layerId, blockedAreaId);
        if (blockedPixels.empty()) {
            throw std::runtime_error("未找到指定的阻塞区域");
        }

        // 4.2 原始像素包围盒（1-based）
        BBoxPx1 originalBBox = ComputeBBoxFromPixels1(blockedPixels);
        // 4.3 像素扩展（有补偿）
        ClipReport clipReport{};
        BBoxPx1 expandedBBoxPx = ExpandBBoxByMPixelsCompensated(originalBBox, expandPixels, grid, &clipReport);

        // === 5) （可选）打印“像素域扩展后的现实四角”（用于对比）
        auto baseCorners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
        //std::vector<std::string> names = { "左上(LT)","右上(RT)","右下(RB)","左下(LB)" };
        //std::cout << "\n=== 像素域扩展后的现实四角(用于对比) ===\n";
        //for (size_t i = 0; i < baseCorners.size(); ++i) {
            //std::cout << names[i] << ": (" << baseCorners[i].first << ", " << baseCorners[i].second << ")\n";
        //}

        // === 6) 最终：现实域 0.5px + k×1px 对齐扩张 + 裁剪到板框 ===
        auto expandedCorners = GetExpandedPCBCornersSnapAndClip(
            grid, layerId, blockedAreaId,
            expandPixels, kWholePixels,
            boardMinX, boardMinY, boardMaxX, boardMaxY
        );

        //std::cout << "\n=== 最终小PCB边界框(0.5px + " << kWholePixels
            //<< "px 对齐 + 板框裁剪后) ===\n";
        //for (size_t i = 0; i < expandedCorners.size(); ++i) {
            //std::cout << names[i] << ": (" << expandedCorners[i].first
                //<< ", " << expandedCorners[i].second << ")\n";
        //}


        // === 6.1) 最终四角 -> 像素顶角坐标(1-based, row,col)，保证落在小板内部 ===


        // === 使用通用 get 函数（基于已得到的 expandedCorners） ===
        auto finalCornerPxVec = GetFinalCornerTopPixelsInsidePx1(grid, expandedCorners);
        //auto finalCornerPxVec = GetLastFinalCornerTopPixelsInsidePx1();

        costUpdater.updateCost(finalCornerPxVec[0].second, finalCornerPxVec[0].first, finalCornerPxVec[2].second, finalCornerPxVec[2].first);
        costUpdater.outputWeightsToFile(finalCornerPxVec[0].second, finalCornerPxVec[0].first, finalCornerPxVec[2].second, finalCornerPxVec[2].first, outputFile);
        std::cout << "\n=== 最终小PCB四角对应的像素顶角坐标(1-based, row,col) [via GetFinalCornerTopPixelsInsidePx1] ===\n";
        const char* nm2[4] = { "LT","RT","RB","LB" };
        for (int i = 0; i < 4; ++i) {
            std::cout << nm2[i] << ": (" << finalCornerPxVec[i].first
                << ", " << finalCornerPxVec[i].second << ")\n";
        }
    }
    catch (const std::exception& e) {
        std::cerr << "�������г���: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "����ִ�����!" << std::endl;
    return 0;
}
*/