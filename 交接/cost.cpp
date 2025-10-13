#include "cost.h"
#include <queue>
#include <fstream>
#include <iostream>
#include "blocked_area_analyzer.h"  

CostUpdater::CostUpdater(Grid& grid) : grid(grid) {}

void CostUpdater::updateCost() {
    int directions[4][2] = { {-1, 0}, {1, 0}, {0, -1}, {0, 1} };//����

    int currentCost = 5;
    //for (int layer = 0; layer < grid.Layers; ++layer) {
    for (int j = 111; j < 274; ++j) {//�߽�����
        for (int k = 83; k < 246; ++k) {
            if (grid.grid[1][j][k].isblocked == 1) {
                grid.grid[1][j][k].cost = 5;
            }
        }
    }
    //}
    for (int iteration = 0; iteration < 5; ++iteration) {
        //for (int layer = 0; layer < grid.Layers; ++layer) {
        for (int j = 111; j < 274; ++j) {//�߽�����
            for (int k = 83; k < 246; ++k) {

                if (grid.grid[1][j][k].cost == currentCost) {
                    for (const auto& dir : directions) {
                        int newJ = j + dir[0];
                        int newK = k + dir[1];

                        if (newJ >= 0 && newJ < 275 && newK >= 0 && newK < 247) {// �ǵøı߽�����
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
//    const int unitStep = 30; // ÿ10�񽵵�һ������
//    const int maxDistance = maxCost * unitStep;
//
//    int maxJ = grid.grid[1].size();
//    int maxK = grid.grid[1][0].size();
//
//    std::queue<std::pair<int, int>> q;
//
//    // Step 1����ʼ���ϰ���Ϊ5������BFS����
//    for (int j = 111; j < 274; ++j) {
//        for (int k = 83; k < 246; ++k) {
//            if (grid.grid[1][j][k].isblocked == 1) {
//                grid.grid[1][j][k].cost = maxCost;
//                q.emplace(j, k); // ���ϰ���������
//            }
//            else {
//                grid.grid[1][j][k].cost = 0; // ��շ��ϰ���
//            }
//        }
//    }
//
//    // Step 2��BFS ��ɢ����¼ÿ���㵽�ϰ������С�����پ���
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
//                // �����û���ʹ��Ҳ����ϰ���
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
//    int maxJ = grid.grid[1].size();            // ����
//    int maxK = grid.grid[1][0].size();         // ����
//
//    // Step 1: ��ʼ������ cost Ϊ 0
//    for (int j = 0; j < maxJ; ++j) {
//        for (int k = 0; k < maxK; ++k) {
//            grid.grid[1][j][k].cost = 0;
//        }
//    }
//
//    // Step 2: �������и��ӣ��������㣬����������60��
//    for (int j = 0; j < maxJ; ++j) {
//        for (int k = 0; k < maxK; ++k) {
//            if (grid.grid[1][j][k].isblocked == 1) {
//                // ��������Ϊ cost = 5
//                grid.grid[1][j][k].cost = cost_blocked;
//
//                // ������չ��� 60 ��
//                for (int offset = 1; offset <= horizontal_distance; ++offset) {
//                    int left = k - offset;
//                    if (left < 0) break;
//
//                    if (grid.grid[1][j][left].isblocked != 1 && grid.grid[1][j][left].cost < cost_horizontal) {
//                        grid.grid[1][j][left].cost = cost_horizontal;
//                    }
//                }
//
//                // ������չ��� 60 ��
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
        std::cerr << "�޷����ļ� " << outputFile << " ����д��\n";
        return;
    }

    //for (int i = 0; i <grid.Layers; i++){
    for (int j = 111; j < 274; j++) { //��������ֺ�����ĸ������vector,j��y��k��x
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

    // �ر��ļ�
    outFile.close();
    std::cout << "Ȩ���ѳɹ������ " << outputFile << "\n";
}



//int main() {
//    Grid grid;
//    grid.setClearance(0.7);
//    grid.setInputScale(10);
//    grid.SetUp("D:\\����ʫ����Ŀ\\Liang�Ŀ���\\���ӵİ汾\\testcase.kicad_pcb");
//    CostUpdater costUpdater(grid);
//    costUpdater.updateCost();
//    std::string outputFile = "D:\\����ʫ����Ŀ\\Liang�Ŀ���\\���ӵİ汾\\weights_output.txt";
//    costUpdater.outputWeightsToFile(outputFile);
//    
//    
//    ////�������ռ�����
//    //for (int i = 0; i <grid.Layers; i++){
//    //for (int j = 36; j < 72; j++) { //��������ֺ�����ĸ������vector,j��y��k��x
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
//    ////������Ӷ������
//    ////for (int i = 0; i <grid.Layers; i++){
//    //for (int j = 9; j < 25; j++) { //��������ֺ�����ĸ������vector,j��y��k��x
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
//    ////������ӳɱ�
//    ////for (int i = 0; i <grid.Layers; i++){
//    //for (int j = 9; j < 25; j++) { //��������ֺ�����ĸ������vector,j��y��k��x
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
//    for (int j = 111; j < 274; j++) { //��������ֺ�����ĸ������vector,j��y��k��x
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
