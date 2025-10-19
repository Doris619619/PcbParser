//#include <vector>
//#include <fstream>
//#include <iostream>
//#include <map>
//#include <set>
//#include "grid.h"
//#include "overlap_analyzer.h"
//#include "blocked_area_analyzer.h"
//#include <stack>
//#include <limits>
//#include <algorithm>
//#include <iomanip>
//
//// 函数声明
//void analyzeRegion5Elements(const Grid& grid);
//
//int main() {
//    Grid grid;
//    grid.SetUp("testcase.kicad_pcb");
//
//    // 创建输出文件
//    std::ofstream outFile("pcb_analyzer1.txt");
//    if (!outFile.is_open()) {
//        std::cerr << "无法创建输出文件 pcb_analyzer.txt" << std::endl;
//        return 1;
//    }
//
//    // 只输出前两层
//    for (int layer = 0; layer < 2 && layer < grid.Layers; layer++) {
//        outFile << "=== Layer " << layer << " ===" << std::endl;
//
//        for (int j = 0; j < grid.height; j++) {
//            for (int k = 0; k < grid.width; k++) {
//                outFile << grid.grid[layer][j][k].isoccupied << " ";
//            }
//            outFile << std::endl;
//        }
//        outFile << "\n" << std::endl;
//    }
//
//    // === 重叠区域分析和打孔判断 ===
//    outFile << "=== 重叠区域分析和打孔建议 ===" << std::endl;
//
//    if (grid.Layers >= 2) {
//        // 使用新的重叠分析器
//        auto overlapResult = OverlapAnalyzer::analyzeLayerOverlap(grid);
//        std::string report = OverlapAnalyzer::generateReport(overlapResult, grid);
//        outFile << report;
//    }
//    else {
//        outFile << "错误：需要至少2层来进行重叠分析" << std::endl;
//    }
//
//    outFile.close();
//    std::cout << "PCB分析结果已保存到 pcb_analyzer.txt" << std::endl;
//    std::cout << "包含重叠区域分析和打孔建议" << std::endl;
//
//    // 新增：检索区域5的所有元器件元素
//    std::cout << "\n=== 开始检索区域5的元器件元素 ===" << std::endl;
//    analyzeRegion5Elements(grid);
//
//    // ===== 新增：阻塞区域分析（带调试信息）=====
//    std::cout << "\n=== 阻塞区域调试信息 ===" << std::endl;
//
//    // 统计每层的阻塞单元格数量和分布
//    for (int layer = 0; layer < grid.Layers && layer < 2; layer++) {
//        int blockedCount = 0;
//        int totalCells = grid.height * grid.width;
//        int boundaryBlocked = 0;
//        int interiorBlocked = 0;
//
//        // 统计边界阻塞和内部阻塞
//        for (int j = 0; j < grid.height; j++) {
//            for (int k = 0; k < grid.width; k++) {
//                if (grid.grid[layer][j][k].isblocked) {
//                    blockedCount++;
//                    // 检查是否是边界
//                    if (j == 0 || j == grid.height - 1 || k == 0 || k == grid.width - 1) {
//                        boundaryBlocked++;
//                    }
//                    else {
//                        interiorBlocked++;
//                    }
//                }
//            }
//        }
//
//        std::cout << "Layer " << layer << " 阻塞统计:" << std::endl;
//        std::cout << "  网格大小: " << grid.height << "x" << grid.width << " = " << totalCells << " 单元格" << std::endl;
//        std::cout << "  总阻塞单元格: " << blockedCount << "/" << totalCells
//            << " (" << std::fixed << std::setprecision(1)
//            << (blockedCount * 100.0 / totalCells) << "%)" << std::endl;
//        std::cout << "  边界阻塞: " << boundaryBlocked << " 单元格" << std::endl;
//        std::cout << "  内部阻塞: " << interiorBlocked << " 单元格" << std::endl;
//
//        // 采样显示一些阻塞单元格的信息
//        std::cout << "  阻塞单元格采样:" << std::endl;
//        int sampleCount = 0;
//        for (int j = 0; j < grid.height && sampleCount < 5; j += grid.height / 10) {
//            for (int k = 0; k < grid.width && sampleCount < 5; k += grid.width / 10) {
//                if (grid.grid[layer][j][k].isblocked) {
//                    std::cout << "    位置(" << j << "," << k << "): regionID="
//                        << grid.grid[layer][j][k].regionID
//                        << ", isoccupied=" << grid.grid[layer][j][k].isoccupied
//                        << ", nodes=" << grid.grid[layer][j][k].nodes.size() << std::endl;
//                    sampleCount++;
//                }
//            }
//        }
//        std::cout << std::endl;
//    }
//
//    // 进行阻塞区域分析
//    std::cout << "\n=== 开始阻塞区域分析 ===" << std::endl;
//    auto blockedAnalysis = BlockedAreaAnalyzer::analyzeAllLayers(grid);
//
//    // 生成阻塞区域分析报告（使用新的完整报告函数）
//    std::string blockedReport = BlockedAreaAnalyzer::generateFullReport(blockedAnalysis, grid);
//
//    // 输出阻塞区域分析结果到单独文件
//    std::ofstream blockedOutFile("blocked_area_analysis.txt");
//    if (blockedOutFile.is_open()) {
//        blockedOutFile << blockedReport;
//        blockedOutFile.close();
//        std::cout << "阻塞区域分析结果已保存到 blocked_area_analysis.txt" << std::endl;
//    }
//    else {
//        std::cerr << "无法创建阻塞区域分析输出文件" << std::endl;
//    }
//
//    // 同时在控制台输出阻塞区域分析摘要
//    std::cout << "\n" << blockedReport << std::endl;
//
//    // ===== 将阻塞区域分析结果也追加到主输出文件中 =====
//    std::ofstream appendFile("pcb_analyzer.txt", std::ios::app);
//    if (appendFile.is_open()) {
//        appendFile << "\n\n=== 阻塞区域分析 ===" << std::endl;
//
//        // 先写入调试信息
//        appendFile << "=== 阻塞区域调试信息 ===" << std::endl;
//        for (int layer = 0; layer < grid.Layers && layer < 2; layer++) {
//            int blockedCount = 0;
//            int totalCells = grid.height * grid.width;
//
//            for (int j = 0; j < grid.height; j++) {
//                for (int k = 0; k < grid.width; k++) {
//                    if (grid.grid[layer][j][k].isblocked) {
//                        blockedCount++;
//                    }
//                }
//            }
//
//            appendFile << "Layer " << layer << ": " << blockedCount << "/" << totalCells
//                << " 个阻塞单元格 (" << std::fixed << std::setprecision(1)
//                << (blockedCount * 100.0 / totalCells) << "%)" << std::endl;
//        }
//
//        // 再写入分析报告
//        appendFile << blockedReport;
//        appendFile.close();
//        std::cout << "阻塞区域分析结果已追加到 pcb_analyzer.txt" << std::endl;
//    }
//
//    return 0;
//}
//
//// 检索区域5所有元器件元素的函数（保持不变）
//void analyzeRegion5Elements(const Grid& grid) {
//    // ... 这个函数保持不变 ...
//}