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
//// ��������
//void analyzeRegion5Elements(const Grid& grid);
//
//int main() {
//    Grid grid;
//    grid.SetUp("testcase.kicad_pcb");
//
//    // ��������ļ�
//    std::ofstream outFile("pcb_analyzer1.txt");
//    if (!outFile.is_open()) {
//        std::cerr << "�޷���������ļ� pcb_analyzer.txt" << std::endl;
//        return 1;
//    }
//
//    // ֻ���ǰ����
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
//    // === �ص���������ʹ���ж� ===
//    outFile << "=== �ص���������ʹ�׽��� ===" << std::endl;
//
//    if (grid.Layers >= 2) {
//        // ʹ���µ��ص�������
//        auto overlapResult = OverlapAnalyzer::analyzeLayerOverlap(grid);
//        std::string report = OverlapAnalyzer::generateReport(overlapResult, grid);
//        outFile << report;
//    }
//    else {
//        outFile << "������Ҫ����2���������ص�����" << std::endl;
//    }
//
//    outFile.close();
//    std::cout << "PCB��������ѱ��浽 pcb_analyzer.txt" << std::endl;
//    std::cout << "�����ص���������ʹ�׽���" << std::endl;
//
//    // ��������������5������Ԫ����Ԫ��
//    std::cout << "\n=== ��ʼ��������5��Ԫ����Ԫ�� ===" << std::endl;
//    analyzeRegion5Elements(grid);
//
//    // ===== ���������������������������Ϣ��=====
//    std::cout << "\n=== �������������Ϣ ===" << std::endl;
//
//    // ͳ��ÿ���������Ԫ�������ͷֲ�
//    for (int layer = 0; layer < grid.Layers && layer < 2; layer++) {
//        int blockedCount = 0;
//        int totalCells = grid.height * grid.width;
//        int boundaryBlocked = 0;
//        int interiorBlocked = 0;
//
//        // ͳ�Ʊ߽��������ڲ�����
//        for (int j = 0; j < grid.height; j++) {
//            for (int k = 0; k < grid.width; k++) {
//                if (grid.grid[layer][j][k].isblocked) {
//                    blockedCount++;
//                    // ����Ƿ��Ǳ߽�
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
//        std::cout << "Layer " << layer << " ����ͳ��:" << std::endl;
//        std::cout << "  �����С: " << grid.height << "x" << grid.width << " = " << totalCells << " ��Ԫ��" << std::endl;
//        std::cout << "  ��������Ԫ��: " << blockedCount << "/" << totalCells
//            << " (" << std::fixed << std::setprecision(1)
//            << (blockedCount * 100.0 / totalCells) << "%)" << std::endl;
//        std::cout << "  �߽�����: " << boundaryBlocked << " ��Ԫ��" << std::endl;
//        std::cout << "  �ڲ�����: " << interiorBlocked << " ��Ԫ��" << std::endl;
//
//        // ������ʾһЩ������Ԫ�����Ϣ
//        std::cout << "  ������Ԫ�����:" << std::endl;
//        int sampleCount = 0;
//        for (int j = 0; j < grid.height && sampleCount < 5; j += grid.height / 10) {
//            for (int k = 0; k < grid.width && sampleCount < 5; k += grid.width / 10) {
//                if (grid.grid[layer][j][k].isblocked) {
//                    std::cout << "    λ��(" << j << "," << k << "): regionID="
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
//    // ���������������
//    std::cout << "\n=== ��ʼ����������� ===" << std::endl;
//    auto blockedAnalysis = BlockedAreaAnalyzer::analyzeAllLayers(grid);
//
//    // ������������������棨ʹ���µ��������溯����
//    std::string blockedReport = BlockedAreaAnalyzer::generateFullReport(blockedAnalysis, grid);
//
//    // ������������������������ļ�
//    std::ofstream blockedOutFile("blocked_area_analysis.txt");
//    if (blockedOutFile.is_open()) {
//        blockedOutFile << blockedReport;
//        blockedOutFile.close();
//        std::cout << "���������������ѱ��浽 blocked_area_analysis.txt" << std::endl;
//    }
//    else {
//        std::cerr << "�޷��������������������ļ�" << std::endl;
//    }
//
//    // ͬʱ�ڿ���̨��������������ժҪ
//    std::cout << "\n" << blockedReport << std::endl;
//
//    // ===== ����������������Ҳ׷�ӵ�������ļ��� =====
//    std::ofstream appendFile("pcb_analyzer.txt", std::ios::app);
//    if (appendFile.is_open()) {
//        appendFile << "\n\n=== ����������� ===" << std::endl;
//
//        // ��д�������Ϣ
//        appendFile << "=== �������������Ϣ ===" << std::endl;
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
//                << " ��������Ԫ�� (" << std::fixed << std::setprecision(1)
//                << (blockedCount * 100.0 / totalCells) << "%)" << std::endl;
//        }
//
//        // ��д���������
//        appendFile << blockedReport;
//        appendFile.close();
//        std::cout << "����������������׷�ӵ� pcb_analyzer.txt" << std::endl;
//    }
//
//    return 0;
//}
//
//// ��������5����Ԫ����Ԫ�صĺ��������ֲ��䣩
//void analyzeRegion5Elements(const Grid& grid) {
//    // ... ����������ֲ��� ...
//}