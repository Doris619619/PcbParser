#pragma once
#ifndef OVERLAP_ANALYZER_H
#define OVERLAP_ANALYZER_H

#include "grid.h"
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <string>

struct OverlapAnalysisResult {
    struct ViaRecommendation {
        int layer0Region;
        int layer1Region;
        int overlapCellCount;
        double overlapRatioToBoard;
        std::vector<std::pair<int, int>> viaPositions; // �������� (row, col)
    };

    struct GndRegionInfo {
        int regionID;
        int layer;
        int size;
        int gndCellCount; // ������GND��Ԫ������
        double gndRatio;  // ������GND��Ԫ�����
        std::vector<std::string> gndElements; // ������GNDԪ������
    };

    int totalBoardArea;
    int minOverlapCells; // �Ƴ� thresholdByArea
    std::vector<ViaRecommendation> recommendations;
    std::map<std::pair<int, int>, int> allOverlaps;
    std::vector<GndRegionInfo> gndRegions; // ����GNDԪ�ص�����
};

class OverlapAnalyzer {
public:
    // ��������֮����ص�����
    static OverlapAnalysisResult analyzeLayerOverlap(
        const Grid& grid,
        int layer0 = 0,
        int layer1 = 1,
        int minOverlapCells = 10  // �޸ģ��Ƴ�overlapThresholdRatio����
    );

    // ��������Ƿ����GNDԪ�أ���ǿ�棺���������ں��ܱߣ�
    static bool hasGndElements(const Grid& grid, int layer, int regionID, double gndThreshold = 0.05);

    // ��ȡ�����е�GNDԪ������
    static std::vector<std::string> getGndElementTypes(const Grid& grid, int layer, int regionID);

    // ͳ�������ڵ�GND��Ϣ
    static void getGndRegionInfo(const Grid& grid, int layer, int regionID, int& gndCellCount, double& gndRatio);

    // ��鵥Ԫ���Ƿ�ΪGND
    static bool isCellGnd(const Grid& grid, int layer, int row, int col);

    // ��������ܱ��Ƿ����GND�����ڼ�ⱻGND��Χ������
    static bool hasSurroundingGnd(const Grid& grid, int layer, int regionID, int searchRadius = 3, double gndThreshold = 0.3);

    // ���ɷ�������
    static std::string generateReport(const OverlapAnalysisResult& result, const Grid& grid);

    // ��ȡ�ض�����Ե��ص�λ��
    static std::vector<std::pair<double, double>> getOverlapPositions(
        const Grid& grid,
        int layer0, int region0,
        int layer1, int region1,
        int maxPositions = 100
    );
};

#endif // OVERLAP_ANALYZER_H