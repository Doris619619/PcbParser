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
        std::vector<std::pair<int, int>> viaPositions; // 网格坐标 (row, col)
    };

    struct GndRegionInfo {
        int regionID;
        int layer;
        int size;
        int gndCellCount; // 新增：GND单元格数量
        double gndRatio;  // 新增：GND单元格比例
        std::vector<std::string> gndElements; // 包含的GND元素类型
    };

    int totalBoardArea;
    int minOverlapCells; // 移除 thresholdByArea
    std::vector<ViaRecommendation> recommendations;
    std::map<std::pair<int, int>, int> allOverlaps;
    std::vector<GndRegionInfo> gndRegions; // 包含GND元素的区域
};

class OverlapAnalyzer {
public:
    // 分析两层之间的重叠区域
    static OverlapAnalysisResult analyzeLayerOverlap(
        const Grid& grid,
        int layer0 = 0,
        int layer1 = 1,
        int minOverlapCells = 10  // 修改：移除overlapThresholdRatio参数
    );

    // 检查区域是否包含GND元素（增强版：包括区域内和周边）
    static bool hasGndElements(const Grid& grid, int layer, int regionID, double gndThreshold = 0.05);

    // 获取区域中的GND元素类型
    static std::vector<std::string> getGndElementTypes(const Grid& grid, int layer, int regionID);

    // 统计区域内的GND信息
    static void getGndRegionInfo(const Grid& grid, int layer, int regionID, int& gndCellCount, double& gndRatio);

    // 检查单元格是否为GND
    static bool isCellGnd(const Grid& grid, int layer, int row, int col);

    // 检查区域周边是否存在GND（用于检测被GND包围的区域）
    static bool hasSurroundingGnd(const Grid& grid, int layer, int regionID, int searchRadius = 3, double gndThreshold = 0.3);

    // 生成分析报告
    static std::string generateReport(const OverlapAnalysisResult& result, const Grid& grid);

    // 获取特定区域对的重叠位置
    static std::vector<std::pair<double, double>> getOverlapPositions(
        const Grid& grid,
        int layer0, int region0,
        int layer1, int region1,
        int maxPositions = 100
    );
};

#endif // OVERLAP_ANALYZER_H