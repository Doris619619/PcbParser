
#ifndef BLOCKED_AREA_ANALYZER_H
#define BLOCKED_AREA_ANALYZER_H

#include "grid.h"
#include "overlap_analyzer.h"
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <string>

struct BlockedAreaInfo {
    int areaID;
    int size;
    std::vector<std::pair<int, int>> cells; // 所有阻塞单元格的位置
    std::set<int> adjacentRegions; // 相邻的region ID
    std::map<int, int> regionProximity; // region ID -> 最小距离
};

struct LayerBlockedAnalysis {
    int layerID;
    std::vector<BlockedAreaInfo> blockedAreas;
    std::map<std::pair<int, int>, std::string> separationInfo; // (region1, region2) -> 描述信息
    std::set<std::pair<int, int>> purpleDashedPairs; // 新增：被紫色虚线分隔的区域对
};

class BlockedAreaAnalyzer {
public:
    // 分析单层的阻塞区域
    static LayerBlockedAnalysis analyzeLayerBlockedAreas(const Grid& grid, int layer);

    // 分析所有层的阻塞区域
    static std::map<int, LayerBlockedAnalysis> analyzeAllLayers(const Grid& grid);

    // 生成分析报告
    static std::string generateReport(const std::map<int, LayerBlockedAnalysis>& analysis);

    // 新增函数：生成连接关系
    static std::vector<std::string> generateConnectionRelations(
        const std::map<int, LayerBlockedAnalysis>& blockedAnalysis,
        const Grid& grid);

    // 新增函数：生成完整报告（包含连接关系）
    static std::string generateFullReport(
        const std::map<int, LayerBlockedAnalysis>& blockedAnalysis,
        const Grid& grid);

    // 检查两个region是否被某个阻塞区域分开
    static bool areRegionsSeparated(const LayerBlockedAnalysis& layerAnalysis, int region1, int region2);

private:
    // 泛洪算法寻找阻塞区域
    static BlockedAreaInfo findBlockedArea(const Grid& grid, int layer, int startRow, int startCol,
        std::vector<std::vector<bool>>& visited, int areaID);

    // 寻找阻塞区域相邻的regions
    static void findAdjacentRegions(const Grid& grid, int layer, BlockedAreaInfo& blockedArea);

    // 计算阻塞区域与相邻region的最小距离
    static void calculateRegionProximity(const Grid& grid, int layer, BlockedAreaInfo& blockedArea);

    // 检查单元格是否在网格范围内
    static bool isValidCell(const Grid& grid, int row, int col);
};

#endif // BLOCKED_AREA_ANALYZER_H#pragma once
