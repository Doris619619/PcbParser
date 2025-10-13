
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
    std::vector<std::pair<int, int>> cells; // ����������Ԫ���λ��
    std::set<int> adjacentRegions; // ���ڵ�region ID
    std::map<int, int> regionProximity; // region ID -> ��С����
};

struct LayerBlockedAnalysis {
    int layerID;
    std::vector<BlockedAreaInfo> blockedAreas;
    std::map<std::pair<int, int>, std::string> separationInfo; // (region1, region2) -> ������Ϣ
    std::set<std::pair<int, int>> purpleDashedPairs; // ����������ɫ���߷ָ��������
};

class BlockedAreaAnalyzer {
public:
    // �����������������
    static LayerBlockedAnalysis analyzeLayerBlockedAreas(const Grid& grid, int layer);

    // �������в����������
    static std::map<int, LayerBlockedAnalysis> analyzeAllLayers(const Grid& grid);

    // ���ɷ�������
    static std::string generateReport(const std::map<int, LayerBlockedAnalysis>& analysis);

    // �����������������ӹ�ϵ
    static std::vector<std::string> generateConnectionRelations(
        const std::map<int, LayerBlockedAnalysis>& blockedAnalysis,
        const Grid& grid);

    // ���������������������棨�������ӹ�ϵ��
    static std::string generateFullReport(
        const std::map<int, LayerBlockedAnalysis>& blockedAnalysis,
        const Grid& grid);

    // �������region�Ƿ�ĳ����������ֿ�
    static bool areRegionsSeparated(const LayerBlockedAnalysis& layerAnalysis, int region1, int region2);

private:
    // �����㷨Ѱ����������
    static BlockedAreaInfo findBlockedArea(const Grid& grid, int layer, int startRow, int startCol,
        std::vector<std::vector<bool>>& visited, int areaID);

    // Ѱ�������������ڵ�regions
    static void findAdjacentRegions(const Grid& grid, int layer, BlockedAreaInfo& blockedArea);

    // ������������������region����С����
    static void calculateRegionProximity(const Grid& grid, int layer, BlockedAreaInfo& blockedArea);

    // ��鵥Ԫ���Ƿ�������Χ��
    static bool isValidCell(const Grid& grid, int row, int col);
};

#endif // BLOCKED_AREA_ANALYZER_H#pragma once
