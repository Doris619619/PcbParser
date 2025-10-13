#include "blocked_area_analyzer.h"
#include <queue>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>

bool BlockedAreaAnalyzer::isValidCell(const Grid& grid, int row, int col) {
    return row >= 0 && row < grid.height && col >= 0 && col < grid.width;
}

BlockedAreaInfo BlockedAreaAnalyzer::findBlockedArea(const Grid& grid, int layer,
    int startRow, int startCol,
    std::vector<std::vector<bool>>& visited,
    int areaID) {
    BlockedAreaInfo area;
    area.areaID = areaID;
    area.size = 0;  // ��ȷ��ʼ����С

    // �����ʼ��Ԫ���Ƿ���Ч������������
    if (!isValidCell(grid, startRow, startCol) || visited[startRow][startCol] ||
        !grid.grid[layer][startRow][startCol].isblocked) {
        return area;  // ���ؿյ�����������Ϣ
    }

    std::queue<std::pair<int, int>> q;
    q.push({ startRow, startCol });
    visited[startRow][startCol] = true;

    // �ĸ������ϡ��¡�����
    const int dr[4] = { -1, 1, 0, 0 };
    const int dc[4] = { 0, 0, -1, 1 };

    while (!q.empty()) {
        auto [row, col] = q.front();
        q.pop();

        // ��ӵ���������
        area.cells.push_back({ row, col });
        area.size++;

        // ����ĸ�����
        for (int i = 0; i < 4; i++) {
            int newRow = row + dr[i];
            int newCol = col + dc[i];

            // �ؼ��������������µ�Ԫ������������ż���
            if (isValidCell(grid, newRow, newCol) &&
                !visited[newRow][newCol] &&
                grid.grid[layer][newRow][newCol].isblocked) {
                visited[newRow][newCol] = true;
                q.push({ newRow, newCol });
            }
        }
    }

    return area;
}

void BlockedAreaAnalyzer::findAdjacentRegions(const Grid& grid, int layer, BlockedAreaInfo& blockedArea) {
    std::set<int> adjacent;

    for (const auto& cell : blockedArea.cells) {
        int row = cell.first;
        int col = cell.second;

        // ����ĸ����򣨲������Խ��ߣ�
        const int dr[4] = { -1, 1, 0, 0 };
        const int dc[4] = { 0, 0, -1, 1 };

        for (int i = 0; i < 4; i++) {
            int newRow = row + dr[i];
            int newCol = col + dc[i];

            if (isValidCell(grid, newRow, newCol)) {
                int regionID = grid.grid[layer][newRow][newCol].regionID;
                if (regionID > 0) {
                    adjacent.insert(regionID);
                }
            }
        }
    }

    blockedArea.adjacentRegions = adjacent;
}

void BlockedAreaAnalyzer::calculateRegionProximity(const Grid& grid, int layer, BlockedAreaInfo& blockedArea) {
    std::map<int, int> minDistances;

    for (int regionID : blockedArea.adjacentRegions) {
        minDistances[regionID] = std::numeric_limits<int>::max();
    }

    // �����������ÿ����Ԫ�񣬼��㵽����region����С����
    for (const auto& blockedCell : blockedArea.cells) {
        int blockedRow = blockedCell.first;
        int blockedCol = blockedCell.second;

        // ������Χ��region��Ԫ��
        for (int searchRadius = 1; searchRadius <= 10; searchRadius++) {
            bool foundAll = true;

            for (int regionID : blockedArea.adjacentRegions) {
                if (minDistances[regionID] <= searchRadius - 1) continue;

                bool foundRegion = false;

                // �����뾶ΪsearchRadius�ķ�������
                for (int dr = -searchRadius; dr <= searchRadius && !foundRegion; dr++) {
                    for (int dc = -searchRadius; dc <= searchRadius && !foundRegion; dc++) {
                        int newRow = blockedRow + dr;
                        int newCol = blockedCol + dc;

                        if (isValidCell(grid, newRow, newCol) &&
                            grid.grid[layer][newRow][newCol].regionID == regionID) {
                            // ʹ�������پ���
                            int distance = std::abs(dr) + std::abs(dc);
                            if (distance < minDistances[regionID]) {
                                minDistances[regionID] = distance;
                            }
                            foundRegion = true;
                        }
                    }
                }

                if (!foundRegion) {
                    foundAll = false;
                }
            }

            if (foundAll) break;
        }
    }

    blockedArea.regionProximity = minDistances;
}

LayerBlockedAnalysis BlockedAreaAnalyzer::analyzeLayerBlockedAreas(const Grid& grid, int layer) {
    LayerBlockedAnalysis analysis;
    analysis.layerID = layer;

    std::vector<std::vector<bool>> visited(grid.height, std::vector<bool>(grid.width, false));
    int areaID = 1;

    // ���ԣ�ͳ��ʵ�ʵ�������Ԫ��
    int actualBlockedCount = 0;
    for (int row = 0; row < grid.height; row++) {
        for (int col = 0; col < grid.width; col++) {
            if (grid.grid[layer][row][col].isblocked) {
                actualBlockedCount++;
            }
        }
    }
    std::cout << "Layer " << layer << " ʵ��������Ԫ��: " << actualBlockedCount << std::endl;

    // Ѱ��������������
    for (int row = 0; row < grid.height; row++) {
        for (int col = 0; col < grid.width; col++) {
            if (!visited[row][col] && grid.grid[layer][row][col].isblocked) {
                BlockedAreaInfo area = findBlockedArea(grid, layer, row, col, visited, areaID);

                // �ؼ�������ֻ�е��ҵ���Ч����������ʱ����ӵ������
                if (area.size > 0 && area.size <= actualBlockedCount) {
                    std::cout << "�ҵ��������� " << areaID << ": ��С=" << area.size << " ��Ԫ��" << std::endl;

                    // ��������regions
                    findAdjacentRegions(grid, layer, area);
                    calculateRegionProximity(grid, layer, area);
                    analysis.blockedAreas.push_back(area);
                    areaID++;
                }
                else if (area.size > actualBlockedCount) {
                    std::cout << "������Ч���������� " << areaID << ": ��С=" << area.size
                        << " (����ʵ��������Ԫ���� " << actualBlockedCount << ")" << std::endl;
                }
            }
        }
    }

    // ������Щregions����������ֿ�
    for (const auto& blockedArea : analysis.blockedAreas) {
        if (blockedArea.adjacentRegions.size() >= 2) {
            std::vector<int> regions(blockedArea.adjacentRegions.begin(),
                blockedArea.adjacentRegions.end());

            for (size_t i = 0; i < regions.size(); i++) {
                for (size_t j = i + 1; j < regions.size(); j++) {
                    std::pair<int, int> regionPair = std::make_pair(regions[i], regions[j]);
                    std::stringstream ss;
                    ss << "����" << regions[i] << "������" << regions[j]
                        << "����������" << blockedArea.areaID << "�ֿ�";
                    analysis.separationInfo[regionPair] = ss.str();

                    // ��¼����ɫ���߷ָ��������
                    analysis.purpleDashedPairs.insert(regionPair);
                }
            }
        }
    }

    return analysis;
}

std::map<int, LayerBlockedAnalysis> BlockedAreaAnalyzer::analyzeAllLayers(const Grid& grid) {
    std::map<int, LayerBlockedAnalysis> allAnalysis;

    // ֻ����ǰ���㣬�������������Ϣ
    int maxLayersToAnalyze = std::min(5, grid.Layers);
    for (int layer = 0; layer < maxLayersToAnalyze; layer++) {
        std::cout << "���� Layer " << layer << "..." << std::endl;
        allAnalysis[layer] = analyzeLayerBlockedAreas(grid, layer);
    }

    return allAnalysis;
}

bool BlockedAreaAnalyzer::areRegionsSeparated(const LayerBlockedAnalysis& layerAnalysis,
    int region1, int region2) {
    std::pair<int, int> pair1 = std::make_pair(region1, region2);
    std::pair<int, int> pair2 = std::make_pair(region2, region1);

    return layerAnalysis.separationInfo.find(pair1) != layerAnalysis.separationInfo.end() ||
        layerAnalysis.separationInfo.find(pair2) != layerAnalysis.separationInfo.end();
}

std::string BlockedAreaAnalyzer::generateReport(const std::map<int, LayerBlockedAnalysis>& analysis) {
    std::ostringstream report;

    report << "=== ��������������� ===" << std::endl;
    report << "PCB�ܲ���: " << analysis.size() << std::endl << std::endl;

    for (const auto& layerEntry : analysis) {
        int layerID = layerEntry.first;
        const LayerBlockedAnalysis& layerAnalysis = layerEntry.second;

        report << "Layer " << layerID << ":" << std::endl;
        report << "  �����������: " << layerAnalysis.blockedAreas.size() << std::endl;

        if (!layerAnalysis.blockedAreas.empty()) {
            report << "  ������������:" << std::endl;
            for (const auto& blockedArea : layerAnalysis.blockedAreas) {
                report << "    �������� " << blockedArea.areaID
                    << ": ��С=" << blockedArea.size
                    << " ��Ԫ��, ����regions=";

                if (!blockedArea.adjacentRegions.empty()) {
                    for (auto it = blockedArea.adjacentRegions.begin();
                        it != blockedArea.adjacentRegions.end(); ++it) {
                        if (it != blockedArea.adjacentRegions.begin()) report << ", ";
                        report << *it;
                    }

                    // ��ʾ��ÿ��region����С����
                    report << std::endl << "      ������region����С����: ";
                    for (const auto& proximity : blockedArea.regionProximity) {
                        report << "����" << proximity.first << "->" << proximity.second << "��Ԫ�� ";
                    }
                }
                else {
                    report << "��";
                }
                report << std::endl;
            }
        }
        else {
            report << "  ����û����������" << std::endl;
        }

        // ��ʾ���ֿ���regions
        if (!layerAnalysis.separationInfo.empty()) {
            report << "  ����������ֿ���regions:" << std::endl;
            for (const auto& separation : layerAnalysis.separationInfo) {
                report << "    " << separation.second << std::endl;
            }
        }
        else {
            report << "  û�б���������ֿ���regions" << std::endl;
        }

        report << std::endl;
    }

    // ͳ����Ϣ
    int totalBlockedAreas = 0;
    int totalSeparations = 0;
    int totalBlockedCells = 0;

    for (const auto& layerEntry : analysis) {
        totalBlockedAreas += layerEntry.second.blockedAreas.size();
        totalSeparations += layerEntry.second.separationInfo.size();
        for (const auto& blockedArea : layerEntry.second.blockedAreas) {
            totalBlockedCells += blockedArea.size;
        }
    }

    report << "=== ͳ��ժҪ ===" << std::endl;
    report << "������������: " << totalBlockedAreas << std::endl;
    report << "��������Ԫ����: " << totalBlockedCells << std::endl;
    report << "��region�������: " << totalSeparations << std::endl;

    return report.str();
}

// ������������ȡ��������Ĵ�С
static int getBlockedAreaSize(const LayerBlockedAnalysis& layerAnalysis, int region1, int region2) {
    for (const auto& blockedArea : layerAnalysis.blockedAreas) {
        if (blockedArea.adjacentRegions.find(region1) != blockedArea.adjacentRegions.end() &&
            blockedArea.adjacentRegions.find(region2) != blockedArea.adjacentRegions.end()) {
            return blockedArea.size;
        }
    }
    return 0;
}

// �޸ĺ�����ӹ�ϵ���ɺ���
std::vector<std::string> BlockedAreaAnalyzer::generateConnectionRelations(
    const std::map<int, LayerBlockedAnalysis>& blockedAnalysis,
    const Grid& grid) {

    std::vector<std::string> connections;
    std::set<std::pair<int, int>> processedPairs; // �����ظ�

    // 1. ����ͬһ���ڵ����ӹ�ϵ
    for (const auto& layerEntry : blockedAnalysis) {
        int layer = layerEntry.first;
        const LayerBlockedAnalysis& analysis = layerEntry.second;

        // ��ȡ�ò������region�����С
        std::map<int, int> regionSizes; // regionID -> ��С
        for (int j = 0; j < grid.height; j++) {
            for (int k = 0; k < grid.width; k++) {
                int regionID = grid.grid[layer][j][k].regionID;
                if (regionID > 0) {
                    regionSizes[regionID]++;
                }
            }
        }

        // ���ͬһ��������region��
        std::vector<int> regions;
        for (const auto& entry : regionSizes) {
            regions.push_back(entry.first);
        }

        for (size_t i = 0; i < regions.size(); i++) {
            for (size_t j = i + 1; j < regions.size(); j++) {
                int region1 = regions[i];
                int region2 = regions[j];
                std::pair<int, int> pair = std::make_pair(region1, region2);
                std::pair<int, int> reversePair = std::make_pair(region2, region1);

                if (processedPairs.find(pair) != processedPairs.end() ||
                    processedPairs.find(reversePair) != processedPairs.end()) {
                    continue;
                }

                int region1Size = regionSizes[region1];
                int region2Size = regionSizes[region2];

                // ����Ƿ���������ָ�
                if (analysis.purpleDashedPairs.find(pair) != analysis.purpleDashedPairs.end() ||
                    analysis.purpleDashedPairs.find(reversePair) != analysis.purpleDashedPairs.end()) {
                    // ��ɫ���ߣ�����������ָ�
                    int blockedSize = getBlockedAreaSize(analysis, region1, region2);
                    std::stringstream ss;
                    ss << "    R" << region1 << "--��ɫ����--R" << region2
                        << "  value:R" << region1 << "(" << region1Size << " cells), R"
                        << region2 << "(" << region2Size << " cells), cost:" << blockedSize << " cells";
                    connections.push_back(ss.str());
                }
                else {
                    // ��ɫ���ߣ�û�б���������ָ�
                    std::stringstream ss;
                    ss << "    R" << region1 << "--��ɫ����--R" << region2
                        << "  value:R" << region1 << "(" << region1Size << " cells), R"
                        << region2 << "(" << region2Size << " cells), cost:0 cells";
                    connections.push_back(ss.str());
                }

                processedPairs.insert(pair);
            }
        }
    }

    // 2. �������²�֮������ӹ�ϵ����׽��飩
    if (blockedAnalysis.size() >= 2) {
        // ʹ���ص�����������ȡ��׽���
        auto overlapResult = OverlapAnalyzer::analyzeLayerOverlap(grid);

        // ��ȡ����������С
        std::map<int, std::map<int, int>> layerRegionSizes;
        for (int layer = 0; layer < 2 && layer < grid.Layers; layer++) {
            for (int j = 0; j < grid.height; j++) {
                for (int k = 0; k < grid.width; k++) {
                    int regionID = grid.grid[layer][j][k].regionID;
                    if (regionID > 0) {
                        layerRegionSizes[layer][regionID]++;
                    }
                }
            }
        }

        for (const auto& recommendation : overlapResult.recommendations) {
            int region0 = recommendation.layer0Region;
            int region1 = recommendation.layer1Region;
            std::pair<int, int> pair = std::make_pair(region0, region1);
            std::pair<int, int> reversePair = std::make_pair(region1, region0);

            if (processedPairs.find(pair) == processedPairs.end() &&
                processedPairs.find(reversePair) == processedPairs.end()) {

                // ��ȡ�����С
                int region0Size = 0, region1Size = 0;
                if (layerRegionSizes.find(0) != layerRegionSizes.end() &&
                    layerRegionSizes[0].find(region0) != layerRegionSizes[0].end()) {
                    region0Size = layerRegionSizes[0][region0];
                }
                if (layerRegionSizes.find(1) != layerRegionSizes.end() &&
                    layerRegionSizes[1].find(region1) != layerRegionSizes[1].end()) {
                    region1Size = layerRegionSizes[1][region1];
                }

                // ��ɫ���ߣ����Դ������
                std::stringstream ss;
                ss << "    R" << region0 << "--��ɫ����--R" << region1
                    << "  value:R" << region0 << "(" << region0Size << " cells), R"
                    << region1 << "(" << region1Size << " cells), cost:"
                    << recommendation.overlapCellCount << " cells";
                connections.push_back(ss.str());
                processedPairs.insert(pair);
            }
        }
    }

    return connections;
}

// �����������������������ӹ�ϵ����
std::string BlockedAreaAnalyzer::generateFullReport(
    const std::map<int, LayerBlockedAnalysis>& blockedAnalysis,
    const Grid& grid) {

    std::ostringstream report;

    // 1. ԭ�е����������������
    report << generateReport(blockedAnalysis);

    // 2. ���������ӹ�ϵ����
    report << "\n=== �������ӹ�ϵ ===" << std::endl;
    report << "ͼ��˵��:" << std::endl;
    report << "  ��ɫ����: ͬһ���������������������" << std::endl;
    report << "  ��ɫ����: ͬһ����������û�б������������" << std::endl;
    report << "  ��ɫ����: ���²������������ͨ���������" << std::endl;
    report << std::endl;

    auto connections = generateConnectionRelations(blockedAnalysis, grid);

    if (!connections.empty()) {
        report << "���ӹ�ϵ�б�:" << std::endl;
        for (const auto& connection : connections) {
            report << connection << std::endl;
        }
    }
    else {
        report << "û���ҵ��������ӹ�ϵ" << std::endl;
    }

    return report.str();
}