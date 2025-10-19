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
    area.size = 0;  // 明确初始化大小

    // 检查起始单元格是否有效且是阻塞区域
    if (!isValidCell(grid, startRow, startCol) || visited[startRow][startCol] ||
        !grid.grid[layer][startRow][startCol].isblocked) {
        return area;  // 返回空的阻塞区域信息
    }

    std::queue<std::pair<int, int>> q;
    q.push({ startRow, startCol });
    visited[startRow][startCol] = true;

    // 四个方向：上、下、左、右
    const int dr[4] = { -1, 1, 0, 0 };
    const int dc[4] = { 0, 0, -1, 1 };

    while (!q.empty()) {
        auto [row, col] = q.front();
        q.pop();

        // 添加到阻塞区域
        area.cells.push_back({ row, col });
        area.size++;

        // 检查四个方向
        for (int i = 0; i < 4; i++) {
            int newRow = row + dr[i];
            int newCol = col + dc[i];

            // 关键修正：必须检查新单元格是阻塞区域才加入
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

        // 检查四个方向（不包括对角线）
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

    // 对阻塞区域的每个单元格，计算到相邻region的最小距离
    for (const auto& blockedCell : blockedArea.cells) {
        int blockedRow = blockedCell.first;
        int blockedCol = blockedCell.second;

        // 搜索周围的region单元格
        for (int searchRadius = 1; searchRadius <= 10; searchRadius++) {
            bool foundAll = true;

            for (int regionID : blockedArea.adjacentRegions) {
                if (minDistances[regionID] <= searchRadius - 1) continue;

                bool foundRegion = false;

                // 搜索半径为searchRadius的方形区域
                for (int dr = -searchRadius; dr <= searchRadius && !foundRegion; dr++) {
                    for (int dc = -searchRadius; dc <= searchRadius && !foundRegion; dc++) {
                        int newRow = blockedRow + dr;
                        int newCol = blockedCol + dc;

                        if (isValidCell(grid, newRow, newCol) &&
                            grid.grid[layer][newRow][newCol].regionID == regionID) {
                            // 使用曼哈顿距离
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

    // 调试：统计实际的阻塞单元格
    int actualBlockedCount = 0;
    for (int row = 0; row < grid.height; row++) {
        for (int col = 0; col < grid.width; col++) {
            if (grid.grid[layer][row][col].isblocked) {
                actualBlockedCount++;
            }
        }
    }
    std::cout << "Layer " << layer << " 实际阻塞单元格: " << actualBlockedCount << std::endl;

    // 寻找所有阻塞区域
    for (int row = 0; row < grid.height; row++) {
        for (int col = 0; col < grid.width; col++) {
            if (!visited[row][col] && grid.grid[layer][row][col].isblocked) {
                BlockedAreaInfo area = findBlockedArea(grid, layer, row, col, visited, areaID);

                // 关键修正：只有当找到有效的阻塞区域时才添加到结果中
                if (area.size > 0 && area.size <= actualBlockedCount) {
                    std::cout << "找到阻塞区域 " << areaID << ": 大小=" << area.size << " 单元格" << std::endl;

                    // 分析相邻regions
                    findAdjacentRegions(grid, layer, area);
                    calculateRegionProximity(grid, layer, area);
                    analysis.blockedAreas.push_back(area);
                    areaID++;
                }
                else if (area.size > actualBlockedCount) {
                    std::cout << "忽略无效的阻塞区域 " << areaID << ": 大小=" << area.size
                        << " (大于实际阻塞单元格数 " << actualBlockedCount << ")" << std::endl;
                }
            }
        }
    }

    // 分析哪些regions被阻塞区域分开
    for (const auto& blockedArea : analysis.blockedAreas) {
        if (blockedArea.adjacentRegions.size() >= 2) {
            std::vector<int> regions(blockedArea.adjacentRegions.begin(),
                blockedArea.adjacentRegions.end());

            for (size_t i = 0; i < regions.size(); i++) {
                for (size_t j = i + 1; j < regions.size(); j++) {
                    std::pair<int, int> regionPair = std::make_pair(regions[i], regions[j]);
                    std::stringstream ss;
                    ss << "区域" << regions[i] << "和区域" << regions[j]
                        << "被阻塞区域" << blockedArea.areaID << "分开";
                    analysis.separationInfo[regionPair] = ss.str();

                    // 记录被紫色虚线分隔的区域对
                    analysis.purpleDashedPairs.insert(regionPair);
                }
            }
        }
    }

    return analysis;
}

std::map<int, LayerBlockedAnalysis> BlockedAreaAnalyzer::analyzeAllLayers(const Grid& grid) {
    std::map<int, LayerBlockedAnalysis> allAnalysis;

    // 只分析前几层，避免输出过多信息
    int maxLayersToAnalyze = std::min(5, grid.Layers);
    for (int layer = 0; layer < maxLayersToAnalyze; layer++) {
        std::cout << "分析 Layer " << layer << "..." << std::endl;
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

    report << "=== 阻塞区域分析报告 ===" << std::endl;
    report << "PCB总层数: " << analysis.size() << std::endl << std::endl;

    for (const auto& layerEntry : analysis) {
        int layerID = layerEntry.first;
        const LayerBlockedAnalysis& layerAnalysis = layerEntry.second;

        report << "Layer " << layerID << ":" << std::endl;
        report << "  阻塞区域个数: " << layerAnalysis.blockedAreas.size() << std::endl;

        if (!layerAnalysis.blockedAreas.empty()) {
            report << "  阻塞区域详情:" << std::endl;
            for (const auto& blockedArea : layerAnalysis.blockedAreas) {
                report << "    阻塞区域 " << blockedArea.areaID
                    << ": 大小=" << blockedArea.size
                    << " 单元格, 相邻regions=";

                if (!blockedArea.adjacentRegions.empty()) {
                    for (auto it = blockedArea.adjacentRegions.begin();
                        it != blockedArea.adjacentRegions.end(); ++it) {
                        if (it != blockedArea.adjacentRegions.begin()) report << ", ";
                        report << *it;
                    }

                    // 显示与每个region的最小距离
                    report << std::endl << "      与相邻region的最小距离: ";
                    for (const auto& proximity : blockedArea.regionProximity) {
                        report << "区域" << proximity.first << "->" << proximity.second << "单元格 ";
                    }
                }
                else {
                    report << "无";
                }
                report << std::endl;
            }
        }
        else {
            report << "  本层没有阻塞区域" << std::endl;
        }

        // 显示被分开的regions
        if (!layerAnalysis.separationInfo.empty()) {
            report << "  被阻塞区域分开的regions:" << std::endl;
            for (const auto& separation : layerAnalysis.separationInfo) {
                report << "    " << separation.second << std::endl;
            }
        }
        else {
            report << "  没有被阻塞区域分开的regions" << std::endl;
        }

        report << std::endl;
    }

    // 统计信息
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

    report << "=== 统计摘要 ===" << std::endl;
    report << "总阻塞区域数: " << totalBlockedAreas << std::endl;
    report << "总阻塞单元格数: " << totalBlockedCells << std::endl;
    report << "总region分离对数: " << totalSeparations << std::endl;

    return report.str();
}

// 新增函数：获取阻塞区域的大小
static int getBlockedAreaSize(const LayerBlockedAnalysis& layerAnalysis, int region1, int region2) {
    for (const auto& blockedArea : layerAnalysis.blockedAreas) {
        if (blockedArea.adjacentRegions.find(region1) != blockedArea.adjacentRegions.end() &&
            blockedArea.adjacentRegions.find(region2) != blockedArea.adjacentRegions.end()) {
            return blockedArea.size;
        }
    }
    return 0;
}

// 修改后的连接关系生成函数
std::vector<std::string> BlockedAreaAnalyzer::generateConnectionRelations(
    const std::map<int, LayerBlockedAnalysis>& blockedAnalysis,
    const Grid& grid) {

    std::vector<std::string> connections;
    std::set<std::pair<int, int>> processedPairs; // 避免重复

    // 1. 分析同一层内的连接关系
    for (const auto& layerEntry : blockedAnalysis) {
        int layer = layerEntry.first;
        const LayerBlockedAnalysis& analysis = layerEntry.second;

        // 获取该层的所有region及其大小
        std::map<int, int> regionSizes; // regionID -> 大小
        for (int j = 0; j < grid.height; j++) {
            for (int k = 0; k < grid.width; k++) {
                int regionID = grid.grid[layer][j][k].regionID;
                if (regionID > 0) {
                    regionSizes[regionID]++;
                }
            }
        }

        // 检查同一层内所有region对
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

                // 检查是否被阻塞区域分隔
                if (analysis.purpleDashedPairs.find(pair) != analysis.purpleDashedPairs.end() ||
                    analysis.purpleDashedPairs.find(reversePair) != analysis.purpleDashedPairs.end()) {
                    // 紫色虚线：被阻塞区域分隔
                    int blockedSize = getBlockedAreaSize(analysis, region1, region2);
                    std::stringstream ss;
                    ss << "    R" << region1 << "--紫色虚线--R" << region2
                        << "  value:R" << region1 << "(" << region1Size << " cells), R"
                        << region2 << "(" << region2Size << " cells), cost:" << blockedSize << " cells";
                    connections.push_back(ss.str());
                }
                else {
                    // 黑色虚线：没有被阻塞区域分隔
                    std::stringstream ss;
                    ss << "    R" << region1 << "--黑色虚线--R" << region2
                        << "  value:R" << region1 << "(" << region1Size << " cells), R"
                        << region2 << "(" << region2Size << " cells), cost:0 cells";
                    connections.push_back(ss.str());
                }

                processedPairs.insert(pair);
            }
        }
    }

    // 2. 分析上下层之间的连接关系（打孔建议）
    if (blockedAnalysis.size() >= 2) {
        // 使用重叠分析器来获取打孔建议
        auto overlapResult = OverlapAnalyzer::analyzeLayerOverlap(grid);

        // 获取各层的区域大小
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

                // 获取区域大小
                int region0Size = 0, region1Size = 0;
                if (layerRegionSizes.find(0) != layerRegionSizes.end() &&
                    layerRegionSizes[0].find(region0) != layerRegionSizes[0].end()) {
                    region0Size = layerRegionSizes[0][region0];
                }
                if (layerRegionSizes.find(1) != layerRegionSizes.end() &&
                    layerRegionSizes[1].find(region1) != layerRegionSizes[1].end()) {
                    region1Size = layerRegionSizes[1][region1];
                }

                // 黄色虚线：可以打孔连接
                std::stringstream ss;
                ss << "    R" << region0 << "--黄色虚线--R" << region1
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

// 新增函数：生成完整的连接关系报告
std::string BlockedAreaAnalyzer::generateFullReport(
    const std::map<int, LayerBlockedAnalysis>& blockedAnalysis,
    const Grid& grid) {

    std::ostringstream report;

    // 1. 原有的阻塞区域分析报告
    report << generateReport(blockedAnalysis);

    // 2. 新增的连接关系部分
    report << "\n=== 区域连接关系 ===" << std::endl;
    report << "图例说明:" << std::endl;
    report << "  紫色虚线: 同一层两个区域被阻塞区域隔开" << std::endl;
    report << "  黑色虚线: 同一层两个区域没有被阻塞区域隔开" << std::endl;
    report << "  黄色虚线: 上下层两个区域可以通过打孔连接" << std::endl;
    report << std::endl;

    auto connections = generateConnectionRelations(blockedAnalysis, grid);

    if (!connections.empty()) {
        report << "连接关系列表:" << std::endl;
        for (const auto& connection : connections) {
            report << connection << std::endl;
        }
    }
    else {
        report << "没有找到区域连接关系" << std::endl;
    }

    return report.str();
}