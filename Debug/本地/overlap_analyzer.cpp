#include "overlap_analyzer.h"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <queue>

// 检查单元格是否为GND
bool OverlapAnalyzer::isCellGnd(const Grid& grid, int layer, int row, int col) {
    if (row < 0 || row >= grid.height || col < 0 || col >= grid.width) {
        return false;
    }

    for (auto node : grid.grid[layer][row][col].nodes) {
        if (node) {
            // 检查segment的net名称
            if (node->name == "segment") {
                for (const auto& child : node->children) {
                    if (child->name == "net") {
                        int netID = std::stoi(child->parameters[0]);
                        auto it = grid.NetNamewithID.find(netID);
                        if (it != grid.NetNamewithID.end() &&
                            it->second.find("GND") != std::string::npos) {
                            return true;
                        }
                    }
                }
            }
            // 检查via的net名称
            else if (node->name == "via") {
                for (const auto& child : node->children) {
                    if (child->name == "net") {
                        int netID = std::stoi(child->parameters[0]);
                        auto it = grid.NetNamewithID.find(netID);
                        if (it != grid.NetNamewithID.end() &&
                            it->second.find("GND") != std::string::npos) {
                            return true;
                        }
                    }
                }
            }
            // 检查pad的net名称
            else if (node->name == "pad") {
                for (const auto& child : node->children) {
                    if (child->name == "net") {
                        std::string netName = child->parameters[0];
                        if (netName.find("GND") != std::string::npos) {
                            return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

// 统计区域内的GND信息
void OverlapAnalyzer::getGndRegionInfo(const Grid& grid, int layer, int regionID, int& gndCellCount, double& gndRatio) {
    gndCellCount = 0;
    int regionCellCount = 0;

    for (int j = 0; j < grid.height; j++) {
        for (int k = 0; k < grid.width; k++) {
            if (grid.grid[layer][j][k].regionID == regionID) {
                regionCellCount++;
                if (isCellGnd(grid, layer, j, k)) {
                    gndCellCount++;
                }
            }
        }
    }

    gndRatio = (regionCellCount > 0) ? static_cast<double>(gndCellCount) / regionCellCount : 0.0;
}

// 检查区域周边是否存在GND
bool OverlapAnalyzer::hasSurroundingGnd(const Grid& grid, int layer, int regionID, int searchRadius, double gndThreshold) {
    std::set<std::pair<int, int>> regionCells;
    std::set<std::pair<int, int>> borderCells;

    // 收集区域内的所有单元格
    for (int j = 0; j < grid.height; j++) {
        for (int k = 0; k < grid.width; k++) {
            if (grid.grid[layer][j][k].regionID == regionID) {
                regionCells.insert({ j, k });
            }
        }
    }

    // 收集区域边界周围的单元格
    for (const auto& cell : regionCells) {
        int j = cell.first;
        int k = cell.second;

        // 检查四个方向
        const int dr[4] = { -1, 1, 0, 0 };
        const int dc[4] = { 0, 0, -1, 1 };

        for (int d = 0; d < 4; d++) {
            int nj = j + dr[d];
            int nk = k + dc[d];

            if (nj >= 0 && nj < grid.height && nk >= 0 && nk < grid.width) {
                if (regionCells.find({ nj, nk }) == regionCells.end()) {
                    borderCells.insert({ nj, nk });
                }
            }
        }
    }

    // 统计边界周围的GND单元格
    int gndCount = 0;
    int totalBorderCells = borderCells.size();

    for (const auto& cell : borderCells) {
        if (isCellGnd(grid, layer, cell.first, cell.second)) {
            gndCount++;
        }
    }

    double gndRatio = (totalBorderCells > 0) ? static_cast<double>(gndCount) / totalBorderCells : 0.0;
    return gndRatio >= gndThreshold;
}

// 增强版GND检测：包括区域内和周边
bool OverlapAnalyzer::hasGndElements(const Grid& grid, int layer, int regionID, double gndThreshold) {
    // 方法1：检查区域内是否有GND
    int gndCellCount;
    double gndRatio;
    getGndRegionInfo(grid, layer, regionID, gndCellCount, gndRatio);

    // 如果区域内GND比例达到阈值，直接返回true
    if (gndRatio >= gndThreshold) {
        return true;
    }

    // 方法2：检查区域周边是否有足够多的GND
    if (hasSurroundingGnd(grid, layer, regionID, 3, 0.3)) {
        return true;
    }

    return false;
}

// 获取区域中的GND元素类型
std::vector<std::string> OverlapAnalyzer::getGndElementTypes(const Grid& grid, int layer, int regionID) {
    std::vector<std::string> elementTypes;
    std::set<std::string> uniqueTypes;

    for (int j = 0; j < grid.height; j++) {
        for (int k = 0; k < grid.width; k++) {
            if (grid.grid[layer][j][k].regionID == regionID) {
                for (auto node : grid.grid[layer][j][k].nodes) {
                    if (node) {
                        std::string elementType = node->name;
                        bool isGnd = false;

                        if (elementType == "segment") {
                            for (const auto& child : node->children) {
                                if (child->name == "net") {
                                    int netID = std::stoi(child->parameters[0]);
                                    auto it = grid.NetNamewithID.find(netID);
                                    if (it != grid.NetNamewithID.end() &&
                                        it->second.find("GND") != std::string::npos) {
                                        isGnd = true;
                                        break;
                                    }
                                }
                            }
                        }
                        else if (elementType == "via") {
                            for (const auto& child : node->children) {
                                if (child->name == "net") {
                                    int netID = std::stoi(child->parameters[0]);
                                    auto it = grid.NetNamewithID.find(netID);
                                    if (it != grid.NetNamewithID.end() &&
                                        it->second.find("GND") != std::string::npos) {
                                        isGnd = true;
                                        break;
                                    }
                                }
                            }
                        }
                        else if (elementType == "pad") {
                            for (const auto& child : node->children) {
                                if (child->name == "net") {
                                    std::string netName = child->parameters[0];
                                    if (netName.find("GND") != std::string::npos) {
                                        isGnd = true;
                                        break;
                                    }
                                }
                            }
                        }

                        if (isGnd && uniqueTypes.find(elementType) == uniqueTypes.end()) {
                            uniqueTypes.insert(elementType);
                            elementTypes.push_back(elementType);
                        }
                    }
                }
            }
        }
    }

    return elementTypes;
}

OverlapAnalysisResult OverlapAnalyzer::analyzeLayerOverlap(
    const Grid& grid,
    int layer0,
    int layer1,
    int minOverlapCells) {  // 修改：移除overlapThresholdRatio参数

    OverlapAnalysisResult result;
    result.totalBoardArea = grid.width * grid.height;
    result.minOverlapCells = minOverlapCells;  // 修改：只保留最小重叠单元格数

    // 统计每层的区域集合和区域大小
    std::set<int> layer0Regions, layer1Regions;
    std::map<int, int> layer0RegionSize, layer1RegionSize;

    for (int j = 0; j < grid.height; j++) {
        for (int k = 0; k < grid.width; k++) {
            if (grid.grid[layer0][j][k].regionID > 0) {
                layer0Regions.insert(grid.grid[layer0][j][k].regionID);
                layer0RegionSize[grid.grid[layer0][j][k].regionID]++;
            }
            if (grid.grid[layer1][j][k].regionID > 0) {
                layer1Regions.insert(grid.grid[layer1][j][k].regionID);
                layer1RegionSize[grid.grid[layer1][j][k].regionID]++;
            }
        }
    }

    // 检查第二层（Layer 1）哪些区域包含GND元素（增强检测）
    for (int region1 : layer1Regions) {
        if (hasGndElements(grid, layer1, region1, 0.05)) {
            OverlapAnalysisResult::GndRegionInfo gndInfo;
            gndInfo.regionID = region1;
            gndInfo.layer = layer1;
            gndInfo.size = layer1RegionSize[region1];
            getGndRegionInfo(grid, layer1, region1, gndInfo.gndCellCount, gndInfo.gndRatio);
            gndInfo.gndElements = getGndElementTypes(grid, layer1, region1);
            result.gndRegions.push_back(gndInfo);
        }
    }

    // 统计所有区域对的重叠情况
    std::map<std::pair<int, int>, int> overlapCount;
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> overlapPositions;

    for (int j = 0; j < grid.height; j++) {
        for (int k = 0; k < grid.width; k++) {
            int region0 = grid.grid[layer0][j][k].regionID;
            int region1 = grid.grid[layer1][j][k].regionID;

            if (region0 > 0 && region1 > 0 &&
                !grid.grid[layer0][j][k].isblocked &&
                !grid.grid[layer1][j][k].isblocked) {

                std::pair<int, int> regionPair = std::make_pair(region0, region1);
                overlapCount[regionPair]++;
                overlapPositions[regionPair].push_back({ j, k });
            }
        }
    }

    result.allOverlaps = overlapCount;

    // 找出符合打孔条件的区域对
    // 新条件：只有第二层区域包含GND元素且重叠单元格数≥10
    for (const auto& entry : overlapCount) {
        int region0 = entry.first.first;
        int region1 = entry.first.second;
        int overlapCountValue = entry.second;

        // 检查第二层区域是否包含GND元素（增强检测）
        bool hasGnd = hasGndElements(grid, layer1, region1, 0.05);

        double overlapRatioToBoard = (result.totalBoardArea > 0) ?
            static_cast<double>(overlapCountValue) / result.totalBoardArea : 0.0;

        bool meetsCellCountCriteria = (overlapCountValue >= result.minOverlapCells);

        // 新条件：必须包含GND元素且满足单元格数量条件
        if (hasGnd && meetsCellCountCriteria) {
            OverlapAnalysisResult::ViaRecommendation recommendation;
            recommendation.layer0Region = region0;
            recommendation.layer1Region = region1;
            recommendation.overlapCellCount = overlapCountValue;
            recommendation.overlapRatioToBoard = overlapRatioToBoard;

            // 收集前10个重叠位置用于显示
            const auto& positions = overlapPositions[entry.first];
            int maxToCollect = std::min(10, static_cast<int>(positions.size()));
            for (int i = 0; i < maxToCollect; i++) {
                recommendation.viaPositions.push_back(positions[i]);
            }

            result.recommendations.push_back(recommendation);
        }
    }

    return result;
}

std::string OverlapAnalyzer::generateReport(const OverlapAnalysisResult& result, const Grid& grid) {
    std::ostringstream report;

    report << "=== 重叠区域分析和打孔建议 ===" << std::endl;
    report << "PCB板网格尺寸: " << grid.width << " × " << grid.height
        << " = " << result.totalBoardArea << " 个单元格" << std::endl;

    // 统计每层的区域集合
    std::set<int> layer0Regions, layer1Regions;
    std::map<int, int> layer0RegionSize, layer1RegionSize;

    for (int j = 0; j < grid.height; j++) {
        for (int k = 0; k < grid.width; k++) {
            if (grid.grid[0][j][k].regionID > 0) {
                layer0Regions.insert(grid.grid[0][j][k].regionID);
                layer0RegionSize[grid.grid[0][j][k].regionID]++;
            }
            if (grid.grid[1][j][k].regionID > 0) {
                layer1Regions.insert(grid.grid[1][j][k].regionID);
                layer1RegionSize[grid.grid[1][j][k].regionID]++;
            }
        }
    }

    report << "Layer 0 区域: ";
    for (int region : layer0Regions) report << region << " ";
    report << std::endl;

    report << "Layer 1 区域: ";
    for (int region : layer1Regions) report << region << " ";
    report << std::endl;

    // 显示包含GND元素的区域（增强信息）
    if (!result.gndRegions.empty()) {
        report << "\n=== 包含GND元素的区域 (增强检测) ===" << std::endl;
        report << "检测方法：1. 区域内GND单元格比例 ≥ 5%  2. 区域周边GND包围比例 ≥ 30%" << std::endl;
        for (const auto& gndRegion : result.gndRegions) {
            report << "Layer " << gndRegion.layer << " 区域 " << gndRegion.regionID
                << " (大小:" << gndRegion.size << ", GND单元格:" << gndRegion.gndCellCount
                << ", GND比例:" << std::fixed << std::setprecision(1) << gndRegion.gndRatio * 100 << "%)";
            if (!gndRegion.gndElements.empty()) {
                report << " [";
                for (size_t i = 0; i < gndRegion.gndElements.size(); i++) {
                    if (i > 0) report << ", ";
                    report << gndRegion.gndElements[i];
                }
                report << "]";
            }
            report << std::endl;
        }
    }
    else {
        report << "\n=== 包含GND元素的区域 ===" << std::endl;
        report << "Layer 1 中没有发现包含GND元素的区域" << std::endl;
        report << "（使用增强检测：区域内GND比例≥5% 或 周边GND包围比例≥30%）" << std::endl;
    }

    report << "\n打孔阈值条件:" << std::endl;
    report << "  - PCB板总面积: " << result.totalBoardArea << " 个单元格" << std::endl;
    report << "  - 最小重叠单元格: " << result.minOverlapCells << " 个" << std::endl;
    report << "  - GND检测: 区域内GND比例≥5% 或 周边GND包围比例≥30%" << std::endl;

    report << "\n区域重叠统计（所有区域对）：" << std::endl;

    bool hasViaRecommendations = false;

    for (const auto& entry : result.allOverlaps) {
        int region0 = entry.first.first;
        int region1 = entry.first.second;
        int overlapCountValue = entry.second;

        double overlapRatioToBoard = (result.totalBoardArea > 0) ?
            static_cast<double>(overlapCountValue) / result.totalBoardArea : 0.0;

        report << "区域对 (" << region0 << " -> " << region1
            << "): " << overlapCountValue << " 个重叠单元格"
            << " (区域" << region0 << "大小: " << layer0RegionSize[region0]
            << ", 区域" << region1 << "大小: " << layer1RegionSize[region1]
            << ", 占PCB总面积: " << std::fixed << std::setprecision(2)
            << overlapRatioToBoard * 100 << "%)";

        // 检查第二层区域是否包含GND元素
        bool hasGnd = OverlapAnalyzer::hasGndElements(grid, 1, region1);

        // 双重条件判断：包含GND元素、重叠单元格数量≥10
        bool meetsCellCountCriteria = (overlapCountValue >= result.minOverlapCells);
        bool meetsGndCriteria = hasGnd;

        if (meetsGndCriteria && meetsCellCountCriteria) {
            report << " ? 建议打孔连接 (区域包含GND且重叠单元格≥10个)";
            hasViaRecommendations = true;
        }
        else {
            report << " ? 不满足打孔条件 (";
            std::vector<std::string> reasons;
            if (!meetsGndCriteria) {
                reasons.push_back("区域不包含GND元素");
            }
            if (!meetsCellCountCriteria) {
                reasons.push_back("单元格数<10");
            }

            for (size_t i = 0; i < reasons.size(); i++) {
                if (i > 0) report << "且";
                report << reasons[i];
            }
            report << ")";
        }
        report << std::endl;
    }

    // 输出具体的打孔位置建议
    if (hasViaRecommendations) {
        report << "\n=== 具体打孔位置建议 (区域包含GND元素且重叠单元格≥10个) ===" << std::endl;

        for (const auto& rec : result.recommendations) {
            report << "区域 " << rec.layer0Region << " (Layer 0, 大小:" << layer0RegionSize[rec.layer0Region]
                << ") -> 区域 " << rec.layer1Region << " (Layer 1, 大小:" << layer1RegionSize[rec.layer1Region]
                << "): " << rec.overlapCellCount << " 个重叠单元格"
                << " (占PCB总面积: " << std::fixed << std::setprecision(2)
                << rec.overlapRatioToBoard * 100 << "%)" << std::endl;

            int viaCount = 0;
            const int MAX_VIAS_TO_SHOW = 10;

            for (size_t i = 0; i < rec.viaPositions.size() && viaCount < MAX_VIAS_TO_SHOW; i++) {
                int row = rec.viaPositions[i].first;
                int col = rec.viaPositions[i].second;
                double actualX = grid.min_x + (col + 0.5) / grid.inputScale;
                double actualY = grid.min_y + (row + 0.5) / grid.inputScale;

                report << "  位置 " << viaCount + 1 << ": 坐标("
                    << std::fixed << std::setprecision(3) << actualX << ", " << actualY
                    << "), 网格(" << row << ", " << col << ")" << std::endl;
                viaCount++;
            }

            if (viaCount == MAX_VIAS_TO_SHOW && rec.overlapCellCount > MAX_VIAS_TO_SHOW) {
                report << "  ... (还有 " << (rec.overlapCellCount - MAX_VIAS_TO_SHOW) << " 个位置)" << std::endl;
            }
            report << std::endl;
        }
    }
    else {
        report << "\n没有找到需要打孔连接的显著重叠区域。" << std::endl;
        report << "原因：没有区域对同时满足以下两个条件：" << std::endl;
        report << "  1. Layer 1区域包含GND元素" << std::endl;
        report << "  2. 重叠单元格数量 ≥ " << result.minOverlapCells << " 个" << std::endl;
    }

    // 输出区域统计信息
    report << "\n=== 区域统计 ===" << std::endl;
    report << "PCB板网格尺寸: " << grid.width << " × " << grid.height << " = " << result.totalBoardArea << " 个单元格" << std::endl;
    report << "Layer 0 实际区域数量: " << layer0Regions.size() << std::endl;
    report << "Layer 1 实际区域数量: " << layer1Regions.size() << std::endl;
    report << "Layer 1 包含GND元素的区域数量: " << result.gndRegions.size() << std::endl;
    report << "发现的重叠区域对数量: " << result.allOverlaps.size() << std::endl;
    report << "打孔阈值条件:" << std::endl;
    report << "  1. Layer 1区域必须包含GND元素" << std::endl;
    report << "  2. 重叠单元格数量 ≥ " << result.minOverlapCells << " 个" << std::endl;

    // 按区域分组显示重叠情况
    report << "\n=== 按区域分组的重叠情况 ===" << std::endl;
    for (int region0 : layer0Regions) {
        report << "Layer 0 区域 " << region0 << " (大小:" << layer0RegionSize[region0]
            << ") 与 Layer 1 的重叠情况:" << std::endl;
        bool hasOverlap = false;

        for (int region1 : layer1Regions) {
            std::pair<int, int> regionPair = std::make_pair(region0, region1);
            if (result.allOverlaps.find(regionPair) != result.allOverlaps.end()) {
                int overlapCountValue = result.allOverlaps.at(regionPair);
                double ratioToBoard = (result.totalBoardArea > 0) ?
                    static_cast<double>(overlapCountValue) / result.totalBoardArea : 0.0;

                bool hasGnd = OverlapAnalyzer::hasGndElements(grid, 1, region1);
                bool meetsCellCountCriteria = (overlapCountValue >= result.minOverlapCells);
                bool meetsAllCriteria = hasGnd && meetsCellCountCriteria;

                report << "  -> 区域 " << region1 << " (大小:" << layer1RegionSize[region1]
                    << ", 包含GND:" << (hasGnd ? "是" : "否")
                    << "): " << overlapCountValue << " 个重叠单元格"
                    << " (占PCB总面积: " << std::fixed << std::setprecision(2) << ratioToBoard * 100 << "%)";
                if (meetsAllCriteria) {
                    report << " ? 符合打孔条件";
                }
                else {
                    report << " ? 不符合";
                }
                report << std::endl;
                hasOverlap = true;
            }
        }

        if (!hasOverlap) {
            report << "  无重叠区域" << std::endl;
        }
        report << std::endl;
    }

    return report.str();
}

std::vector<std::pair<double, double>> OverlapAnalyzer::getOverlapPositions(
    const Grid& grid,
    int layer0, int region0,
    int layer1, int region1,
    int maxPositions) {

    std::vector<std::pair<double, double>> positions;

    for (int j = 0; j < grid.height && positions.size() < static_cast<size_t>(maxPositions); j++) {
        for (int k = 0; k < grid.width && positions.size() < static_cast<size_t>(maxPositions); k++) {
            if (grid.grid[layer0][j][k].regionID == region0 &&
                grid.grid[layer1][j][k].regionID == region1 &&
                !grid.grid[layer0][j][k].isblocked &&
                !grid.grid[layer1][j][k].isblocked) {

                double actualX = grid.min_x + (k + 0.5) / grid.inputScale;
                double actualY = grid.min_y + (j + 0.5) / grid.inputScale;
                positions.push_back({ actualX, actualY });
            }
        }
    }

    return positions;
}