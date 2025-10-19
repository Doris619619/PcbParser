#include "overlap_analyzer.h"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <queue>

// ��鵥Ԫ���Ƿ�ΪGND
bool OverlapAnalyzer::isCellGnd(const Grid& grid, int layer, int row, int col) {
    if (row < 0 || row >= grid.height || col < 0 || col >= grid.width) {
        return false;
    }

    for (auto node : grid.grid[layer][row][col].nodes) {
        if (node) {
            // ���segment��net����
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
            // ���via��net����
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
            // ���pad��net����
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

// ͳ�������ڵ�GND��Ϣ
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

// ��������ܱ��Ƿ����GND
bool OverlapAnalyzer::hasSurroundingGnd(const Grid& grid, int layer, int regionID, int searchRadius, double gndThreshold) {
    std::set<std::pair<int, int>> regionCells;
    std::set<std::pair<int, int>> borderCells;

    // �ռ������ڵ����е�Ԫ��
    for (int j = 0; j < grid.height; j++) {
        for (int k = 0; k < grid.width; k++) {
            if (grid.grid[layer][j][k].regionID == regionID) {
                regionCells.insert({ j, k });
            }
        }
    }

    // �ռ�����߽���Χ�ĵ�Ԫ��
    for (const auto& cell : regionCells) {
        int j = cell.first;
        int k = cell.second;

        // ����ĸ�����
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

    // ͳ�Ʊ߽���Χ��GND��Ԫ��
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

// ��ǿ��GND��⣺���������ں��ܱ�
bool OverlapAnalyzer::hasGndElements(const Grid& grid, int layer, int regionID, double gndThreshold) {
    // ����1������������Ƿ���GND
    int gndCellCount;
    double gndRatio;
    getGndRegionInfo(grid, layer, regionID, gndCellCount, gndRatio);

    // ���������GND�����ﵽ��ֵ��ֱ�ӷ���true
    if (gndRatio >= gndThreshold) {
        return true;
    }

    // ����2����������ܱ��Ƿ����㹻���GND
    if (hasSurroundingGnd(grid, layer, regionID, 3, 0.3)) {
        return true;
    }

    return false;
}

// ��ȡ�����е�GNDԪ������
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
    int minOverlapCells) {  // �޸ģ��Ƴ�overlapThresholdRatio����

    OverlapAnalysisResult result;
    result.totalBoardArea = grid.width * grid.height;
    result.minOverlapCells = minOverlapCells;  // �޸ģ�ֻ������С�ص���Ԫ����

    // ͳ��ÿ������򼯺Ϻ������С
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

    // ���ڶ��㣨Layer 1����Щ�������GNDԪ�أ���ǿ��⣩
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

    // ͳ����������Ե��ص����
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

    // �ҳ����ϴ�������������
    // ��������ֻ�еڶ����������GNDԪ�����ص���Ԫ������10
    for (const auto& entry : overlapCount) {
        int region0 = entry.first.first;
        int region1 = entry.first.second;
        int overlapCountValue = entry.second;

        // ���ڶ��������Ƿ����GNDԪ�أ���ǿ��⣩
        bool hasGnd = hasGndElements(grid, layer1, region1, 0.05);

        double overlapRatioToBoard = (result.totalBoardArea > 0) ?
            static_cast<double>(overlapCountValue) / result.totalBoardArea : 0.0;

        bool meetsCellCountCriteria = (overlapCountValue >= result.minOverlapCells);

        // ���������������GNDԪ�������㵥Ԫ����������
        if (hasGnd && meetsCellCountCriteria) {
            OverlapAnalysisResult::ViaRecommendation recommendation;
            recommendation.layer0Region = region0;
            recommendation.layer1Region = region1;
            recommendation.overlapCellCount = overlapCountValue;
            recommendation.overlapRatioToBoard = overlapRatioToBoard;

            // �ռ�ǰ10���ص�λ��������ʾ
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

    report << "=== �ص���������ʹ�׽��� ===" << std::endl;
    report << "PCB������ߴ�: " << grid.width << " �� " << grid.height
        << " = " << result.totalBoardArea << " ����Ԫ��" << std::endl;

    // ͳ��ÿ������򼯺�
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

    report << "Layer 0 ����: ";
    for (int region : layer0Regions) report << region << " ";
    report << std::endl;

    report << "Layer 1 ����: ";
    for (int region : layer1Regions) report << region << " ";
    report << std::endl;

    // ��ʾ����GNDԪ�ص�������ǿ��Ϣ��
    if (!result.gndRegions.empty()) {
        report << "\n=== ����GNDԪ�ص����� (��ǿ���) ===" << std::endl;
        report << "��ⷽ����1. ������GND��Ԫ����� �� 5%  2. �����ܱ�GND��Χ���� �� 30%" << std::endl;
        for (const auto& gndRegion : result.gndRegions) {
            report << "Layer " << gndRegion.layer << " ���� " << gndRegion.regionID
                << " (��С:" << gndRegion.size << ", GND��Ԫ��:" << gndRegion.gndCellCount
                << ", GND����:" << std::fixed << std::setprecision(1) << gndRegion.gndRatio * 100 << "%)";
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
        report << "\n=== ����GNDԪ�ص����� ===" << std::endl;
        report << "Layer 1 ��û�з��ְ���GNDԪ�ص�����" << std::endl;
        report << "��ʹ����ǿ��⣺������GND������5% �� �ܱ�GND��Χ������30%��" << std::endl;
    }

    report << "\n�����ֵ����:" << std::endl;
    report << "  - PCB�������: " << result.totalBoardArea << " ����Ԫ��" << std::endl;
    report << "  - ��С�ص���Ԫ��: " << result.minOverlapCells << " ��" << std::endl;
    report << "  - GND���: ������GND������5% �� �ܱ�GND��Χ������30%" << std::endl;

    report << "\n�����ص�ͳ�ƣ���������ԣ���" << std::endl;

    bool hasViaRecommendations = false;

    for (const auto& entry : result.allOverlaps) {
        int region0 = entry.first.first;
        int region1 = entry.first.second;
        int overlapCountValue = entry.second;

        double overlapRatioToBoard = (result.totalBoardArea > 0) ?
            static_cast<double>(overlapCountValue) / result.totalBoardArea : 0.0;

        report << "����� (" << region0 << " -> " << region1
            << "): " << overlapCountValue << " ���ص���Ԫ��"
            << " (����" << region0 << "��С: " << layer0RegionSize[region0]
            << ", ����" << region1 << "��С: " << layer1RegionSize[region1]
            << ", ռPCB�����: " << std::fixed << std::setprecision(2)
            << overlapRatioToBoard * 100 << "%)";

        // ���ڶ��������Ƿ����GNDԪ��
        bool hasGnd = OverlapAnalyzer::hasGndElements(grid, 1, region1);

        // ˫�������жϣ�����GNDԪ�ء��ص���Ԫ��������10
        bool meetsCellCountCriteria = (overlapCountValue >= result.minOverlapCells);
        bool meetsGndCriteria = hasGnd;

        if (meetsGndCriteria && meetsCellCountCriteria) {
            report << " ? ���������� (�������GND���ص���Ԫ���10��)";
            hasViaRecommendations = true;
        }
        else {
            report << " ? ������������ (";
            std::vector<std::string> reasons;
            if (!meetsGndCriteria) {
                reasons.push_back("���򲻰���GNDԪ��");
            }
            if (!meetsCellCountCriteria) {
                reasons.push_back("��Ԫ����<10");
            }

            for (size_t i = 0; i < reasons.size(); i++) {
                if (i > 0) report << "��";
                report << reasons[i];
            }
            report << ")";
        }
        report << std::endl;
    }

    // �������Ĵ��λ�ý���
    if (hasViaRecommendations) {
        report << "\n=== ������λ�ý��� (�������GNDԪ�����ص���Ԫ���10��) ===" << std::endl;

        for (const auto& rec : result.recommendations) {
            report << "���� " << rec.layer0Region << " (Layer 0, ��С:" << layer0RegionSize[rec.layer0Region]
                << ") -> ���� " << rec.layer1Region << " (Layer 1, ��С:" << layer1RegionSize[rec.layer1Region]
                << "): " << rec.overlapCellCount << " ���ص���Ԫ��"
                << " (ռPCB�����: " << std::fixed << std::setprecision(2)
                << rec.overlapRatioToBoard * 100 << "%)" << std::endl;

            int viaCount = 0;
            const int MAX_VIAS_TO_SHOW = 10;

            for (size_t i = 0; i < rec.viaPositions.size() && viaCount < MAX_VIAS_TO_SHOW; i++) {
                int row = rec.viaPositions[i].first;
                int col = rec.viaPositions[i].second;
                double actualX = grid.min_x + (col + 0.5) / grid.inputScale;
                double actualY = grid.min_y + (row + 0.5) / grid.inputScale;

                report << "  λ�� " << viaCount + 1 << ": ����("
                    << std::fixed << std::setprecision(3) << actualX << ", " << actualY
                    << "), ����(" << row << ", " << col << ")" << std::endl;
                viaCount++;
            }

            if (viaCount == MAX_VIAS_TO_SHOW && rec.overlapCellCount > MAX_VIAS_TO_SHOW) {
                report << "  ... (���� " << (rec.overlapCellCount - MAX_VIAS_TO_SHOW) << " ��λ��)" << std::endl;
            }
            report << std::endl;
        }
    }
    else {
        report << "\nû���ҵ���Ҫ������ӵ������ص�����" << std::endl;
        report << "ԭ��û�������ͬʱ������������������" << std::endl;
        report << "  1. Layer 1�������GNDԪ��" << std::endl;
        report << "  2. �ص���Ԫ������ �� " << result.minOverlapCells << " ��" << std::endl;
    }

    // �������ͳ����Ϣ
    report << "\n=== ����ͳ�� ===" << std::endl;
    report << "PCB������ߴ�: " << grid.width << " �� " << grid.height << " = " << result.totalBoardArea << " ����Ԫ��" << std::endl;
    report << "Layer 0 ʵ����������: " << layer0Regions.size() << std::endl;
    report << "Layer 1 ʵ����������: " << layer1Regions.size() << std::endl;
    report << "Layer 1 ����GNDԪ�ص���������: " << result.gndRegions.size() << std::endl;
    report << "���ֵ��ص����������: " << result.allOverlaps.size() << std::endl;
    report << "�����ֵ����:" << std::endl;
    report << "  1. Layer 1����������GNDԪ��" << std::endl;
    report << "  2. �ص���Ԫ������ �� " << result.minOverlapCells << " ��" << std::endl;

    // �����������ʾ�ص����
    report << "\n=== �����������ص���� ===" << std::endl;
    for (int region0 : layer0Regions) {
        report << "Layer 0 ���� " << region0 << " (��С:" << layer0RegionSize[region0]
            << ") �� Layer 1 ���ص����:" << std::endl;
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

                report << "  -> ���� " << region1 << " (��С:" << layer1RegionSize[region1]
                    << ", ����GND:" << (hasGnd ? "��" : "��")
                    << "): " << overlapCountValue << " ���ص���Ԫ��"
                    << " (ռPCB�����: " << std::fixed << std::setprecision(2) << ratioToBoard * 100 << "%)";
                if (meetsAllCriteria) {
                    report << " ? ���ϴ������";
                }
                else {
                    report << " ? ������";
                }
                report << std::endl;
                hasOverlap = true;
            }
        }

        if (!hasOverlap) {
            report << "  ���ص�����" << std::endl;
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