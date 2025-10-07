// Separate2.cpp
// 坐标约定：像素 1-based，左上为(1,1)，行向下增大、列向右增大。
// 现实坐标：mm，y 向下增大（与屏幕一致）。

#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>
#include <array>
#include <climits>

#include "blocked_area_analyzer.h"   // 内含 grid.h

// ========== 像素中心 <-> 现实坐标 ==========
static inline std::pair<double, double>
PixelCenterPx1ToRealMM(const Grid& grid, int row1, int col1)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    row1 = std::max(1, std::min(row1, grid.height));
    col1 = std::max(1, std::min(col1, grid.width));
    const double s = (double)grid.inputScale;
    const double x = grid.min_x + ((double)col1 - 0.5) / s;
    const double y = grid.min_y + ((double)row1 - 0.5) / s; // 向下增大
    return { x, y };
}

static inline std::pair<int, int>
RealMMToPixelPx1(const Grid& grid, double x_mm, double y_mm)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    const double s = (double)grid.inputScale;
    int col1 = (int)std::floor((x_mm - grid.min_x) * s) + 1;
    int row1 = (int)std::floor((y_mm - grid.min_y) * s) + 1;
    col1 = std::max(1, std::min(col1, grid.width));
    row1 = std::max(1, std::min(row1, grid.height));
    return { row1, col1 };
}

// ========== 像素边界 -> 现实外沿矩形 ==========
struct RealRectMM { double x_left, y_top, x_right, y_bot; };

static inline RealRectMM PixelBoxEdgesPx1ToRealRect(const Grid& grid,
    int row_min, int col_min,
    int row_max, int col_max)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    const double s = (double)grid.inputScale;
    RealRectMM rr;
    rr.x_left = grid.min_x + ((double)col_min - 1.0) / s;
    rr.x_right = grid.min_x + ((double)col_max) / s;
    rr.y_top = grid.min_y + ((double)row_min - 1.0) / s;
    rr.y_bot = grid.min_y + ((double)row_max) / s;
    return rr;
}

// ========== 获取阻塞像素（1-based） ==========
static std::vector<std::pair<int, int>>
GetBlockedAreaPixels_1Based(const Grid& grid, int layer, int blockedAreaID)
{
    if (layer < 0 || layer >= grid.Layers)
        throw std::runtime_error("层号超出范围");

    LayerBlockedAnalysis la = BlockedAreaAnalyzer::analyzeLayerBlockedAreas(grid, layer);

    const BlockedAreaInfo* target = nullptr;
    for (const auto& area : la.blockedAreas) {
        if (area.areaID == blockedAreaID) { target = &area; break; }
    }
    if (!target) return {};

    std::vector<std::pair<int, int>> out;
    out.reserve(target->cells.size());
    for (const auto& rc0 : target->cells) {
        out.emplace_back(rc0.first + 1, rc0.second + 1); // 0-based -> 1-based
    }
    return out;
}

// ========== 包围盒（1-based） ==========
struct BBoxPx1 {
    int row_min, col_min; // top-left (含)
    int row_max, col_max; // bottom-right (含)
};

static inline BBoxPx1 ComputeBBoxFromPixels1(const std::vector<std::pair<int, int>>& pxs1)
{
    if (pxs1.empty()) throw std::runtime_error("ComputeBBoxFromPixels1: 输入像素为空");
    int rmin = INT_MAX, cmin = INT_MAX, rmax = INT_MIN, cmax = INT_MIN;
    for (auto [r1, c1] : pxs1) {
        rmin = std::min(rmin, r1);
        cmin = std::min(cmin, c1);
        rmax = std::max(rmax, r1);
        cmax = std::max(cmax, c1);
    }
    return { rmin, cmin, rmax, cmax };
}

// 像素四角（LB, RB, RT, LT）
static inline std::array<std::pair<int, int>, 4> BBoxCorners_LB_RB_RT_LT(const BBoxPx1& b)
{
    return {
        std::pair<int,int>{ b.row_max, b.col_min }, // LB
        std::pair<int,int>{ b.row_max, b.col_max }, // RB
        std::pair<int,int>{ b.row_min, b.col_max }, // RT
        std::pair<int,int>{ b.row_min, b.col_min }  // LT
    };
}

// 现实四角（严格 LT, LB, RT, RB，用于切板/绘制边框）
static inline std::array<std::pair<double, double>, 4>
RealCorners_LT_LB_RT_RB(const Grid& grid, const BBoxPx1& b)
{
    const RealRectMM rr = PixelBoxEdgesPx1ToRealRect(grid, b.row_min, b.col_min, b.row_max, b.col_max);
    return {
        std::pair<double,double>{ rr.x_left , rr.y_top }, // LT
        std::pair<double,double>{ rr.x_left , rr.y_bot }, // LB
        std::pair<double,double>{ rr.x_right, rr.y_top }, // RT
        std::pair<double,double>{ rr.x_right, rr.y_bot }  // RB
    };
}

// ========== m 像素扩张（带“边界补偿”与剪裁报告） ==========
struct ClipReport {
    bool left = false, right = false, top = false, bottom = false;
    int  left_deficit_px = 0, right_deficit_px = 0, top_deficit_px = 0, bottom_deficit_px = 0; // 该侧缺了多少像素
};

// 逻辑：希望左右各扩 m、上下各扩 m；若某侧不够扩，把缺的补到对侧。
static inline BBoxPx1 ExpandBBoxByMPixelsCompensated(const BBoxPx1& b, int m,
    const Grid& g, ClipReport* rep = nullptr)
{
    const int W = g.width, H = g.height;

    // 该区域到边界的“可扩”空间
    const int avail_left = b.col_min - 1;
    const int avail_right = W - b.col_max;
    const int avail_top = b.row_min - 1;
    const int avail_bottom = H - b.row_max;

    const int left_def = std::max(0, m - avail_left);
    const int right_def = std::max(0, m - avail_right);
    const int top_def = std::max(0, m - avail_top);
    const int bot_def = std::max(0, m - avail_bottom);

    if (rep) {
        rep->left = (left_def > 0);
        rep->right = (right_def > 0);
        rep->top = (top_def > 0);
        rep->bottom = (bot_def > 0);
        rep->left_deficit_px = left_def;
        rep->right_deficit_px = right_def;
        rep->top_deficit_px = top_def;
        rep->bottom_deficit_px = bot_def;
    }

    // 对侧补偿：希望总扩张量 ≈ 2m
    const int extend_left = m + right_def; // 右侧不够，把缺的补到左侧
    const int extend_right = m + left_def;  // 左侧不够，把缺的补到右侧
    const int extend_top = m + bot_def;   // 下侧不够，补到上侧
    const int extend_bottom = m + top_def;   // 上侧不够，补到下侧

    BBoxPx1 e;
    e.col_min = std::max(1, b.col_min - extend_left);
    e.col_max = std::min(W, b.col_max + extend_right);
    e.row_min = std::max(1, b.row_min - extend_top);
    e.row_max = std::min(H, b.row_max + extend_bottom);

    // 防翻转（极端窄板+超大 m）
    if (e.col_min > e.col_max) { e.col_min = 1; e.col_max = W; }
    if (e.row_min > e.row_max) { e.row_min = 1; e.row_max = H; }

    return e;
}


// 根据层编号获取层名称
std::string getLayerNameById(const Grid& grid, int layerId) {
    if (layerId >= 0 && layerId < grid.layerName.size()) {
        return grid.layerName[layerId]; // 返回对应层编号的层名称
    }
    else {
        return "无效的层编号"; // 如果层编号超出范围，返回错误信息
    }
}


// 根据层名称获取层编号
int getLayerIdByName(const Grid& grid, const std::string& layerName) {
    for (int i = 0; i < grid.layerName.size(); i++) {
        if (grid.layerName[i] == layerName) {
            return i; // 返回对应的层编号
        }
    }
    return -1; // 如果找不到，返回-1
}




// ========== 获取扩展后小PCB的四个现实角点坐标 ==========
static std::vector<std::pair<double, double>>
GetExpandedPCBRealCorners(const Grid& grid, int layer, int blockedAreaID, int expandPixels = 5)
{
    // 1. 获取阻塞区域像素
    auto blockedPixels = GetBlockedAreaPixels_1Based(grid, layer, blockedAreaID);
    if (blockedPixels.empty()) {
        throw std::runtime_error("未找到指定的阻塞区域");
    }

    // 2. 计算原始包围盒
    BBoxPx1 originalBBox = ComputeBBoxFromPixels1(blockedPixels);

    // 3. 扩展包围盒
    ClipReport clipReport;
    BBoxPx1 expandedBBox = ExpandBBoxByMPixelsCompensated(originalBBox, expandPixels, grid, &clipReport);

    // 4. 获取现实坐标的四个角点（顺时针顺序：LT, RT, RB, LB）
    auto realCorners = RealCorners_LT_LB_RT_RB(grid, expandedBBox);

    // 5. 重新排列为顺时针顺序：LT -> RT -> RB -> LB
    std::vector<std::pair<double, double>> result;
    result.reserve(4);

    // LT (左上)
    result.push_back(realCorners[0]);
    // RT (右上)  
    result.push_back(realCorners[2]);
    // RB (右下)
    result.push_back(realCorners[3]);
    // LB (左下)
    result.push_back(realCorners[1]);

    return result;
}

// 使用示例：
void exampleUsage(const Grid& grid) {
    try {
        int layer = 0;          // 层号
        int blockedAreaID = 1;  // 阻塞区域ID
        int expandPixels = 5;   // 扩展像素数

        auto corners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);

        // 打印四个角点坐标
        std::cout << "小PCB的四个角点坐标（顺时针顺序）：" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        for (size_t i = 0; i < corners.size(); ++i) {
            std::cout << "角点 " << i + 1 << ": ("
                << corners[i].first << " mm, "
                << corners[i].second << " mm)" << std::endl;
        }

        // 如果需要直接使用这个vector
        // corners 现在包含了四个点的坐标，按顺时针顺序：LT->RT->RB->LB

    }
    catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }
}

// ========== 计算线段与矩形边界的交点 ==========
static std::vector<std::pair<double, double>>
ComputeLineRectIntersections(double x1, double y1, double x2, double y2, const RealRectMM& rect)
{
    std::vector<std::pair<double, double>> intersections;

    // 检查与四条边的交点
    // 左边界 (x = rect.x_left)
    if (x1 != x2) {
        double t = (rect.x_left - x1) / (x2 - x1);
        if (t >= 0 && t <= 1) {
            double y = y1 + t * (y2 - y1);
            if (y >= rect.y_top && y <= rect.y_bot) {
                intersections.emplace_back(rect.x_left, y);
            }
        }
    }

    // 右边界 (x = rect.x_right)
    if (x1 != x2) {
        double t = (rect.x_right - x1) / (x2 - x1);
        if (t >= 0 && t <= 1) {
            double y = y1 + t * (y2 - y1);
            if (y >= rect.y_top && y <= rect.y_bot) {
                intersections.emplace_back(rect.x_right, y);
            }
        }
    }

    // 上边界 (y = rect.y_top)
    if (y1 != y2) {
        double t = (rect.y_top - y1) / (y2 - y1);
        if (t >= 0 && t <= 1) {
            double x = x1 + t * (x2 - x1);
            if (x >= rect.x_left && x <= rect.x_right) {
                intersections.emplace_back(x, rect.y_top);
            }
        }
    }

    // 下边界 (y = rect.y_bot)
    if (y1 != y2) {
        double t = (rect.y_bot - y1) / (y2 - y1);
        if (t >= 0 && t <= 1) {
            double x = x1 + t * (x2 - x1);
            if (x >= rect.x_left && x <= rect.x_right) {
                intersections.emplace_back(x, rect.y_bot);
            }
        }
    }

    return intersections;
}

// ========== 判断点是否在矩形内部 ==========
static bool IsPointInsideRect(double x, double y, const RealRectMM& rect)
{
    return (x >= rect.x_left && x <= rect.x_right &&
        y >= rect.y_top && y <= rect.y_bot);
}

// ========== 修正后的获取线段信息函数 ==========
static bool GetSegmentInfo(const std::shared_ptr<Node>& segment,
    double& x1, double& y1, double& x2, double& y2, std::string& layer)
{
    x1 = y1 = x2 = y2 = 0.0;
    layer.clear();

    bool hasStart = false, hasEnd = false, hasLayer = false;

    for (const auto& child : segment->children) {
        if (child->name == "start") {
            if (child->parameters.size() >= 2) {
                x1 = std::stod(child->parameters[0]);
                y1 = std::stod(child->parameters[1]);
                hasStart = true;
            }
        }
        else if (child->name == "end") {
            if (child->parameters.size() >= 2) {
                x2 = std::stod(child->parameters[0]);
                y2 = std::stod(child->parameters[1]);
                hasEnd = true;
            }
        }
        else if (child->name == "layer") {
            if (!child->parameters.empty()) {
                layer = child->parameters[0];
                hasLayer = true;
            }
        }
    }

    return (hasStart && hasEnd && hasLayer);
}

// ========== 判断点是否为连接点 ==========
static bool IsConnectionPoint(double x, double y, const std::vector<std::shared_ptr<Node>>& segments,
    const std::string& targetLayer, double tolerance = 0.001)
{
    int connectionCount = 0;

    for (const auto& segment : segments) {
        double x1, y1, x2, y2;
        std::string layer;

        if (!GetSegmentInfo(segment, x1, y1, x2, y2, layer)) {
            continue;
        }

        // 检查是否在目标层
        if (layer != targetLayer) {
            continue;
        }

        // 检查点是否与线段起点匹配
        if (std::abs(x1 - x) < tolerance && std::abs(y1 - y) < tolerance) {
            connectionCount++;
        }

        // 检查点是否与线段终点匹配
        if (std::abs(x2 - x) < tolerance && std::abs(y2 - y) < tolerance) {
            connectionCount++;
        }

        // 如果连接数超过1，说明是连接点
        if (connectionCount > 1) {
            return true;
        }
    }

    return false;
}

// ========== 修正后的主要处理函数 ==========
static std::vector<std::pair<double, double>>
FindIntersectionsAndEndpointsInSmallPCB(const Grid& grid, int layerId, int blockedAreaId, int expandPixels = 5)
{
    std::vector<std::pair<double, double>> result;

    try {
        // 1. 获取小PCB的边界框
        auto corners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
        RealRectMM pcbBox;
        pcbBox.x_left = corners[0].first;
        pcbBox.y_top = corners[0].second;
        pcbBox.x_right = corners[1].first;
        pcbBox.y_bot = corners[2].second;

        std::cout << "小PCB边界框: 左=" << pcbBox.x_left << ", 上=" << pcbBox.y_top
            << ", 右=" << pcbBox.x_right << ", 下=" << pcbBox.y_bot << std::endl;

        // 2. 获取当前层的名称
        std::string targetLayerName = getLayerNameById(grid, layerId);
        std::cout << "目标层: " << targetLayerName << std::endl;

        // 3. 直接使用 grid.segments
        const auto& segments = grid.segments;
        std::cout << "找到 " << segments.size() << " 个线段" << std::endl;

        // 4. 遍历所有线段
        int processedCount = 0;
        for (const auto& segment : segments) {
            double x1, y1, x2, y2;
            std::string layer;

            if (!GetSegmentInfo(segment, x1, y1, x2, y2, layer)) {
                continue;
            }

            // 检查是否在目标层
            if (layer != targetLayerName) {
                continue;
            }

            processedCount++;
            std::cout << "处理线段 " << processedCount << ": (" << x1 << "," << y1 << ") -> (" << x2 << "," << y2
                << "), 层: " << layer << std::endl;

            // 5. 计算线段与小PCB边界的交点
            auto intersections = ComputeLineRectIntersections(x1, y1, x2, y2, pcbBox);

            // 6. 检查线段端点是否在小PCB内部
            bool startInside = IsPointInsideRect(x1, y1, pcbBox);
            bool endInside = IsPointInsideRect(x2, y2, pcbBox);

            std::cout << "  起点在内部: " << (startInside ? "是" : "否")
                << ", 终点在内部: " << (endInside ? "是" : "否")
                << ", 交点数量: " << intersections.size() << std::endl;

            // 7. 收集交点（不检查是否为连接点，因为交点在边界上）
            for (const auto& intersection : intersections) {
                // 交点通常在边界上，不需要检查是否为连接点
                result.push_back(intersection);
                std::cout << "  交点: (" << intersection.first << ", " << intersection.second << ")" << std::endl;
            }

            // 8. 检查终止点，排除连接点
            if (startInside && !endInside) {
                // 起点在内部，终点在外部 → 起点是终止点
                // 检查起点是否为连接点
                if (!IsConnectionPoint(x1, y1, segments, targetLayerName)) {
                    result.emplace_back(x1, y1);
                    std::cout << "  终止点(起点在内部): (" << x1 << ", " << y1 << ")" << std::endl;
                }
                else {
                    std::cout << "  排除连接点(起点): (" << x1 << ", " << y1 << ")" << std::endl;
                }
            }
            else if (!startInside && endInside) {
                // 起点在外部，终点在内部 → 终点是终止点
                // 检查终点是否为连接点
                if (!IsConnectionPoint(x2, y2, segments, targetLayerName)) {
                    result.emplace_back(x2, y2);
                    std::cout << "  终止点(终点在内部): (" << x2 << ", " << y2 << ")" << std::endl;
                }
                else {
                    std::cout << "  排除连接点(终点): (" << x2 << ", " << y2 << ")" << std::endl;
                }
            }
            else if (startInside && endInside) {
                // 线段完全在内部，两个端点都在小PCB内
                // 检查两个端点是否为连接点
                if (!IsConnectionPoint(x1, y1, segments, targetLayerName)) {
                    result.emplace_back(x1, y1);
                    std::cout << "  终止点(起点在内部): (" << x1 << ", " << y1 << ")" << std::endl;
                }
                else {
                    std::cout << "  排除连接点(起点): (" << x1 << ", " << y1 << ")" << std::endl;
                }

                if (!IsConnectionPoint(x2, y2, segments, targetLayerName)) {
                    result.emplace_back(x2, y2);
                    std::cout << "  终止点(终点在内部): (" << x2 << ", " << y2 << ")" << std::endl;
                }
                else {
                    std::cout << "  排除连接点(终点): (" << x2 << ", " << y2 << ")" << std::endl;
                }
            }
        }

        std::cout << "共处理了 " << processedCount << " 个在目标层上的线段" << std::endl;

        // 9. 去重（可选，确保结果中没有重复点）
        auto last = std::unique(result.begin(), result.end(),
            [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                return std::abs(a.first - b.first) < 0.001 && std::abs(a.second - b.second) < 0.001;
            });
        result.erase(last, result.end());

    }
    catch (const std::exception& e) {
        std::cerr << "在FindIntersectionsAndEndpointsInSmallPCB中出错: " << e.what() << std::endl;
    }

    return result;
}

// ========== 主函数 ==========
int main() {
    try {
        // 1. 初始化Grid，读取PCB文件
        Grid grid;
        std::string filename = "testcase.kicad_pcb";  // 请替换为实际文件路径
        std::cout << "正在读取PCB文件: " << filename << std::endl;
        grid.SetUp(filename);

        std::cout << "Grid初始化成功!" << std::endl;
        std::cout << "网格尺寸: " << grid.width << " x " << grid.height << std::endl;
        std::cout << "层数: " << grid.Layers << std::endl;
        std::cout << "输入比例: " << grid.inputScale << std::endl;
        std::cout << "坐标范围: x[" << grid.min_x << ", " << grid.max_x << "], y[" << grid.min_y << ", " << grid.max_y << "]" << std::endl;

        // 2. 显示可用的层信息
        std::cout << "\n=== 可用层信息 ===" << std::endl;
        for (int i = 0; i < grid.layerName.size(); ++i) {
            std::cout << "层 " << i << ": " << grid.layerName[i] << std::endl;
        }

        // 3. 分析每层的阻塞区域
        std::cout << "\n=== 各层阻塞区域分析 ===" << std::endl;
        for (int layer = 0; layer < grid.Layers; ++layer) {
            try {
                LayerBlockedAnalysis la = BlockedAreaAnalyzer::analyzeLayerBlockedAreas(grid, layer);
                std::cout << "层 " << layer << " (" << getLayerNameById(grid, layer)
                    << "): 找到 " << la.blockedAreas.size() << " 个阻塞区域" << std::endl;

                for (const auto& area : la.blockedAreas) {
                    std::cout << "  阻塞区域ID: " << area.areaID
                        << ", 像素数: " << area.cells.size() << std::endl;
                }
            }
            catch (const std::exception& e) {
                std::cout << "层 " << layer << " 分析失败: " << e.what() << std::endl;
            }
        }

        // 4. 测试获取扩展后小PCB的角点坐标
        std::cout << "\n=== 测试获取扩展后小PCB角点坐标 ===" << std::endl;

        // 选择测试参数（可以根据上面的分析结果调整）
        int testLayer = 1;          // 测试层号
        int testBlockedAreaID = 1;  // 测试阻塞区域ID
        int expandPixels = 4;       // 扩展像素数

        std::cout << "测试参数: 层=" << testLayer
            << "(" << getLayerNameById(grid, testLayer) << ")"
            << ", 阻塞区域ID=" << testBlockedAreaID
            << ", 扩展像素=" << expandPixels << std::endl;

        // 获取阻塞区域像素
        auto blockedPixels = GetBlockedAreaPixels_1Based(grid, testLayer, testBlockedAreaID);
        if (blockedPixels.empty()) {
            std::cout << "警告: 未找到阻塞区域 " << testBlockedAreaID << "，将尝试其他区域..." << std::endl;

            // 尝试第一个可用的阻塞区域
            LayerBlockedAnalysis la = BlockedAreaAnalyzer::analyzeLayerBlockedAreas(grid, testLayer);
            if (!la.blockedAreas.empty()) {
                testBlockedAreaID = la.blockedAreas[0].areaID;
                blockedPixels = GetBlockedAreaPixels_1Based(grid, testLayer, testBlockedAreaID);
                std::cout << "改为使用阻塞区域ID: " << testBlockedAreaID << std::endl;
            }
        }

        if (!blockedPixels.empty()) {
            std::cout << "找到阻塞区域像素数: " << blockedPixels.size() << std::endl;

            // 计算原始包围盒
            BBoxPx1 originalBBox = ComputeBBoxFromPixels1(blockedPixels);
            std::cout << "原始包围盒: 行[" << originalBBox.row_min << "-" << originalBBox.row_max
                << "], 列[" << originalBBox.col_min << "-" << originalBBox.col_max << "]" << std::endl;

            // 获取扩展后的小PCB角点坐标
            auto corners = GetExpandedPCBRealCorners(grid, testLayer, testBlockedAreaID, expandPixels);

            // 输出结果
            std::cout << "\n=== 扩展后小PCB的四个角点坐标 ===" << std::endl;
            std::cout << std::fixed << std::setprecision(6);
            std::vector<std::string> cornerNames = { "左上(LT)", "右上(RT)", "右下(RB)", "左下(LB)" };

            for (size_t i = 0; i < corners.size(); ++i) {
                std::cout << cornerNames[i] << ": ("
                    << corners[i].first << " mm, "
                    << corners[i].second << " mm)" << std::endl;
            }

            // 5. 将角点坐标存储到vector中
            std::vector<std::pair<double, double>> pcbCorners = corners;
            std::cout << "\n角点坐标已存储到vector中，大小: " << pcbCorners.size() << std::endl;

            // 6. 创建对应的RealRectMM
            RealRectMM pcbBox;
            pcbBox.x_left = corners[0].first;   // 左上x
            pcbBox.y_top = corners[0].second;   // 左上y  
            pcbBox.x_right = corners[1].first;  // 右上x
            pcbBox.y_bot = corners[2].second;   // 右下y

            std::cout << "\n=== 小PCB边界框 ===" << std::endl;
            std::cout << "左边界 x_left: " << pcbBox.x_left << " mm" << std::endl;
            std::cout << "上边界 y_top: " << pcbBox.y_top << " mm" << std::endl;
            std::cout << "右边界 x_right: " << pcbBox.x_right << " mm" << std::endl;
            std::cout << "下边界 y_bot: " << pcbBox.y_bot << " mm" << std::endl;

            // 7. 将边界框存储到vector中
            std::vector<double> pcbBounds;
            pcbBounds.reserve(4);
            pcbBounds.push_back(pcbBox.x_left);
            pcbBounds.push_back(pcbBox.y_top);
            pcbBounds.push_back(pcbBox.x_right);
            pcbBounds.push_back(pcbBox.y_bot);

            std::cout << "\n边界框已存储到vector中，大小: " << pcbBounds.size() << std::endl;
            std::vector<std::string> boundNames = { "左边界", "上边界", "右边界", "下边界" };
            for (size_t i = 0; i < pcbBounds.size(); ++i) {
                std::cout << boundNames[i] << ": " << pcbBounds[i] << " mm" << std::endl;
            }

            // 8. 查找小PCB内的所有交点和终止点
            std::cout << "\n=== 开始查找小PCB内的交点和终止点 ===" << std::endl;
            auto intersectionPoints = FindIntersectionsAndEndpointsInSmallPCB(grid, testLayer, testBlockedAreaID, expandPixels);

            std::cout << "\n=== 最终结果汇总 ===" << std::endl;
            std::cout << "共找到 " << intersectionPoints.size() << " 个点:" << std::endl;
            for (size_t i = 0; i < intersectionPoints.size(); ++i) {
                std::cout << "点 " << i + 1 << ": (" << intersectionPoints[i].first
                    << ", " << intersectionPoints[i].second << ")" << std::endl;
            }

            // 9. 可选：将交点坐标保存到文件或其他处理
            std::cout << "\n所有处理完成!" << std::endl;

        }
        else {
            std::cout << "错误: 在层 " << testLayer << " 上未找到任何阻塞区域" << std::endl;
        }

    }
    catch (const std::exception& e) {
        std::cerr << "程序运行出错: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\n程序执行完成!" << std::endl;
    return 0;
}
