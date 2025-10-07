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

// ========== 修改后的点信息结构 ==========
struct PointWithNet {
    double x;
    double y;
    std::string net;

    PointWithNet(double x_val, double y_val, const std::string& net_val)
        : x(x_val), y(y_val), net(net_val) {
    }
};

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

// ========== 修改后的获取线段信息函数 ==========
static bool GetSegmentInfo(const std::shared_ptr<Node>& segment,
    double& x1, double& y1, double& x2, double& y2,
    std::string& layer, std::string& net)
{
    x1 = y1 = x2 = y2 = 0.0;
    layer.clear();
    net.clear();

    bool hasStart = false, hasEnd = false, hasLayer = false, hasNet = false;

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
        else if (child->name == "net") {
            if (!child->parameters.empty()) {
                net = child->parameters[0];
                hasNet = true;
            }
        }
    }

    return (hasStart && hasEnd && hasLayer && hasNet);
}

// ========== 修改后的判断连接点函数 ==========
static bool IsConnectionPoint(double x, double y, const std::vector<std::shared_ptr<Node>>& segments,
    const std::string& targetLayer, double tolerance = 0.001)
{
    int connectionCount = 0;

    for (const auto& segment : segments) {
        double x1, y1, x2, y2;
        std::string layer, net;

        if (!GetSegmentInfo(segment, x1, y1, x2, y2, layer, net)) {
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





static std::vector<PointWithNet>
FindIntersectionsAndEndpointsInSmallPCB(const Grid& grid, int layerId, int blockedAreaId, int expandPixels = 5)
{
    std::vector<PointWithNet> result;

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
            std::string layer, net;

            if (!GetSegmentInfo(segment, x1, y1, x2, y2, layer, net)) {
                continue;
            }

            // 检查是否在目标层
            if (layer != targetLayerName) {
                continue;
            }

            processedCount++;
            std::cout << "处理线段 " << processedCount << ": (" << x1 << "," << y1 << ") -> (" << x2 << "," << y2
                << "), 层: " << layer << ", 网络: " << net << std::endl;

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
                result.emplace_back(intersection.first, intersection.second, net);
                std::cout << "  交点: (" << intersection.first << ", " << intersection.second
                    << "), 网络: " << net << std::endl;
            }

            // 8. 检查终止点，排除连接点
            if (startInside && !endInside) {
                // 起点在内部，终点在外部 → 起点是终止点
                // 检查起点是否为连接点
                if (!IsConnectionPoint(x1, y1, segments, targetLayerName)) {
                    result.emplace_back(x1, y1, net);
                    std::cout << "  终止点(起点在内部): (" << x1 << ", " << y1
                        << "), 网络: " << net << std::endl;
                }
                else {
                    std::cout << "  排除连接点(起点): (" << x1 << ", " << y1
                        << "), 网络: " << net << std::endl;
                }
            }
            else if (!startInside && endInside) {
                // 起点在外部，终点在内部 → 终点是终止点
                // 检查终点是否为连接点
                if (!IsConnectionPoint(x2, y2, segments, targetLayerName)) {
                    result.emplace_back(x2, y2, net);
                    std::cout << "  终止点(终点在内部): (" << x2 << ", " << y2
                        << "), 网络: " << net << std::endl;
                }
                else {
                    std::cout << "  排除连接点(终点): (" << x2 << ", " << y2
                        << "), 网络: " << net << std::endl;
                }
            }
            else if (startInside && endInside) {
                // 线段完全在内部，两个端点都在小PCB内
                // 检查两个端点是否为连接点
                if (!IsConnectionPoint(x1, y1, segments, targetLayerName)) {
                    result.emplace_back(x1, y1, net);
                    std::cout << "  终止点(起点在内部): (" << x1 << ", " << y1
                        << "), 网络: " << net << std::endl;
                }
                else {
                    std::cout << "  排除连接点(起点): (" << x1 << ", " << y1
                        << "), 网络: " << net << std::endl;
                }

                if (!IsConnectionPoint(x2, y2, segments, targetLayerName)) {
                    result.emplace_back(x2, y2, net);
                    std::cout << "  终止点(终点在内部): (" << x2 << ", " << y2
                        << "), 网络: " << net << std::endl;
                }
                else {
                    std::cout << "  排除连接点(终点): (" << x2 << ", " << y2
                        << "), 网络: " << net << std::endl;
                }
            }
        }

        std::cout << "共处理了 " << processedCount << " 个在目标层上的线段" << std::endl;

        // 9. 去重（确保结果中没有重复点）
        auto last = std::unique(result.begin(), result.end(),
            [](const PointWithNet& a, const PointWithNet& b) {
                return std::abs(a.x - b.x) < 0.001 &&
                    std::abs(a.y - b.y) < 0.001 &&
                    a.net == b.net;
            });
        result.erase(last, result.end());

    }
    catch (const std::exception& e) {
        std::cerr << "在FindIntersectionsAndEndpointsInSmallPCB中出错: " << e.what() << std::endl;
    }

    return result;
}


// ========== 扩张四个顶点 ==========
static std::vector<std::pair<double, double>>
ExpandCornersByOffset(const std::vector<std::pair<double, double>>& originalCorners, double offsetMM)
{
    if (originalCorners.size() != 4) {
        throw std::runtime_error("扩张顶点函数需要4个输入点");
    }

    std::vector<std::pair<double, double>> expandedCorners;
    expandedCorners.reserve(4);

    // 计算矩形的中心点
    double center_x = 0.0, center_y = 0.0;
    for (const auto& corner : originalCorners) {
        center_x += corner.first;
        center_y += corner.second;
    }
    center_x /= 4.0;
    center_y /= 4.0;

    // 对每个顶点进行扩张
    for (const auto& corner : originalCorners) {
        // 计算从中心指向顶点的方向向量
        double dx = corner.first - center_x;
        double dy = corner.second - center_y;

        // 计算向量的长度
        double length = std::sqrt(dx * dx + dy * dy);

        if (length > 0) {
            // 单位化方向向量
            dx /= length;
            dy /= length;

            // 沿着方向向量移动offset距离
            double new_x = corner.first + dx * offsetMM;
            double new_y = corner.second + dy * offsetMM;

            expandedCorners.emplace_back(new_x, new_y);
        }
        else {
            // 如果长度为零，直接使用原坐标
            expandedCorners.push_back(corner);
        }
    }

    return expandedCorners;
}

// ========== 修正精度的简约主函数 ==========
int main() {
    try {
        // 1. 初始化Grid
        Grid grid;
        std::string filename = "testcase.kicad_pcb";
        grid.SetUp(filename);
        std::cout << "Grid初始化成功! 网格尺寸: " << grid.width << " x " << grid.height << std::endl;

        // 2. 设置参数
        int layerId = 1;           // 层序号
        int blockedAreaId = 2;     // 阻塞区域编号
        int expandPixels = 3;      // 扩展像素数
        double cornerOffset = 3.3; // 顶点扩张距离(mm)

        std::cout << "参数: 层=" << layerId << "(" << getLayerNameById(grid, layerId)
            << "), 阻塞区域=" << blockedAreaId << ", 扩展=" << expandPixels << "像素"
            << ", 顶点扩张=" << cornerOffset << "mm" << std::endl;

        // 3. 获取小PCB边界框
        auto corners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);

        // 设置高精度输出
        std::cout << std::fixed << std::setprecision(5);
        std::cout << "\n=== 原始小PCB边界框 ===" << std::endl;
        std::vector<std::string> cornerNames = { "左上(LT)", "右上(RT)", "右下(RB)", "左下(LB)" };
        for (size_t i = 0; i < corners.size(); ++i) {
            std::cout << cornerNames[i] << ": (" << corners[i].first << ", " << corners[i].second << ")" << std::endl;
        }

        // 4. 扩张顶点
        auto expandedCorners = ExpandCornersByOffset(corners, cornerOffset);

        std::cout << "\n=== 扩张后的小PCB边界框 ===" << std::endl;
        for (size_t i = 0; i < expandedCorners.size(); ++i) {
            std::cout << cornerNames[i] << ": (" << expandedCorners[i].first << ", " << expandedCorners[i].second << ")" << std::endl;
        }

        // 5. 将扩张后的顶点存储到vector中
        std::vector<std::pair<double, double>> finalExpandedCorners = expandedCorners;
        std::cout << "\n扩张后的顶点已存储到finalExpandedCorners vector中，大小: " << finalExpandedCorners.size() << std::endl;

        // 6. 查找交点和终止点（使用原始边界框）
        auto points = FindIntersectionsAndEndpointsInSmallPCB(grid, layerId, blockedAreaId, expandPixels);

        // 7. 输出最终结果
        std::cout << "\n=== 线段交点和终止点 ===" << std::endl;
        std::cout << "共找到 " << points.size() << " 个点:" << std::endl;
        for (size_t i = 0; i < points.size(); ++i) {
            std::cout << "点 " << i + 1 << ": (" << points[i].x << ", " << points[i].y
                << "), 网络: " << points[i].net << std::endl;
        }

    }
    catch (const std::exception& e) {
        std::cerr << "程序运行出错: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "程序执行完成!" << std::endl;
    return 0;
}
