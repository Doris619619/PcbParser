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
#include "separate3.h"
#include <fstream>
#include <sstream>

#include <map>

#include "First_Part12.h"   // KiCadParser 声明
#include <algorithm>
#include <cmath>
#include <stdexcept>


#include "blocked_area_analyzer.h"   // 内含 grid.h

// ========== 像素中心 <-> 现实坐标 ==========
std::pair<double, double>
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

std::pair<int, int>
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

RealRectMM PixelBoxEdgesPx1ToRealRect(const Grid& grid,
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
std::vector<std::pair<int, int>>
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


BBoxPx1 ComputeBBoxFromPixels1(const std::vector<std::pair<int, int>>& pxs1)
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
std::array<std::pair<int, int>, 4> BBoxCorners_LB_RB_RT_LT(const BBoxPx1& b)
{
    return {
        std::pair<int,int>{ b.row_max, b.col_min }, // LB
        std::pair<int,int>{ b.row_max, b.col_max }, // RB
        std::pair<int,int>{ b.row_min, b.col_max }, // RT
        std::pair<int,int>{ b.row_min, b.col_min }  // LT
    };
}

// 现实四角（严格 LT, LB, RT, RB，用于切板/绘制边框）
std::array<std::pair<double, double>, 4>
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



// 逻辑：希望左右各扩 m、上下各扩 m；若某侧不够扩，把缺的补到对侧。
 BBoxPx1 ExpandBBoxByMPixelsCompensated(const BBoxPx1& b, int m,
    const Grid& g, ClipReport* rep)
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






static inline RealRectMM ExpandRect_HalfPlusK(const RealRectMM& in, const Grid& grid, int kWholePixels) {
    if (grid.inputScale <= 0) throw std::runtime_error("inputScale 必须为正");
    const double px_mm = 1.0 / static_cast<double>(grid.inputScale); // 1px = 1/inputScale mm
    const double d = (0.5 + static_cast<double>(kWholePixels)) * px_mm;
    RealRectMM r = in;
    r.x_left -= d;  r.x_right += d;
    r.y_top -= d;  r.y_bot += d;
    return r;
}

static inline RealRectMM ClipRectToBoard(const RealRectMM& in,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY) {
    RealRectMM r = in;
    if (r.x_left < boardMinX) r.x_left = boardMinX;
    if (r.x_right > boardMaxX) r.x_right = boardMaxX;
    if (r.y_top < boardMinY) r.y_top = boardMinY;
    if (r.y_bot > boardMaxY) r.y_bot = boardMaxY;
    // 极端防翻转
    if (r.x_left > r.x_right) { r.x_left = boardMinX; r.x_right = boardMaxX; }
    if (r.y_top > r.y_bot) { r.y_top = boardMinY; r.y_bot = boardMaxY; }
    return r;
}

// 方案A主函数（建议用这个）：像素域扩张 -> 现实域 0.5px+kpx -> 裁剪板框 -> 返回 LT,RT,RB,LB
std::vector<std::pair<double, double>>
GetExpandedPCBCornersSnapAndClip(const Grid& grid, int layer, int blockedAreaID,
    int expandPixels, int kWholePixels,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY)
{
    // 1) 先得到“像素域扩张后的现实四角”（LT,RT,RB,LB）
    auto baseCorners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);
    if (baseCorners.size() != 4) {
        throw std::runtime_error("GetExpandedPCBCornersSnapAndClip: 角点数 != 4");
    }

    // 2) 现实域对称外扩 0.5px + k×1px
    RealRectMM baseRect = CornersToRect_LT_RT_RB_LB(baseCorners);
    RealRectMM grown = ExpandRect_HalfPlusK(baseRect, grid, kWholePixels);

    // 3) 裁剪到真实板框
    RealRectMM clipped = ClipRectToBoard(grown, boardMinX, boardMinY, boardMaxX, boardMaxY);

    // 4) 回四角（LT,RT,RB,LB）
    return RectToCorners_LT_RT_RB_LB(clipped);
}


// ========== 获取扩展后小PCB的四个现实角点坐标 ==========
std::vector<std::pair<double, double>>
GetExpandedPCBRealCorners(const Grid& grid, int layer, int blockedAreaID, int expandPixels)
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
std::vector<std::pair<double, double>>
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
bool IsPointInsideRect(double x, double y, const RealRectMM& rect)
{
    return (x >= rect.x_left && x <= rect.x_right &&
        y >= rect.y_top && y <= rect.y_bot);
}

// ========== 修改后的获取线段信息函数 ==========
bool GetSegmentInfo(const std::shared_ptr<Node>& segment,
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
bool IsConnectionPoint(double x, double y, const std::vector<std::shared_ptr<Node>>& segments,
    const std::string& targetLayer, double tolerance)
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





std::vector<PointWithNet>
FindIntersectionsAndEndpointsInSmallPCB(const Grid& grid, int layerId, int blockedAreaId, int expandPixels)
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


// === 像素大小（mm/px） ===
static inline double PixelSizeMM(const Grid& g) {
    if (g.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    return 1.0 / static_cast<double>(g.inputScale);
}

// === 用四角（LT,RT,RB,LB）组装矩形 ===
static inline RealRectMM CornersToRect_LT_RT_RB_LB(const std::vector<std::pair<double, double>>& c) {
    if (c.size() != 4) throw std::runtime_error("CornersToRect_LT_RT_RB_LB: 需要4个角点");
    RealRectMM r;
    // 约定：c[0]=LT, c[1]=RT, c[2]=RB, c[3]=LB
    r.x_left = c[0].first;
    r.y_top = c[0].second;
    r.x_right = c[1].first;
    r.y_bot = c[2].second;
    return r;
}

// === 把矩形转回四角（LT,RT,RB,LB） ===
static inline std::vector<std::pair<double, double>> RectToCorners_LT_RT_RB_LB(const RealRectMM& r) {
    return {
        { r.x_left , r.y_top },   // LT
        { r.x_right, r.y_top },   // RT
        { r.x_right, r.y_bot },   // RB
        { r.x_left , r.y_bot }    // LB
    };
}






// —— 辅助：把“在边界上的现实坐标”稳定地映射到小板内部的像素 ——
// 注意：右边界/下边界要往里收一个极小量，否则会落到外侧像素。
static inline std::pair<int, int> RealToPixel_Inside_LT(const Grid& g, double x, double y) {
    return RealMMToPixelPx1(g, x, y); // 左/上边界直接 floor 落在内侧
}
static inline std::pair<int, int> RealToPixel_Inside_RT(const Grid& g, double x, double y) {
    const double eps = 1.0 / (g.inputScale * 1024.0); // < 1px 的很小量
    return RealMMToPixelPx1(g, x - eps, y);
}
static inline std::pair<int, int> RealToPixel_Inside_RB(const Grid& g, double x, double y) {
    const double eps = 1.0 / (g.inputScale * 1024.0);
    return RealMMToPixelPx1(g, x - eps, y - eps);
}
static inline std::pair<int, int> RealToPixel_Inside_LB(const Grid& g, double x, double y) {
    const double eps = 1.0 / (g.inputScale * 1024.0);
    return RealMMToPixelPx1(g, x, y - eps);
}

// —— 对外1：四个角的“锚点像素”（LT,RT,RB,LB） ——
// 先拿现实角点(LT,RT,RB,LB)，再做“内侧偏置”的像素映射。
std::array<std::pair<int, int>, 4>
GetExpandedPCBCornerAnchorPixels(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels)
{
    auto corners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);
    if (corners.size() != 4) throw std::runtime_error("角点数量应为4");

    // corners[0]=LT, [1]=RT, [2]=RB, [3]=LB （当前实现返回顺序即如此）
    auto lt = RealToPixel_Inside_LT(grid, corners[0].first, corners[0].second);
    auto rt = RealToPixel_Inside_RT(grid, corners[1].first, corners[1].second);
    auto rb = RealToPixel_Inside_RB(grid, corners[2].first, corners[2].second);
    auto lb = RealToPixel_Inside_LB(grid, corners[3].first, corners[3].second);
    return { lt, rt, rb, lb };
}

// —— 辅助：根据“角的方向”生成 N×N 的紫色块，并做 1..H/1..W 裁剪 ——
// dr_sign: 行方向(+1向下, -1向上)；dc_sign: 列方向(+1向右, -1向左)
static inline void AppendCornerPatch(std::vector<std::pair<int, int>>& out,
    int r0, int c0, int N,
    int dr_sign, int dc_sign,
    const Grid& g)
{
    for (int di = 0; di < N; ++di) {
        int r = r0 + dr_sign * di;
        if (r < 1 || r > g.height) continue;
        for (int dj = 0; dj < N; ++dj) {
            int c = c0 + dc_sign * dj;
            if (c < 1 || c > g.width) continue;
            out.emplace_back(r, c);
        }
    }
}

// —— 对外2：四个角的紫色块像素，按 LT→RT→RB→LB 顺序拼接到一个 vector ——
// 例：patchN=3 就是每个角 3×3；若只要一个像素，传 patchN=1 即可。
std::vector<std::pair<int, int>>
GetCornerPixelsPx1(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels,
    int patchN)
{
    if (patchN < 1) patchN = 1;

    auto anchors = GetExpandedPCBCornerAnchorPixels(grid, layer, blockedAreaID, expandPixels);
    std::vector<std::pair<int, int>> out;
    out.reserve(4 * patchN * patchN);

    // LT：向“下、右”扩 N×N
    AppendCornerPatch(out, anchors[0].first, anchors[0].second, patchN, +1, +1, grid);
    // RT：向“下、左”扩
    AppendCornerPatch(out, anchors[1].first, anchors[1].second, patchN, +1, -1, grid);
    // RB：向“上、左”扩
    AppendCornerPatch(out, anchors[2].first, anchors[2].second, patchN, -1, -1, grid);
    // LB：向“上、右”扩
    AppendCornerPatch(out, anchors[3].first, anchors[3].second, patchN, -1, +1, grid);

    return out;
}

// ========== 生成完整PCB模板文件 ==========
void GenerateCompletePCBTemplate(
    const std::vector<PointWithNet>& modulePoints,
    const std::vector<std::pair<double, double>>& expandedCorners,
    const std::string& outputFilename)
{
    // 四个module的模板
    std::string template_four_modules = R"((kicad_pcb (version 20171130) (host pcbnew "(5.1.2)-2")

  (general
    (thickness 1.6)
  )

  (page A4)
  (layers
    (0 Top signal)
    (31 Bottom signal)
    (32 B.Adhes user)
    (33 F.Adhes user)
    (34 B.Paste user)
    (35 F.Paste user)
    (36 B.SilkS user)
    (37 F.SilkS user)
    (38 B.Mask user)
    (39 F.Mask user)
    (40 Dwgs.User user)
    (41 Cmts.User user)
    (42 Eco1.User user)
    (43 Eco2.User user)
    (44 Edge.Cuts user)
    (45 Margin user)
    (46 B.CrtYd user)
    (47 F.CrtYd user)
    (48 B.Fab user)
    (49 F.Fab user)
  )

  (setup
    (pad_to_mask_clearance 0.1016)
    (aux_axis_origin -120.62617 153.22756)
    (grid_origin -120.62617 153.22756)
    (pcbplotparams
      (layerselection 0x00010fc_ffffffff)
      (plot_on_all_layers_selection 0x0000000_00000000)
      (disableapertmacros false)
      (usegerberextensions false)
      (usegerberattributes true)
      (usegerberadvancedattributes true)
      (creategerberjobfile true)
      (dashed_line_dash_ratio 12)
      (dashed_line_gap_ratio 3)
      (svgprecision 4)
      (plotframeref false)
      (viasonmask false)
      (mode 1)
      (useauxorigin false)
      (hpglpennumber 1)
      (hpglpenspeed 20)
      (hpglpendiameter 15)
      (dxfpolygonmode true)
      (dxfimperialunits true)
      (dxfusepcbnewfont true)
      (psnegative false)
      (psa4output false)
      (plotreference true)
      (plotvalue true)
      (plotinvisibletext false)
      (sketchpadsonfab false)
      (subtractmaskfromsilk false)
      (outputformat 1)
      (mirror false)
      (drillshape 1)
      (scaleselection 1)
      (outputdirectory "")
    )
  )

  (net 0 "")

  (net 1 Y+)
  (net 2 Y-)
  (net 3 NetR3_2)
  (net 4 W-)
  (net 5 W+)
  (net 6 VIN)
  (net 7 VCC_5V)
  (net 8 TK6)
  (net 9 TK5)
  (net 10 TK4)
  (net 11 TK3)
  (net 12 TK2)
  (net 13 TK1)
  (net 14 SK4)
  (net 15 SK3)
  (net 16 SK2)
  (net 17 SK1)
  (net 18 S_LED7)
  (net 19 S_LED6)
  (net 20 S_LED5)
  (net 21 S_LED4)
  (net 22 S_LED3)
  (net 23 S_LED2)
  (net 24 S_LED1)
  (net 25 PWM2)
  (net 26 PWM1)
  (net 27 NetR19_2)
  (net 28 NetR16_2)
  (net 29 NetR15_2)
  (net 30 NetR14_1)
  (net 31 NetR12_1)
  (net 32 NetD2_1)
  (net 33 NetD1_1)
  (net 34 NetC5_2)
  (net 35 K_LED4)
  (net 36 K_LED3)
  (net 37 K_LED2)
  (net 38 K_LED1)
  (net 39 GND)
  (net 40 ALS_IN)

  (net_class Default "This is the default net class."
    (add_net ALS_IN)
    (add_net K_LED1)
    (add_net K_LED2)
    (add_net K_LED3)
    (add_net K_LED4)
    (add_net NetC5_2)
    (add_net NetD1_1)
    (add_net NetD2_1)
    (add_net NetR12_1)
    (add_net NetR14_1)
    (add_net NetR15_2)
    (add_net NetR16_2)
    (add_net NetR19_2)
    (add_net NetR3_2)
    (add_net PWM1)
    (add_net PWM2)
    (add_net SK1)
    (add_net SK2)
    (add_net SK3)
    (add_net SK4)
    (add_net S_LED1)
    (add_net S_LED2)
    (add_net S_LED3)
    (add_net S_LED4)
    (add_net S_LED5)
    (add_net S_LED6)
    (add_net S_LED7)
    (add_net TK1)
    (add_net TK2)
    (add_net TK3)
    (add_net TK4)
    (add_net TK5)
    (add_net TK6)
    (add_net VIN)
    (add_net W+)
    (add_net W-)
    (add_net Y+)
    (add_net Y-)
    (trace_width 0.254)
  )
  (net_class Power "This is the default net class."
    (add_net GND)
    (add_net VCC_5V)
    (trace_width 0.8)
  )

 
  
  (module "MIJI_ADPcbLib.PcbLib:R1206-JP_4" (layer Bottom) (tedit 0) (tstamp 0)
    (at 127.600901 96.36789899999977 270)
    (fp_text reference 3 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
    
    (pad 1 smd rect (at 0 0 0) (size 1.5 1) (layers Bottom B.Paste B.Mask)
(net 3 NetR3_2))
     )

  (module "MIJI_ADPcbLib.PcbLib:SMD-EC6.3X6.3X10_5" (layer Bottom) (tedit 0) (tstamp 0)
    (at 114.551001 88.57789899999983 270)
    (fp_text reference 4 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
       (pad 1 smd rect (at 0 0 0) (size 3.3 2.132) (layers Bottom B.Paste B.Mask)
(net 3 NetR3_2))
      )

  
  (module "Miscellaneous Devices LC.PcbLib:0805_C_7" (layer Bottom) (tedit 0) (tstamp 0)
    (at 109.04720099999992 105.11449900000001 270)
    (fp_text reference 6 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
      (pad 1 smd rect (at 0 0 0) (size 1.16 1.47) (layers Bottom B.Paste B.Mask)
(net 6 VIN))

  )

  (module "Miscellaneous Devices LC.PcbLib:1206_C_8" (layer Bottom) (tedit 0) (tstamp 0)
    (at 138.290501 96.60859899999969 180)
    (fp_text reference 7 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
       (pad 1 smd rect (at 0 0 0) (size 1.13 1.8) (layers Bottom B.Paste B.Mask)
(net 6 VIN))
     )

 
  (gr_line (start 160.5011 75.003598) (end 98.501099 75.003598) (layer Bottom) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 75.003598) (end 98.501099 120.003599) (layer Bottom) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 120.003599) (end 160.5011 120.003599) (layer Bottom) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 160.5011 120.003599) (end 160.5011 75.003598) (layer Bottom) (width 0.05) (tstamp 42A38F4))
))";

    // 两个module的模板
    std::string template_two_modules = R"((kicad_pcb (version 20171130) (host pcbnew "(5.1.2)-2")

  (general
    (thickness 1.6)
  )

  (page A4)
  (layers
    (0 Top signal)
    (31 Bottom signal)
    (32 B.Adhes user)
    (33 F.Adhes user)
    (34 B.Paste user)
    (35 F.Paste user)
    (36 B.SilkS user)
    (37 F.SilkS user)
    (38 B.Mask user)
    (39 F.Mask user)
    (40 Dwgs.User user)
    (41 Cmts.User user)
    (42 Eco1.User user)
    (43 Eco2.User user)
    (44 Edge.Cuts user)
    (45 Margin user)
    (46 B.CrtYd user)
    (47 F.CrtYd user)
    (48 B.Fab user)
    (49 F.Fab user)
  )

  (setup
    (pad_to_mask_clearance 0.1016)
    (aux_axis_origin -120.62617 153.22756)
    (grid_origin -120.62617 153.22756)
    (pcbplotparams
      (layerselection 0x00010fc_ffffffff)
      (plot_on_all_layers_selection 0x0000000_00000000)
      (disableapertmacros false)
      (usegerberextensions false)
      (usegerberattributes true)
      (usegerberadvancedattributes true)
      (creategerberjobfile true)
      (dashed_line_dash_ratio 12)
      (dashed_line_gap_ratio 3)
      (svgprecision 4)
      (plotframeref false)
      (viasonmask false)
      (mode 1)
      (useauxorigin false)
      (hpglpennumber 1)
      (hpglpenspeed 20)
      (hpglpendiameter 15)
      (dxfpolygonmode true)
      (dxfimperialunits true)
      (dxfusepcbnewfont true)
      (psnegative false)
      (psa4output false)
      (plotreference true)
      (plotvalue true)
      (plotinvisibletext false)
      (sketchpadsonfab false)
      (subtractmaskfromsilk false)
      (outputformat 1)
      (mirror false)
      (drillshape 1)
      (scaleselection 1)
      (outputdirectory "")
    )
  )

  (net 0 "")

  (net 1 Y+)
  (net 2 Y-)
  (net 3 NetR3_2)
  (net 4 W-)
  (net 5 W+)
  (net 6 VIN)
  (net 7 VCC_5V)
  (net 8 TK6)
  (net 9 TK5)
  (net 10 TK4)
  (net 11 TK3)
  (net 12 TK2)
  (net 13 TK1)
  (net 14 SK4)
  (net 15 SK3)
  (net 16 SK2)
  (net 17 SK1)
  (net 18 S_LED7)
  (net 19 S_LED6)
  (net 20 S_LED5)
  (net 21 S_LED4)
  (net 22 S_LED3)
  (net 23 S_LED2)
  (net 24 S_LED1)
  (net 25 PWM2)
  (net 26 PWM1)
  (net 27 NetR19_2)
  (net 28 NetR16_2)
  (net 29 NetR15_2)
  (net 30 NetR14_1)
  (net 31 NetR12_1)
  (net 32 NetD2_1)
  (net 33 NetD1_1)
  (net 34 NetC5_2)
  (net 35 K_LED4)
  (net 36 K_LED3)
  (net 37 K_LED2)
  (net 38 K_LED1)
  (net 39 GND)
  (net 40 ALS_IN)

  (net_class Default "This is the default net class."
    (add_net ALS_IN)
    (add_net K_LED1)
    (add_net K_LED2)
    (add_net K_LED3)
    (add_net K_LED4)
    (add_net NetC5_2)
    (add_net NetD1_1)
    (add_net NetD2_1)
    (add_net NetR12_1)
    (add_net NetR14_1)
    (add_net NetR15_2)
    (add_net NetR16_2)
    (add_net NetR19_2)
    (add_net NetR3_2)
    (add_net PWM1)
    (add_net PWM2)
    (add_net SK1)
    (add_net SK2)
    (add_net SK3)
    (add_net SK4)
    (add_net S_LED1)
    (add_net S_LED2)
    (add_net S_LED3)
    (add_net S_LED4)
    (add_net S_LED5)
    (add_net S_LED6)
    (add_net S_LED7)
    (add_net TK1)
    (add_net TK2)
    (add_net TK3)
    (add_net TK4)
    (add_net TK5)
    (add_net TK6)
    (add_net VIN)
    (add_net W+)
    (add_net W-)
    (add_net Y+)
    (add_net Y-)
    (trace_width 0.254)
  )
  (net_class Power "This is the default net class."
    (add_net GND)
    (add_net VCC_5V)
    (trace_width 0.8)
  )


 
  
  (module "MIJI_ADPcbLib.PcbLib:R1206-JP_4" (layer Bottom) (tedit 0) (tstamp 0)
    (at 127.600901 96.36789899999977 270)
    (fp_text reference 3 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
    
    (pad 1 smd rect (at 0 0 0) (size 1.5 1) (layers Bottom B.Paste B.Mask)
(net 3 NetR3_2))
     )

  (module "MIJI_ADPcbLib.PcbLib:SMD-EC6.3X6.3X10_5" (layer Bottom) (tedit 0) (tstamp 0)
    (at 114.551001 88.57789899999983 270)
    (fp_text reference 4 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
       (pad 1 smd rect (at 0 0 0) (size 3.3 2.132) (layers Bottom B.Paste B.Mask)
(net 3 NetR3_2))
      )

 
  (gr_line (start 160.5011 75.003598) (end 98.501099 75.003598) (layer Bottom) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 75.003598) (end 98.501099 120.003599) (layer Bottom) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 120.003599) (end 160.5011 120.003599) (layer Bottom) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 160.5011 120.003599) (end 160.5011 75.003598) (layer Bottom) (width 0.05) (tstamp 42A38F4))
))";

    // 处理module坐标替换
    std::vector<PointWithNet> replacementPoints;

    if (modulePoints.size() == 2) {
        // 情况1: 只有2个点，直接使用这两个点替换前两个module
        replacementPoints = modulePoints;
        std::cout << "使用2个点替换前两个module" << std::endl;
    }
    else if (modulePoints.size() >= 4) {
        // 情况2: 有4个或更多点，按网络分组处理
        // 找到有两个相同网络的点
        std::map<std::string, std::vector<PointWithNet>> pointsByNet;
        for (const auto& point : modulePoints) {
            pointsByNet[point.net].push_back(point);
        }

        bool foundMatchingNet = false;
        for (const auto& [net, netPoints] : pointsByNet) {
            if (netPoints.size() >= 2) {
                // 使用前两个相同网络的点
                replacementPoints.push_back(netPoints[0]);
                replacementPoints.push_back(netPoints[1]);
                foundMatchingNet = true;
                std::cout << "找到网络 '" << net << "' 的2个点用于前两个module" << std::endl;
                break;
            }
        }

        // 如果没有找到相同网络的点，使用前4个点
        if (!foundMatchingNet) {
            for (size_t i = 0; i < std::min(modulePoints.size(), size_t(4)); i++) {
                replacementPoints.push_back(modulePoints[i]);
            }
            std::cout << "没有找到相同网络的点，使用前4个点" << std::endl;
        }

        // 添加剩余的点
        size_t addedCount = replacementPoints.size();
        for (const auto& point : modulePoints) {
            if (addedCount >= 4) break;

            // 检查是否已经添加过
            bool alreadyAdded = false;
            for (const auto& addedPoint : replacementPoints) {
                if (point.x == addedPoint.x && point.y == addedPoint.y && point.net == addedPoint.net) {
                    alreadyAdded = true;
                    break;
                }
            }

            if (!alreadyAdded) {
                replacementPoints.push_back(point);
                addedCount++;
            }
        }
    }
    else {
        // 其他情况: 使用所有点
        replacementPoints = modulePoints;
        std::cout << "使用所有 " << modulePoints.size() << " 个点" << std::endl;
    }

    // 选择模板
    std::string template_str;
    if (replacementPoints.size() == 2) {
        template_str = template_two_modules;
        std::cout << "使用两个module的模板" << std::endl;
    }
    else {
        template_str = template_four_modules;
        std::cout << "使用四个module的模板" << std::endl;
    }

    // 处理模板替换
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);

    // 分割模板并按行处理
    std::istringstream templateStream(template_str);
    std::string line;
    int moduleCount = 0;
    int grLineCount = 0;

    while (std::getline(templateStream, line)) {
        // 处理module坐标替换
        if (line.find("(module") != std::string::npos) {
            // 输出module行
            oss << line << "\n";

            // 读取下一行(应该包含at坐标)
            if (std::getline(templateStream, line)) {
                if (line.find("(at") != std::string::npos) {
                    // 替换坐标
                    if (moduleCount < replacementPoints.size()) {
                        const auto& point = replacementPoints[moduleCount];
                        oss << "    (at " << point.x << " " << point.y << " 270)\n";
                        std::cout << "替换module " << (moduleCount + 1) << " 坐标为: ("
                            << point.x << ", " << point.y << "), 网络: " << point.net << std::endl;
                    }
                    else {
                        // 如果没有足够的点，保留原坐标
                        oss << line << "\n";
                    }
                    moduleCount++;
                }
                else {
                    oss << line << "\n";
                }
            }
        }
        // 处理gr_line坐标替换
        else if (line.find("(gr_line") != std::string::npos && expandedCorners.size() == 4) {
            if (grLineCount == 0) {
                // 第一条线: 上边 (LT -> RT)
                oss << "  (gr_line (start " << expandedCorners[0].first << " " << expandedCorners[0].second
                    << ") (end " << expandedCorners[1].first << " " << expandedCorners[1].second
                    << ") (layer Bottom) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 1) {
                // 第二条线: 右边 (RT -> RB)
                oss << "  (gr_line (start " << expandedCorners[1].first << " " << expandedCorners[1].second
                    << ") (end " << expandedCorners[2].first << " " << expandedCorners[2].second
                    << ") (layer Bottom) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 2) {
                // 第三条线: 下边 (RB -> LB)
                oss << "  (gr_line (start " << expandedCorners[2].first << " " << expandedCorners[2].second
                    << ") (end " << expandedCorners[3].first << " " << expandedCorners[3].second
                    << ") (layer Bottom) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 3) {
                // 第四条线: 左边 (LB -> LT)
                oss << "  (gr_line (start " << expandedCorners[3].first << " " << expandedCorners[3].second
                    << ") (end " << expandedCorners[0].first << " " << expandedCorners[0].second
                    << ") (layer Bottom) (width 0.05) (tstamp 42A38F4))\n";
            }
            grLineCount++;
        }
        else {
            oss << line << "\n";
        }


    }

    // 写入文件
    std::ofstream file(outputFilename);
    if (!file.is_open()) {
        throw std::runtime_error("无法创建输出文件: " + outputFilename);
    }

    file << oss.str();
    file.close();

    std::cout << "完整PCB模板文件已生成: " << outputFilename << std::endl;
    std::cout << "共输出 " << moduleCount << " 个module" << std::endl;
}




// ========== 全局变量存储最新生成的文件名 ==========
std::string g_latestGeneratedFile = "";





//

int main() {
    try {
        // === 1) 初始化 Grid ===
        std::string filename = "testcase.kicad_pcb";   // 按需修改
        Grid grid;
        grid.SetUp(filename);
        std::cout << "Grid初始化成功! 网格尺寸: " << grid.width << " x " << grid.height << std::endl;
        
        
        std::cout<< std::fixed << std::setprecision(6)<<grid.min_x<<' '<<grid.min_y<<' '<<grid.max_x<<' '<<grid.max_y;


       
        // === 2) 参数 ===
        int layerId = 0;   // 层序号
        int blockedAreaId = 1;   // 阻塞区域编号
        int expandPixels = 4;   // 像素域扩展像素数（有补偿）
        int kWholePixels = 2;   // 现实域再扩 k 个整像素（与 0.5px 一起）
        double pxSizeMM = 1.0 / static_cast<double>(grid.inputScale);

        std::cout << "参数: 层=" << layerId << " (" << getLayerNameById(grid, layerId)
            << "), 阻塞区域=" << blockedAreaId
            << ", 像素域扩展=" << expandPixels << " px"
            << ", 现实域最终再扩 = 0.5px + " << kWholePixels << "px"
            << " ≈ " << std::fixed << std::setprecision(5)
            << (0.5 + kWholePixels) * pxSizeMM << " mm\n";

        // === 3) 读取真实板框（Edge.Cuts 外接矩形） ===
        KiCadParser parser2;
        if (!parser2.parseFile(filename)) {
            throw std::runtime_error("解析 KiCad 文件失败，无法获取板框");
        }
        KiCadParser::Point2D bl, tl, tr, br; // 左下、左上、右上、右下
        if (!parser2.getBoardCorners(bl, tl, tr, br)) {
            throw std::runtime_error("未找到 Edge.Cuts 板框");
        }
        const double boardMinX = bl.x;
        const double boardMinY = bl.y;
        const double boardMaxX = tr.x;
        const double boardMaxY = tr.y;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[Board] min=(" << boardMinX << "," << boardMinY
            << "), max=(" << boardMaxX << "," << boardMaxY << ")\n";

        // === 4) 打印“像素扩展阶段”的像素坐标（有补偿） ===
        // 4.1 阻塞区域像素（1-based）
        auto blockedPixels = GetBlockedAreaPixels_1Based(grid, layerId, blockedAreaId);
        if (blockedPixels.empty()) {
            throw std::runtime_error("未找到指定的阻塞区域");
        }

        // 4.2 原始像素包围盒（1-based）
        BBoxPx1 originalBBox = ComputeBBoxFromPixels1(blockedPixels);
        std::cout << "\n=== 原始像素包围盒(1-based, col,x  row,y) ===\n";
        std::cout << "LT(col,row)=(" << originalBBox.col_min << "," << originalBBox.row_min << ")\n";
        std::cout << "RT(col,row)=(" << originalBBox.col_max << "," << originalBBox.row_min << ")\n";
        std::cout << "RB(col,row)=(" << originalBBox.col_max << "," << originalBBox.row_max << ")\n";
        std::cout << "LB(col,row)=(" << originalBBox.col_min << "," << originalBBox.row_max << ")\n";
        std::cout << "W x H (px) = " << (originalBBox.col_max - originalBBox.col_min + 1)
            << " x " << (originalBBox.row_max - originalBBox.row_min + 1) << "\n";

        // 4.3 像素扩展（有补偿）
        ClipReport clipReport{};
        BBoxPx1 expandedBBoxPx = ExpandBBoxByMPixelsCompensated(originalBBox, expandPixels, grid, &clipReport);

        std::cout << "\n=== 像素扩展后包围盒(有补偿, 1-based) ===\n";
        std::cout << "LT(col,row)=(" << expandedBBoxPx.col_min << "," << expandedBBoxPx.row_min << ")\n";
        std::cout << "RT(col,row)=(" << expandedBBoxPx.col_max << "," << expandedBBoxPx.row_min << ")\n";
        std::cout << "RB(col,row)=(" << expandedBBoxPx.col_max << "," << expandedBBoxPx.row_max << ")\n";
        std::cout << "LB(col,row)=(" << expandedBBoxPx.col_min << "," << expandedBBoxPx.row_max << ")\n";
        std::cout << "W x H (px) = " << (expandedBBoxPx.col_max - expandedBBoxPx.col_min + 1)
            << " x " << (expandedBBoxPx.row_max - expandedBBoxPx.row_min + 1) << "\n";
        std::cout << "触边: left=" << clipReport.left
            << ", right=" << clipReport.right
            << ", top=" << clipReport.top
            << ", bottom=" << clipReport.bottom << "\n";

        // === 5) （可选）打印“像素域扩展后的现实四角”（用于对比）
        auto baseCorners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
        std::vector<std::string> names = { "左上(LT)","右上(RT)","右下(RB)","左下(LB)" };
        std::cout << "\n=== 像素域扩展后的现实四角(用于对比) ===\n";
        for (size_t i = 0; i < baseCorners.size(); ++i) {
            std::cout << names[i] << ": (" << baseCorners[i].first << ", " << baseCorners[i].second << ")\n";
        }

        // === 6) 最终：现实域 0.5px + k×1px 对齐扩张 + 裁剪到板框 ===
        auto expandedCorners = GetExpandedPCBCornersSnapAndClip(
            grid, layerId, blockedAreaId,
            expandPixels, kWholePixels,
            boardMinX, boardMinY, boardMaxX, boardMaxY
        );

        std::cout << "\n=== 最终小PCB边界框(0.5px + " << kWholePixels
            << "px 对齐 + 板框裁剪后) ===\n";
        for (size_t i = 0; i < expandedCorners.size(); ++i) {
            std::cout << names[i] << ": (" << expandedCorners[i].first
                << ", " << expandedCorners[i].second << ")\n";
        }


        // === 6.1) 最终四角 -> 像素顶角坐标(1-based, row,col)，保证落在小板内部 ===
        auto toPxInside = [&](double x_mm, double y_mm, int cornerIdx) -> std::pair<int, int> {
            // cornerIdx: 0=LT, 1=RT, 2=RB, 3=LB
            const double eps = 1.0 / (grid.inputScale * 1024.0); // 远小于1px
            double xx = x_mm, yy = y_mm;
            if (cornerIdx == 1 || cornerIdx == 2) xx -= eps; // RT/RB：向左收一点
            if (cornerIdx == 2 || cornerIdx == 3) yy -= eps; // RB/LB：向上收一点
            return RealMMToPixelPx1(grid, xx, yy); // 返回( row , col )，1-based
            };

        std::array<std::pair<int, int>, 4> finalCornerPx;
        for (int i = 0; i < 4; ++i) {
            finalCornerPx[i] = toPxInside(expandedCorners[i].first, expandedCorners[i].second, i);
        }

        std::cout << "\n=== 最终小PCB四角对应的像素顶角坐标(1-based, row,col) ===\n";
        const char* nm[4] = { "LT","RT","RB","LB" };
        for (int i = 0; i < 4; ++i) {
            std::cout << nm[i] << ": (" << finalCornerPx[i].first
                << ", " << finalCornerPx[i].second << ")\n";
        }

        // === 7) 交点/终止点（仍按像素域扩展的矩形做，保持你原逻辑） ===
        auto points = FindIntersectionsAndEndpointsInSmallPCB(grid, layerId, blockedAreaId, expandPixels);

        std::cout << "\n=== 线段交点和终止点 ===\n";
        std::cout << "共找到 " << points.size() << " 个点\n";

        // === 8) 生成模板（用最终四角） ===
        GenerateCompletePCBTemplate(points, expandedCorners, "complete_pcb_template.kicad_pcb");

    }
    catch (const std::exception& e) {
        std::cerr << "程序运行出错: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "程序执行完成!" << std::endl;
    return 0;
}



