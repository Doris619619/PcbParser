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


#include "First_Part12.h"
#include "grid.h"

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



#include <iostream>
#include <vector>
#include <string>
#include "grid.h"
#include "First_Part12.h"

#include <iostream>
#include <vector>
#include <string>
#include "grid.h"
#include "First_Part12.h"

// 计算两点间的距离
double distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}

// 判断线段是否与矩形相交
bool isIntersecting(double x1, double y1, double x2, double y2, const RealRectMM& rect) {
    // 判断线段是否与矩形边界相交
    return !(x2 < rect.x_left || x1 > rect.x_right || y2 < rect.y_top || y1 > rect.y_bot);
}

// 判断点是否在矩形内部
bool isPointInsideRect(double px, double py, const RealRectMM& rect) {
    return px >= rect.x_left && px <= rect.x_right && py >= rect.y_top && py <= rect.y_bot;
}

// 输出交点的现实坐标
void outputIntersection(double x, double y) {
    std::cout << "交点坐标: (" << x << ", " << y << ")" << std::endl;
}


// 计算并输出所有交点和终点
void findIntersectionAndEndPoints(const Grid& grid, int layerId, int blockedAreaId, const RealRectMM& pcbBox) {
    std::string layerName = getLayerNameById(grid, layerId); // 通过层编号获取层名称

    // 获取阻塞区域的像素坐标
    std::vector<std::pair<int, int>> blockedAreaPixels = GetBlockedAreaPixels_1Based(grid, layerId, blockedAreaId);

    if (blockedAreaPixels.empty()) {
        std::cerr << "未找到阻塞区域 " << blockedAreaId << "!" << std::endl;
        return;
    }

    std::cout << "正在检查阻塞区域 " << blockedAreaId << " 在层 " << layerName << " 上的线段..." << std::endl;

    // 遍历所有线段
    for (const auto& segment : grid.segments) {
        double x1, y1, x2, y2;
        std::string layer;

        // 获取线段的起始和结束坐标，以及所在层
        for (std::shared_ptr<Node>& child : segment->children) {
            if (child->name == "start") {
                x1 = std::stod(child->parameters[0]);
                y1 = std::stod(child->parameters[1]);
            }
            if (child->name == "end") {
                x2 = std::stod(child->parameters[0]);
                y2 = std::stod(child->parameters[1]);
            }
            if (child->name == "layer") {
                layer = child->parameters[0];  // 获取层信息
            }
        }

        // 判断该线段是否在指定的层中
        if (layer != layerName) continue;

        // 检查线段是否与小PCB的边界相交
        if (isIntersecting(x1, y1, x2, y2, pcbBox)) {
            // 输出交点
            outputIntersection(x1, y1);
            outputIntersection(x2, y2);

            // 判断线段的交点是否在小PCB内部或边界上
            if (isPointInsideRect(x1, y1, pcbBox)) {
                std::cout << "交点 (" << x1 << ", " << y1 << ") 在小PCB区域内或边界上" << std::endl;
            }

            if (isPointInsideRect(x2, y2, pcbBox)) {
                std::cout << "交点 (" << x2 << ", " << y2 << ") 在小PCB区域内或边界上" << std::endl;
            }

            // 判断线段终点是否在小PCB内部，如果是则输出该终点
            if (distance(x1, y1, pcbBox.x_left, pcbBox.y_top) < distance(x2, y2, pcbBox.x_left, pcbBox.y_top)) {
                if (isPointInsideRect(x2, y2, pcbBox)) {
                    std::cout << "终止点坐标: (" << x2 << ", " << y2 << ") 在小PCB区域内或边界上" << std::endl;
                }
            }
            else {
                if (isPointInsideRect(x1, y1, pcbBox)) {
                    std::cout << "终止点坐标: (" << x1 << ", " << y1 << ") 在小PCB区域内或边界上" << std::endl;
                }
            }
        }
    }
}

int main() {
    // 假设已成功解析KiCad文件并初始化Grid
    Grid grid;
    std::string filename = "testcase.kicad_pcb";  // 你的KiCad PCB文件路径
    grid.SetUp(filename);  // 初始化 Grid，读取 PCB 文件并解析

    // 获取小PCB区域的边界（现实坐标）
    RealRectMM pcbBox = PixelBoxEdgesPx1ToRealRect(grid, 1, 1, grid.height, grid.width); // 这里使用的只是一个示例，实际坐标要根据需要来调整

    // 输入层的序号和阻塞区域编号
    int layerId = 1;  // 假设选择第0层
    int blockedAreaId = 1;  // 假设选择第1个阻塞区域

    // 输出交点和终点
    findIntersectionAndEndPoints(grid, layerId, blockedAreaId, pcbBox);

    return 0;
}
