// Separate.cpp
// 约定：对外像素坐标均为 1-based，左上角像素为 (1,1)；列向右增大、行向下增大。
// 现实坐标单位：mm（KiCad 坐标系），并与屏幕一致：向下 y 增大。

#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>
#include <array>     // for std::array
#include <climits>   // for INT_MAX, INT_MIN

#include "blocked_area_analyzer.h"   // 内含 grid.h（而 grid.h 内含 First_Part12.h）

// ====== 可调参数（你可按项目实际修改） ======
static constexpr double kDefaultClearanceMM = 0.20; // 缺省工艺间隙
static constexpr double kRouterExtraMM = 1.00; // 额外走线/操作空间
static constexpr int    kSafetyPx = 2;    // 安全像素冗余

// —— 数据结构：把像素与现实坐标一起打包 —— //
struct BlockedPixel {
    int row1;      // 1-based，左上为 (1,1)
    int col1;      // 1-based
    double x_mm;   // 该像素正中心的 KiCad X（mm）
    double y_mm;   // 该像素正中心的 KiCad Y（mm），向下增大
};

// —— 像素中心 <-> 现实坐标（mm） —— //
// 左上 (1,1) 像素中心： x = min_x + (c - 0.5)/scale,  y = min_y + (r - 0.5)/scale
inline std::pair<double, double>
PixelCenterPx1ToRealMM(const Grid& grid, int row1, int col1)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    row1 = std::max(1, std::min(row1, grid.height));
    col1 = std::max(1, std::min(col1, grid.width));
    const double scale = static_cast<double>(grid.inputScale);
    const double x_mm = grid.min_x + (static_cast<double>(col1) - 0.5) / scale;
    const double y_mm = grid.min_y + (static_cast<double>(row1) - 0.5) / scale; // 向下增大
    return { x_mm, y_mm };
}

// 现实坐标 -> 1-based 像素
inline std::pair<int, int>
RealMMToPixelPx1(const Grid& grid, double x_mm, double y_mm)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    const double scale = static_cast<double>(grid.inputScale);
    int col1 = static_cast<int>(std::floor((x_mm - grid.min_x) * scale)) + 1;
    int row1 = static_cast<int>(std::floor((y_mm - grid.min_y) * scale)) + 1;
    col1 = std::max(1, std::min(col1, grid.width));
    row1 = std::max(1, std::min(row1, grid.height));
    return { row1, col1 };
}

// —— 现实“像素边界”矩形（mm） —— //
// 由 1-based 像素包围盒（含端点像素）得到外沿矩形的 mm 顶点：
// x_left  = min_x + (col_min - 1)/scale,  x_right = min_x + col_max/scale
// y_top   = min_y + (row_min - 1)/scale,  y_bot   = min_y + row_max/scale
struct RealRectMM {
    double x_left, y_top, x_right, y_bot;
};
inline RealRectMM PixelBoxEdgesPx1ToRealRect(const Grid& grid, int row_min, int col_min,
    int row_max, int col_max) {
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    const double s = static_cast<double>(grid.inputScale);
    RealRectMM rr;
    rr.x_left = grid.min_x + (static_cast<double>(col_min) - 1.0) / s;
    rr.x_right = grid.min_x + (static_cast<double>(col_max)) / s;
    rr.y_top = grid.min_y + (static_cast<double>(row_min) - 1.0) / s;
    rr.y_bot = grid.min_y + (static_cast<double>(row_max)) / s;
    return rr;
}

// —— 输入层号 + 阻塞区域编号 → 返回像素（1-based）—— //
std::vector<std::pair<int, int>>
GetBlockedAreaPixels_1Based(const Grid& grid, int layer, int blockedAreaID)
{
    if (layer < 0 || layer >= grid.Layers) {
        throw std::runtime_error("层号超出范围");
    }
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

// —— 像素 + 现实中心打包 —— //
std::vector<BlockedPixel>
CollectBlockedPixelsWithReal(const Grid& grid, int layer, int blockedAreaID)
{
    std::vector<BlockedPixel> result;
    auto pxs = GetBlockedAreaPixels_1Based(grid, layer, blockedAreaID);
    result.reserve(pxs.size());
    for (auto [r1, c1] : pxs) {
        auto [x_mm, y_mm] = PixelCenterPx1ToRealMM(grid, r1, c1);
        result.push_back({ r1, c1, x_mm, y_mm });
    }
    return result;
}

// —— 包围盒（1-based） —— //
struct BBoxPx1 {
    int row_min, col_min; // top-left (含)
    int row_max, col_max; // bottom-right (含)
};
static BBoxPx1 ComputeBBoxFromPixels1(const std::vector<std::pair<int, int>>& pxs1)
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

static std::array<std::pair<int, int>, 4>
BBoxCorners_LB_RB_RT_LT(const BBoxPx1& b)
{
    // 左上为 (1,1)，因此 “左下” 是 (row_max, col_min)
    std::array<std::pair<int, int>, 4> corners;
    corners[0] = { b.row_max, b.col_min }; // 左下 LB
    corners[1] = { b.row_max, b.col_max }; // 右下 RB
    corners[2] = { b.row_min, b.col_max }; // 右上 RT
    corners[3] = { b.row_min, b.col_min }; // 左上 LT
    return corners;
}

// ====== 新增：判断 pad 是否作用于某层（保守策略） ======
static bool PadCoversLayer(const Grid& grid, int layer, const std::shared_ptr<Node>& padNode)
{
    // 在 pad 子节点里找 (layer ...) 或 (layers ...)
    bool foundLayerSpec = false;
    bool covers = false;

    for (const auto& ch : padNode->children) {
        if (!ch) continue;
        if (ch->name != "layer" && ch->name != "layers") continue;
        foundLayerSpec = true;

        for (const auto& param : ch->parameters) {
            // 显式层名
            int id = grid.getLayerId(param);
            if (id == layer) { covers = true; break; }

            // 通配：*.Cu / F&B 通用（保守：认为覆盖所有铜层）
            if (param == "*.Cu" || param.find(".Cu") != std::string::npos) {
                covers = true; break;
            }
        }
        if (covers) break;
    }

    // 没写层也按“覆盖所有层”处理（保守）
    if (!foundLayerSpec) covers = true;
    return covers;
}

// ====== 新增：计算某层“最大 pad 外形尺寸”（mm） ======
static double GetLayerMaxPadSizeMM(const Grid& grid, int layer)
{
    double maxPad = 0.0; // mm
    // 通过 grid.footprints 遍历所有 pad
    for (const auto& fp : grid.footprints) {
        if (!fp) continue;
        for (const auto& ch : fp->children) {
            if (!ch || ch->name != "pad") continue;
            if (!PadCoversLayer(grid, layer, ch)) continue;

            double w = 0.0, h = 0.0;
            // 找 (size w h)
            for (const auto& sub : ch->children) {
                if (!sub) continue;
                if (sub->name == "size" && sub->parameters.size() >= 2) {
                    try {
                        w = std::stod(sub->parameters[0]);
                        h = std::stod(sub->parameters[1]);
                    }
                    catch (...) { /* 忽略异常 */ }
                    break;
                }
            }
            // 取主轴长度（矩形/圆/圆角矩形都用 max(w,h)）
            maxPad = std::max(maxPad, std::max(w, h));
        }
    }

    // 如果什么都没找到（比如解析不到 pad），为了安全用一个保守最小值
    if (maxPad <= 0) maxPad = 1.0; // 1mm 兜底，避免 0 扩张
    return maxPad;
}

// ====== 新增：按“能放下最大 pad + 走线余量”扩张包围盒 ======
static BBoxPx1 ExpandBBoxForPadAndRouting(const Grid& grid,
    const BBoxPx1& inBox,
    int layer,
    double clearance_mm = kDefaultClearanceMM,
    double router_extra = kRouterExtraMM,
    int safety_px = kSafetyPx)
{
    const double max_pad_mm = GetLayerMaxPadSizeMM(grid, layer);
    const double edge_strip_mm = max_pad_mm + 2.0 * std::max(0.0, clearance_mm);
    const double total_margin_m = edge_strip_mm + std::max(0.0, router_extra);

    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    int margin_px = static_cast<int>(std::ceil(total_margin_m * grid.inputScale)) + std::max(0, safety_px);

    BBoxPx1 out = inBox;
    out.row_min = std::max(1, inBox.row_min - margin_px);
    out.col_min = std::max(1, inBox.col_min - margin_px);
    out.row_max = std::min(grid.height, inBox.row_max + margin_px);
    out.col_max = std::min(grid.width, inBox.col_max + margin_px);

    // 若尺寸过窄（极小区域），保证最小厚度 >= edge_strip（两边各 edge_strip）
    int min_half = std::max(1, margin_px);
    if (out.row_max - out.row_min + 1 < 2 * min_half)
    {
        int mid = (inBox.row_min + inBox.row_max) / 2;
        out.row_min = std::max(1, mid - min_half);
        out.row_max = std::min(grid.height, mid + min_half);
    }
    if (out.col_max - out.col_min + 1 < 2 * min_half)
    {
        int mid = (inBox.col_min + inBox.col_max) / 2;
        out.col_min = std::max(1, mid - min_half);
        out.col_max = std::min(grid.width, mid + min_half);
    }

    return out;
}

// ====== 新增：返回扩张后四角（像素 & 现实） ======
struct ExpandedBoxOut {
    std::array<std::pair<int, int>, 4> corners_px1_LB_RB_RT_LT;
    std::array<std::pair<double, double>, 4> corners_mm_LB_RB_RT_LT;
    BBoxPx1 box_px1;
    RealRectMM rect_mm;
};

ExpandedBoxOut
GetExpandedBoxForArea(const Grid& grid, int layer, int areaID,
    double clearance_mm = kDefaultClearanceMM,
    double router_extra = kRouterExtraMM,
    int safety_px = kSafetyPx)
{
    auto cells1 = GetBlockedAreaPixels_1Based(grid, layer, areaID);
    if (cells1.empty()) {
        throw std::runtime_error("指定区域不存在或没有像素");
    }
    BBoxPx1 bb0 = ComputeBBoxFromPixels1(cells1);
    BBoxPx1 bb = ExpandBBoxForPadAndRouting(grid, bb0, layer, clearance_mm, router_extra, safety_px);

    // 像素四角（与先前函数一致顺序：LB, RB, RT, LT）
    auto px4 = BBoxCorners_LB_RB_RT_LT(bb);

    // 现实矩形（mm）：像素边界
    RealRectMM rr = PixelBoxEdgesPx1ToRealRect(grid, bb.row_min, bb.col_min, bb.row_max, bb.col_max);

    // 转换为四个 mm 顶点（LB, RB, RT, LT）
    std::array<std::pair<double, double>, 4> mm4 = {
        std::pair<double,double>{ rr.x_left,  rr.y_bot }, // LB
        std::pair<double,double>{ rr.x_right, rr.y_bot }, // RB
        std::pair<double,double>{ rr.x_right, rr.y_top }, // RT
        std::pair<double,double>{ rr.x_left,  rr.y_top }  // LT
    };

    return { px4, mm4, bb, rr };
}

// ====== 旧：获取原始包围盒四角（像素） ======
std::array<std::pair<int, int>, 4>
GetBlockedAreaBBoxCorners_1Based(const Grid& grid, int layer, int areaID)
{
    auto cells1 = GetBlockedAreaPixels_1Based(grid, layer, areaID);
    if (cells1.empty()) {
        return std::array<std::pair<int, int>, 4>{
            std::pair<int, int>{0, 0},
                std::pair<int, int>{0, 0},
                std::pair<int, int>{0, 0},
                std::pair<int, int>{0, 0}
        };
    }
    auto bbox = ComputeBBoxFromPixels1(cells1);
    return BBoxCorners_LB_RB_RT_LT(bbox);
}

// ====== 演示 main ======
// 用法：./your_app [kicad_text_file]
// 若未给文件名，默认读取 "testcase.txt"
// 演示：layer=0, areaID=1
int main(int argc, char** argv)
{
    try {
        std::string filename = (argc >= 2) ? argv[1] : "testcase.kicad_pcb";

        Grid grid;
        // grid.setInputScale(1); // 需要的话可以改
        grid.SetUp(filename);

        int layer = 0;
        int areaID = 1;

        // 1) 打印该区域内所有像素及其中心的现实坐标
        auto items = CollectBlockedPixelsWithReal(grid, layer, areaID);
        std::cout << std::fixed << std::setprecision(5);
        std::cout << "像素总数: " << items.size() << "\n";
        for (const auto& bp : items) {
            std::cout << "Pixel(row=" << bp.row1 << ", col=" << bp.col1 << ")  ->  "
                << "Real(mm): (x=" << bp.x_mm << ", y=" << bp.y_mm << ")\n";
        }

        // 2) 计算并打印扩张后的四角（像素 & 现实）
        auto out = GetExpandedBoxForArea(grid, layer, areaID);
        const char* names[4] = { "LB","RB","RT","LT" };

        std::cout << "\n[扩张边框] 四角像素 (1-based):\n";
        for (int i = 0; i < 4; i++) {
            std::cout << names[i] << ": (row=" << out.corners_px1_LB_RB_RT_LT[i].first
                << ", col=" << out.corners_px1_LB_RB_RT_LT[i].second << ")\n";
        }

        std::cout << "\n[扩张边框] 四角现实 (mm):\n";
        for (int i = 0; i < 4; i++) {
            std::cout << names[i] << ": (x=" << out.corners_mm_LB_RB_RT_LT[i].first
                << ", y=" << out.corners_mm_LB_RB_RT_LT[i].second << ")\n";
        }

        // 附：打印用于扩张的几个关键数值，便于调参
        double maxPad = GetLayerMaxPadSizeMM(grid, layer);
        std::cout << "\n本层最大 pad 尺寸（主轴 mm）: " << maxPad << "\n";
        std::cout << "inputScale: " << grid.inputScale << " px/mm\n";
        std::cout << "扩张后像素包围盒: rows[" << out.box_px1.row_min << "," << out.box_px1.row_max
            << "] cols[" << out.box_px1.col_min << "," << out.box_px1.col_max << "]\n";

    }
    catch (const std::exception& ex) {
        std::cerr << "[错误] " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
