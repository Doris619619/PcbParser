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

#include "blocked_area_analyzer.h"   // 内含 grid.h

// —— 数据结构：把像素与现实坐标一起打包 —— //
struct BlockedPixel {
    int row1;      // 1-based，左上为 (1,1)
    int col1;      // 1-based
    double x_mm;   // 该像素正中心的 KiCad X（mm）
    double y_mm;   // 该像素正中心的 KiCad Y（mm），向下增大
};

// —— 基础换算（grid.inputScale = 每毫米多少格，>0）—— //
// 屏幕一致的 y 方向：向下增大
// 左上 (1,1) 的像素中心：
//   x_mm = min_x + (col1 - 0.5)/scale
//   y_mm = min_y + (row1 - 0.5)/scale      <<< 修正点
inline std::pair<double, double>
PixelCenterPx1ToRealMM(const Grid& grid, int row1, int col1)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    // 夹紧到有效像素范围
    row1 = std::max(1, std::min(row1, grid.height));
    col1 = std::max(1, std::min(col1, grid.width));

    const double scale = static_cast<double>(grid.inputScale);
    const double x_mm = grid.min_x + (static_cast<double>(col1) - 0.5) / scale;
    const double y_mm = grid.min_y + (static_cast<double>(row1) - 0.5) / scale; // 向下增大
    return { x_mm, y_mm };
}

// 现实坐标 (x_mm, y_mm) -> 1-based 像素 (row1, col1)
// 对应像素包含判断采用 floor；越界将被夹紧到 [1..H/W].
inline std::pair<int, int>
RealMMToPixelPx1(const Grid& grid, double x_mm, double y_mm)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    const double scale = static_cast<double>(grid.inputScale);

    // 列（x 向右）：从 min_x 开始
    int col1 = static_cast<int>(std::floor((x_mm - grid.min_x) * scale)) + 1;
    // 行（y 向下）：从 min_y 开始（与屏幕/鼠标一致）
    int row1 = static_cast<int>(std::floor((y_mm - grid.min_y) * scale)) + 1;

    // 夹紧
    col1 = std::max(1, std::min(col1, grid.width));
    row1 = std::max(1, std::min(row1, grid.height));
    return { row1, col1 };
}

// —— 主功能 1：输入层号 + 阻塞区域编号 → 返回该区域所有像素（1-based, 左上为原点）—— //
std::vector<std::pair<int, int>>
GetBlockedAreaPixels_1Based(const Grid& grid, int layer, int blockedAreaID)
{
    if (layer < 0 || layer >= grid.Layers) {
        throw std::runtime_error("层号超出范围");
    }

    // 用现有分析器找该层所有阻塞区域
    LayerBlockedAnalysis la = BlockedAreaAnalyzer::analyzeLayerBlockedAreas(grid, layer);

    // 找到目标 areaID
    const BlockedAreaInfo* target = nullptr;
    for (const auto& area : la.blockedAreas) {
        if (area.areaID == blockedAreaID) {
            target = &area;
            break;
        }
    }
    if (!target) {
        return {};
    }

    std::vector<std::pair<int, int>> out;
    out.reserve(target->cells.size());

    // BlockedAreaAnalyzer::cells 里是 0-based 的 (row0, col0)
    // 转成 1-based，且符合左上原点：(row1=row0+1, col1=col0+1)
    for (const auto& rc0 : target->cells) {
        int row1 = rc0.first + 1;  // 0-based -> 1-based
        int col1 = rc0.second + 1;
        out.emplace_back(row1, col1);
    }
    return out;
}

// —— 主功能 2：把像素与现实坐标一并打包返回 —— //
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

// —— 计算包围盒（1-based）与角点 —— //
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
    // 注意：左上为 (1,1)，因此“左下”是 (row_max, col_min)
    std::array<std::pair<int, int>, 4> corners;
    corners[0] = { b.row_max, b.col_min }; // 左下
    corners[1] = { b.row_max, b.col_max }; // 右下
    corners[2] = { b.row_min, b.col_max }; // 右上
    corners[3] = { b.row_min, b.col_min }; // 左上
    return corners;
}

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

// —— 演示 main —— //
// 用法：./your_app [kicad_text_file]
// 若未给文件名，默认读取 "testcase.txt"
// 固定演示：layer=0, areaID=1
int main(int argc, char** argv)
{
    try {
        std::string filename = (argc >= 2) ? argv[1] : "testcase.kicad_pcb";

        Grid grid;
        // 如需修改 inputScale，请在 SetUp 之前调用：
        // grid.setInputScale(1);  // 例如：1 像素 = 1 mm
        grid.SetUp(filename);

        int layer = 0;
        int areaID = 1;

        auto items = CollectBlockedPixelsWithReal(grid, layer, areaID);

        if (items.empty()) {
            std::cout << "[提示] 未找到层 " << layer << " 的阻塞区域 " << areaID << " 或该区域为空。\n";
            return 0;
        }

        std::cout << "文件: " << filename << "\n"
            << "层: " << layer << "，阻塞区域: " << areaID
            << "，像素总数: " << items.size() << "\n\n";

        std::cout << std::fixed << std::setprecision(5);
        for (const auto& bp : items) {
            std::cout << "Pixel( row=" << bp.row1 << ", col=" << bp.col1 << " )  ->  "
                << "Real(mm): (x=" << bp.x_mm << ", y=" << bp.y_mm << ")\n";
        }

        // （可选）打印包围盒四角
        auto corners = GetBlockedAreaBBoxCorners_1Based(grid, layer, areaID);
        std::cout << "\nBBox corners (LB, RB, RT, LT) 1-based:\n";
        const char* names[4] = { "LB","RB","RT","LT" };
        for (int i = 0; i < 4; i++) {
            std::cout << names[i] << ": (row=" << corners[i].first
                << ", col=" << corners[i].second << ")\n";
        }

    }
    catch (const std::exception& ex) {
        std::cerr << "[错误] " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
