// Divide.cpp
// 对外统一坐标标准：左下角像素 = (1,1)，所有对外坐标一律 1-based。
// 内部仍使用 Grid 的 0-based (row,col) 存储，通过边界转换适配。
// 功能：
// 1) 输入（层号、阻塞区域编号）→ 返回该阻塞区域内所有像素坐标（1-based）
// 2) 像素坐标(1-based) ↔ 现实坐标(KiCad, mm) 双向转换
// 3) 示例 main：默认读取 A30.txt，演示 layer=0、areaID=1 的输出

#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <string>
#include "blocked_area_analyzer.h"  // 已包含 grid.h

// ---------------- 类型别名 ----------------
using Pixel1 = std::pair<int, int>;            // (row1, col1) —— 1-based，对外可见
using PixelList1 = std::vector<Pixel1>;
using RealPoint = std::pair<double, double>;   // (x, y) in KiCad units (mm)

// ---------------- 工具：内部<->外部坐标转换 ----------------
// 外 -> 内：1-based (r1,c1) → 0-based (r0,c0)
static inline std::pair<int, int> ToZeroBased(int r1, int c1) {
    return { r1 - 1, c1 - 1 };
}

// 内 -> 外：0-based (r0,c0) → 1-based (r1,c1)
static inline Pixel1 ToOneBased(int r0, int c0) {
    return { r0 + 1, c0 + 1 };
}

// ---------------- 阻塞区域查询（对外 1-based） ----------------
// 输入：层号、阻塞区域编号；输出：该阻塞区域内所有像素坐标（1-based）
PixelList1 GetBlockedAreaPixels_1Based(const Grid& grid, int layer, int blockedAreaID) {
    if (layer < 0 || layer >= grid.Layers)
        throw std::invalid_argument("GetBlockedAreaPixels_1Based: layer 越界");
    if (blockedAreaID <= 0)
        throw std::invalid_argument("GetBlockedAreaPixels_1Based: blockedAreaID 必须从 1 开始");

    // 调用同学写的阻塞区域分析（内部是 0-based 行列）
    LayerBlockedAnalysis ana = BlockedAreaAnalyzer::analyzeLayerBlockedAreas(grid, layer);
    //（该分析返回的结构里，阻塞像素坐标存放在 cells 向量中）:contentReference[oaicite:0]{index=0}

    for (const auto& area : ana.blockedAreas) {
        if (area.areaID == blockedAreaID) {
            PixelList1 out;
            out.reserve(area.cells.size());
            for (const auto& rc0 : area.cells) {
                // 将内部 0-based (row,col) 统一转换成 1-based
                out.push_back(ToOneBased(rc0.first, rc0.second));
            }
            return out;
        }
    }
    return {}; // 未找到该阻塞区域
}

// ---------------- 像素(1-based) ↔ 现实坐标 ----------------
// 约定：原点在板框左下角，x 向右、y 向上；grid.inputScale = (格/毫米)
// 一个像素(cell)的几何范围为 [col0, col0+1] × [row0, row0+1]（格）；中心在 +0.5

// (1) 1-based 像素中心 → 现实坐标(mm)
RealPoint PixelCenterToReal_1Based(const Grid& grid, int row1, int col1) {
    auto [row0, col0] = ToZeroBased(row1, col1);
    if (row0 < 0 || row0 >= grid.height || col0 < 0 || col0 >= grid.width)
        throw std::invalid_argument("PixelCenterToReal_1Based: 像素索引越界");

    // 注意：grid.min_x/min_y 是板框左下角；row/col 越大，物理 y/x 越大
    double x = grid.min_x + (static_cast<double>(col0) + 0.5) / static_cast<double>(grid.inputScale);
    double y = grid.min_y + (static_cast<double>(row0) + 0.5) / static_cast<double>(grid.inputScale);
    return { x, y };
}

// (2) 现实坐标 → 所在像素(1-based)
// col0 = floor( (x - min_x) * inputScale )，row0 同理；最后统一 +1
Pixel1 RealToPixel_1Based(const Grid& grid, double x, double y) {
    double gx = (x - grid.min_x) * static_cast<double>(grid.inputScale);
    double gy = (y - grid.min_y) * static_cast<double>(grid.inputScale);
    int col0 = static_cast<int>(std::floor(gx));
    int row0 = static_cast<int>(std::floor(gy));

    // 边界 clamp，避免浮点误差落到 [width,height] 之外
    if (col0 < 0) col0 = 0;
    if (row0 < 0) row0 = 0;
    if (col0 >= grid.width)  col0 = grid.width - 1;
    if (row0 >= grid.height) row0 = grid.height - 1;

    return ToOneBased(row0, col0);
}

// ---------------- 示例 main ----------------
// 用法：./your_app [kicad_text_file]
// 不传参时默认读取 "A30.txt"（你也可以改成 L92.txt / testcase.txt 等）
int main(int argc, char** argv) {
    try {
        std::string filename = (argc >= 2) ? argv[1] : "testcase.kicad_pcb";

        Grid grid;
        grid.SetUp(filename);  // 解析 KiCad 文本并构建网格、离散化等（你们已有流程）:contentReference[oaicite:1]{index=1}

        int layer = 0;
        int areaID = 1;

        PixelList1 cells1 = GetBlockedAreaPixels_1Based(grid, layer, areaID);

        std::cout << "Layer " << layer << ", AreaID " << areaID
            << " 的像素数量 = " << cells1.size() << std::endl;

        // 打印所有像素（1-based）
        for (const auto& rc1 : cells1) {
            std::cout << "(row=" << rc1.first << ", col=" << rc1.second << ")";
            // 同时演示：该像素中心对应的现实坐标（mm）
            auto p = PixelCenterToReal_1Based(grid, rc1.first, rc1.second);
            std::cout << " -> center(x=" << p.first << ", y=" << p.second << ")\n";
        }

        // 演示：现实坐标 → 像素(1-based)
        {
            double sx = 100.02183, sy = 91.97056;
            auto rc1 = RealToPixel_1Based(grid, sx, sy);
            std::cout << "示例坐标 (" << sx << ", " << sy
                << ") 所在像素(1-based) = (row=" << rc1.first
                << ", col=" << rc1.second << ")\n";
        }
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "[Divide.cpp] 运行出错: " << e.what() << std::endl;
        return 1;
    }
}
