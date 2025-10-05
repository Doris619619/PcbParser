// Divide4.cpp  —— 统一像素坐标为“左上角(1,1)”的正式版
// 对外：像素坐标 1-based，左上原点；对内：沿用 Grid 的 0-based、左下原点。
// 功能：
//   1) 输入(layer, areaID) → 返回区域内像素集合 (1-based, 左上)
//   2) 计算区域最小包围盒 & 最终扩张包围盒（像素四角）
//   3) 像素(1-based, 左上) ↔ 现实坐标(mm) 双向转换
//   4) main: 读取A30.txt，layer=0, areaID=1，打印两个包围盒四角（像素+物理中心）

#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <string>
#include <array>
#include <limits>
#include <algorithm>

#include "blocked_area_analyzer.h"  // 内含 grid.h

// ---------------- 类型别名 ----------------
using Pixel1 = std::pair<int, int>;          // (row1, col1) —— 1-based, 左上原点
using PixelList1 = std::vector<Pixel1>;
using RealPoint = std::pair<double, double>;    // (x, y) in mm (KiCad)

// ---------------- 外(左上1基) ↔ 内(左下0基) 的唯一换算口 ----------------
// ！！！行用 grid.height 翻转，列仅做 ±1，严禁把 width 用到行上！！！
static inline std::pair<int, int> ToZeroBased(const Grid& grid, int row1, int col1) {
    int row0 = grid.height - row1; // 上→下（外） 变成 下→上（内）
    int col0 = col1 - 1;
    return { row0, col0 };
}
static inline std::pair<int, int> ToOneBased(const Grid& grid, int row0, int col0) {
    int row1 = grid.height - row0; // 下→上（内） 变成 上→下（外）
    int col1 = col0 + 1;
    return { row1, col1 };
}

// ---------------- mm <-> px 小工具 ----------------
static inline int    mmToPxCeil(const Grid& grid, double mm) { return (int)std::ceil(mm * (double)grid.inputScale); }
static inline int    mmToPxFloor(const Grid& grid, double mm) { return (int)std::floor(mm * (double)grid.inputScale); }
static inline double pxToMm(const Grid& grid, int px) { return (double)px / (double)grid.inputScale; }

// ---------------- 1-based像素最小包围盒 ----------------
struct BBoxPx1 {
    // 左上原点：row_min = top(数值小)，row_max = bottom(数值大)
    int row_min, col_min; // top-left  (含)
    int row_max, col_max; // bottom-right(含)
};

// 按 “左下, 右下, 右上, 左上” 顺序返回四个顶点（1-based, 左上原点语义）
static inline std::array<Pixel1, 4> BBoxCornersLB_RB_RT_LT(const BBoxPx1& b) {
    return {
        Pixel1{b.row_max, b.col_min}, // LB  bottom-left
        Pixel1{b.row_max, b.col_max}, // RB  bottom-right
        Pixel1{b.row_min, b.col_max}, // RT  top-right
        Pixel1{b.row_min, b.col_min}  // LT  top-left
    };
}

// 由像素集合(1-based, 左上)计算包围盒
static inline BBoxPx1 ComputeBBoxFromPixels1(const PixelList1& pxs1) {
    if (pxs1.empty()) throw std::runtime_error("ComputeBBoxFromPixels1: 输入像素为空");
    int rmin = std::numeric_limits<int>::max();
    int cmin = std::numeric_limits<int>::max();
    int rmax = std::numeric_limits<int>::min();
    int cmax = std::numeric_limits<int>::min();
    for (auto [r1, c1] : pxs1) {
        rmin = std::min(rmin, r1);
        cmin = std::min(cmin, c1);
        rmax = std::max(rmax, r1);
        cmax = std::max(cmax, c1);
    }
    return { rmin, cmin, rmax, cmax };
}

// ---------------- 获取阻塞区域像素 (1-based, 左上) ----------------
PixelList1 GetBlockedAreaPixels_1Based(const Grid& grid, int layer, int blockedAreaID) {
    if (layer < 0 || layer >= grid.Layers)
        throw std::invalid_argument("GetBlockedAreaPixels_1Based: layer 越界");
    if (blockedAreaID <= 0)
        throw std::invalid_argument("GetBlockedAreaPixels_1Based: blockedAreaID 必须从 1 开始");

    LayerBlockedAnalysis ana = BlockedAreaAnalyzer::analyzeLayerBlockedAreas(grid, layer);
    for (const auto& area : ana.blockedAreas) {
        if (area.areaID == blockedAreaID) {
            PixelList1 out;
            out.reserve(area.cells.size());
            for (const auto& rc0 : area.cells) {
                // area.cells 为内部 0-based(左下)；统一转为对外 1-based(左上)
                auto rc1 = ToOneBased(grid, rc0.first, rc0.second);
                out.push_back(rc1);
            }
            return out;
        }
    }
    return {};
}

// ---------------- 像素(1-based, 左上) ↔ 现实坐标(mm) ----------------
// 像素中心 -> 现实坐标
RealPoint PixelCenterToReal_1Based(const Grid& grid, int row1, int col1) {
    auto [row0, col0] = ToZeroBased(grid, row1, col1);
    if (row0 < 0 || row0 >= grid.height || col0 < 0 || col0 >= grid.width)
        throw std::invalid_argument("PixelCenterToReal_1Based: 像素索引越界");
    double x = grid.min_x + ((double)col0 + 0.5) / (double)grid.inputScale;
    double y = grid.min_y + ((double)row0 + 0.5) / (double)grid.inputScale;
    return { x, y };
}
// 现实坐标 -> 所在像素(1-based, 左上)
Pixel1 RealToPixel_1Based(const Grid& grid, double x, double y) {
    int col0 = (int)std::floor((x - grid.min_x) * (double)grid.inputScale);
    int row0 = (int)std::floor((y - grid.min_y) * (double)grid.inputScale);
    col0 = std::clamp(col0, 0, grid.width - 1);
    row0 = std::clamp(row0, 0, grid.height - 1);
    return ToOneBased(grid, row0, col0);
}

// ---------------- 该层最大 pad 外接尺寸（mm） ----------------
double GetMaxPadExtentMM_OnLayer(const Grid& grid, int layer) {
    double max_mm = 0.0;
    for (const auto& fp : grid.footprints) {
        if (!fp) continue;
        for (const auto& node : fp->children) {
            if (!node || node->name != "pad") continue;

            double w = 0.0, h = 0.0;
            std::vector<int> layersID;
            for (const auto& child : node->children) {
                if (!child) continue;
                if (child->name == "size" && child->parameters.size() >= 2) {
                    w = std::stod(child->parameters[0]);
                    h = std::stod(child->parameters[1]);
                }
                else if ((child->name == "layer" || child->name == "layers") && !child->parameters.empty()) {
                    for (const auto& lname : child->parameters) {
                        layersID.push_back(const_cast<Grid&>(grid).getLayerId(lname));
                    }
                }
            }
            if (!layersID.empty() &&
                std::find(layersID.begin(), layersID.end(), layer) != layersID.end()) {
                max_mm = std::max(max_mm, std::max(w, h));
            }
        }
    }
    return max_mm; // 若该层无pad，返回0
}

// ---------------- 扩张参数与规划 ----------------
enum class EdgePref { Any, Bottom, Right, Top, Left };
struct BoxPlanParams {
    double pad_clearance_mm = 0.25;
    double track_width_mm = 0.25;
    double extra_mm = 0.20;
    double pad_max_fallback_mm = 1.50;
    EdgePref preferred_edge = EdgePref::Any; // 若倾向某一侧放pad可改
};

// 按既定标准规划最终包围盒（1-based, 左上）
BBoxPx1 PlanFinalBox_1Based(const Grid& grid, int layer, const BBoxPx1& inBox,
    const BoxPlanParams& p) {
    // 1) pad 尺寸
    double pad_mm = GetMaxPadExtentMM_OnLayer(grid, layer);
    if (pad_mm <= 0.0) pad_mm = p.pad_max_fallback_mm;

    // 2) 量化
    int pad_px = mmToPxCeil(grid, pad_mm);
    int clear_px = mmToPxCeil(grid, p.pad_clearance_mm);
    int router_px = mmToPxCeil(grid, std::max(p.track_width_mm, p.pad_clearance_mm));
    int extra_px = mmToPxCeil(grid, p.extra_mm);

    // 3) 四边统一边距（保证任一侧可放一个最大pad并起线）
    int M_side_px = (pad_px + 1) / 2 + clear_px + router_px + extra_px;

    // 4) 内部最小宽/高
    int NeedDim_px = pad_px + 2 * (clear_px + router_px) + 2 * extra_px;

    // 5) 当前尺寸
    int cur_w = inBox.col_max - inBox.col_min + 1;
    int cur_h = inBox.row_max - inBox.row_min + 1;

    // 6) 初步统一外扩
    BBoxPx1 out = {
        std::max(1,           inBox.row_min - M_side_px),           // 向上：row减小
        std::max(1,           inBox.col_min - M_side_px),           // 向左：col减小
        std::min(grid.height, inBox.row_max + M_side_px),           // 向下：row增大
        std::min(grid.width , inBox.col_max + M_side_px)            // 向右：col增大
    };

    // 7) 对称补偿，确保达到 NeedDim
    auto widen_sym = [&]() {
        int w = out.col_max - out.col_min + 1;
        int h = out.row_max - out.row_min + 1;
        int add_w = std::max(0, NeedDim_px - w);
        int add_h = std::max(0, NeedDim_px - h);
        int left_add = (add_w + 1) / 2, right_add = add_w / 2;
        int up_add = (add_h + 1) / 2, down_add = add_h / 2;

        out.col_min = std::max(1, out.col_min - left_add);
        out.col_max = std::min(grid.width, out.col_max + right_add);
        // 左上原点：向“上”减小 row，向“下”增大 row
        out.row_min = std::max(1, out.row_min - up_add);
        out.row_max = std::min(grid.height, out.row_max + down_add);
        };
    widen_sym();

    // 8) 偏好侧 bonus（给该侧再多留一点“pad沟槽”）
    auto add_side_bonus = [&](EdgePref side) {
        int bonus = (pad_px + 1) / 2 + router_px;
        switch (side) {
        case EdgePref::Bottom: out.row_max = std::min(grid.height, out.row_max + bonus); break;
        case EdgePref::Top:    out.row_min = std::max(1, out.row_min - bonus); break;
        case EdgePref::Left:   out.col_min = std::max(1, out.col_min - bonus); break;
        case EdgePref::Right:  out.col_max = std::min(grid.width, out.col_max + bonus); break;
        default: break;
        }
        };
    add_side_bonus(p.preferred_edge);

    // 9) 最终 clamp
    out.row_min = std::max(1, out.row_min);
    out.col_min = std::max(1, out.col_min);
    out.row_max = std::min(grid.height, out.row_max);
    out.col_max = std::min(grid.width, out.col_max);

    // 10) 仍不足时给提示（通常是靠板边）
    int final_w = out.col_max - out.col_min + 1;
    int final_h = out.row_max - out.row_min + 1;
    if (final_w < NeedDim_px || final_h < NeedDim_px) {
        std::cerr << "[Warn] 扩张后仍不足以容纳 pad/走线/DRC（邻近板边） "
            << "need=" << NeedDim_px << "px, got=(" << final_w << "x" << final_h << ")\n";
    }
    return out;
}

// ---------------- 便捷API：最小/最终包围盒四角 ----------------
std::array<Pixel1, 4> GetBlockedAreaBBoxCorners_1Based(const Grid& grid, int layer, int areaID) {
    auto cells1 = GetBlockedAreaPixels_1Based(grid, layer, areaID);
    auto bbox = ComputeBBoxFromPixels1(cells1);
    return BBoxCornersLB_RB_RT_LT(bbox);
}
std::array<Pixel1, 4> GetFinalBBoxCorners_1Based(const Grid& grid, int layer, int areaID, const BoxPlanParams& plan) {
    auto cells1 = GetBlockedAreaPixels_1Based(grid, layer, areaID);
    auto bbox_min = ComputeBBoxFromPixels1(cells1);
    auto bbox_ex = PlanFinalBox_1Based(grid, layer, bbox_min, plan);
    return BBoxCornersLB_RB_RT_LT(bbox_ex);
}

// ---------------- main ----------------
int main(int argc, char** argv) {
    try {
        std::string filename = (argc >= 2) ? argv[1] : "testcase.kicad_pcb";

        Grid grid;
        grid.SetUp(filename);

        int layer = 0;
        int areaID = 1;

        // ① 最小包围盒
        auto corners_bbox = GetBlockedAreaBBoxCorners_1Based(grid, layer, areaID);
        std::cout << "[BoundingBox 1-based, Top-Left origin]\n";
        const char* labels[4] = { "LB", "RB", "RT", "LT" };
        for (int i = 0; i < 4; ++i) {
            auto [r, c] = corners_bbox[i];
            auto [x, y] = PixelCenterToReal_1Based(grid, r, c);
            std::cout << "  " << labels[i] << ": (row=" << r << ", col=" << c << ") -> center(x=" << x << ", y=" << y << ")\n";
        }

        // ② 扩张后的最终包围盒
        BoxPlanParams plan;
        plan.pad_clearance_mm = 0.25;
        plan.track_width_mm = 0.25;
        plan.extra_mm = 0.20;
        plan.pad_max_fallback_mm = 1.50;
        plan.preferred_edge = EdgePref::Any;

        auto corners_final = GetFinalBBoxCorners_1Based(grid, layer, areaID, plan);
        std::cout << "\n[Final Expanded Box 1-based, Top-Left origin]\n";
        for (int i = 0; i < 4; ++i) {
            auto [r, c] = corners_final[i];
            auto [x, y] = PixelCenterToReal_1Based(grid, r, c);
            std::cout << "  " << labels[i] << ": (row=" << r << ", col=" << c << ") -> center(x=" << x << ", y=" << y << ")\n";
        }

        // —— 自检：把一个像素中心反查回像素 —— （可留作调试）
        auto [rx, ry] = PixelCenterToReal_1Based(grid, 1, 1);
        auto rc_back = RealToPixel_1Based(grid, rx, ry);
        // std::cerr << "[SelfCheck] (1,1) -> ("<<rx<<","<<ry<<") -> ("
        //           << rc_back.first << "," << rc_back.second << ")\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] Divide4 main: " << e.what() << std::endl;
        return 1;
    }
}
