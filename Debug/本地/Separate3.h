#ifndef SEPARATE6_H
#define SEPARATE6_H

#include <utility>
#include <array>
#include <vector>
#include <string>
#include <memory>
#pragma once
// SegmentOps.h — 公共声明（实现放在 Separate3.cpp 里，或单独 SegmentOps.cpp 里）

#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <map>

namespace SegOps {

    // ========== 基础结构 ==========
    struct Segment {
        double startX = 0, startY = 0;
        double endX = 0, endY = 0;
        double width = 0;
        std::string layer;
        int net = 0;
    };

    struct SegmentSpan { size_t begin = 0, end = 0; };

    struct MatchParams {
        double coordTol = 0.05;
        bool   requireLayerEqual = true;
    };

    struct ClusterNode { double x = 0, y = 0; std::vector<int> segs; };

    // ========== 读取/解析 ==========
    /** 从 .kicad_pcb 文本中抽取全部 (segment ...) 片段，解析为 Segment，并记录其在原文中的 [begin,end) 范围。
     *  @param pcb   路径
     *  @param outS  输出：所有段
     *  @param outSp 输出：每段在原文中的起止区间
     *  @param raw   可选：返回整个文件原文
     */
    void LoadSegmentsWithSpans(const std::string& pcb,
        std::vector<Segment>& outS,
        std::vector<SegmentSpan>& outSp,
        std::string* raw = nullptr);

    // ========== 匹配/投票 ==========
    /** 统计小段两端点在原始段集合上的“命中”，并给出推荐的 net/width（按多数投票与条件约束）。*/
    void pickOrigNetAndWidthForSmallSegment(const Segment& s,
        const std::vector<Segment>& orig,
        const MatchParams& mp,
        int& net, double& w,
        std::vector<std::pair<int, double>>* outHits = nullptr);

    // ========== 端点聚类成路径 ==========
    /** 以共享端点（容差 tol）为边，构建邻接表（按 layer 分组后构图）。*/
    void buildAdjacencyBySharedEndpoints(const std::vector<Segment>& s,
        double tol,
        std::vector<std::vector<int>>& adj);

    /** 求连通分量（每个分量代表一条（或一簇）路径）。*/
    std::vector<std::vector<int>> connectedComponents(const std::vector<std::vector<int>>& adj);

    // ========== 文本重写/导出 ==========
    /** 按 sp 的切片顺序，将 raw 中对应 (segment ...) 片段替换为 newNet/newW（若你需要真正改写 net/width，可在实现里调用 replaceNum）。*/
    std::string rewriteSmall(const std::string& raw,
        const std::vector<SegmentSpan>& sp,
        const std::vector<int>& newNet,
        const std::vector<double>& newW);

    /** 把对比结果导出成 CSV。*/
    void SaveCSV(const std::string& path,
        const std::vector<Segment>& s,
        const std::vector<int>& n,
        const std::vector<double>& w);

} // namespace SegOps
//===========================================================================================================


// ========== 数据结构（供外部复用） ==========
struct RealRectMM { double x_left, y_top, x_right, y_bot; };

struct BBoxPx1 {
    int row_min, col_min; // top-left (含)
    int row_max, col_max; // bottom-right (含)
};

struct ClipReport {
    bool left = false, right = false, top = false, bottom = false;
    int  left_deficit_px = 0, right_deficit_px = 0, top_deficit_px = 0, bottom_deficit_px = 0;
};

struct PointWithNet {
    double x;
    double y;
    std::string net;
    PointWithNet(double x_val, double y_val, const std::string& net_val)
        : x(x_val), y(y_val), net(net_val) {
    }
};


// 只需要 Grid / Node 等类型声明，包含 grid.h 即可
#include "grid.h"

// —— 角点像素（锚点）与角落紫色块 ——
// 返回四个角的“内侧锚点像素”（顺序：LT, RT, RB, LB），像素坐标 1-based。
std::array<std::pair<int, int>, 4>
GetExpandedPCBCornerAnchorPixels(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels = 5);

// 返回四个角各自 N×N 的“紫色块”像素，并按 LT→RT→RB→LB 顺序拼成一个向量。
// N>=1，奇数更好（3、5…）。超界会自动裁剪。
std::vector<std::pair<int, int>>
GetCornerPixelsPx1(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels = 5,
    int patchN = 1);





std::vector<std::pair<int, int>>
GetFinalCornerTopPixelsInsidePx1(
    const Grid& grid,
    const std::vector<std::pair<double, double>>& finalCornersLT_RT_RB_LB);

std::vector<std::pair<int, int>>
GetFinalCornerTopPixelsInsidePx1_Auto(
    const Grid& grid,
    int layerId, int blockedAreaId,
    int expandPixels, int kWholePixels,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY);

const std::vector<std::pair<int, int>>&
GetLastFinalCornerTopPixelsInsidePx1();












// ========== 坐标换算 ==========
std::pair<double, double> PixelCenterPx1ToRealMM(const Grid& grid, int row1, int col1);
std::pair<int, int>       RealMMToPixelPx1(const Grid& grid, double x_mm, double y_mm);

// ========== 像素框/矩形 ==========
RealRectMM PixelBoxEdgesPx1ToRealRect(const Grid& grid,
    int row_min, int col_min,
    int row_max, int col_max);

std::vector<std::pair<int, int>>
GetBlockedAreaPixels_1Based(const Grid& grid, int layer, int blockedAreaID);

BBoxPx1 ComputeBBoxFromPixels1(const std::vector<std::pair<int, int>>& pxs1);

std::array<std::pair<int, int>, 4>
BBoxCorners_LB_RB_RT_LT(const BBoxPx1& b);

std::array<std::pair<double, double>, 4>
RealCorners_LT_LB_RT_RB(const Grid& grid, const BBoxPx1& b);

BBoxPx1 ExpandBBoxByMPixelsCompensated(const BBoxPx1& b, int m,
    const Grid& g, ClipReport* rep = nullptr);

// ========== 层名/层号 ==========
std::string getLayerNameById(const Grid& grid, int layerId);
int         getLayerIdByName(const Grid& grid, const std::string& layerName);

// ========== 小PCB边界/顶点 ==========
std::vector<std::pair<double, double>>
GetExpandedPCBRealCorners(const Grid& grid, int layer, int blockedAreaID,
    int expandPixels = 5);





// —— 新增：矩形/角点互转与对齐扩张、板框裁剪 ——

// 角点(LT,RT,RB,LB) -> 矩形
RealRectMM CornersToRect_LT_RT_RB_LB(const std::vector<std::pair<double, double>>& c);
// 矩形 -> 角点(LT,RT,RB,LB)
std::vector<std::pair<double, double>> RectToCorners_LT_RT_RB_LB(const RealRectMM& r);

// 现实域：在四边各扩 (0.5 + k)*pixel_len_mm
RealRectMM ExpandRect_HalfPlusK(const RealRectMM& in, const Grid& grid, int kWholePixels);

// 现实域：按板框外接矩形做裁剪
RealRectMM ClipRectToBoard(const RealRectMM& in,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY);

// 方案A主函数：像素域扩张 -> 半像素+整像素对齐 -> 板框裁剪 -> 返回角点(LT,RT,RB,LB)
std::vector<std::pair<double, double>>
GetExpandedPCBCornersSnapAndClip(const Grid& grid, int layer, int blockedAreaID,
    int expandPixels, int kWholePixels,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY);



// 示例函数（可选）
void exampleUsage(const Grid& grid);

// ========== 线段/几何 ==========
std::vector<std::pair<double, double>>
ComputeLineRectIntersections(double x1, double y1, double x2, double y2, const RealRectMM& rect);

bool IsPointInsideRect(double x, double y, const RealRectMM& rect);

bool GetSegmentInfo(const std::shared_ptr<Node>& segment,
    double& x1, double& y1, double& x2, double& y2,
    std::string& layer, std::string& net);

bool IsConnectionPoint(double x, double y,
    const std::vector<std::shared_ptr<Node>>& segments,
    const std::string& targetLayer, double tolerance = 0.001);

// ========== 交点/终止点收集 ==========
std::vector<PointWithNet>
FindIntersectionsAndEndpointsInSmallPCB(const Grid& grid, int layerId, int blockedAreaId,
    int expandPixels = 5);

// ========== 输出模板 ==========
void GenerateCompletePCBTemplate(
    const std::vector<PointWithNet>& modulePoints,
    const std::vector<std::pair<double, double>>& expandedCorners,
    const std::string& outputFilename);

#endif // SEPARATE6_H
#pragma once



// 收集：返回在“第一个矩形”内（含边界）的所有 segment SEID（去重并升序）
std::vector<int>
CollectSegmentsTouchingFirstRect(const Grid& grid, int layerId, int blockedAreaId, int expandPixels = 5);

// 删除：基于上面的收集，调用 KiCadParser::removeSegmentById 批量删除，返回删除条数
size_t
DeleteSegmentsTouchingFirstRect(KiCadParser& parser, const Grid& grid, int layerId, int blockedAreaId, int expandPixels = 5);


