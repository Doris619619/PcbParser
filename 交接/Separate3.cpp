// Separate2.cpp
// ����Լ�������� 1-based������Ϊ(1,1)����������������������
// ��ʵ���꣺mm��y ������������Ļһ�£���

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

#include "First_Part12.h"   // KiCadParser ����
#include <algorithm>
#include <cmath>
#include <stdexcept>


#include "blocked_area_analyzer.h"   // �ں� grid.h



#include <unordered_set>

// ���� С���ߣ��� [[SEID:x]] ������� x��ʧ�ܷ��� -1
static int ExtractSEID(const std::shared_ptr<Node>& segment) {
    if (!segment || segment->name != "segment" || segment->parameters.empty())
        return -1;
    const std::string& p0 = segment->parameters.front(); // ���� "[[SEID:x]]"
    const std::string prefix = "[[SEID:";
    const std::string suffix = "]]";
    if (p0.rfind(prefix, 0) != 0) return -1;
    size_t colon = p0.find(':');
    size_t end = p0.rfind(suffix);
    if (colon == std::string::npos || end == std::string::npos || end <= colon + 1) return -1;
    try {
        return std::stoi(p0.substr(colon + 1, end - colon - 1));
    }
    catch (...) { return -1; }
}


//--------------------------------------------ɾ��end����start�ھ���1����������߶�-------------------------------------
// ���� �ռ� SEID��ֻ������һ�����Ρ������������� �� ��ʵ�������أ������߽�
std::vector<int>
CollectSegmentsTouchingFirstRect(const Grid& grid, int layerId, int blockedAreaId, int expandPixels) {
    std::vector<int> out;

    // ��һ�����Σ����һ����������չ����õ�����ʵ���ؾ���
    // ���к��� GetExpandedPCBRealCorners(layerId, blockedAreaId, expandPixels)
    // ����˳�� LT, RT, RB, LB������ת�ɾ��κ���
    auto corners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
    if (corners.size() != 4) return out;
    RealRectMM firstRect = {
        /*x_left*/  corners[0].first,
        /*y_top*/   corners[0].second,
        /*x_right*/ corners[1].first,
        /*y_bot*/   corners[2].second
    };

    const std::string targetLayer = getLayerNameById(grid, layerId);

    std::unordered_set<int> uniq;//����һ�������set�����洢segment����ţ�set��������uniq
    uniq.reserve(grid.segments.size());

    for (const auto& seg : grid.segments) {
        double x1, y1, x2, y2;
        std::string layer, net;
        if (!GetSegmentInfo(seg, x1, y1, x2, y2, layer, net)) continue;
        if (layer != targetLayer) continue;

        // �����������ϡ� �� IsPointInsideRect() ������Ǳ������жϣ�>= && <=��
        if (IsPointInsideRect(x1, y1, firstRect) || IsPointInsideRect(x2, y2, firstRect)) {
            int id = ExtractSEID(seg);
            if (id > 0) uniq.insert(id);//��ID����unordered_set�������Ĵ洢��
        }
    }

    out.assign(uniq.begin(), uniq.end());// ��set���Ƶ�vector������ת�ƣ�
    std::sort(out.begin(), out.end());
    return out;
}

// ���� ����ɾ������ First_Part12.cpp ���ɾ���ӿ�
size_t
DeleteSegmentsTouchingFirstRect(KiCadParser& parser, const Grid& grid, int layerId, int blockedAreaId, int expandPixels) {
    auto ids = CollectSegmentsTouchingFirstRect(grid, layerId, blockedAreaId, expandPixels);
    size_t removedTotal = 0;
    for (int id : ids) {
        // KiCadParser::removeSegmentById �����ڲ��ؽ����棨�� First_Part12.cpp��
        removedTotal += parser.removeSegmentById(id);
    }
    return removedTotal;
}

//--------------------------------------------ɾ��end����start�ھ���1����������߶�-------------------------------------








// ========== �������� <-> ��ʵ���� ==========
std::pair<double, double>
PixelCenterPx1ToRealMM(const Grid& grid, int row1, int col1)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale ����Ϊ��");
    row1 = std::max(1, std::min(row1, grid.height));
    col1 = std::max(1, std::min(col1, grid.width));
    const double s = (double)grid.inputScale;
    const double x = grid.min_x + ((double)col1 - 0.5) / s;
    const double y = grid.min_y + ((double)row1 - 0.5) / s; // ��������
    return { x, y };
}

std::pair<int, int>
RealMMToPixelPx1(const Grid& grid, double x_mm, double y_mm)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale ����Ϊ��");
    const double s = (double)grid.inputScale;
    int col1 = (int)std::floor((x_mm - grid.min_x) * s) + 1;
    int row1 = (int)std::floor((y_mm - grid.min_y) * s) + 1;
    col1 = std::max(1, std::min(col1, grid.width));
    row1 = std::max(1, std::min(row1, grid.height));
    return { row1, col1 };
}
// ========== ���ر߽� -> ��ʵ���ؾ��� ==========
//�����߽��������꣬���������ʵ����
RealRectMM PixelBoxEdgesPx1ToRealRect(const Grid& grid,
    int row_min, int col_min,
    int row_max, int col_max)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale ����Ϊ��");
    const double s = (double)grid.inputScale;
    RealRectMM rr;
    rr.x_left = grid.min_x + ((double)col_min - 1.0) / s;
    rr.x_right = grid.min_x + ((double)col_max) / s;
    rr.y_top = grid.min_y + ((double)row_min - 1.0) / s;
    rr.y_bot = grid.min_y + ((double)row_max) / s;
    return rr;
}

// ========== ��ȡ�������أ�1-based�� ==========
std::vector<std::pair<int, int>>
GetBlockedAreaPixels_1Based(const Grid& grid, int layer, int blockedAreaID)
{
    if (layer < 0 || layer >= grid.Layers)
        throw std::runtime_error("��ų�����Χ");

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


// ========== �������ر߽��1-based�� ==========
BBoxPx1 ComputeBBoxFromPixels1(const std::vector<std::pair<int, int>>& pxs1)
{
    if (pxs1.empty()) throw std::runtime_error("ComputeBBoxFromPixels1: ��������Ϊ��");
    int rmin = INT_MAX, cmin = INT_MAX, rmax = INT_MIN, cmax = INT_MIN;
    for (auto [r1, c1] : pxs1) {
        rmin = std::min(rmin, r1);
        cmin = std::min(cmin, c1);
        rmax = std::max(rmax, r1);
        cmax = std::max(cmax, c1);
    }
    return { rmin, cmin, rmax, cmax };
}

// �����Ľǣ�LB, RB, RT, LT��
std::array<std::pair<int, int>, 4> BBoxCorners_LB_RB_RT_LT(const BBoxPx1& b)
{
    return {
        std::pair<int,int>{ b.row_max, b.col_min }, // LB
        std::pair<int,int>{ b.row_max, b.col_max }, // RB
        std::pair<int,int>{ b.row_min, b.col_max }, // RT
        std::pair<int,int>{ b.row_min, b.col_min }  // LT
    };
}





// ��ʵ�Ľǣ��ϸ� LT, LB, RT, RB�������а�/���Ʊ߿�
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

//// ��ʵ�Ľǣ�ʹ�������������꣬�ϸ� LT, LB, RT, RB�������а�/���Ʊ߿�
//std::array<std::pair<double, double>, 4>
//RealCorners_LT_LB_RT_RB(const Grid& grid, const BBoxPx1& b)
//{
//    // ֱ��ʹ�õ�һ�����������ĸ��ǵ�������������
//    return {
//        PixelCenterPx1ToRealMM(grid, b.row_min, b.col_min), // ���Ͻ� (LT)
//        PixelCenterPx1ToRealMM(grid, b.row_max, b.col_min), // ���½� (LB)
//        PixelCenterPx1ToRealMM(grid, b.row_min, b.col_max), // ���Ͻ� (RT)
//        PixelCenterPx1ToRealMM(grid, b.row_max, b.col_max)  // ���½� (RB)
//    };
//}



// �߼���ϣ�����Ҹ��� m�����¸��� m����ĳ�಻��������ȱ�Ĳ����Բࡣ
 BBoxPx1 ExpandBBoxByMPixelsCompensated(const BBoxPx1& b, int m,
    const Grid& g, ClipReport* rep)
{
    const int W = g.width, H = g.height;

    // �����򵽱߽�ġ��������ռ�
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

    // �Բಹ����ϣ���������� �� 2m
    const int extend_left = m + right_def; // �Ҳ಻������ȱ�Ĳ������
    const int extend_right = m + left_def;  // ��಻������ȱ�Ĳ����Ҳ�
    const int extend_top = m + bot_def;   // �²಻���������ϲ�
    const int extend_bottom = m + top_def;   // �ϲ಻���������²�

    BBoxPx1 e;
    e.col_min = std::max(1, b.col_min - extend_left);
    e.col_max = std::min(W, b.col_max + extend_right);
    e.row_min = std::max(1, b.row_min - extend_top);
    e.row_max = std::min(H, b.row_max + extend_bottom);

    // ����ת������խ��+���� m��
    if (e.col_min > e.col_max) { e.col_min = 1; e.col_max = W; }
    if (e.row_min > e.row_max) { e.row_min = 1; e.row_max = H; }

    return e;
}


// ���ݲ��Ż�ȡ������
std::string getLayerNameById(const Grid& grid, int layerId) {
    if (layerId >= 0 && layerId < grid.layerName.size()) {
        return grid.layerName[layerId]; // ���ض�Ӧ���ŵĲ�����
    }
    else {
        return "��Ч�Ĳ���"; // ������ų�����Χ�����ش�����Ϣ
    }
}


// ���ݲ����ƻ�ȡ����
int getLayerIdByName(const Grid& grid, const std::string& layerName) {
    for (int i = 0; i < grid.layerName.size(); i++) {
        if (grid.layerName[i] == layerName) {
            return i; // ���ض�Ӧ�Ĳ���
        }
    }
    return -1; // ����Ҳ���������-1
}






static inline RealRectMM ExpandRect_HalfPlusK(const RealRectMM& in, const Grid& grid, int kWholePixels) {
    if (grid.inputScale <= 0) throw std::runtime_error("inputScale ����Ϊ��");
    const double px_mm = 1.0 / static_cast<double>(grid.inputScale); // 1px = 1/inputScale mm
    const double d = (static_cast<double>(kWholePixels)) * px_mm;
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
    // ���˷���ת
    if (r.x_left > r.x_right) { r.x_left = boardMinX; r.x_right = boardMaxX; }
    if (r.y_top > r.y_bot) { r.y_top = boardMinY; r.y_bot = boardMaxY; }
    return r;
}

// ����A������������������������������� -> ��ʵ�� 0.5px+kpx -> �ü���� -> ���� LT,RT,RB,LB
std::vector<std::pair<double, double>>
GetExpandedPCBCornersSnapAndClip(const Grid& grid, int layer, int blockedAreaID,
    int expandPixels, int kWholePixels,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY)
{
    // 1) �ȵõ������������ź����ʵ�Ľǡ���LT,RT,RB,LB��
    auto baseCorners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);
    if (baseCorners.size() != 4) {
        throw std::runtime_error("GetExpandedPCBCornersSnapAndClip: �ǵ��� != 4");
    }

    // 2) ��ʵ��Գ����� 0.5px + k��1px
    RealRectMM baseRect = CornersToRect_LT_RT_RB_LB(baseCorners);
    RealRectMM grown = ExpandRect_HalfPlusK(baseRect, grid, kWholePixels);

    // 3) �ü�����ʵ���
    RealRectMM clipped = ClipRectToBoard(grown, boardMinX, boardMinY, boardMaxX, boardMaxY);

    // 4) ���Ľǣ�LT,RT,RB,LB��
    return RectToCorners_LT_RT_RB_LB(clipped);
}


// ========== ��ȡ��չ��СPCB���ĸ���ʵ�ǵ����� ==========
std::vector<std::pair<double, double>>
GetExpandedPCBRealCorners(const Grid& grid, int layer, int blockedAreaID, int expandPixels)
{    



    // 1. ��ȡ������������
    auto blockedPixels = GetBlockedAreaPixels_1Based(grid, layer, blockedAreaID);
    if (blockedPixels.empty()) {
        throw std::runtime_error("δ�ҵ�ָ������������");
    }

    // 2. ����ԭʼ��Χ��
    BBoxPx1 originalBBox = ComputeBBoxFromPixels1(blockedPixels);

    // 3. ��չ��Χ��
    ClipReport clipReport;
    BBoxPx1 expandedBBox = ExpandBBoxByMPixelsCompensated(originalBBox, expandPixels, grid, &clipReport);

    // 4. ��ȡ��ʵ������ĸ��ǵ㣨˳ʱ��˳��LT, RT, RB, LB��
    auto realCorners = RealCorners_LT_LB_RT_RB(grid, expandedBBox);

    // 5. ��������Ϊ˳ʱ��˳��LT -> RT -> RB -> LB
    std::vector<std::pair<double, double>> result;
    result.reserve(4);

    // LT (����)
    result.push_back(realCorners[0]);
    // RT (����)  
    result.push_back(realCorners[2]);
    // RB (����)
    result.push_back(realCorners[3]);
    // LB (����)
    result.push_back(realCorners[1]);

    return result;
}

// ʹ��ʾ����
void exampleUsage(const Grid& grid) {
    try {
        int layer = 0;          // ���
        int blockedAreaID = 1;  // ��������ID
        int expandPixels = 5;   // ��չ������

        auto corners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);

        // ��ӡ�ĸ��ǵ�����
        std::cout << "СPCB���ĸ��ǵ����꣨˳ʱ��˳�򣩣�" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        for (size_t i = 0; i < corners.size(); ++i) {
            std::cout << "�ǵ� " << i + 1 << ": ("
                << corners[i].first << " mm, "
                << corners[i].second << " mm)" << std::endl;
        }

        // �����Ҫֱ��ʹ�����vector
        // corners ���ڰ������ĸ�������꣬��˳ʱ��˳��LT->RT->RB->LB

    }
    catch (const std::exception& e) {
        std::cerr << "����: " << e.what() << std::endl;
    }
}



// ========== �����߶�����α߽�Ľ��� ==========
std::vector<std::pair<double, double>>
ComputeLineRectIntersections(double x1, double y1, double x2, double y2, const RealRectMM& rect)
{
    std::vector<std::pair<double, double>> intersections;

    // ����������ߵĽ���
    // ��߽� (x = rect.x_left)
    if (x1 != x2) {
        double t = (rect.x_left - x1) / (x2 - x1);
        if (t >= 0 && t <= 1) {
            double y = y1 + t * (y2 - y1);
            if (y >= rect.y_top && y <= rect.y_bot) {
                intersections.emplace_back(rect.x_left, y);
            }
        }
    }

    // �ұ߽� (x = rect.x_right)
    if (x1 != x2) {
        double t = (rect.x_right - x1) / (x2 - x1);
        if (t >= 0 && t <= 1) {
            double y = y1 + t * (y2 - y1);
            if (y >= rect.y_top && y <= rect.y_bot) {
                intersections.emplace_back(rect.x_right, y);
            }
        }
    }

    // �ϱ߽� (y = rect.y_top)
    if (y1 != y2) {
        double t = (rect.y_top - y1) / (y2 - y1);
        if (t >= 0 && t <= 1) {
            double x = x1 + t * (x2 - x1);
            if (x >= rect.x_left && x <= rect.x_right) {
                intersections.emplace_back(x, rect.y_top);
            }
        }
    }

    // �±߽� (y = rect.y_bot)
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

// ========== �жϵ��Ƿ��ھ����ڲ� ==========
bool IsPointInsideRect(double x, double y, const RealRectMM& rect)
{
    return (x >= rect.x_left && x <= rect.x_right &&
        y >= rect.y_top && y <= rect.y_bot);
}

// ========== �޸ĺ�Ļ�ȡ�߶���Ϣ���� ==========
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

// ========== �޸ĺ���ж����ӵ㺯�� ==========
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

        // ����Ƿ���Ŀ���
        if (layer != targetLayer) {
            continue;
        }

        // �����Ƿ����߶����ƥ��
        if (std::abs(x1 - x) < tolerance && std::abs(y1 - y) < tolerance) {
            connectionCount++;
        }

        // �����Ƿ����߶��յ�ƥ��
        if (std::abs(x2 - x) < tolerance && std::abs(y2 - y) < tolerance) {
            connectionCount++;
        }

        // �������������1��˵�������ӵ�
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
        // 1. ��ȡСPCB�ı߽��
        auto corners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
        RealRectMM pcbBox;
        pcbBox.x_left = corners[0].first;
        pcbBox.y_top = corners[0].second;
        pcbBox.x_right = corners[1].first;
        pcbBox.y_bot = corners[2].second;

        std::cout << "СPCB�߽��: ��=" << pcbBox.x_left << ", ��=" << pcbBox.y_top
            << ", ��=" << pcbBox.x_right << ", ��=" << pcbBox.y_bot << std::endl;

        // 2. ��ȡ��ǰ�������
        std::string targetLayerName = getLayerNameById(grid, layerId);
        std::cout << "Ŀ���: " << targetLayerName << std::endl;

        // 3. ֱ��ʹ�� grid.segments
        const auto& segments = grid.segments;
        std::cout << "�ҵ� " << segments.size() << " ���߶�" << std::endl;

        // 4. ���������߶�
        int processedCount = 0;
        for (const auto& segment : segments) {
            double x1, y1, x2, y2;
            std::string layer, net;

            if (!GetSegmentInfo(segment, x1, y1, x2, y2, layer, net)) {
                continue;
            }

            // ����Ƿ���Ŀ���
            if (layer != targetLayerName) {
                continue;
            }

            processedCount++;
            std::cout << "�����߶� " << processedCount << ": (" << x1 << "," << y1 << ") -> (" << x2 << "," << y2
                << "), ��: " << layer << ", ����: " << net << std::endl;

            // 5. �����߶���СPCB�߽�Ľ���
            auto intersections = ComputeLineRectIntersections(x1, y1, x2, y2, pcbBox);

            // 6. ����߶ζ˵��Ƿ���СPCB�ڲ�
            bool startInside = IsPointInsideRect(x1, y1, pcbBox);
            bool endInside = IsPointInsideRect(x2, y2, pcbBox);

            std::cout << "  ������ڲ�: " << (startInside ? "��" : "��")
                << ", �յ����ڲ�: " << (endInside ? "��" : "��")
                << ", ��������: " << intersections.size() << std::endl;

            // 7. �ռ����㣨������Ƿ�Ϊ���ӵ㣬��Ϊ�����ڱ߽��ϣ�
            for (const auto& intersection : intersections) {
                // ����ͨ���ڱ߽��ϣ�����Ҫ����Ƿ�Ϊ���ӵ�
                result.emplace_back(intersection.first, intersection.second, net);
                std::cout << "  ����: (" << intersection.first << ", " << intersection.second
                    << "), ����: " << net << std::endl;
            }

            // 8. �����ֹ�㣬�ų����ӵ�
            if (startInside && !endInside) {
                // ������ڲ����յ����ⲿ �� �������ֹ��
                // �������Ƿ�Ϊ���ӵ�
                if (!IsConnectionPoint(x1, y1, segments, targetLayerName)) {
                    result.emplace_back(x1, y1, net);
                    std::cout << "  ��ֹ��(������ڲ�): (" << x1 << ", " << y1
                        << "), ����: " << net << std::endl;
                }
                else {
                    std::cout << "  �ų����ӵ�(���): (" << x1 << ", " << y1
                        << "), ����: " << net << std::endl;
                }
            }
            else if (!startInside && endInside) {
                // ������ⲿ���յ����ڲ� �� �յ�����ֹ��
                // ����յ��Ƿ�Ϊ���ӵ�
                if (!IsConnectionPoint(x2, y2, segments, targetLayerName)) {
                    result.emplace_back(x2, y2, net);
                    std::cout << "  ��ֹ��(�յ����ڲ�): (" << x2 << ", " << y2
                        << "), ����: " << net << std::endl;
                }
                else {
                    std::cout << "  �ų����ӵ�(�յ�): (" << x2 << ", " << y2
                        << "), ����: " << net << std::endl;
                }
            }
            else if (startInside && endInside) {
                // �߶���ȫ���ڲ��������˵㶼��СPCB��
                // ��������˵��Ƿ�Ϊ���ӵ�
                if (!IsConnectionPoint(x1, y1, segments, targetLayerName)) {
                    result.emplace_back(x1, y1, net);
                    std::cout << "  ��ֹ��(������ڲ�): (" << x1 << ", " << y1
                        << "), ����: " << net << std::endl;
                }
                else {
                    std::cout << "  �ų����ӵ�(���): (" << x1 << ", " << y1
                        << "), ����: " << net << std::endl;
                }

                if (!IsConnectionPoint(x2, y2, segments, targetLayerName)) {
                    result.emplace_back(x2, y2, net);
                    std::cout << "  ��ֹ��(�յ����ڲ�): (" << x2 << ", " << y2
                        << "), ����: " << net << std::endl;
                }
                else {
                    std::cout << "  �ų����ӵ�(�յ�): (" << x2 << ", " << y2
                        << "), ����: " << net << std::endl;
                }
            }
        }

        std::cout << "�������� " << processedCount << " ����Ŀ����ϵ��߶�" << std::endl;

        // 9. ȥ�أ�ȷ�������û���ظ��㣩
        auto last = std::unique(result.begin(), result.end(),
            [](const PointWithNet& a, const PointWithNet& b) {
                return std::abs(a.x - b.x) < 0.001 &&
                    std::abs(a.y - b.y) < 0.001 &&
                    a.net == b.net;
            });
        result.erase(last, result.end());

    }
    catch (const std::exception& e) {
        std::cerr << "��FindIntersectionsAndEndpointsInSmallPCB�г���: " << e.what() << std::endl;
    }

    return result;
}


// === ���ش�С��mm/px�� ===
static inline double PixelSizeMM(const Grid& g) {
    if (g.inputScale <= 0) throw std::runtime_error("grid.inputScale ����Ϊ��");
    return 1.0 / static_cast<double>(g.inputScale);
}

// === ���Ľǣ�LT,RT,RB,LB����װ���� ===
static inline RealRectMM CornersToRect_LT_RT_RB_LB(const std::vector<std::pair<double, double>>& c) {
    if (c.size() != 4) throw std::runtime_error("CornersToRect_LT_RT_RB_LB: ��Ҫ4���ǵ�");
    RealRectMM r;
    // Լ����c[0]=LT, c[1]=RT, c[2]=RB, c[3]=LB
    r.x_left = c[0].first;
    r.y_top = c[0].second;
    r.x_right = c[1].first;
    r.y_bot = c[2].second;
    return r;
}

// === �Ѿ���ת���Ľǣ�LT,RT,RB,LB�� ===
static inline std::vector<std::pair<double, double>> RectToCorners_LT_RT_RB_LB(const RealRectMM& r) {
    return {
        { r.x_left , r.y_top },   // LT
        { r.x_right, r.y_top },   // RT
        { r.x_right, r.y_bot },   // RB
        { r.x_left , r.y_bot }    // LB
    };
}






// ���� �������ѡ��ڱ߽��ϵ���ʵ���ꡱ�ȶ���ӳ�䵽С���ڲ������� ����
// ע�⣺�ұ߽�/�±߽�Ҫ������һ����С����������䵽������ء�
static inline std::pair<int, int> RealToPixel_Inside_LT(const Grid& g, double x, double y) {
    return RealMMToPixelPx1(g, x, y); // ��/�ϱ߽�ֱ�� floor �����ڲ�
}
static inline std::pair<int, int> RealToPixel_Inside_RT(const Grid& g, double x, double y) {
    const double eps = 1.0 / (g.inputScale * 1024.0); // < 1px �ĺ�С��
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

// ���� ����1���ĸ��ǵġ�ê�����ء���LT,RT,RB,LB�� ����
// ������ʵ�ǵ�(LT,RT,RB,LB)���������ڲ�ƫ�á�������ӳ�䡣
std::array<std::pair<int, int>, 4>
GetExpandedPCBCornerAnchorPixels(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels)
{
    auto corners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);
    if (corners.size() != 4) throw std::runtime_error("�ǵ�����ӦΪ4");

    // corners[0]=LT, [1]=RT, [2]=RB, [3]=LB ����ǰʵ�ַ���˳����ˣ�
    auto lt = RealToPixel_Inside_LT(grid, corners[0].first, corners[0].second);
    auto rt = RealToPixel_Inside_RT(grid, corners[1].first, corners[1].second);
    auto rb = RealToPixel_Inside_RB(grid, corners[2].first, corners[2].second);
    auto lb = RealToPixel_Inside_LB(grid, corners[3].first, corners[3].second);
    return { lt, rt, rb, lb };
}

// ���� ���������ݡ��ǵķ������� N��N ����ɫ�飬���� 1..H/1..W �ü� ����
// dr_sign: �з���(+1����, -1����)��dc_sign: �з���(+1����, -1����)
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

// ���� ����2���ĸ��ǵ���ɫ�����أ��� LT��RT��RB��LB ˳��ƴ�ӵ�һ�� vector ����
// ����patchN=3 ����ÿ���� 3��3����ֻҪһ�����أ��� patchN=1 ���ɡ�
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

    // LT�����¡��ҡ��� N��N
    AppendCornerPatch(out, anchors[0].first, anchors[0].second, patchN, +1, +1, grid);
    // RT�����¡�����
    AppendCornerPatch(out, anchors[1].first, anchors[1].second, patchN, +1, -1, grid);
    // RB�����ϡ�����
    AppendCornerPatch(out, anchors[2].first, anchors[2].second, patchN, -1, -1, grid);
    // LB�����ϡ��ҡ���
    AppendCornerPatch(out, anchors[3].first, anchors[3].second, patchN, -1, +1, grid);

    return out;
}

// ========== ��������PCBģ���ļ� ==========
void GenerateCompletePCBTemplate(
    const std::vector<PointWithNet>& modulePoints,
    const std::vector<std::pair<double, double>>& expandedCorners,
    const std::string& outputFilename)
{
    // �ĸ�module��ģ��
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

 
  (gr_line (start 160.5011 75.003598) (end 98.501099 75.003598) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 75.003598) (end 98.501099 120.003599) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 120.003599) (end 160.5011 120.003599) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 160.5011 120.003599) (end 160.5011 75.003598) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
))";

    // ����module��ģ��
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

    // ����module�����滻
    std::vector<PointWithNet> replacementPoints;

    if (modulePoints.size() == 2) {
        // ���1: ֻ��2���㣬ֱ��ʹ�����������滻ǰ����module
        replacementPoints = modulePoints;
        std::cout << "ʹ��2�����滻ǰ����module" << std::endl;
    }
    else if (modulePoints.size() >= 4) {
        // ���2: ��4�������㣬��������鴦��
        // �ҵ���������ͬ����ĵ�
        std::map<std::string, std::vector<PointWithNet>> pointsByNet;
        for (const auto& point : modulePoints) {
            pointsByNet[point.net].push_back(point);
        }

        bool foundMatchingNet = false;
        for (const auto& [net, netPoints] : pointsByNet) {
            if (netPoints.size() >= 2) {
                // ʹ��ǰ������ͬ����ĵ�
                replacementPoints.push_back(netPoints[0]);
                replacementPoints.push_back(netPoints[1]);
                foundMatchingNet = true;
                std::cout << "�ҵ����� '" << net << "' ��2��������ǰ����module" << std::endl;
                break;
            }
        }

        // ���û���ҵ���ͬ����ĵ㣬ʹ��ǰ4����
        if (!foundMatchingNet) {
            for (size_t i = 0; i < std::min(modulePoints.size(), size_t(4)); i++) {
                replacementPoints.push_back(modulePoints[i]);
            }
            std::cout << "û���ҵ���ͬ����ĵ㣬ʹ��ǰ4����" << std::endl;
        }

        // ���ʣ��ĵ�
        size_t addedCount = replacementPoints.size();
        for (const auto& point : modulePoints) {
            if (addedCount >= 4) break;

            // ����Ƿ��Ѿ���ӹ�
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
        // �������: ʹ�����е�
        replacementPoints = modulePoints;
        std::cout << "ʹ������ " << modulePoints.size() << " ����" << std::endl;
    }

    // ѡ��ģ��
    std::string template_str;
    if (replacementPoints.size() == 2) {
        template_str = template_two_modules;
        std::cout << "ʹ������module��ģ��" << std::endl;
    }
    else {
        template_str = template_four_modules;
        std::cout << "ʹ���ĸ�module��ģ��" << std::endl;
    }

    // ����ģ���滻
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);

    // �ָ�ģ�岢���д���
    std::istringstream templateStream(template_str);
    std::string line;
    int moduleCount = 0;
    int grLineCount = 0;

    while (std::getline(templateStream, line)) {
        // ����module�����滻
        if (line.find("(module") != std::string::npos) {
            // ���module��
            oss << line << "\n";

            // ��ȡ��һ��(Ӧ�ð���at����)
            if (std::getline(templateStream, line)) {
                if (line.find("(at") != std::string::npos) {
                    // �滻����
                    if (moduleCount < replacementPoints.size()) {
                        const auto& point = replacementPoints[moduleCount];
                        oss << "    (at " << point.x << " " << point.y << " 270)\n";
                        std::cout << "�滻module " << (moduleCount + 1) << " ����Ϊ: ("
                            << point.x << ", " << point.y << "), ����: " << point.net << std::endl;
                    }
                    else {
                        // ���û���㹻�ĵ㣬����ԭ����
                        oss << line << "\n";
                    }
                    moduleCount++;
                }
                else {
                    oss << line << "\n";
                }
            }
        }
        // ����gr_line�����滻
        else if (line.find("(gr_line") != std::string::npos && expandedCorners.size() == 4) {
            if (grLineCount == 0) {
                // ��һ����: �ϱ� (LT -> RT)
                oss << "  (gr_line (start " << expandedCorners[0].first << " " << expandedCorners[0].second
                    << ") (end " << expandedCorners[1].first << " " << expandedCorners[1].second
                    << ") (layer Bottom) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 1) {
                // �ڶ�����: �ұ� (RT -> RB)
                oss << "  (gr_line (start " << expandedCorners[1].first << " " << expandedCorners[1].second
                    << ") (end " << expandedCorners[2].first << " " << expandedCorners[2].second
                    << ") (layer Bottom) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 2) {
                // ��������: �±� (RB -> LB)
                oss << "  (gr_line (start " << expandedCorners[2].first << " " << expandedCorners[2].second
                    << ") (end " << expandedCorners[3].first << " " << expandedCorners[3].second
                    << ") (layer Bottom) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 3) {
                // ��������: ��� (LB -> LT)
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

    // д���ļ�
    std::ofstream file(outputFilename);
    if (!file.is_open()) {
        throw std::runtime_error("�޷���������ļ�: " + outputFilename);
    }

    file << oss.str();
    file.close();

    std::cout << "����PCBģ���ļ�������: " << outputFilename << std::endl;
    std::cout << "����� " << moduleCount << " ��module" << std::endl;
}




// ========== ȫ�ֱ����洢�������ɵ��ļ��� ==========
std::string g_latestGeneratedFile = "";





//

int main() {
    try {
        // === 1) ��ʼ�� Grid ===
        std::string filename = "testcase.kicad_pcb";   // �����޸�
        Grid grid;
        grid.SetUp(filename);
        std::cout << "Grid��ʼ���ɹ�! ����ߴ�: " << grid.width << " x " << grid.height << std::endl;
        
        
        std::cout<< std::fixed << std::setprecision(6)<<grid.min_x<<' '<<grid.min_y<<' '<<grid.max_x<<' '<<grid.max_y;


       
        // === 2) ���� ===
        int layerId = 1;   // �����
        int blockedAreaId = 2;   // ����������
        int expandPixels = 15;   // ��������չ���������в�����
        int kWholePixels = 20;   // ��ʵ������ k �������أ��� 0.5px һ��
        double pxSizeMM = 1.0 / static_cast<double>(grid.inputScale);

        std::cout << "����: ��=" << layerId << " (" << getLayerNameById(grid, layerId)
            << "), ��������=" << blockedAreaId
            << ", ��������չ=" << expandPixels << " px"
            << ", ��ʵ���������� = 0.5px + " << kWholePixels << "px"
            << " �� " << std::fixed << std::setprecision(5)
            << (0.5 + kWholePixels) * pxSizeMM << " mm\n";

        // === 3) ��ȡ��ʵ���Edge.Cuts ��Ӿ��Σ� ===
        KiCadParser parser2;
        if (!parser2.parseFile(filename)) {
            throw std::runtime_error("���� KiCad �ļ�ʧ�ܣ��޷���ȡ���");
        }
        KiCadParser::Point2D bl, tl, tr, br; // ���¡����ϡ����ϡ�����
        if (!parser2.getBoardCorners(bl, tl, tr, br)) {
            throw std::runtime_error("δ�ҵ� Edge.Cuts ���");
        }
        const double boardMinX = bl.x;
        const double boardMinY = bl.y;
        const double boardMaxX = tr.x;
        const double boardMaxY = tr.y;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[Board] min=(" << boardMinX << "," << boardMinY
            << "), max=(" << boardMaxX << "," << boardMaxY << ")\n";

        // === 4) ��ӡ��������չ�׶Ρ����������꣨�в����� ===
        // 4.1 �����������أ�1-based��
        auto blockedPixels = GetBlockedAreaPixels_1Based(grid, layerId, blockedAreaId);
        if (blockedPixels.empty()) {
            throw std::runtime_error("δ�ҵ�ָ������������");
        }

        // 4.2 ԭʼ���ذ�Χ�У�1-based��
        BBoxPx1 originalBBox = ComputeBBoxFromPixels1(blockedPixels);
        std::cout << "\n=== ԭʼ���ذ�Χ��(1-based, col,x  row,y) ===\n";
        std::cout << "LT(col,row)=(" << originalBBox.col_min << "," << originalBBox.row_min << ")\n";
        std::cout << "RT(col,row)=(" << originalBBox.col_max << "," << originalBBox.row_min << ")\n";
        std::cout << "RB(col,row)=(" << originalBBox.col_max << "," << originalBBox.row_max << ")\n";
        std::cout << "LB(col,row)=(" << originalBBox.col_min << "," << originalBBox.row_max << ")\n";
        std::cout << "W x H (px) = " << (originalBBox.col_max - originalBBox.col_min + 1)
            << " x " << (originalBBox.row_max - originalBBox.row_min + 1) << "\n";

        // 4.3 ������չ���в�����
        ClipReport clipReport{};
        BBoxPx1 expandedBBoxPx = ExpandBBoxByMPixelsCompensated(originalBBox, expandPixels, grid, &clipReport);

        std::cout << "\n=== ������չ���Χ��(�в���, 1-based) ===\n";
        std::cout << "LT(col,row)=(" << expandedBBoxPx.col_min << "," << expandedBBoxPx.row_min << ")\n";
        std::cout << "RT(col,row)=(" << expandedBBoxPx.col_max << "," << expandedBBoxPx.row_min << ")\n";
        std::cout << "RB(col,row)=(" << expandedBBoxPx.col_max << "," << expandedBBoxPx.row_max << ")\n";
        std::cout << "LB(col,row)=(" << expandedBBoxPx.col_min << "," << expandedBBoxPx.row_max << ")\n";
        std::cout << "W x H (px) = " << (expandedBBoxPx.col_max - expandedBBoxPx.col_min + 1)
            << " x " << (expandedBBoxPx.row_max - expandedBBoxPx.row_min + 1) << "\n";
        std::cout << "����: left=" << clipReport.left
            << ", right=" << clipReport.right
            << ", top=" << clipReport.top
            << ", bottom=" << clipReport.bottom << "\n";

        // === 5) ����ѡ����ӡ����������չ�����ʵ�Ľǡ������ڶԱȣ�
        auto baseCorners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
        std::vector<std::string> names = { "����(LT)","����(RT)","����(RB)","����(LB)" };
        std::cout << "\n=== ��������չ�����ʵ�Ľ�(���ڶԱ�) ===\n";
        for (size_t i = 0; i < baseCorners.size(); ++i) {
            std::cout << names[i] << ": (" << baseCorners[i].first << ", " << baseCorners[i].second << ")\n";
        }

        // === 6) ���գ���ʵ�� 0.5px + k��1px �������� + �ü������ ===
        auto expandedCorners = GetExpandedPCBCornersSnapAndClip(
            grid, layerId, blockedAreaId,
            expandPixels, kWholePixels,
            boardMinX, boardMinY, boardMaxX, boardMaxY
        );

        std::cout << "\n=== ����СPCB�߽��(0.5px + " << kWholePixels
            << "px ���� + ���ü���) ===\n";
        for (size_t i = 0; i < expandedCorners.size(); ++i) {
            std::cout << names[i] << ": (" << expandedCorners[i].first
                << ", " << expandedCorners[i].second << ")\n";
        }


        // === 6.1) �����Ľ� -> ���ض�������(1-based, row,col)����֤����С���ڲ� ===
        auto toPxInside = [&](double x_mm, double y_mm, int cornerIdx) -> std::pair<int, int> {
            // cornerIdx: 0=LT, 1=RT, 2=RB, 3=LB
            const double eps = 1.0 / (grid.inputScale * 1024.0); // ԶС��1px
            double xx = x_mm, yy = y_mm;
            if (cornerIdx == 1 || cornerIdx == 2) xx -= eps; // RT/RB��������һ��
            if (cornerIdx == 2 || cornerIdx == 3) yy -= eps; // RB/LB��������һ��
            return RealMMToPixelPx1(grid, xx, yy); // ����( row , col )��1-based
            };

        std::array<std::pair<int, int>, 4> finalCornerPx;
        for (int i = 0; i < 4; ++i) {
            finalCornerPx[i] = toPxInside(expandedCorners[i].first, expandedCorners[i].second, i);
        }

        std::cout << "\n=== ����СPCB�ĽǶ�Ӧ�����ض�������(1-based, row,col) ===\n";
        const char* nm[4] = { "LT","RT","RB","LB" };
        for (int i = 0; i < 4; ++i) {
            std::cout << nm[i] << ": (" << finalCornerPx[i].first
                << ", " << finalCornerPx[i].second << ")\n";
        }

        // === 7) ����/��ֹ�㣨�԰���������չ�ľ�������������ԭ�߼��� ===
        auto points = FindIntersectionsAndEndpointsInSmallPCB(grid, layerId, blockedAreaId, expandPixels);

        std::cout << "\n=== �߶ν������ֹ�� ===\n";
        std::cout << "���ҵ� " << points.size() << " ����\n";

        // === 8) ����ģ�壨�������Ľǣ� ===
        GenerateCompletePCBTemplate(points, expandedCorners, "complete_pcb_template.kicad_pcb");


        // === [NEW] �ռ���ɾ������һ�����������γɵ���ʵ���ؾ��Ρ���/���ϵ��߶� ===
        {
            // 1) �ռ���ѡ SEID��ע�⣺expandPixels �������һ��������չ����һ�£�
            std::vector<int> candIds = CollectSegmentsTouchingFirstRect(
                grid,           // ��ǰ Grid���� segments / layerName �ȣ�
                layerId,        // ���
                blockedAreaId,  // �����������
                expandPixels    // ��һ�����������õ�������
            );
            std::sort(candIds.begin(), candIds.end());
            candIds.erase(std::unique(candIds.begin(), candIds.end()), candIds.end());

            std::cout << "\n=== ��ɾ���߶Σ���һ������/���ϣ�===\n";
            if (candIds.empty()) {
                std::cout << "��\n";
            }
            else {
                std::cout << "�� " << candIds.size() << " ����";
                for (size_t i = 0; i < candIds.size(); ++i) {
                    if (i) std::cout << ", ";
                    std::cout << candIds[i];
                }
                std::cout << "\n";

                // ���� ��ϸ��ӡ��SEID + ���յ� + �� + net
                auto extractSeid = [](const std::shared_ptr<Node>& seg) -> int {
                    if (!seg || seg->parameters.empty()) return -1;
                    const std::string& p0 = seg->parameters.front(); // ���� [[SEID:123]]
                    const std::string prefix = "[[SEID:";
                    const std::string suffix = "]]";
                    if (p0.rfind(prefix, 0) != 0) return -1;
                    size_t end = p0.rfind(suffix);
                    if (end == std::string::npos) return -1;
                    try { return std::stoi(p0.substr(prefix.size(), end - prefix.size())); }
                    catch (...) { return -1; }
                    };

                std::cout << "���� ��ϸ ����\n";
                for (const auto& seg : grid.segments) {
                    int sid = extractSeid(seg);
                    if (sid <= 0) continue;
                    if (!std::binary_search(candIds.begin(), candIds.end(), sid)) continue;

                    double x1, y1, x2, y2; std::string lyr, net;
                    if (GetSegmentInfo(seg, x1, y1, x2, y2, lyr, net)) {
                        std::cout << "  SEID " << sid << " [" << lyr << "]: "
                            << "(" << x1 << ", " << y1 << ") -> ("
                            << x2 << ", " << y2 << ")  net=" << net << "\n";
                    }
                }
            }

            // 2) ����ɾ����������ļ�
            if (!candIds.empty()) {
                KiCadParser parser;
                if (!parser.parseFile(filename)) {
                    std::cerr << "����ԭʼ Kicad �ļ�ʧ��: " << filename << "\n";
                }
                else {
                    size_t removed = 0;
                    for (int id : candIds) removed += parser.removeSegmentById(id);
                    std::cout << "��ɾ���߶�����: " << removed << "\n";

                    // ����ļ�����<ԭ��ȥ��׺>.firstrect_pruned.kicad_pcb
                    std::string outPath = filename;
                    auto dot = outPath.rfind('.');
                    if (dot != std::string::npos) outPath = outPath.substr(0, dot);
                    outPath += ".firstrect_pruned.kicad_pcb";

                    if (parser.saveAsKicadPcb(outPath, /*stripTempIds=*/true)) {
                        std::cout << "��д��ɾ������ļ�: " << outPath << "\n";
                    }
                    else {
                        std::cerr << "д��ʧ��\n";
                    }
                }
            }
        }


    }
    catch (const std::exception& e) {
        std::cerr << "�������г���: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "����ִ�����!" << std::endl;
    return 0;
}



