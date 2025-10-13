#ifndef SEPARATE6_H
#define SEPARATE6_H

#include <utility>
#include <array>
#include <vector>
#include <string>
#include <memory>

// ֻ��Ҫ Grid / Node ���������������� grid.h ����
#include "grid.h"

// ���� �ǵ����أ�ê�㣩�������ɫ�� ����
// �����ĸ��ǵġ��ڲ�ê�����ء���˳��LT, RT, RB, LB������������ 1-based��
std::array<std::pair<int, int>, 4>
GetExpandedPCBCornerAnchorPixels(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels = 5);

// �����ĸ��Ǹ��� N��N �ġ���ɫ�顱���أ����� LT��RT��RB��LB ˳��ƴ��һ��������
// N>=1���������ã�3��5������������Զ��ü���
std::vector<std::pair<int, int>>
GetCornerPixelsPx1(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels = 5,
    int patchN = 1);


// ========== ���ݽṹ�����ⲿ���ã� ==========
struct RealRectMM { double x_left, y_top, x_right, y_bot; };

struct BBoxPx1 {
    int row_min, col_min; // top-left (��)
    int row_max, col_max; // bottom-right (��)
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

// ========== ���껻�� ==========
std::pair<double, double> PixelCenterPx1ToRealMM(const Grid& grid, int row1, int col1);
std::pair<int, int>       RealMMToPixelPx1(const Grid& grid, double x_mm, double y_mm);

// ========== ���ؿ�/���� ==========
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

// ========== ����/��� ==========
std::string getLayerNameById(const Grid& grid, int layerId);
int         getLayerIdByName(const Grid& grid, const std::string& layerName);

// ========== СPCB�߽�/���� ==========
std::vector<std::pair<double, double>>
GetExpandedPCBRealCorners(const Grid& grid, int layer, int blockedAreaID,
    int expandPixels = 5);





// ���� ����������/�ǵ㻥ת��������š����ü� ����

// �ǵ�(LT,RT,RB,LB) -> ����
RealRectMM CornersToRect_LT_RT_RB_LB(const std::vector<std::pair<double, double>>& c);
// ���� -> �ǵ�(LT,RT,RB,LB)
std::vector<std::pair<double, double>> RectToCorners_LT_RT_RB_LB(const RealRectMM& r);

// ��ʵ�����ı߸��� (0.5 + k)*pixel_len_mm
RealRectMM ExpandRect_HalfPlusK(const RealRectMM& in, const Grid& grid, int kWholePixels);

// ��ʵ�򣺰������Ӿ������ü�
RealRectMM ClipRectToBoard(const RealRectMM& in,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY);

// ����A������������������ -> ������+�����ض��� -> ���ü� -> ���ؽǵ�(LT,RT,RB,LB)
std::vector<std::pair<double, double>>
GetExpandedPCBCornersSnapAndClip(const Grid& grid, int layer, int blockedAreaID,
    int expandPixels, int kWholePixels,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY);



// ʾ����������ѡ��
void exampleUsage(const Grid& grid);

// ========== �߶�/���� ==========
std::vector<std::pair<double, double>>
ComputeLineRectIntersections(double x1, double y1, double x2, double y2, const RealRectMM& rect);

bool IsPointInsideRect(double x, double y, const RealRectMM& rect);

bool GetSegmentInfo(const std::shared_ptr<Node>& segment,
    double& x1, double& y1, double& x2, double& y2,
    std::string& layer, std::string& net);

bool IsConnectionPoint(double x, double y,
    const std::vector<std::shared_ptr<Node>>& segments,
    const std::string& targetLayer, double tolerance = 0.001);

// ========== ����/��ֹ���ռ� ==========
std::vector<PointWithNet>
FindIntersectionsAndEndpointsInSmallPCB(const Grid& grid, int layerId, int blockedAreaId,
    int expandPixels = 5);

// ========== ���ģ�� ==========
void GenerateCompletePCBTemplate(
    const std::vector<PointWithNet>& modulePoints,
    const std::vector<std::pair<double, double>>& expandedCorners,
    const std::string& outputFilename);

#endif // SEPARATE6_H
#pragma once



// �ռ��������ڡ���һ�����Ρ��ڣ����߽磩������ segment SEID��ȥ�ز�����
std::vector<int>
CollectSegmentsTouchingFirstRect(const Grid& grid, int layerId, int blockedAreaId, int expandPixels = 5);

// ɾ��������������ռ������� KiCadParser::removeSegmentById ����ɾ��������ɾ������
size_t
DeleteSegmentsTouchingFirstRect(KiCadParser& parser, const Grid& grid, int layerId, int blockedAreaId, int expandPixels = 5);


